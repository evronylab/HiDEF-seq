// Build with the HTSlib already supplied by the pinned HiDEF-seq container:
// g++ -O2 -std=c++17 splitBamByZmw.cpp -o splitBamByZmw -lhts
// A single reader dispatches unchanged BAM records to the exact legacy ID
// partitions. All compressed streams share one bounded HTSlib thread pool.
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
struct Options {
    std::string input, ids, prefix;
    size_t chunks = 0;
    int threads = 1;
    size_t maxWriters = 128;
};

long long integer(const std::string& text, const std::string& label) {
    size_t used = 0;
    long long result;
    try { result = std::stoll(text, &used); }
    catch (const std::exception&) { throw std::runtime_error("invalid " + label + ": " + text); }
    if (used != text.size()) throw std::runtime_error("invalid " + label + ": " + text);
    return result;
}

Options parse(int argc, char** argv) {
    Options result;
    for (int i = 1; i < argc; ++i) {
        const std::string name = argv[i];
        if (name == "--help") {
            std::cout << "splitBamByZmw --input input.bam --ids ordered.zmwIDs.txt "
                         "--chunks N --output-prefix sample [--threads 1] [--max-open-writers 128]\n";
            std::exit(0);
        }
        if (++i == argc) throw std::runtime_error("missing value for " + name);
        const std::string value = argv[i];
        if (name == "--input") result.input = value;
        else if (name == "--ids") result.ids = value;
        else if (name == "--output-prefix") result.prefix = value;
        else if (name == "--chunks" || name == "--threads" || name == "--max-open-writers") {
            const auto number = integer(value, name);
            if (number < 1) throw std::runtime_error(name + " must be positive");
            if (name != "--chunks" && number > 4096)
                throw std::runtime_error(name + " must not exceed 4096");
            if (name == "--chunks") result.chunks = static_cast<size_t>(number);
            else if (name == "--threads") result.threads = static_cast<int>(number);
            else result.maxWriters = static_cast<size_t>(number);
        } else throw std::runtime_error("unknown option " + name);
    }
    if (result.input.empty() || result.ids.empty() || result.prefix.empty() || !result.chunks)
        throw std::runtime_error("--input, --ids, --chunks and --output-prefix are required");
    return result;
}

using Destinations = std::unordered_map<int32_t, std::vector<size_t>>;

Destinations assignments(const Options& options) {
    std::ifstream input(options.ids);
    if (!input) throw std::runtime_error("cannot open ID enumeration " + options.ids);
    std::vector<int32_t> ids;
    std::string line;
    while (std::getline(input, line)) {
        const auto begin = line.find_first_not_of(" \t\r");
        const auto end = line.find_last_not_of(" \t\r");
        if (begin == std::string::npos) throw std::runtime_error("empty ID enumeration line");
        const auto value = integer(line.substr(begin, end - begin + 1), "ZMW ID");
        if (value < 0 || value > std::numeric_limits<int32_t>::max())
            throw std::runtime_error("ZMW ID outside nonnegative int32 range");
        ids.push_back(static_cast<int32_t>(value));
    }
    if (!input.eof()) throw std::runtime_error("error reading ID enumeration");
    if (ids.size() < options.chunks)
        throw std::runtime_error("ID count is smaller than chunk count; legacy chunks would be empty");
    const size_t quotient = ids.size() / options.chunks;
    const size_t remainder = ids.size() % options.chunks;
    Destinations result;
    result.reserve(ids.size());
    size_t offset = 0;
    for (size_t chunk = 0; chunk < options.chunks; ++chunk) {
        const size_t size = quotient + (chunk < remainder);
        for (size_t index = 0; index < size; ++index) {
            auto& destinations = result[ids[offset++]];
            // Legacy include treats the file as a set within each chunk. The
            // same numeric ID in different chunks intentionally duplicates all
            // matching records across those chunks, including run collisions.
            if (destinations.empty() || destinations.back() != chunk)
                destinations.push_back(chunk);
        }
    }
    return result;
}

int32_t holeNumber(const bam1_t* record) {
    const uint8_t* tag = bam_aux_get(record, "zm");
    // The pinned pbindex/zmwfilter pair indexes a missing zm tag as zero,
    // including when the QNAME contains a different hole number. Match that
    // observed include/show-all behavior rather than parsing a QNAME fallback.
    if (!tag) return 0;
    switch (*tag) {
        case 'c': case 'C': case 's': case 'S': case 'i': case 'I': break;
        default: throw std::runtime_error("zm tag must have an integer BAM type");
    }
    errno = 0;
    const int64_t value = bam_aux2i(tag);
    if (errno || value < 0 || value > std::numeric_limits<int32_t>::max())
        throw std::runtime_error("invalid numeric zm tag");
    return static_cast<int32_t>(value);
}

struct Streams {
    htsThreadPool pool {nullptr, 0};
    samFile* reader = nullptr;
    std::vector<samFile*> writers;
    ~Streams() {
        for (auto* writer : writers) if (writer) sam_close(writer);
        if (reader) sam_close(reader);
        if (pool.pool) hts_tpool_destroy(pool.pool);
    }
};

void splitGroup(const Options& options, const Destinations& destinations,
                const std::vector<std::string>& paths, size_t first, size_t last) {
    Streams streams;
    streams.pool.pool = hts_tpool_init(options.threads);
    streams.pool.qsize = 2;  // Bounded outstanding jobs for each stream.
    if (!streams.pool.pool) throw std::runtime_error("cannot create shared HTSlib thread pool");
    streams.reader = sam_open(options.input.c_str(), "rb");
    if (!streams.reader) throw std::runtime_error("cannot open input BAM");
    if (hts_get_format(streams.reader)->format != bam)
        throw std::runtime_error("input must be BAM");
    const int eofStatus = hts_check_EOF(streams.reader);
    if (eofStatus == 0) throw std::runtime_error("input BAM is missing its BGZF EOF marker");
    if (eofStatus < 0) throw std::runtime_error("cannot verify input BAM EOF marker");
    if (hts_set_thread_pool(streams.reader, &streams.pool) < 0)
        throw std::runtime_error("cannot assign reader thread pool");
    const std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> header(sam_hdr_read(streams.reader), sam_hdr_destroy);
    if (!header) throw std::runtime_error("cannot read BAM header");
    for (size_t chunk = first; chunk < last; ++chunk) {
        const auto& path = paths[chunk];
        samFile* writer = sam_open(path.c_str(), "wb");
        if (!writer) throw std::runtime_error("cannot open output " + path);
        streams.writers.push_back(writer);
        if (hts_set_thread_pool(writer, &streams.pool) < 0 || sam_hdr_write(writer, header.get()) < 0)
            throw std::runtime_error("cannot initialize output " + path);
    }
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    if (!record) throw std::runtime_error("cannot allocate BAM record");
    std::vector<uint64_t> counts(last - first, 0);
    int status;
    while ((status = sam_read1(streams.reader, header.get(), record.get())) >= 0) {
        const auto match = destinations.find(holeNumber(record.get()));
        if (match == destinations.end())
            throw std::runtime_error(std::string("BAM ZMW is missing from ID enumeration: ") + bam_get_qname(record.get()));
        for (const auto chunk : match->second) {
            if (chunk < first || chunk >= last) continue;
            if (sam_write1(streams.writers[chunk - first], header.get(), record.get()) < 0)
                throw std::runtime_error("error writing " + paths[chunk]);
            ++counts[chunk - first];
        }
    }
    if (status < -1) throw std::runtime_error("error reading BAM records");
    for (size_t chunk = first; chunk < last; ++chunk) {
        const int result = sam_close(streams.writers[chunk - first]);
        streams.writers[chunk - first] = nullptr;
        if (result < 0) throw std::runtime_error("error finalizing " + paths[chunk]);
        if (!counts[chunk - first]) throw std::runtime_error("output chunk contains no BAM records");
        std::cout << (chunk + 1) << '\t' << counts[chunk - first] << '\t' << paths[chunk] << '\n';
    }
}

void split(const Options& options) {
    const auto destinations = assignments(options);
    std::vector<std::string> paths;
    for (size_t chunk = 0; chunk < options.chunks; ++chunk) {
        paths.push_back(options.prefix + ".chunk" + std::to_string(chunk + 1) + ".bam");
        if (std::ifstream(paths.back()).good())
            throw std::runtime_error("refusing to overwrite existing output " + paths.back());
    }
    // Large user-requested chunk counts remain supported. Each bounded group
    // requires one sequential BAM pass, without changing global chunk IDs,
    // enumeration partitions, or the membership/order of any output chunk.
    for (size_t first = 0; first < options.chunks; first += options.maxWriters) {
        splitGroup(options, destinations, paths, first,
                   std::min(options.chunks, first + options.maxWriters));
    }
}
}  // namespace

int main(int argc, char** argv) {
    try { split(parse(argc, argv)); }
    catch (const std::exception& error) {
        std::cerr << "splitBamByZmw: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
