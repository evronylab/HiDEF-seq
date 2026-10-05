// Experimental component only; not wired into calculateBurdens.R.
// g++ -O3 -std=c++17 annotateCoverage.cpp -o annotateCoverage -lhts
// Read one sorted, non-overlapping, nonzero four-column coverage BED. Cache
// only its current reference chromosome, emit optional five-column per-base
// BED, and accumulate the same observed contexts as the legacy BED/awk path.
#include <htslib/faidx.h>

#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
struct Options {
    std::string bed, fasta, fai, row, counts, bedOutput;
};

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string name = argv[i];
        if (name == "--help") {
            std::cout << "annotateCoverage --bed runs.bed --fasta reference.fa --fai reference.fa.fai "
                         "--row-id 1 --counts 1.reftnc_plus_strand.tsv [--bed-output -|output.bed]\n";
            std::exit(0);
        }
        if (++i == argc) throw std::runtime_error("missing value for " + name);
        const std::string value = argv[i];
        if (name == "--bed") options.bed = value;
        else if (name == "--fasta") options.fasta = value;
        else if (name == "--fai") options.fai = value;
        else if (name == "--row-id") options.row = value;
        else if (name == "--counts") options.counts = value;
        else if (name == "--bed-output") options.bedOutput = value;
        else throw std::runtime_error("unknown option " + name);
    }
    if (options.bed.empty() || options.fasta.empty() || options.fai.empty() ||
        options.row.empty() || options.counts.empty())
        throw std::runtime_error("--bed, --fasta, --fai, --row-id and --counts are required");
    if (options.row.find_first_of("\t\r\n") != std::string::npos)
        throw std::runtime_error("row ID must not contain delimiters");
    return options;
}

int64_t coordinate(const std::string& text) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("BED coordinates must be nonnegative integers");
    size_t used;
    const auto value = std::stoll(text, &used);
    if (used != text.size()) throw std::runtime_error("invalid BED coordinate");
    return value;
}

std::array<std::string, 4> fields(const std::string& line) {
    std::array<std::string, 4> result;
    size_t begin = 0;
    for (size_t index = 0; index < 3; ++index) {
        const auto end = line.find('\t', begin);
        if (end == std::string::npos) throw std::runtime_error("coverage BED must have exactly four tab-separated fields");
        result[index] = line.substr(begin, end - begin);
        begin = end + 1;
    }
    result[3] = line.substr(begin);
    if (result[3].empty() || result[3].find_first_of("\t\r\n") != std::string::npos)
        throw std::runtime_error("coverage BED must have exactly four tab-separated fields");
    return result;
}

unsigned char normalizedBase(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // N/n and every other FASTA symbol become N.
    }
}

std::array<std::string, 126> contexts() {
    std::array<std::string, 126> result;
    const std::string alphabet = "ACGTN";
    for (size_t a = 0; a < 5; ++a)
        for (size_t b = 0; b < 5; ++b)
            for (size_t c = 0; c < 5; ++c)
                result[a * 25 + b * 5 + c] = std::string{alphabet[a], alphabet[b], alphabet[c]};
    result[125] = ".";
    return result;
}

void printAwkNumber(std::ostream& output, double value) {
    // awk prints exactly representable integers in full; its default OFMT for
    // other numeric values is %.6g. Keep that legacy pre-R formatting policy.
    if (std::floor(value) == value && std::fabs(value) < 9223372036854775808.0)
        output << std::fixed << std::setprecision(0) << value;
    else output << std::defaultfloat << std::setprecision(6) << value;
}

void annotate(const Options& options) {
    using FaiPtr = std::unique_ptr<faidx_t, decltype(&fai_destroy)>;
    FaiPtr reference(fai_load3(options.fasta.c_str(), options.fai.c_str(), nullptr, 0), fai_destroy);
    if (!reference) throw std::runtime_error("cannot load the existing FASTA index");
    std::unordered_map<std::string, int> order;
    for (int index = 0; index < faidx_nseq(reference.get()); ++index) {
        const std::string name = faidx_iseq(reference.get(), index);
        // The existing seqkit/awk reference BED splits names at ':' and '-'.
        // Do not silently fix that independent legacy behavior in this helper.
        if (name.find_first_of(":-\t\r\n ") != std::string::npos)
            throw std::runtime_error("legacy-fallback-required: reference contig name contains a legacy awk delimiter: " + name);
        order.emplace(name, index);
    }
    std::ifstream input(options.bed);
    if (!input) throw std::runtime_error("cannot open coverage BED");
    std::ofstream bedFile;
    std::ostream* bedOutput = nullptr;
    if (options.bedOutput == "-") bedOutput = &std::cout;
    else if (!options.bedOutput.empty()) {
        bedFile.open(options.bedOutput);
        if (!bedFile) throw std::runtime_error("cannot create BED output");
        bedOutput = &bedFile;
    }
    const auto names = contexts();
    std::array<double, 126> sums {};
    std::array<bool, 126> seen {};
    using SequencePtr = std::unique_ptr<char, decltype(&std::free)>;
    SequencePtr sequence(nullptr, std::free);
    std::string chromosome;
    hts_pos_t sequenceLength = 0;
    int chromosomeRank = -1;
    int64_t previousEnd = 0;
    bool any = false;
    std::string line;
    while (std::getline(input, line)) {
        const auto row = fields(line);
        const auto found = order.find(row[0]);
        if (found == order.end()) throw std::runtime_error("coverage contig absent from FASTA: " + row[0]);
        const int64_t start = coordinate(row[1]);
        const int64_t end = coordinate(row[2]);
        char* depthEnd = nullptr;
        errno = 0;
        const double depth = std::strtod(row[3].c_str(), &depthEnd);
        if (errno || depthEnd == row[3].c_str() || *depthEnd || !std::isfinite(depth) || depth <= 0)
            throw std::runtime_error("coverage must be a finite positive number; zero rows are not accepted");
        if (row[0] != chromosome) {
            if (found->second <= chromosomeRank)
                throw std::runtime_error("coverage BED must follow FASTA contig order");
            chromosome = row[0];
            chromosomeRank = found->second;
            previousEnd = 0;
            const hts_pos_t expectedLength = faidx_seq_len64(reference.get(), chromosome.c_str());
            if (expectedLength < 1) throw std::runtime_error("reference contig is empty");
            // Release the previous chromosome before allocating the next.
            sequence.reset();
            sequence.reset(faidx_fetch_seq64(reference.get(), chromosome.c_str(), 0,
                                             expectedLength - 1, &sequenceLength));
            if (!sequence || sequenceLength != expectedLength)
                throw std::runtime_error("reference sequence length does not match its index");
            for (hts_pos_t index = 0; index < sequenceLength; ++index)
                sequence.get()[index] = static_cast<char>(normalizedBase(sequence.get()[index]));
        }
        if (start < previousEnd || end <= start || end > sequenceLength)
            throw std::runtime_error("coverage intervals must be non-overlapping and inside their reference contig");
        previousEnd = end;
        for (int64_t position = start; position < end; ++position) {
            size_t context = 125;
            if (position > 0 && position + 1 < sequenceLength) {
                const auto* bases = reinterpret_cast<const unsigned char*>(sequence.get());
                context = bases[position - 1] * 25 + bases[position] * 5 + bases[position + 1];
            }
            seen[context] = true;
            sums[context] += depth;  // Same per-base addition order as legacy awk.
            if (bedOutput)
                *bedOutput << chromosome << '\t' << position << '\t' << position + 1
                           << '\t' << row[3] << '\t' << names[context] << '\n';
        }
        any = true;
    }
    if (!input.eof()) throw std::runtime_error("error reading coverage BED");
    if (bedOutput) {
        bedOutput->flush();
        if (!*bedOutput) throw std::runtime_error("error writing annotated BED");
    }
    std::ofstream counts(options.counts);
    if (!counts) throw std::runtime_error("cannot create context count output");
    if (!any) counts << options.row << "\tNA\t0\n";
    else {
        for (size_t context = 0; context < seen.size(); ++context) {
            if (!seen[context]) continue;
            counts << options.row << '\t' << names[context] << '\t';
            printAwkNumber(counts, sums[context]);
            counts << '\n';
        }
    }
    counts.flush();
    if (!counts) throw std::runtime_error("error writing context count output");
}
}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try { annotate(parse(argc, argv)); }
    catch (const std::exception& error) {
        std::cerr << "annotateCoverage: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
