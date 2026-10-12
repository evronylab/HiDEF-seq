#!/usr/bin/env python3
"""Check production demultiplex channel matching without scientific workers.

Generate a fixture, run fixture.nf with Nextflow on an allocated node, then
verify. The fixture uses the production matching/grouping fragments, with
synthetic demultiplex outputs and merged paths; it does not execute Lima.
"""
import argparse
import hashlib
import json
from pathlib import Path


def generate(directory, repo):
    directory.mkdir(parents=True, exist_ok=False)
    source = (repo / "main.nf").read_text()
    first_start = source.index("  mergeDemuxBams_round1_input_ch = ")
    first_end = source.index("\n\n", first_start)
    second_start = source.index("  limaDemux_round2_samples_ch = ")
    second_end = source.index("\n  limaDemux_round2 = ", second_start)
    fragments = source[first_start:first_end] + "\n" + r'''
  mergeDemuxBams_round1 = mergeDemuxBams_round1_input_ch.map { run, individual, sample, outer, bams ->
    tuple(run, individual, sample, outer, file("${run}.${sample}.${outer}.merged.bam"), file("${run}.${sample}.${outer}.merged.bam.pbi"))
  }
''' + source[second_start:second_end]
    configs, observed_bams, expected_first, expected_second = [], [], [], []
    def case(run, sample, outer, inner=None, orientations=None):
        outer_ids = outer.split("-")
        canonical = "-".join(sorted(outer_ids))
        configs.append(dict(run_id=run, individual_id="ind", sample_id=sample,
                            barcode_ids=outer, barcode_ids_parsed=dict(canonical=canonical),
                            barcode_ids_round2=inner, round2_enabled=inner is not None,
                            barcode_ids_round2_parsed=None if inner is None else dict(mode="same", canonical=inner)))
        names = orientations or [outer_ids[0] + "--" + outer_ids[-1]]
        paths = sorted(run + ".demux." + name + ".bam" for name in names)
        for path in reversed(paths):
            row = [run, canonical, path]
            if row not in observed_bams:
                observed_bams.append(row)
        row = [run, "ind", sample, outer, paths]
        if row not in expected_first:
            expected_first.append(row)
        if inner is not None:
            expected_second.append([run, "ind", sample, outer, inner, inner,
                                    run + "." + sample + "." + outer + ".merged.bam",
                                    run + "." + sample + "." + outer + ".merged.bam.pbi",
                                    run + ".round2.fasta", "same", ""])
    # Two different samples legitimately share their first-round barcode.
    case("shared", "sampleA", "outer", "innerA")
    case("shared", "sampleB", "outer", "innerB")
    # One sample can occur twice in a run, with distinct second-round barcodes.
    case("repeated", "sampleC", "outer", "innerA")
    case("repeated", "sampleC", "outer", "innerB")
    # --different can produce both orientations for one canonical pair.
    case("orientations", "sampleD", "right-left", "innerA", ["left--right", "right--left"])
    # Existing single-round same and different modes retain all input BAMs.
    case("single", "sampleE", "outer")
    case("single_pair", "sampleF", "right-left", orientations=["left--right", "right--left"])
    spec = dict(configs=configs, bams=observed_bams,
                expected_first=sorted(expected_first), expected_second=sorted(expected_second))
    (directory / "spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    prefix = r'''nextflow.enable.dsl=2
params.fixture_spec = null
params.fixture_output = null
workflow {
  def spec = new groovy.json.JsonSlurper().parseText(file(params.fixture_spec).text)
  run_sample_configs = spec.configs
  params.lima_round2_supplemental_settings = ''
  limaDemux_round1_map_ch = channel.fromList(spec.bams).map { run, key, bam -> tuple(run, key, file(bam)) }
  round2_barcodes_fasta_ch = channel.fromList(spec.configs.findAll { it.round2_enabled }.collect { it.run_id }.unique())
    .map { run -> tuple(run, file("${run}.round2.fasta")) }
'''
    suffix = r'''
  first = mergeDemuxBams_round1_input_ch.map { run, individual, sample, outer, bams ->
    tuple('first', [run, individual, sample, outer, bams.collect { it.name }])
  }
  second = limaDemux_round2_input_ch.map { run, individual, sample, outer, inner, innerkey, bam, pbi, fasta, mode, flags ->
    tuple('second', [run, individual, sample, outer, inner, innerkey, bam.name, pbi.name, fasta.name, mode, flags])
  }
  first.mix(second).collect(flat: false).view { rows ->
    file(params.fixture_output).text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(rows)) + '\n'
    "Recorded ${rows.size()} demultiplex channel groups"
  }
}
'''
    (directory / "fixture.nf").write_text(prefix + fragments + suffix)
    (directory / "source.json").write_text(json.dumps(dict(main_sha256=hashlib.sha256(source.encode()).hexdigest()), indent=2) + "\n")


def verify(directory):
    spec = json.loads((directory / "spec.json").read_text())
    rows = json.loads((directory / "observed.json").read_text())
    for kind in ("first", "second"):
        actual = sorted(row[1] for row in rows if row[0] == kind)
        expected = spec["expected_" + kind]
        assert actual == expected, json.dumps(dict(kind=kind, expected=expected, actual=actual), indent=2)
    print("Exact first-round groups and second-round fan-out passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("generate", "verify"))
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    if args.mode == "generate":
        generate(args.directory.resolve(), args.repo.resolve())
    else:
        verify(args.directory.resolve())


if __name__ == "__main__":
    main()
