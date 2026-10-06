#!/usr/bin/env python3
"""Accept pinned full BAM diagnostics under an explicitly approved tie-order rule.

This reads existing evidence; it neither rescans BAMs nor changes the original
strict-order comparison. Its plan must pin every report/source and current input
identity. The diagnostic's raw-record comparison is stricter than logical-record
equality: only permutations within identical (reference, position) groups differ.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re


RULE = ('BAMs must contain exactly the same logical records with the same '
        'multiplicities and coordinates, and coordinate ordering must be '
        'preserved. Records whose genomic sort coordinates are identical may '
        'occur in a different relative order.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path):
    stat = path.stat()
    return {'resolved': str(path.resolve()), 'inode': stat.st_ino,
            'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
            'ctime_ns': stat.st_ctime_ns}


def validate_reports(original, record, indices, completion):
    """Validate scientific evidence independently of filesystem bindings."""
    require(original['kind'] == 'bam' and original['status'] == 'failure',
            'Expected original strict-order BAM failure')
    require(original['records']['status'] == 'failure' and
            original['records']['reason'] == 'ordered record bytes differ',
            'Original failure is not the known ordered-record difference')
    require(not original.get('input_changed_during_validation'), 'Original input changed')
    require(original['headers']['status'] == 'pass' and
            original['reference_dictionary']['status'] == 'pass',
            'Original header or dictionary comparison did not pass')
    require(len(original['quickcheck']) == 2 and all(
        x['returncode'] == 0 and not x['stderr'] for x in original['quickcheck']),
        'Original BAM quickchecks did not pass')
    require(record['headers'] == original['headers'] and
            record['reference_dictionary_exact'] is True,
            'Diagnostic header/dictionary differs from original evidence')
    totals = record['totals']
    for side in ('reference', 'candidate'):
        key = side + '_records'
        require(totals[key] == original['records'][key] and totals[key] > 0,
                'Diagnostic record count mismatch')
    require(totals['reference_records'] == totals['candidate_records'],
            'Reference/candidate record counts differ')
    require(totals['different_group_coordinates'] == 0,
            'Coordinates or coordinate-group ordering differ')
    require(totals['different_record_multisets'] == 0,
            'Record contents or duplicate multiplicities differ')
    require(0 < totals['groups_with_order_only_difference'] <= totals['groups'] <=
            totals['reference_records'], 'Invalid coordinate-group totals')
    require(completion['diagnostic_complete'] is True and
            completion['inputs_unchanged'] is True, 'Incomplete or changed-input diagnostic')
    require(set(indices) == {'reference', 'candidate'} and
            set(completion['index_statuses']) == set(indices), 'Missing index side')
    for side, checks in indices.items():
        require(len(checks) == 2 and {x['index'] for x in checks} == {'.bai', '.pbi'},
                'Missing or duplicate BAI/PBI check')
        require(completion['index_statuses'][side] == [x['status'] for x in checks],
                'Completion and index statuses disagree')
        for check in checks:
            require(check['status'] == 'pass' and check['execution']['returncode'] == 0
                    and not check['execution']['stderr'], 'Index rebuild did not pass')
            expected = check['original_payload_sha256']
            require(re.fullmatch('[0-9a-f]{64}', expected) is not None and
                    check['rebuilt_payload_sha256'] == expected,
                    'Index payload differs from its own-file rebuild')
    return totals


def accept(plan_path):
    plan_path = Path(plan_path).resolve()
    plan = json.loads(plan_path.read_text())
    require(plan['approval']['rule'] == RULE and
            plan['approval']['source'] == 'explicit user message' and
            plan['approval']['date'] == '2026-10-06', 'Missing approved rule')
    root = Path(plan['evidence_root'])
    pins = plan['evidence_sha256']
    for path, expected in pins.items():
        require(sha(root / path) == expected, 'Evidence checksum mismatch: ' + path)

    def read(path):
        require(path in pins, 'Unpinned report: ' + path)
        return json.loads((root / path).read_text())

    results = read(plan['original_results'])
    summary = read(plan['original_summary'])
    manifest = read(plan['original_manifest'])
    audit = read(plan['publication_source_audit'])
    require(audit['status'] == 'pass', 'Original producer audit did not pass')
    audit_rows = {x['path']: x for x in audit['checks']}
    require(len(audit_rows) == len(audit['checks']), 'Duplicate producer audit path')
    require(summary['scope'] == 'full-publication' and summary['status'] == 'failure',
            'Original failed full-publication report must remain preserved')
    require(summary['counts'] == dict(Counter(x['status'] for x in results)),
            'Original result counts disagree')
    rows = {x['path']: x for x in results}
    require(len(rows) == len(results), 'Duplicate original result path')
    inventory_paths = {x['path'] for x in manifest['inventory']}
    require(len(inventory_paths) == len(manifest['inventory']) and
            set(rows) == inventory_paths | {'@source-config', '@publication-source-audit'},
            'Original full result path set differs from manifest')
    required = {x['path']: x['kind'] for x in manifest['required_publications']['required']}
    require(len(required) == summary['required_publication_count'] == 711,
            'Required publication inventory changed')
    require(len(results) == len(manifest['inventory']) + 2, 'Original results incomplete')
    bam_paths = {x['path'] for x in plan['diagnoses']}
    require(len(bam_paths) == len(plan['diagnoses']) == 2 and
            bam_paths == {p for p, k in required.items() if k == 'bam'},
            'Supplement must bind exactly the two required BAMs')
    require({x['path'] for x in results if x['status'] in ('failure', 'review')} == bam_paths,
            'Unresolved failure/review outside the two diagnosed BAMs')
    accepted = []
    for item in plan['diagnoses']:
        path = item['path']
        original = rows[path]
        directory = item['directory']
        failure = read(directory + '/failed-comparison.json')
        require(failure == original, 'Diagnostic binds a different original BAM result')
        record = read(directory + '/record-diagnosis.json')
        indices = read(directory + '/index-diagnosis.json')
        completion = read(directory + '/completion.json')
        require(record['job_id'] == item['job_id'], 'Unexpected diagnostic job')
        start_ns = min((root / directory / name).stat().st_mtime_ns for name in
                       ('compare_bam.py', 'diagnose.py', 'failed-comparison.json'))
        require(start_ns == item['diagnostic_start_ns'], 'Diagnostic start snapshot changed')
        for label, filename in (('helper', 'compare_bam.py'), ('diagnostic', 'diagnose.py')):
            source = directory + '/' + filename
            require(source in pins and record['source_sha256'][label] == pins[source],
                    'Diagnostic source pin mismatch')
        totals = validate_reports(original, record, indices, completion)
        for side in ('reference', 'candidate'):
            bam = Path(original[side])
            require(bam == Path(manifest[side]) / path, 'BAM is outside original manifest root')
            link = root / directory / side / 'input.bam'
            require(link.is_symlink() and link.resolve() == bam.resolve(),
                    'Index scratch symlink targets a different BAM')
            for suffix in ('', '.bai', '.pbi'):
                source = Path(str(bam) + suffix)
                require(str(source) in plan['input_identities'], 'Missing input identity')
                current = identity(source)
                require(current == plan['input_identities'][str(source)], 'Input identity changed')
                if side == 'candidate':
                    producer = audit_rows[path + suffix]
                    require(producer['status'] == 'pass' and
                            producer['publication_stat'][1:] == [current[k] for k in
                                ('inode', 'size', 'mtime_ns', 'ctime_ns')],
                            'Candidate identity differs from frozen producer audit')
                require(max(current['mtime_ns'], current['ctime_ns']) <= start_ns,
                        'Input was modified after diagnostic startup')
            for check in indices[side]:
                rebuilt = Path(str(link) + check['index'])
                require(check['rebuilt'] == str(rebuilt) and rebuilt.is_file(),
                        'Unexpected/missing rebuilt index')
                tool = original['tools']['samtools' if check['index'] == '.bai' else 'pbindex']
                expected_command = ([tool, 'index', str(link), str(rebuilt)]
                                    if check['index'] == '.bai' else [tool, str(link)])
                require(check['execution']['command'] == expected_command,
                        'Index command used a different BAM or output')
        accepted.append({'path': path, 'kind': 'bam', 'status': 'pass',
                         'diagnostic': directory, 'job_id': item['job_id'], 'totals': totals})
        for suffix in ('.bai', '.pbi'):
            index_path = path + suffix
            companion = rows[index_path]
            require(required[index_path] == companion['kind'] == 'bam_index' and
                    companion['status'] == 'covered' and companion['owner'] == path,
                    'BAM-index companion ownership differs')
            accepted.append({'path': index_path, 'kind': 'bam_index', 'status': 'pass',
                             'owner': path, 'reason': 'Both own-file index rebuilds passed'})
    supplementary = {x['path'] for x in accepted}
    for path, kind in required.items():
        require(rows[path]['kind'] == kind, 'Required output kind changed')
        require(path in supplementary or rows[path]['status'] == 'pass',
                'Required output lacks accepted evidence')
    require(rows['@source-config']['status'] == rows['@publication-source-audit']['status'] == 'pass',
            'Configuration or producer provenance did not pass')
    return {'schema_version': 1, 'status': 'pass', 'scope': 'supplemental-bam-coordinate-order',
            'approved_rule': RULE, 'approval': plan['approval'], 'accepted': accepted,
            'original_report_status': summary['status'], 'original_report_modified': False,
            'required_publications_accepted': len(required),
            'non_bam_original_passes': len(required) - len(accepted),
            'plan': str(plan_path), 'plan_sha256': sha(plan_path),
            'gate_sha256': sha(Path(__file__)), 'evidence_sha256': pins}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = accept(args.plan)
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'accepted_paths': len(result['accepted']),
                      'original_report_status': result['original_report_status']}))


if __name__ == '__main__':
    main()
