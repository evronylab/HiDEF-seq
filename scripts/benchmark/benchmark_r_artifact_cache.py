#!/usr/bin/env python3
"""Compare frozen Python/R cache CLIs on real products in private scratch entries.

The cold operation includes the identical `cp --reflink=never` product builder,
checksum validation, durable publication, and CLI startup. Warm operations include
verification and hard-link restoration. Preparation algorithms are unchanged and
are deliberately not rerun. Fresh measurement workers include waited-for child CPU and report
maximum process RSS, not simultaneous aggregate process-tree memory.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import resource
import subprocess
import sys
import time


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''): digest.update(block)
    return digest.hexdigest()


def product_inventory(path):
    paths = [path] if path.is_file() else sorted(item for item in path.rglob('*') if item.is_file())
    return {('' if item == path else item.relative_to(path).as_posix()): dict(bytes=item.stat().st_size, sha256=sha256(item)) for item in paths}


def input_stat(path):
    paths = [path] if path.is_file() else [path, *sorted(path.rglob('*'))]
    return {str(item): (item.stat().st_size, item.stat().st_mtime_ns, item.stat().st_ctime_ns) for item in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-helper', required=True)
    parser.add_argument('--r-helper', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--product', action='append', required=True, help='LABEL=ABSOLUTE_PATH')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    root = Path(args.out).resolve(); root.mkdir(parents=True, exist_ok=False)
    helpers = dict(python=[sys.executable, str(Path(args.python_helper).resolve())], R=['Rscript', str(Path(args.r_helper).resolve())])
    inputs = {label: Path(path).resolve() for label,path in (item.split('=',1) for item in args.product)}
    sources = [Path(args.python_helper), Path(args.r_helper), Path(args.r_helper).with_name('sharedFunctions.R'), Path(__file__)]
    provenance = {str(path.resolve()): sha256(path) for path in sources}
    before = {label: input_stat(path) for label,path in inputs.items()}
    sizes = {}
    (root/'provenance.json').write_text(json.dumps(dict(sources=provenance, inputs={label:dict(path=str(path),stat=before[label]) for label,path in inputs.items()}, job=os.environ.get('SLURM_JOB_ID')),indent=2)+'\n')
    records = []
    for label,source in inputs.items():
        inventory = product_inventory(source)
        sizes[label] = sum(item['bytes'] for item in inventory.values())
        source_hash = hashlib.sha256(json.dumps(inventory,sort_keys=True,separators=(',', ':')).encode()).hexdigest()
        expected_products = {('product' if not name else 'product/'+name): entry for name,entry in inventory.items()}
        for repeat in range(args.repeats):
            order = ['python','R'] if repeat % 2 == 0 else ['R','python']
            for language in order:
                base=root/f'{label}-{repeat}-{language}'; base.mkdir()
                cache=base/'cache'; cold=base/'cold'; warm=base/'warm'; cold.mkdir(); warm.mkdir()
                identity=dict(schema=1,namespace='cache-benchmark',settings=dict(scenario=label,source_sha256=source_hash))
                identity_file=base/'identity.json'; identity_file.write_text(json.dumps(identity))
                for phase, work, builder in [('cold',cold,['cp','-a','--reflink=never','--',str(source),'product']),('warm',warm,['sh','-c','echo unexpected-builder >&2; exit 99'])]:
                    command=helpers[language]+['run','--root',str(cache),'--identity',str(identity_file),'--product','product','--',*builder]
                    timing=base/f'{phase}.time.json'; stdout=base/f'{phase}.stdout'; stderr=base/f'{phase}.stderr'
                    started=time.monotonic()
                    with stdout.open('w') as out, stderr.open('w') as err:
                        result=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--measure',str(timing),*command],cwd=work,stdout=out,stderr=err)
                    measured=json.loads(timing.read_text())
                    user,system,elapsed,rss,status=[measured[key] for key in ('user','system','elapsed','rss','status')]
                    record=dict(scenario=label,bytes=sizes[label],repeat=repeat,language=language,phase=phase,user_seconds=float(user),system_seconds=float(system),cpu_seconds=float(user)+float(system),elapsed_seconds=float(elapsed),max_rss_kib=int(rss),status=int(status),wrapper_returncode=result.returncode,command=command,outer_elapsed_seconds=time.monotonic()-started)
                    records.append(record)
                    (root/'measurements.json').write_text(json.dumps(records,indent=2)+'\n')
                    if result.returncode: raise RuntimeError(stderr.read_text())
                    output=work/'product'
                    if product_inventory(output)!=inventory: raise AssertionError('Restored/built content mismatch')
                    destination=Path(stdout.read_text().strip())
                    manifest=json.loads((destination/'manifest.complete.json').read_text())
                    if manifest['products']!=expected_products: raise AssertionError('Manifest mismatch')
                    if phase=='warm':
                        for name in inventory:
                            relative=Path('product')/name if name else Path('product')
                            if (work/relative).stat().st_ino!=(destination/relative).stat().st_ino: raise AssertionError('Expected hard-link restoration')
                    shutil.copy2(destination/'manifest.complete.json',base/f'{phase}.manifest.json')
                # Remove only entries created by this benchmark; never source caches.
                shutil.rmtree(cache); shutil.rmtree(cold); shutil.rmtree(warm)
    if provenance!={str(path.resolve()):sha256(path) for path in sources}: raise AssertionError('Benchmark code changed')
    if before!={label:input_stat(path) for label,path in inputs.items()}: raise AssertionError('Source products changed')
    summary=[]
    for label in inputs:
        for phase in ['cold','warm']:
            row=dict(scenario=label,phase=phase,bytes=sizes[label])
            for language in helpers:
                selected=[item for item in records if item['scenario']==label and item['phase']==phase and item['language']==language]
                row[language]={metric:statistics.median(item[metric] for item in selected) for metric in ['cpu_seconds','elapsed_seconds','max_rss_kib']}
            row['R_over_python']={metric:row['R'][metric]/row['python'][metric] for metric in row['R']}
            summary.append(row)
    (root/'result.json').write_text(json.dumps(dict(status='PASS',scope=__doc__,repeats=args.repeats,summary=summary,measurements=records,source_hashes=provenance),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--measure':
        # One fresh wrapper per command makes ru_maxrss independent between
        # languages/repeats instead of inheriting an earlier child's maximum.
        started=time.perf_counter()
        result=subprocess.run(sys.argv[3:])
        usage=resource.getrusage(resource.RUSAGE_CHILDREN)
        Path(sys.argv[2]).write_text(json.dumps(dict(user=usage.ru_utime,system=usage.ru_stime,elapsed=time.perf_counter()-started,rss=usage.ru_maxrss,status=result.returncode))+'\n')
        sys.exit(result.returncode)
    main()
