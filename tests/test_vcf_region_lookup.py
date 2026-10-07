#!/usr/bin/env python3
"""Exercise real bcftools target selection using small deterministic VCFs."""
import argparse
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('repo', nargs='?', default='.')
args = parser.parse_args()
repo = Path(args.repo).resolve(strict=True)
with tempfile.TemporaryDirectory(prefix='hidef-vcf-targets-') as temporary:
    out = Path(temporary)
    seq='ACGT'*12500
    with (out/'reference.fa').open('x') as f:
        for chrom in ('chrA','chrB'):
            f.write('>'+chrom+'\n')
            for i in range(0,len(seq),60):f.write(seq[i:i+60]+'\n')
    subprocess.run(['samtools','faidx',str(out/'reference.fa')],check=True)
    for gt in (True,False):
        path=out/('with_gt.vcf' if gt else 'without_gt.vcf')
        with path.open('x') as f:
            f.write('##fileformat=VCFv4.2\n')
            for chrom in ('chrA','chrB'):f.write(f'##contig=<ID={chrom},length=50000>\n')
            f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
            if gt:f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
            for chrom in ('chrA','chrB'):
                for pos in (100,101,500,510,700,900,1000,1100,2101,3101,5000,7300):
                    ref=seq[pos-1]; alt=next(b for b in 'ACGT' if b!=ref)
                    ad='20,5';genotype='0/1'
                    if pos==500:ref=seq[pos-1:pos+4];alt=ref[0]
                    if pos==700:alt=ref+'GG'
                    if pos==5000:ref=seq[pos-1:7205];alt=ref[0]
                    if pos==900:
                        alt=','.join(b for b in 'ACGT' if b!=ref);ad='20,5,3,0';genotype='1/2'
                    if pos==1000:
                        ref=seq[pos-1:pos+2]
                        alt=''.join(next(b for b in 'ACGT' if b!=c) for c in ref)
                    fmt='GT:GQ:AD' if gt else 'AD'
                    sample=f'{genotype}:30:{ad}' if gt else ad
                    f.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t12.3456\t.\t.\t{fmt}\t{sample}\n')
                    if pos==100:
                        second=f'{genotype}:30:20,6' if gt else '20,6'
                        f.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t24.1234\t.\t.\t{fmt}\t{second}\n')
        with Path(str(path)+'.gz').open('xb') as f:
            subprocess.run(['bgzip','-c',str(path)],stdout=f,check=True)
        subprocess.run(['tabix','-p','vcf',str(path)+'.gz'],check=True)
    subprocess.run(['Rscript', '--vanilla', str(Path(__file__).with_suffix('.R')),
                    str(repo), str(out)], check=True)
