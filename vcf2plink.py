#!/usr/bin/env python3
"""
   vcf2plink.py - converts a VCF file to plink format
   Copyright (C) 2016 Giulio Genovese

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Written by Giulio Genovese <giulio.genovese@gmail.com>
"""

import argparse, os, sys
from subprocess import call, Popen, PIPE

parser = argparse.ArgumentParser(description = 'vcf2plink.py: converts a VCF file to plink format (Oct 3rd 2016)', add_help = False, usage = 'vcf2plink.py [options]')
parser.add_argument('--vcf', metavar = '<in.vcf.gz>', type = str, default = '-', help = 'Specify a VCF file to be converted')
parser.add_argument('--ref', metavar = '<file>', type = str, required = True, help = 'reference sequence (e.g. human_g1k_v37.fasta)')
parser.add_argument('--filter', metavar = '<expr>', type = str, help = 'exclude sites for which the expression is true (e.g. "FORMAT/DP<10 || FORMAT/GQ<20") (see bcftools manual for details)')
parser.add_argument('--set-GTs', metavar = '<char>', type = str, default = '.', help = 'set genotypes of failed samples to missing (.) or ref (0)')
parser.add_argument('--out', metavar = '[prefix]', type = str, default = 'plink', help = 'Specify prefix for output files')
parser.add_argument('--mem', metavar = '<int>', type = int, help = 'main workspace size, in GB')
parser.add_argument('--impute-sex', metavar = '{female max F} {male min F}', type = float, nargs = '*', help = 'Impute sex of samples from X chromosome inbreeding coefficients [.4 .4]')
parser.add_argument('--pdf', metavar = '[filename]', type = str, help = 'Generate F statistic histogram for sex imputation')

try:
  parser.error = parser.exit
  args = parser.parse_args()
  if args.impute_sex == []:
    args.impute_sex = [.4, .4]
except SystemExit:
  parser.print_help()
  exit(2)

# make sure the ouput directory exists
outdir = os.path.dirname(args.out)
if outdir != '' and not os.path.isdir(outdir):
  os.makedirs(outdir)

# invoke bcftools to preprocess VCF file
sys.stderr.write('vcf2plink.py: Convert VCF file to plink\n')
if args.vcf == '-':
  p1 = Popen(['bcftools', 'norm', '-Ou', '-m', '-any'], stdin = sys.stdin.buffer, stdout = PIPE)
else:
  p1 = Popen(['bcftools', 'norm', '-Ou', '-m', '-any', args.vcf], stdout = PIPE)

p2 = Popen(['bcftools', 'norm', '-Ou', '-f', args.ref], stdin = p1.stdout, stdout = PIPE)

if args.filter:
  p2a = Popen(['bcftools', 'filter', '-Ou', '-e', args.filter, '--set-GTs', args.set_GTs], stdin = p2.stdout, stdout = PIPE)

p3 = Popen(['bcftools', 'annotate', '-Ob', '-x', 'ID', '-I', '+%CHROM:%POS:%REF:%ALT'], stdin = p2a.stdout if args.filter else p2.stdout, stdout = PIPE)


# invoke plink to convert VCF file
if call(['plink', '--make-bed',
  '--keep-allele-order',
  '--bcf', '/dev/stdin',
  '--vcf-idspace-to', '_',
  '--const-fid',
  '--allow-extra-chr', '0',
  '--split-x', 'b37', 'no-fail',
  '--out', args.out] + (['--memory', str(1024 * args.mem)] if args.mem else []), stdin = p3.stdout):
    raise Exception('vcf2plink.py: Problems converting VCF file to plink')

# impute sex using chromosome X calls
if args.impute_sex:
  if call(['plink', '--make-bed',
  '--bfile', args.out,
  '--keep-allele-order',
  '--impute-sex', 'ycount', str(args.impute_sex[0]), str(args.impute_sex[1]), '10000', '0',
  '--out', args.out]):
    raise Exception('vcf2plink.py: Problems imputing sex')
  if call(['rm', args.out + '.bed~', args.out + '.bim~', args.out + '.fam~']):
    raise Exception('vcf2plink.py: Problems removing temporary files')
  if args.pdf:
    import pandas as pd
    from matplotlib.backends.backend_pdf import PdfPages
    df = pd.read_csv(args.out + '.sexcheck', delim_whitespace = True)
    with PdfPages(args.pdf) as pdf:
      ax = df.plot(kind = 'scatter', x = 'F', y = 'YCOUNT', alpha = 1/2)
      ax.set_xlabel('Chromosome X inbreeding coefficient (F)')
      ax.set_ylabel('Chromosome Y non-missing genotypes (YCOUNT)')
      pdf.savefig(ax.get_figure())
