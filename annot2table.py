#!/usr/bin/env python3
"""
   annot2table.py - convert VCF annotation file to table format
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

import argparse, sys, pandas as pd
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description = 'annot2table.py: convert VCF annotation file to table format (Oct 3rd 2016)', add_help = False, usage = 'annot2table.py [options]')
parser.add_argument('-jar', metavar = '[jar]', type = str, required = True, help = 'SnpSift jar file')
parser.add_argument('-mem', metavar = '[int]', type = int, default = 4, help = 'Memory for java virtual machine in GB [4]')
parser.add_argument('-vcf', metavar = 'FILE', type = str, default = '-', help = 'Input VCF file [stdin]')
parser.add_argument('-x', metavar = 'FILE', type = str, help = 'Table with extra columns')
parser.add_argument('-c', metavar = 'str', nargs='+', type = str, help = 'Columns to add')
parser.add_argument('-s', metavar = 'INT', type = int, default = 65536, help = 'Number of lines to process at each iteration [65536]')
parser.add_argument('-fld', metavar = 'str', nargs='+', type = str, default = ['FILTER', 'RS', 'CPG', 'AC', 'AN', 'KGPhase3', 'EXACAC', 'NONPSYCHAC', 'dbNSFP_COSMIC_CNT', 'dbNSFP_COSMIC_ID'], help = 'Fields to retain [FILTER RS CPG AC AN KGPhase3 EXACAC NONPSYCHAC dbNSFP_COSMIC_CNT dbNSFP_COSMIC_ID]')
parser.add_argument('-ann', metavar = 'str', nargs='+', type = str, default = ['EFFECT', 'IMPACT', 'GENE', 'FEATUREID', 'HGVS_C', 'HGVS_P'], help = 'ANN fields to retain [EFFECT IMPACT GENE FEATUREID HGVS_C HGVS_P] (see http://snpeff.sourceforge.net/SnpEff_manual.html#ann)')
parser.add_argument('-i', metavar = '[FILE]', type = str, help = 'File with variants to include')
parser.add_argument('-e', metavar = '[FILE]', type = str, help = 'File with variants to exclude')
parser.add_argument('-f', action = 'store_true', default = False, help = 'Whether to filter on the filter column')
parser.add_argument('-d', metavar = '[FILE]', type = str, help = 'File with variants classified as damaging')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

# create a dictionary with the xdf columns to add to the matrix
if args.x and args.c:
  try:
    xdf = pd.read_csv(args.x, sep = '\t', compression = 'gzip')
  except:
    xdf = pd.read_csv(args.x, sep = '\t')
  if not 'SNP' in xdf:
    xdf['SNP'] = xdf.loc[:, 'CHROM'].astype(str) + ':' + xdf.loc[:, 'POS'].astype(str) + ':' + xdf.loc[:, 'REF'] + ':' + xdf.loc[:, 'ALT']
  xdict = {x: pd.Series(xdf.loc[:, x].values, index = xdf.loc[:, 'SNP']).to_dict() for x in args.c}

if args.i:
  include = set(pd.read_csv(args.i, header = None)[0])

if args.e:
  exclude = set(pd.read_csv(args.e, header = None)[0])

if args.d:
  damaging = set(pd.read_csv(args.d, header = None)[0])

sys.stderr.write('annot2table.py: Extract fields with SnpSift\n')

fields = ['CHROM', 'POS', 'REF', 'ALT'] + (['FILTER'] if args.f and not 'FILTER' in args.fld else []) + args.fld + ['ANN[0].' + x for x in args.ann]
snpsift_args = ['java', '-Xmx' + str(args.mem) + 'g', '-jar', args.jar, 'extractFields', args.vcf] + fields
sys.stderr.write('annot2table.py: ' + ' '.join(snpsift_args) + '\n')
p1 = Popen(snpsift_args, stdin = sys.stdin.buffer if args.vcf == '-' else None, stdout = PIPE)

tp = pd.read_csv(p1.stdout, iterator = True, chunksize = args.s, sep = '\t', low_memory = False)
flag = True
for chunk in tp:
  replace = dict([(x, x.replace('#', '', 1)) for x in chunk.columns if x[0]=='#'] + [(x, x.replace('dbNSFP_', '', 1)) for x in chunk.columns if x[0:7]=='dbNSFP_'] + [(x, x.replace('ANN[0].', '', 1)) for x in chunk.columns if x[0:7]=='ANN[0].'])
  chunk.rename(columns=replace, inplace = True)
  if not 'SNP' in chunk:
    chunk['SNP'] = chunk.loc[:, 'CHROM'].astype(str) + ':' + chunk.loc[:, 'POS'].astype(str) + ':' + chunk.loc[:, 'REF'] + ':' + chunk.loc[:, 'ALT']

  # add columns to table
  if args.x and args.c:
    for x in args.c:
      chunk.loc[:, x] = chunk.loc[:, 'SNP'].map(xdict[x])

  # create a damaging impact annotation
  if args.d and 'IMPACT' in chunk:
    idx = (chunk.loc[:, 'IMPACT'] != 'HIGH') & chunk.loc[:, 'SNP'].apply(lambda x: x in damaging)
    chunk.loc[idx, 'IMPACT'] = 'DAMAGING'
    
  # output columns
  idx = pd.Series(index = chunk.index, dtype = bool)
  if args.f and 'FILTER' in chunk and not chunk.loc[:, 'FILTER'].isnull().all():
    idx &= (chunk.loc[:, 'FILTER'] == 'PASS') | (chunk.loc[:, 'FILTER'] == '.')
  if args.i:
    idx &= [x in include for x in chunk.loc[:, 'SNP']]
  if args.e:
    idx &= [not x in exclude for x in chunk.loc[:, 'SNP']]
  columns = ['CHROM', 'POS', 'REF', 'ALT'] + [x.replace('dbNSFP_', '', 1) if x[0:7]=='dbNSFP_' else x for x in args.fld] + args.ann + (args.c if args.c else [])
  chunk.loc[idx,].to_csv(sys.stdout, sep = '\t', na_rep = '.', columns = columns, header = flag, index = False)
  flag = False
