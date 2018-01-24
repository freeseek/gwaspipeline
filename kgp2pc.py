#!/usr/bin/env python3
"""
   kgp2pc.py - compute PCs from merged KGP datasets
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

import argparse, sys, pandas as pd, numpy as np
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description = 'kgp2pc.py: compute PC components from merged KGP datasets (Oct 3rd 2016)', add_help = False, usage = 'kgp2pc.py [options]')
parser.add_argument("--grm-bin", metavar = "[prefix]", type = str, help = "Prefix for binary GRM file")
parser.add_argument("--pop", metavar = "[filename]", type = str, required = True, help = "File with population information")
parser.add_argument('--fam', metavar = '[filename]', type = str, help = 'Specify full name of .fam file')
parser.add_argument("--pca", metavar = "[int]", type = int, default = 20, help = "number of PCs")
parser.add_argument("--out", metavar = "[prefix]", type = str, default = "plink", help = "Specify prefix for output files")
parser.add_argument("--groups", metavar = "[groups]", type = str, default = 'ALL,AFAM,EUR', help = "Groups to be used for the principal component computations")
parser.add_argument("--remove", metavar = "[filename]", type = str, help = "Exclude all samples named in the file")
parser.add_argument('--xlsx', action = 'store_true', default = False, help = 'Whether the output table is an xlsx file.')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHE']
dtypes = {'FID': str, 'IID': str, 'PAT': str, 'MAT': str, 'SEX': int, 'PHE': int}
df_fam = pd.read_csv(args.fam, delim_whitespace = True, header = None, names = names, dtype = dtypes)
dtypes = {'FID': str, 'IID': str, 'POP': str}
df_pop = pd.read_csv(args.pop, delim_whitespace = True, dtype = dtypes)
df = df_fam[['FID', 'IID', 'SEX']].merge(df_pop, how = 'left')
df.loc[df['POP'].isnull(), 'POP'] = 'SET'

pop = dict()
pop['EAS'] = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV']
pop['EUR'] = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']
pop['AFR'] = ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB']
pop['AMR'] = ['MXL', 'PUR', 'CLM', 'PEL']
pop['SAS'] = ['GIH', 'PJL', 'BEB', 'STU', 'ITU']
for lbl in pop:
  pop[lbl] = pop[lbl] + [lbl + '-' + x for x in pop[lbl]]
pop['SET'] = ['SET']
pop['ALL'] = pop['EAS'] + pop['EUR'] + pop['AFR'] + pop['AMR'] + pop['SAS']
pop['AFAM'] = pop['EUR'] + pop['AFR']

idx = dict()
for lbl in pop:
  idx[lbl] = np.array([x in pop[lbl] for x in df['POP']])

# invoke gcta64 to compute the principal components
for group in args.groups.split(','):
  if not group in idx:
    sys.stderr.write('kgp2pc.py: Group ' + group + ' undefined\n')
    continue
  sys.stderr.write('kgp2pc.py: Computing PC for group ' + group + '\n')
  gcta64_args = ["gcta64",
    '--grm-bin', args.grm_bin,
    '--keep', '/dev/stdin',
    '--pca', str(args.pca),
                 '--out', args.out + '.' + group.lower()]
  if args.remove:
    gcta64_args += ['--remove', args.remove]

  p = Popen(gcta64_args, stdin = PIPE, universal_newlines = True)
  df[np.logical_or(idx[group], idx['SET'])][['FID', 'IID']].to_csv(p.stdin, sep = '\t', index = False)
  p.communicate()
  if p.returncode:
    raise Exception('kgp2pc.py: Problems running gcta64')

  # add principal components information
  eigenval = pd.read_csv(args.out + '.' + group.lower() + '.eigenval', header = None, squeeze = True)
  names = ['FID', 'IID'] + ['PC' + str(i + 1) for i in range(args.pca)]
  dtypes = {'FID': str, 'IID': str}
  dtypes.update({'PC' + str(i + 1): float for i in range(args.pca)})
  eigenvec = pd.read_csv(args.out + '.' + group.lower() + '.eigenvec', delim_whitespace = True, header = None, names = names, dtype = dtypes)
  for i in range(args.pca):
    eigenvec['PC' + str(i + 1)] *= eigenval[i]
  df_pca = df[np.logical_or(idx[group], idx['SET'])].merge(eigenvec, how = 'left')

  # generate output table
  if (args.xlsx):
    writer = pd.ExcelWriter(args.out + '.' + group.lower() + '.xlsx', engine = 'xlsxwriter')
    df_pca.to_excel(writer, sheet_name = 'Sheet1', index = False)
    writer.save()
  else:
    df_pca.to_csv(args.out + '.' + group.lower() + ".pca", sep = "\t", index = False)
