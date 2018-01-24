#!/usr/bin/env python3
"""
   plink2score.py - compute polygenic scores
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

import argparse, os, sys, pandas as pd, numpy as np, functools
from subprocess import call

try:
  import pandas.rpy.common as com, rpy2
except ImportError:
  sys.stderr.write('You need to install the rpy2 module first\n')
  sys.stderr.write('(run this in your terminal: "python3 -m pip install rpy2" or "python3 -m pip install --user rpy2")\n')
  exit(2)

parser = argparse.ArgumentParser(description = 'plink2score.py: Run plink software to compute polygenic scores (Oct 3rd 2016)', add_help = False, usage = 'plink2score.py [options]')
parser.add_argument('--bfile', metavar = '{prefix}', type = str, help = 'Specify .bed + .bim + .fam prefix')
parser.add_argument('--bed', metavar = '[filename]', type = str, help = 'Specify full name of .bed file')
parser.add_argument('--bim', metavar = '[filename]', type = str, help = 'Specify full name of .bim file')
parser.add_argument('--fam', metavar = '[filename]', type = str, help = 'Specify full name of .fam file')
parser.add_argument('--extract', metavar = '<range>', type = str, help = 'Exclude all variants not named in the file.')
parser.add_argument('--exclude', metavar = '<range>', type = str, help = 'Exclude all variants named in the file.')
parser.add_argument('--keep', metavar = '[filename]', type = str, help = 'Exclude all samples not named in the file.')
parser.add_argument('--remove', metavar = '[filename]', type = str, help = 'Exclude all samples named in the file.')
parser.add_argument('--dosage', metavar = '[filename]', type = str, nargs = '+', help = 'a master list with one entry per line (see plink --help --dosage)')
parser.add_argument('--map', metavar = '[filename]', type = str, help = 'Specify full name of .map file')
parser.add_argument('--score', metavar = '[filename]', type = str, required = True, nargs = '+', help = 'Linear scoring system (see plink --help --score)')
parser.add_argument('--q-score-range', metavar = '[range file]', type = str, help = 'p-value ranges (see plink --help --score)')
parser.add_argument('--covar', metavar = '[filename]', type = str, help = 'File with covariates to regress polygenic scores')
parser.add_argument('--covar-name', metavar = '[filename]', type = str, default = 'PC1,PC2,PC3,PC4,PC5' , help = 'Use column in the covariate file with the given name')
parser.add_argument('--pop', metavar = '[filename]', type = str, help = 'Use table to include samples labels')
parser.add_argument('--adj-pop', metavar = '[pop]', type = str, help = 'Label group to use to perform the adjustments')
parser.add_argument('--debug', metavar = '[prefix]', type = str, default = 'plink', help = 'Specify debugging file with polygenic scores before covariate adjusting')
parser.add_argument('--out', metavar = '[prefix]', type = str, default = 'plink', help = 'Specify prefix for output files')
parser.add_argument('--xlsx', action='store_true', default = False, help='Whether the output table is an xlsx file.')
parser.add_argument('--noclean', action='store_true', default = False, help='Clean up temporary files.')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

# make sure the ouput directory exists
if args.out:
  outdir = os.path.dirname(args.out)
  if outdir != '' and not os.path.isdir(outdir):
    try:
      os.makedirs(outdir)
    except FileExistsError:
      if not os.path.isdir(outdir):
        raise Exception('plink2score.py: Could not create output directory')

if args.covar:
  df_covar = pd.read_csv(args.covar, delim_whitespace = True, dtype = {'FID':np.object, 'IID':np.object})
  pc_names = args.covar_name.split(',')
  if any([x not in df_covar.columns for x in pc_names]):
    raise Exception('plink2score.py: missing covariate')

if args.pop:
  df_pop = pd.read_csv(args.pop, delim_whitespace = True, dtype = {'FID':np.object, 'IID':np.object})
  if not any(df_pop['POP'] == args.adj_pop):
    raise Exception('plink2score.py: missing samples label to perform adjustment')

# uses default set of ranges if none given
if not args.q_score_range:
  args.q_score_range = args.out + '.range'
  f = open(args.q_score_range, 'w')
  f.write("S1 0 1e-8\nS2 0 1e-7\nS3 0 1e-6\nS4 0 1e-5\nS5 0 1e-4\nS6 0 0.001\nS7 0 0.05\nS8 0 0.1\nS9 0 0.2\nS10 0 0.3\nS11 0 0.4\nS12 0 0.5\nS13 0 1\n")
  f.close()

ranges = pd.read_csv(args.q_score_range, delim_whitespace = True, header = None, names = ['V1', 'V2', 'V3'])

list_tables = list()
if args.pop:
  list_tables.append(df_pop)
if args.covar:
  if args.pop and 'POP' in df_covar:
    df_covar = df_covar.drop('POP', 1)
  list_tables.append(df_covar)

for score in args.score:
  # remove possible output files as plink is not guaranteed to overwrite it
  for rng in ranges['V1']:
    if os.path.isfile(args.out + '.' + rng + '.profile'):
      os.remove(args.out + '.' + rng + '.profile')
  SCORE = '.'.join(os.path.basename(score).upper().split('.')[0:-1])
  sys.stderr.write('plink2score.py: Run plink to compute polygenic scores\n')
  plink_args = ['plink', '--score', score, '1', '2', '3', '--q-score-range', args.q_score_range, score, '1', '4']
  if (args.bfile):
    plink_args += ['--bfile', args.bfile]
  if (args.bed):
    plink_args += ['--bed', args.bed]
  if (args.bim):
    plink_args += ['--bim', args.bim]
  if (args.fam):
    plink_args += ['--fam', args.fam]
  if (args.extract):
    plink_args += ['--extract', args.extract]
  if (args.exclude):
    plink_args += ['--exclude', args.exclude]
  if (args.keep):
    plink_args += ['--keep', args.keep]
  if (args.remove):
    plink_args += ['--remove', args.remove]
  if (args.dosage):
    plink_args += ['--dosage'] + args.dosage
  if (args.map):
    plink_args += ['--map', args.map]
  if (args.out):
    plink_args += ['--out', args.out]
  print(' '.join(plink_args))
  if call(plink_args):
    raise Exception('plink2score.py: Problems running plink')
  for rng in ranges['V1']:
    if args.bfile or args.bed:
      names = ['FID', 'IID', 'PHENO', SCORE + '_CNT_' + rng, SCORE + '_CNT2_' + rng, SCORE + '_SCORE_' + rng]
      dtypes = {'FID':np.object, 'IID':np.object, 'PHENO':np.float, SCORE + '_CNT_' + rng:np.float, SCORE + '_CNT2_' + rng:np.float32, SCORE + '_SCORE_' + rng:np.float}
    elif args.dosage:
      names = ['FID', 'IID', 'PHENO', SCORE + '_SCORE_' + rng]
      dtypes = {'FID':np.object, 'IID':np.object, 'PHENO':np.float, SCORE + '_SCORE_' + rng:np.float}
    if os.path.isfile(args.out + '.' + rng + '.profile'):
      list_tables.append(pd.read_csv(args.out + '.' + rng + '.profile', delim_whitespace = True,
        names = names, dtype = dtypes, skiprows = 1).drop('PHENO', 1))

df = functools.reduce(lambda x, y: x.merge(y, how = 'left', on = ['FID', 'IID']), list_tables)

# if requested, output intermediate table before polygenic score adjustment
if args.debug:
  df.to_csv(args.debug, sep = "\t", na_rep = 'NA', index = False)

# compute adjusted scores by regressing on principal components and renormalizing
if args.covar:
  r_dataframe = com.convert_to_r_dataframe(df)
  for score in args.score:
    SCORE = '.'.join(os.path.basename(score).upper().split('.')[0:-1])
    for rng in ranges['V1']:
      if SCORE + '_SCORE_' + rng in df.columns:
        formula = SCORE + '_SCORE_' + rng + ' ~ ' + ' + '.join(pc_names)
        fit = rpy2.robjects.r.lm(formula, r_dataframe)
        predscore = rpy2.robjects.r.predict(fit, r_dataframe)
        df[SCORE + '_ADJSCORE_' + rng] = df[SCORE + '_SCORE_' + rng] - predscore
        if args.adj_pop:
          idx = df['POP'] == args.adj_pop
          if not any(idx):
            raise Exception('plink2score.py: missing samples label to perform adjustment')
        else:
          idx = df.index
        df[SCORE + '_ADJSCORE_' + rng] -= np.mean(df[SCORE + '_ADJSCORE_' + rng][idx])
        df[SCORE + '_ADJSCORE_' + rng] /= np.std(df[SCORE + '_ADJSCORE_' + rng][idx])  

# generate output table
if (args.xlsx):
  writer = pd.ExcelWriter(args.out + '.score.xlsx', engine = 'xlsxwriter')
  df.to_excel(writer, sheet_name = 'Sheet1', index = False)
  writer.save()
else:
  df.to_csv(args.out + ".score.tsv", sep = "\t", index = False)

# clean up
if (not args.noclean):
  sys.stderr.write('plink2score.py: Clean up temporary files\n')
  for sfx in ["hh", "range"] + list(ranges['V1'] + '.profile'):
    if os.path.isfile(args.out + "." + sfx):
      os.remove(args.out + "." + sfx)
