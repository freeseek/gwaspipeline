#!/usr/bin/env python3
"""
   kin2ped.py - create a plink ped file from a plink genome file
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

import argparse, pandas as pd

parser = argparse.ArgumentParser(description = 'kin2ped.py: create a plink ped file from a kinship file (Apr 17th 2017)', add_help = False, usage = 'kin2ped.py [options]')
parser.add_argument('--kin', metavar = '<kinship>', type = str, default = '/dev/stdin', help = 'Specify kinship file [stdin]')
parser.add_argument('--zip', action = 'store_true', default = False, help = 'whether input kinship file is compressed [FALSE]')
parser.add_argument('--fam', metavar = '[filename]', type = str, help = 'Specify fam file')
parser.add_argument('--out', metavar = '[filename]', type = str, default = '/dev/stdout', help = 'Specify output filename [stdout]')
parser.add_argument('--pdf', metavar = '[filename]', type = str, help = 'Specify output filename for pdf')
parser.add_argument('--max-ibs0-duos', metavar = 'FLOAT', type = float, default = 0.001, help = 'maximum IBS0 for parent child duos [0.001]')
parser.add_argument('--min-kin-dups', metavar = 'FLOAT', type = float, default = 0.45, help = 'maximum kinship for full siblings [0.45]')
parser.add_argument('--min-kin-full-sibs', metavar = 'FLOAT', type = float, default = 0.17, help = 'minimum kinship for full siblings [0.17]')
parser.add_argument('--max-kin-full-sibs', metavar = 'FLOAT', type = float, default = 0.35, help = 'maximum kinship for full siblings [0.35]')
parser.add_argument('--min-kin-half-sibs', metavar = 'FLOAT', type = float, default = 0.09, help = 'minimum kinship for half siblings [0.09]')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

# download the kinship table
dtypes = {'FID': str, 'ID1': str, 'ID2': str, 'N_SNP': int, 'Z0': float, 'Phi': float, 'HetHet': float, 'IBS0': float, 'Kinship': float, 'Error': float}
if args.zip:
  df = pd.read_csv(args.kin, delim_whitespace = True, dtype = dtypes, compression = 'gzip', low_memory = False)
else:
  df = pd.read_csv(args.kin, delim_whitespace = True, dtype = dtypes, low_memory = False)
if df.empty:
  exit()
  
dfrels1 = df[df['Error'] == 1]
if dfrels1.empty:
  exit()
dfrels2 = dfrels1.copy()
cols = list(dfrels2)
cols[1], cols[2] = cols[2], cols[1]
dfrels2.columns = cols
dfrels = pd.concat([dfrels1, dfrels2])
rels = dfrels.groupby('ID1')['ID2'].apply(lambda x: x.tolist()).to_dict()
idx = (dfrels['IBS0'] <= args.max_ibs0_duos) & (dfrels['Kinship'] > args.min_kin_dups)
dups = dfrels[idx].groupby('ID1')['ID2'].apply(lambda x: x.tolist()).to_dict()
idx = (dfrels['IBS0'] <= args.max_ibs0_duos) & (dfrels['Kinship'] <= args.max_kin_full_sibs)
duos = dfrels[idx].groupby('ID1')['ID2'].apply(lambda x: x.tolist()).to_dict()
idx = (dfrels['IBS0'] > args.max_ibs0_duos) & (dfrels['Kinship'] > args.min_kin_full_sibs) & (dfrels['Kinship'] < args.max_kin_full_sibs)
sibs = dfrels[idx].groupby('ID1')['ID2'].apply(lambda x: x.tolist()).to_dict()
kinship = dfrels.groupby(('ID1','ID2'))['Kinship'].apply(lambda x: x.tolist()[0]).to_dict()

# download the fam table
names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHE']
dtypes = {'FID': str, 'IID': str, 'PAT': str, 'MAT': str, 'SEX': int, 'PHE': int}
fam = pd.read_csv(args.fam, delim_whitespace = True, header = None, names = names, dtype = dtypes, low_memory = False)
fid = {x[0]: x[1] for x in zip(fam['IID'], fam['FID'])}
sex = {x[0]: x[1] for x in zip(fam['IID'], fam['SEX'])}
phe = {x[0]: x[1] for x in zip(fam['IID'], fam['PHE'])}

# identify trios
ped = []
done = set()
pat = dict()
mat = dict()
for c in duos:
  for f in (x for x in duos[c] if (sex[x] == 1)):
    for m in (x for x in duos[c] if (sex[x] == 2)):
      if not m in rels[f]:
        ped += [(fid[c], c, f, m, sex[c], phe[c])]
        done.add((c, f))
        done.add((f, c))
        done.add((c, m))
        done.add((m, c))
        pat[c] = f
        mat[c] = m

# identify duos/sibs triplets
for p in duos:
  for c in (x for x in duos[p] if (not (x, p) in done)):
    if any([x for x in duos[p] if (x != c) and (c in sibs) and (x in sibs[c])]) or (p in sibs) and any([x for x in sibs[p] if (x!=c) and (not x in duos[c])]):
      if sex[p] == 1:
        ped += [(fid[c], c, p, '.', sex[c], phe[c])]
        pat[c] = p
      elif sex[p] == 2:
        ped += [(fid[c], c, '.', p, sex[c], phe[c])]
        mat[c] = p
      done.add((c, p))
      done.add((p, c))

# by eclusion if father already found
for p in (x for x in duos if (x in pat)):
  for c in (x for x in duos[p] if (sex[x] == 1) and (not (x, p) in done)):
    if (not c in dups) or (not pat[p] in dups[c]):
      if sex[p] == 1:
        ped += [(fid[c], c, p, '.', sex[c], phe[c])]
        pat[c] = p
      elif sex[p] == 2:
        ped += [(fid[c], c, '.', p, sex[c], phe[c])]
        mat[c] = p
      done.add((c, p))
      done.add((p, c))

# by eclusion if mother already found
for p in (x for x in duos if (x in mat)):
  for c in (x for x in duos[p] if (sex[x] == 2) and (not (x, p) in done)):
    if (not c in dups) or (not mat[p] in dups[c]):
      if sex[p] == 1:
        ped += [(fid[c], c, p, '.', sex[c], phe[c])]
        pat[c] = p
      elif sex[p] == 2:
        ped += [(fid[c], c, '.', p, sex[c], phe[c])]
        mat[c] = p
      done.add((c, p))
      done.add((p, c))

# identify duos using http://apol1.blogspot.com/2015/10/distinguishing-parents-from-children.html
#for a in duos:
#  for b in duos[p]:
#    for b in (x for x in rels[p] if x!=a and not x in duos[p]):
#      if kinship[(p, a)] < 

if ped:
  pd.DataFrame(ped).to_csv(args.out, header = False, index = False, sep = '\t')

if args.pdf:
  from matplotlib.backends.backend_pdf import PdfPages
  with PdfPages(args.pdf) as pdf:
    ax = df.plot(kind='scatter', x='IBS0', y='Kinship', xlim=(0,.006), ylim=(0,1))
    ax.set_xlabel('IBS0')
    ax.set_ylabel('Kinship')
    ax.axhline(args.max_kin_full_sibs, color='r', linestyle='--', lw=2)
    ax.axhline(args.min_kin_half_sibs, color='r', linestyle='--', lw=2)
    ax.axvline(args.max_ibs0_duos, color='r', linestyle='--', lw=2)
    pdf.savefig(ax.get_figure())
