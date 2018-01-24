#!/usr/bin/env python3
"""
   skat.py - wrapper to run SKAT software
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

import argparse, sys, pandas as pd, hashlib, numpy as np, os
from subprocess import call

parser = argparse.ArgumentParser(description = 'skat.py: Run SKAT software (Oct 3rd 2016)', add_help = False, usage = 'skat.py [options]')
parser.add_argument('--bfile', metavar = '{prefix}', type = str, help = 'Specify .bed + .bim + .fam prefix.')
parser.add_argument('--bed', metavar = '[filename]', type = str, help = 'Specify full name of .bed file.')
parser.add_argument('--bim', metavar = '[filename]', type = str, help = 'Specify full name of .bim file.')
parser.add_argument('--fam', metavar = '[filename]', type = str, help = 'Specify full name of .fam file.')
parser.add_argument('--out', metavar = '[prefix]', type = str, default = 'plink', help = 'Specify prefix for output files.')
parser.add_argument('--min-maf', metavar = '[freq]', type = str, help = 'Exclude variants with minor allele frequency lower than a threshold.')
parser.add_argument('--max-maf', metavar = '[freq]', type = str, help = 'Exclude variants with MAF greater than the threshold.')
parser.add_argument('--min-mac', metavar = '[ct]', type = str, help = 'Exclude variants with minor allele count lower than the given threshold')
parser.add_argument('--max-mac', metavar = '[ct]', type = str, help = 'Exclude variants with minor allele count greater than the given threshold.')
parser.add_argument('--extract', metavar = '<range>', type = str, help = 'Exclude all variants not named in the file.')
parser.add_argument('--exclude', metavar = '<range>', type = str, help = 'Exclude all variants named in the file.')
parser.add_argument('--keep', metavar = '[filename]', type = str, help = 'Exclude all samples not named in the file.')
parser.add_argument('--remove', metavar = '[filename]', type = str, help = 'Exclude all samples named in the file.')
parser.add_argument('--pheno', metavar = '[filename]', type = str, help = 'Load phenotype data from the specified file, instead of using the values in the main input fileset.')
parser.add_argument('--pheno-name', metavar = '[c]', type = str, help = 'If --pheno file has a header row, use column with the given name.')
parser.add_argument('--setid', metavar = '[SetID file]', type = str, required = True, help = 'Load file indicating gene sets.')
parser.add_argument('--covar', metavar = '[SetID file]', type = str, help = 'Specify covariate file.')
parser.add_argument('--continuous', action='store_true', default = False, help = 'Whether the phenotype is continuous.')
parser.add_argument('--rcorr', metavar = '[r.corr]', type = str, default = '1.0', help = 'Correlation [1.0].')
parser.add_argument('--weights-beta', metavar = '[beta weights]', type = float, nargs = 2, help = 'Beta weights [1 25].')
parser.add_argument('--weight', metavar = '[weight file]', type = str, help = 'Weights.')
parser.add_argument('--method', metavar = '[method]', type = str, default = 'davies', help = 'Method [davies].')
parser.add_argument('--impute', metavar = '[impute-method]', type = str, default = 'bestguess', help = 'Imputation method [bestguess].')
parser.add_argument('--noplink', action='store_true', default = False, help='Whether to not run the plink extraction part.')
parser.add_argument('--noskat', action='store_true', default = False, help='Whether to not run the SKAT part.')
parser.add_argument('--noclean', action='store_true', default = False, help='Clean up temporary files.')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

if not args.noplink:
    # extract plink subset for SKAT (this needs to be made into a loop)
    sys.stderr.write('skat.py: Extract plink subset\n')
    # notice that --set-hh-missing and --fill-missing-a2 leave calls on the Y chromosome for females as missing
    plink_args = ['plink', '--make-bed', '--set-hh-missing', '--fill-missing-a2']
    if (args.bfile):
      plink_args += ['--bfile', args.bfile]
    if (args.bed):
      plink_args += ['--bed', args.bed]
    if (args.bim):
      plink_args += ['--bim', args.bim]
    if (args.fam):
      plink_args += ['--fam', args.fam]
    if (args.min_maf):
      plink_args += ['--maf', args.min_maf]
    if (args.max_maf):
      plink_args += ['--max-maf', args.max_maf]
    if (args.min_mac):
      plink_args += ['--mac', args.min_mac]
    if (args.max_mac):
      plink_args += ['--max-mac', args.max_mac]
    if (args.extract):
      plink_args += ['--extract', args.extract]
    if (args.exclude):
      plink_args += ['--exclude', args.exclude]
    if (args.keep):
      plink_args += ['--keep', args.keep]
    if (args.remove):
      plink_args += ['--remove', args.remove]
    if (args.pheno):
      plink_args += ['--pheno', args.pheno]
    if (args.pheno_name):
      plink_args += ['--pheno-name', args.pheno_name]
    if (args.out):
      plink_args += ['--out', args.out]
    if call(plink_args):
      raise Exception('skat.py: Problems running plink')
    # notice that --fill-missing-a2 is of extreme importance here
    
    # adjust bim file
    sys.stderr.write('skat.py: Convert bim file to bim.md5 file\n')
    bim = pd.read_csv(args.out + '.bim', delim_whitespace = True, header = None, names = ['CHROM', 'SNP', 'MAP', 'POS', 'A1', 'A2'], low_memory = False)
    bim['SNP'] = bim['SNP'].apply(lambda x: hashlib.md5(x.encode()).hexdigest())
    bim.to_csv(args.out + '.bim.md5', sep = '\t', header = False, index = False)
    snplist = bim['SNP']
    del bim
    
    # adjust SetID file
    sys.stderr.write('skat.py: Convert SetID file to SetID.md5 file\n')
    setid = pd.read_csv(args.setid, delim_whitespace = True, header = None, names = ['SET', 'SNP'], low_memory = False)
    setid['SNP'] = setid['SNP'].apply(lambda x: hashlib.md5(x.encode()).hexdigest())
    setid[setid['SNP'].isin(snplist)].to_csv(args.out + '.SetID.md5', sep = '\t', header = False, index = False)
    del setid

if not args.noskat:
    # call Generate_SSD_SetID SKAT function
    sys.stderr.write('skat.py: Generate SSD and SSD.Info files\n')
    rscript = open(args.out + '.1.R', 'w')
    rscript.write('library(SKAT)\n')
    rscript.write('Generate_SSD_SetID(\'' + args.out + '.bed\', \'' + args.out + '.bim.md5\', \'' + args.out + '.fam\', \'' + args.out + '.SetID.md5\', \'' + args.out + '.SSD\', \'' + args.out + '.SSD.Info\')\n')
    rscript.close()
    if call(['Rscript', args.out + '.1.R']):
      raise Exception('skat.py: Problems running Rscript')
    
    # adjust weight file
    if (args.weight):
      sys.stderr.write('skat.py: Convert weight file to weight.md5 file\n')
      idweight = pd.read_csv(args.weight, delim_whitespace = True, header = None, names = ['SNP', 'WEIGHT'], low_memory = False)
      idweight['SNP'] = idweight['SNP'].apply(lambda x: hashlib.md5(x.encode()).hexdigest())
      idweight[idweight['SNP'].isin(snplist)].to_csv(args.out + '.weight.md5', sep = '\t', na_rep='0', header = False, index = False)
    
    # call SKAT.SSD.All SKAT function
    sys.stderr.write('skat.py: Run SKAT burden test\n')
    rscript = open(args.out + '.2.R', 'w')
    rscript.write('library(SKAT)\n')

    if (args.covar):
      rscript.write('FAM <- Read_Plink_FAM_Cov(\'' + args.out + '.fam\', \'' + args.covar + '\', Is.binary = FALSE, cov_header = TRUE)\n')
      rscript.write('formula <- as.formula(paste(\'y ~\', paste(paste0(\'FAM$\', names(FAM)[7:ncol(FAM)]), collapse = \' + \')))\n')
      formula = 'formula'
    else:
      rscript.write('FAM <- Read_Plink_FAM(\'' + args.out + '.fam\', Is.binary = FALSE)\n')
      formula = 'y ~ 1'

    if (args.continuous):
      rscript.write('y <- FAM$Phenotype\n')
      out_type = 'C'
    else:
      rscript.write('y <- FAM$Phenotype - 1\n')
      rscript.write('y[FAM$Phenotype==-9] <- NA\n')
      out_type = 'D'

    rscript.write('obj <- SKAT_Null_Model(' + formula + ', out_type = \'' + out_type + '\')\n')
    rscript.write('SSD.INFO <- Open_SSD(\'' + args.out + '.SSD\', \'' + args.out + '.SSD.Info\')\n')

    if (args.weight):
      rscript.write('weight <- Read_SNP_WeightFile(\'' + args.out + '.weight.md5\')\n')
      opt = ', obj.SNPWeight = weight'
    elif (args.weights_beta):
      opt = ',  weights.beta=c(' + str(args.weights_beta[0]) + ',' + str(args.weights_beta[1]) + ')'
    else:
      opt = ''

    if args.continuous:
      rscript.write('out <- SKAT.SSD.All(SSD.INFO, obj, method = \'' + args.method + '\', impute.method =\'' + args.impute + '\', r.corr = ' + args.rcorr + opt + ')\n')
    else:
      rscript.write('out <- SKATBinary.SSD.All(SSD.INFO, obj, method = \'' + args.method + '\', impute.method =\'' + args.impute + '\', r.corr = ' + args.rcorr + opt + ')\n')

    rscript.write('Close_SSD()\n')
    rscript.write('write.table(out$results, file = \'' + args.out + '.out\', row.names = FALSE, quote = FALSE)\n')
    rscript.close()
    if call(['Rscript', args.out + '.2.R']):
      raise Exception('skat.py: Problems running Rscript')

# perform case control burden count
if not args.continuous:
  sys.stderr.write('skat.py: Compute burden count\n')
  if call(['plink', '--bed', args.out + '.bed', '--bim', args.out + '.bim.md5', '--fam', args.out + '.fam', '--allow-no-sex', '--assoc', 'counts', '--out', args.out]):
    raise Exception('skat.py: Problems running plink')
  fam = pd.read_csv(args.out + '.fam', delim_whitespace = True, header = None, names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHE'])
  assoc = pd.read_csv(args.out + '.assoc', delim_whitespace = True)
  setid = pd.read_csv(args.out + '.SetID.md5', delim_whitespace = True, header = None, names = ['SET', 'SNP'], low_memory = False)
  counts = pd.merge(setid, assoc[['SNP', 'C_A', 'C_U']][assoc['C_A'] + assoc['C_U']>0])
  df = pd.pivot_table(counts, values=['C_A', 'C_U'], index = 'SET', aggfunc=np.sum)
  df['SetID'] = df.index
  df['OR'] = df['C_A'] / sum(fam['PHE']==2) / df['C_U'] * sum(fam['PHE']==1)
  skat = pd.read_csv(args.out + '.out', delim_whitespace = True)
  table = pd.merge(skat, df)
  table.to_csv(args.out + '.tsv', sep = '\t', na_rep = 'NA', index = False)

# clean up
if not args.noclean:
  sys.stderr.write('skat.py: Clean up temporary files\n')
  for sfx in ['hh', 'nosex', 'weight.md5', 'bed', 'bim', 'fam', 'bim.md5', 'SetID.md5', 'SSD_LOG.txt', 'SSD', 'SSD.Info.TEMP.txt', 'SSD.Info', '1.R', '2.R', 'assoc', 'log']:
    if os.path.isfile(args.out + '.' + sfx):
      os.remove(args.out + '.' + sfx)
