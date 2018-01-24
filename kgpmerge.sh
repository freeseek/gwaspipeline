#!/bin/bash
###
#  kgpmerge.sh - merge plink dataset with 1000 Genomes project phase 1 dataset
#  Copyright (C) 2016 Giulio Genovese
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Written by Giulio Genovese <giulio.genovese@gmail.com>
###

# merges a plink dataset with KGP samples
# requires plink prefix of the input dataset
# requires plink prefix of the output dataset
# optional file with list of markers to be extracted (usually a prune.in file)
# optional file with list of samples to be removed

set -e -o pipefail

if [ $# -lt 3 ]; then
  echo "About:   Merge plink dataset with 1000 Genomes project phase 1 dataset (Oct 3rd 2016)"
  echo "Usage:   kgpmerge.sh <input prefix> <output prefix> <kgp prefix> [extract range] [remove range]"
  exit 1
fi

in="$1"
out="$2"
kgp="$3"

# makes sure the ouput directory exists
dir=$(dirname $out)
mkdir -p $dir

# if list of markers with at least 1% frequency in KGP set does not exist, create it
if [ ! -f $kgp.maf.001 ]; then
  if [ ! -f $kgp.frq.counts.gz ]; then
    plink --bfile $kgp --keep-allele-order --freq counts --out $kgp && gzip $kgp.frq.counts
  fi
  gzip -cd $kgp.frq.counts.gz | tail -n+2 | awk '$5/($5+$6)>.01 || $6/($5+$6)>.01 {print $2}' | sort > $kgp.maf.001
fi

# generate a list of markers to extract from each dataset
if [ $# -ge 4 ]; then
  cat "$4"
else
  awk '{print $2}' $in.bim
fi | sort | join - $kgp.maf.001 > $out.extract

if [ $# -ge 5 ]; then
  opt="--remove $5"
fi

# extract markers from the provided dataset and the KGP dataset and merge the two datasets
plink --bfile $in --keep-allele-order --extract $out.extract $opt --make-bed --out $out.set
plink --bfile $kgp --keep-allele-order --extract $out.extract --make-bed --out $out.kgp
plink --bfile $out.set --keep-allele-order --bmerge $out.kgp --make-bed --out $out

# cleanup
rm $out.extract $out.set.bed $out.set.bim $out.set.fam $out.set.log $out.kgp.bed $out.kgp.bim $out.kgp.fam $out.kgp.log

# create population categorical data
awk 'BEGIN {print "FID IID POP"} NR==FNR {pop[$2]=$3}
  !($2 in pop) {pop[$2]="SET"} NR>FNR {print $1,$2,pop[$2]}' $kgp.pop $out.fam > $out.pop

# compute similarity matrix
plink --bfile $out --keep-allele-order --autosome --make-grm-bin --out $out --thread-num 10
