#!/bin/bash
###
#  markerqc.sh - perform marker quality control using plink
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

# performs basic quality control of markers
# requires plink prefix of the input dataset
# requires plink prefix of the output dataset
# optional file with list of samples to be removed
# optional minimum frequency of missingness allowed (I use .02 for Illumina, and .05 for Affymetrix)

set -e -o pipefail

if [ $# -lt 2 ]; then
  echo "About:   Perform marker quality control using plink (Oct 3rd 2016)"
  echo "Usage:   markerqc.sh <input prefix> <output prefix> [remove range] [lmiss threshold]"
  exit 1
fi

in="$1"
out="$2"

# makes sure the ouput directory exists                                                                                                 
dir=$(dirname $out)
mkdir -p $dir

if [ $# -ge 3 ]; then
  opt="--remove $3"
fi

if [ $# -ge 4 ]; then
  lmiss="$4"
else
  lmiss=".02"
fi

# identifies markers not on autosomes
awk '$1=="0" || $4=="0" {print $2}' $in.bim | sort > $out.noref

# identifies markers with excess missingness
plink --bfile $in --keep-allele-order $opt --freq counts gz --out $out
gzip -cd $out.frq.counts.gz | awk 'NR>1 && ($5=="NA" || $5==0 || $6==0) {print $2}' | sort > $out.nofrq
gzip -cd $out.frq.counts.gz | awk 'NR>1 && $5!="NA" && $5!=0 && $5!=1 && ($5/($5+$6)<.005 || $5/($5+$6)>.995) {print $2}' | sort > $out.frqex

# identifies markers with too elevated levels of missingness
plink --bfile $in --keep-allele-order $opt --exclude $out.nofrq --missing gz --out $out
gzip -cd $out.lmiss.gz | awk -v f=$lmiss 'NR>1 && $5>=f {print $2}' | sort > $out.miss.ex

# generates historgrams of missingness for individuals and for markers
gzip -cd $out.lmiss.gz | awk 'NR>1 {x[$3]++} END {for (i in x) print i"\t"x[i]}' | sort -k1,1n | awk '{SUM+=$2; print $0"\t"SUM}' > $out.lhist
gzip -cd $out.imiss.gz | awk 'NR>1 {x[$4]++} END {for (i in x) print i"\t"x[i]}' | sort -k1,1n | awk '{SUM+=$2; print $0"\t"SUM}' > $out.ihist

# identifies markers failing Hardy Weinberg equilibrium test due to excess heterozygousness
plink --bfile $in --keep-allele-order $opt --hardy gz --out $out
gzip -cd $out.hwe.gz | awk 'NR>1 && $3~"ALL" && ($9<1e-6 && $7>$8) {print $2}' | sort > $out.hwe.ex

# generates file with list of markers to be removed
if [ -f $out.nopass ]; then
  cat $out.{nopass,noref,nofrq,frqex,miss.ex,hwe.ex} | sort | uniq > $out.ex
else
  cat $out.{noref,nofrq,frqex,miss.ex,hwe.ex} | sort | uniq> $out.ex
fi
