#!/bin/bash
###
#  sampleqc.sh - perform sample quality control using plink
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

# performs basic quality control of samples
# requires plink prefix of the input dataset
# requires plink prefix of the output dataset
# optional file with list of markers to be extracted
# optional minimum frequency of miss
# outputs dups file with duplicate pairs in genome format
# outputs rels file with directly related pairs in genome format
# outputs dups.rm file with list of duplicates to remove
# outputs rels.rm file with list of closely related samples to remove

set -e -o pipefail

if [ $# -lt 2 ]; then
  echo "About:   Perform sample quality control using plink (Oct 3rd 2016)"
  echo "Usage:   sampleqc.sh <input prefix> <output prefix> [exclude range] [imiss threshold]"
  exit 1
fi

in="$1"
out="$2"

# makes sure the ouput directory exists
dir=$(dirname $out)
mkdir -p $dir

if [ $# -ge 3 ]; then
  opt="--exclude $3"
fi

if [ $# -ge 4 ]; then
  imiss="$4"
else
  imiss=".12"
fi


# compute set of independent markers
plink --bfile $in --keep-allele-order $opt --maf .01 --indep 50 5 2 --out $out

# compute pairwise relationships
plink --bfile $in --keep-allele-order --extract $out.prune.in --genome gz full --out $out

# write interesting pairs
gzip -cd $out.genome.gz | awk 'NR==1 || $10>=.35 && $10!="nan"' | sort -k10,10n > $out.rels.genome

# write report
awk -v OFS="\t" 'BEGIN {x["FID1:IID1"]="SEX1"; x["FID2:IID2"]="SEX2"}
  NR==FNR {x[$1":"$2]=$5} NR>FNR && FNR==1 {col1="REL"}
  NR>FNR && FNR>1 && $10>=.7 {col1="duplicate"}
  NR>FNR && FNR>1 && $10>=.35 && $10<.7 && $8>.7 {col1="parent-child"}
  NR>FNR && FNR>1 && $10>=.35 && $10<.7 && $8<.7 {col1="sibling"}
  NR>FNR {print col1,$1,$2,x[$1":"$2],$3,$4,x[$3":"$4],$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' \
  $in.fam $out.rels.genome > $out.rels.tsv

if [ ! -f $out.imiss.gz ]; then
  plink --bfile $in --keep-allele-order $opt --missing gz --out $out
fi

if [ ! -f $out.miss.rm ]; then
  gzip -cd $out.imiss.gz | awk -v f=$imiss 'NR>1 && $6>=f {print $1,$2}' | sort > $out.miss.rm
fi

# computes what duplicates to remove based on missingness
read -a iid1 <<< $(awk 'NR>1 && $10>=.7 {print $1":"$2}' $out.rels.genome)
read -a	iid2 <<< $(awk 'NR>1 && $10>=.7 {print $3":"$4}' $out.rels.genome)
n=$((${#iid1[@]}-1))
for i in $(seq 0 $n); do
  gzip -cd $out.imiss.gz | awk -v iid1="${iid1[$i]}" -v iid2="${iid2[$i]}" '$1":"$2==iid1 || $1":"$2==iid2' |
    sort -k6,6n | awk 'NR==2 {print $1,$2}'
done | sort | uniq > $out.dups.rm

# computes what samples to remove based on relatedness
# if the phenotype is discordant, removes the control, otherwise remove based on missingness
(echo "FID IID"; cat $out.dups.rm) |
  awk 'NR==FNR {x[$1":"$2]++} NR>FNR && FNR>1 && !x[$1":"$2] && !x[$3":"$4] {print $1,$2,$3,$4}' - $out.rels.genome |
  awk 'NR==FNR {x[$1":"$2]=$6} NR>FNR {print $1,$2,x[$1":"$2],$3,$4,x[$3":"$4]}' $in.fam - |
  awk 'NR==FNR {x[$1":"$2]=$6; next} $3==$6 || $3!~"[12]" && $6!~"[12]" {if (x[$1":"$2]>x[$4":"$5]) print $1,$2; else print $4,$5}
  $3==1 && $6==2 || $3!~"[12]" && $6~"[12]" {print $1,$2} $3==2 && $6==1 || $3~"[12]" && $6!~"[12]" {print $4,$5}' <(gzip -cd $out.imiss.gz) - |
  sort | uniq > $out.rels.rm

# computes set of duplicates / related samples to be removed
cat $out.{dups,rels}.rm | sort | uniq > $out.genome.rm
