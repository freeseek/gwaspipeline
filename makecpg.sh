#!/bin/bash
###
#  makecpg.sh - generate a VCF file with all CpG mutations in GRCh37
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

set -e -o pipefail

# install basic tools (Debian specific)
sudo apt-get install bedtools gzip samtools tabix wget

# preparation steps
mkdir -p ~/res/

# download human genome reference (GRCh37) if not already available
if [ ! -f ~/res/human_g1k_v37.fasta ]; then
  wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz |
    gzip -d > ~/res/human_g1k_v37.fasta
fi
if [ ! -f ~/res/human_g1k_v37.fasta.fai ]; then
  samtools faidx ~/res/human_g1k_v37.fasta
fi

# create a VCF file with all CpG mutations (it takes several hours)
awk '{for (i=1; i<$2; i++)
  print $1"\t"i-1"\t"i+1"\t"$1" "i}' \
  ~/res/human_g1k_v37.fasta.fai |
  bedtools getfasta -name \
  -fi ~/res/human_g1k_v37.fasta \
  -bed /dev/stdin \
  -fo /dev/stdout |
  tr '\n' ' ' |
  tr '>' '\n' |
  grep "CG $" |
  awk -v OFS="\t" 'BEGIN {print "##fileformat=VCFv4.1";
  print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}
  {print $1"\t"$2"\t.\tC\tT\t.\t.\t.";
  print $1"\t"$2+1"\t.\tG\tA\t.\t.\t."}' |
  bgzip > ~/res/cpg.vcf.gz &&
  tabix -f ~/res/cpg.vcf.gz
