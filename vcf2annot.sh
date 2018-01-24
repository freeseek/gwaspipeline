#!/bin/bash
###
#  vcf2annot.sh - annotate a VCF file
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

# annotates a VCF file using SnpEff, KGP, ExAC, and dbNSFP
# requires input VCF file as first parameter (- for standard input)
# requires output VCF file prefix as second paramter

set -e -o pipefail

if [ $# -lt 4 ]; then
  echo "About:   Annotate a VCF file (Oct 3rd 2016)"
  echo "Usage:   vcf2annot.sh <input VCF> <output prefix> <ref> <respath>"
  exit 1
fi

vcf="$1"
out="$2"
ref="$3"
res="$4"

# makes sure the ouput directory exists                                                                                                 
dir=$(dirname $out)
mkdir -p $dir

dbsnp="$res/dbsnp.vcf.gz"
cpg="$res/cpg.vcf.gz"
exac="$res/ExAC.r0.3.sites.vep.norm.vcf.gz"
hrc="$res/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz"
nonpsych="$res/ExAC.r0.3.nonpsych.sites.norm.vcf.gz"
gnomad="$res/gnomad.genomes.r2.0.1.sites.vcf.gz"
dbnsfp="$res/dbNSFP2.9.txt.gz"

vcf_fields="CHROM POS REF ALT QUAL FILTER AC AN ANN LOF NMD OM KGPhase3 RS CPG EXACAC HRCAC NONPSYCHAC NOPASS ANCLOF"
dbnsfp_fields="$(echo {SIFT,Polyphen2_H{DIV,VAR},LRT,Mutation{Taste,Assesso}r,FATHMM,PROVEAN}_pred clinvar_{clnsig,trait} COSMIC_{ID,CNT})"
header="$(echo "$vcf_fields $dbnsfp_fields" | sed 's/_pred//g' | tr ' ' '\t')"
format="$(echo "$vcf_fields $(echo $dbnsfp_fields | sed 's/\(^\| \)/\1dbNSFP_/g')" | sed 's/\(^\| \)/\1%/g' | tr ' ' '\t')\n"

if [ $vcf == "-" ]; then
  bcftools norm -Ou -m -any <&0
else
  bcftools norm -Ou -m -any $vcf
fi |
  bcftools norm -Ou -f $ref |
  bcftools annotate -Ou -x ID |
  bcftools view -c 0 |
  java -Xmx6g -jar $res/snpEff.jar ann -v -c $res/snpEff.config GRCh37.75 - |
  java -Xmx2g -jar $res/SnpSift.jar annotate -noId -info RS,PM,KGPhase3,OM $dbsnp - |
  java -Xmx2g -jar $res/SnpSift.jar annotate -exists CPG -noInfo $cpg - |
  java -Xmx2g -jar $res/SnpSift.jar annotate -name HRC -info AC $hrc - |
  java -Xmx2g -jar $res/SnpSift.jar annotate -name GNOMAD -info AC $gnomad - |
  java -Xmx2g -jar $res/SnpSift.jar annotate -name EXAC -info AC $exac - |
  java -Xmx2g -jar $res/SnpSift.jar annotate -name NONPSYCH -info AC $nonpsych - |
  java -Xmx2g -jar $res/SnpSift.jar dbnsfp -v -f $(echo $dbnsfp_fields|tr ' ' ',') -db $dbnsfp - |
  grep -v "^##bcftools\|^##SnpEff\|^##SnpSift" |
  bgzip > $out.annot.vcf.gz &&
  tabix -f $out.annot.vcf.gz

# cleanup
if [ -f snpEff_summary.html ]; then rm snpEff_summary.html; fi
if [ -f snpEff_genes.txt ]; then rm snpEff_genes.txt; fi
