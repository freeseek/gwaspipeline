#!/bin/bash
###
#  annot2set.sh - create variant sets from annotations
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

if [ $# -lt 3 ]; then
  echo "About:   Create variant sets from annotations (Oct 3rd 2016)"
  echo "Usage:   annot2set.sh <input VCF> <output prefix> <respath>"
  exit 1
fi

vcf="$1"
out="$2"
res="$3"

# makes sure the ouput directory exists
dir=$(dirname $out)
mkdir -p $dir

snpsift="java -Xmx2g -jar $res/SnpSift.jar"
known="$res/knownGene.txt.gz"

bcftools query -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.all

bcftools query -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq -c | awk '$1>1 {print $2}' > $out.vcfdups

bcftools query -i "REF~'^[ACGT]$' && ALT~'^[ACGT]$'" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.snps

bcftools query -i "REF~'[ACGT][ACGT]' || ALT~'[ACGT][ACGT]'" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.indels

bcftools query -i "FILTER=='PASS'" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.pass

bcftools query -i "FILTER!='PASS'" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.nopass

bcftools query -i "PM==1" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.pm

bcftools query -i "KGPhase3==1" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.kgp3

bcftools query -i "OM==1" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.om

bcftools query -i "CPG==1" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.cpg

bcftools query -i "EXACAC>0" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.exac

bcftools query -i "HRCAC>0" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.hrc

bcftools query -i "NONPSYCHAC>0" -f "%CHROM:%POS:%REF:%ALT\n" $vcf | sort | uniq > $out.nonpsych

$snpsift filter "ANN[0].IMPACT = 'HIGH'" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.high

$snpsift filter "(ANN[0].IMPACT = 'HIGH') & (ANN[0].EFFECT = 'frameshift_variant')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.frameshift

$snpsift filter "(ANN[0].IMPACT = 'HIGH') & (ANN[0].EFFECT = 'splice_acceptor_variant&intron_variant')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.splice_acceptor

$snpsift filter "(ANN[0].IMPACT = 'HIGH') & (ANN[0].EFFECT = 'splice_donor_variant&intron_variant')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.splice_donor

$snpsift filter "(ANN[0].IMPACT = 'HIGH') & (ANN[0].EFFECT = 'stop_gained')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.stop_gained

$snpsift filter "(ANN[0].IMPACT = 'MODERATE') & (ANN[0].EFFECT = 'inframe_insertion' | ANN[0].EFFECT = 'inframe_deletion')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.nondisruptive_inframe

$snpsift filter "(ANN[0].IMPACT = 'MODERATE') & (ANN[0].EFFECT = 'disruptive_inframe_insertion' | ANN[0].EFFECT = 'disruptive_inframe_deletion')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.disruptive_inframe

$snpsift filter "(ANN[0].IMPACT = 'MODERATE') & (ANN[0].EFFECT =~ 'inframe')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.inframe

$snpsift filter "(ANN[0].IMPACT = 'HIGH') & (ANN[0].EFFECT =~ 'protein_protein_contact')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.protein_protein_contact

$snpsift filter "ANN[0].IMPACT = 'MODERATE'" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.moderate

$snpsift filter "ANN[0].IMPACT = 'LOW'" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.low

$snpsift filter "ANN[0].IMPACT = 'MODIFIER'" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.modifier

$snpsift filter "(exists LOF[*].PERC) & (LOF[*].PERC >= 0.35)" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.lof

$snpsift filter "(exists NMD[*].PERC) & (NMD[*].PERC >= 0.35)" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.nmd

for pred in Polyphen2_H{DIV,VAR} SIFT LRT FATHMM PROVEAN; do
  $snpsift filter "(dbNSFP_${pred}_pred =~ 'D')" $vcf | bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.${pred,,}
done
pred=MutationTaster
$snpsift filter "(dbNSFP_${pred}_pred =~ '[AD]')" $vcf | bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.${pred,,}
pred=MutationAssessor
$snpsift filter "(dbNSFP_${pred}_pred =~ '[HM]')" $vcf | bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq > $out.${pred,,}

###########################################################################
## KNOWN GENES RELATED DEFINITIONS                                       ##
###########################################################################

gzip -cd $known | awk -F"\t" '{split($9, a, ","); split($10, b, ","); for (i=1; i<length(a); i++) {from=$6>a[i]?$6:a[i]; 
  to=$7<b[i]?$7:b[i]; if (to>from) print $2"\t"from-2"\t"to+2}}' | sed 's/^chr//' |
  bedtools sort | bedtools merge | bedtools intersect -a <(bcftools query -f "%CHROM\t%POS\t%POS\t%CHROM:%POS:%REF:%ALT\n" $vcf) -b - |
  cut -f4 | sort | uniq | join - <(cat $out.{high,moderate,low} | sort | uniq) > $out.known

$snpsift filter "(ANN[0].EFFECT = 'synonymous_variant')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq | join - $out.known > $out.synonymous

$snpsift filter "(ANN[0].EFFECT = 'missense_variant') | (ANN[0].EFFECT = 'protein_protein_contact')" $vcf |
  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" | sort | uniq | join - $out.known > $out.missense

join $out.missense $out.polyphen2_hdiv | join - $out.polyphen2_hvar | join - $out.sift | join - $out.lrt |
  join - $out.mutationtaster | join - $out.mutationassessor | join - $out.provean > $out.misdamaging

join -a1 -v1 $out.missense $out.misdamaging | join -a1 -v1 - $out.protein_protein_contact | join - $out.known > $out.misbenign

join -a1 -a2 $out.misdamaging $out.inframe | join -a1 -a2 - $out.protein_protein_contact | join - $out.known > $out.damaging

join -a1 -v1 $out.high $out.protein_protein_contact | join - $out.known > $out.disruptive

join -a1 -a2 $out.damaging $out.disruptive > $out.deleterious
