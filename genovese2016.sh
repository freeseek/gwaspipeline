#!/bin/bash
###
#  genovese2016.sh - script for the analyses on the Sweden exome cohort
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

# variables to be set by the user
vcf="swedish_scz_exomes.vcf.gz"
dict="dict.tsv"
sex="sex.tsv"
batch="batch.tsv"
birth="birth.tsv"
pheno="pheno.tsv"
kit="kit.tsv"
pfx="swe"

# install latest version of GWAS scripts
git clone https://github.com/freeseek/gwaspipeline
chmod a+x gwaspipeline/*.{py,sh,R}
/bin/mv gwaspipeline/*.{py,sh,R} ~/bin/

# install required software and resources
install.sh

# generate files containing gene sets
makegenes.sh

# convert VCF to plink format (rename samples on the fly)
awk '{print $0"\n0"$0"\n00"$0"\n000"$0}' $dict |
  bcftools reheader -s /dev/stdin $vcf |
  vcf2plink.py --ref ~/res/human_g1k_v37.fasta \
    --filter "FORMAT/DP<10" \
    --out $pfx \
    --impute-sex .3 .3 \
    --pdf $pfx.sexcheck.pdf

# create lite version of the dataset
plink2lite.sh $pfx $pfx.lite ~/res

# write sex report
sort $pfx.sexcheck | join -1 2 - $sex |
  awk -v OFS="\t" 'BEGIN {print "FID","IID","PEDSEX","SNPSEX","STATUS","F"}
  $4!=$7 {print $2,$1,$7,$4,"PROBLEM",$6}' > out/$pfx.sex.tsv

# generate list of individuals with mismatching sex nad perform quality control analyses
sort $pfx.fam | join -1 2 - $sex | awk '$5!=$7 {print $2,$1}' > out/$pfx.sex.rm
markerqc.sh $pfx out/$pfx
sampleqc.sh $pfx.lite out/$pfx.lite out/$pfx.ex

# write final sample report
csv2xlsx.py -d tab -b -f 1 0 -t Sex Duos Trios -i out/$pfx.sex.tsv out/$pfx.lite.rels.tsv out/$pfx.ped -o $pfx.report.xlsx

# annotate VCF file (remove samples on the fly)
bcftools view -Ou -c 0 -G $vcf |
  vcf2annot.sh - $pfx ~/res/human_g1k_v37.fasta ~/res

# generate variant sets for further analysis
annot2set.sh $pfx.annot.vcf.gz out/$pfx ~/res

# compute principal components with 1000 Genome project samples
kgpmerge.sh $pfx.lite kgp/$pfx.lite ~/res out/$pfx.lite.prune.in
kgp2pc.py --grm-bin kgp/$pfx.lite --fam kgp/$pfx.lite.fam --out kgp/$pfx.lite --pop ~/res/kgp.pop

# generate list of gene membership for all variants
awk -F":" '{print $1"\t"$2-1"\t"$2-1+length($3)"\t"$0}' out/$pfx.all |
  bedtools intersect -wb -a ~/res/ensGene.coding.bed -b - |
  cut -f4,8 | sort | uniq > out/$pfx.genes.SetID

# compute frequency for all variants
plink --bfile $pfx \
  --remove out/$pfx.lite.genome.rm \
  --set-hh-missing \
  --freq counts gz \
  --out out/$pfx.skat

# compute singleton count
plink --bfile $pfx \
  --mac 1 --max-mac 1 \
  --extract out/$pfx.pass \
  --exclude out/$pfx.nonpsych \
  --remove out/$pfx.lite.genome.rm \
  --set-hh-missing \
  --fill-missing-a2 \
  --make-bed \
  --out $pfx.urv

# generate list of singleton variants
plink --bfile $pfx.urv \
  --recode rlist omit-nonmale-y \
  --out $pfx.urv

# generate covariate file
awk 'NR==FNR {x[$6]++} NR>FNR {print $2"\t"x[$2]+0}' $pfx.urv.rlist $pfx.urv.fam |
  sort -t$'\t' -k1,1 > out/$pfx.urv.count.tsv
cut -f2,3,5- kgp/$pfx.lite.all.pca | sort -t$'\t' -k1,1 |
  join -t$'\t' - <(echo -e "IID\tURV"; cat out/$pfx.urv.count.tsv) |
  join -a1 -t$'\t' - <(echo -e "IID\tKIT"; cat $kit) |
  join -a1 -t$'\t' - <(echo -e "IID\tBIRTH"; cat $birth) |
  awk -F"\t" -v OFS="\t" 'NR==1 {print "FID\t"$0}
  NF==23 {$24="NA"} NF==24 {$25="NA"} NR>1 {print "0\t"$0}' > out/$pfx.skat.cov

# generate list of individuals that will be excluded from burden analyses
(echo -e "0 NA12878\n0 NA12891\n0 NA12892"
  cat out/$pfx.{lite.genome,sex}.rm
  cut -d" " -f5,6 $pfx.urv.rlist |
  sort | uniq -c | sort -k1,1n |
  awk '$1>100 {print $2,$3}') |
  sort | uniq > out/$pfx.rm

# generate table with URV variants and phenotypes
(echo -e "SNP\tPHENOTYPE\tGENDER";
cut -d" " -f1,6 $pfx.urv.rlist | sort -k2,2 |
  join -a1 -v1 -1 2 -2 2 - out/$pfx.rm |
  join - $pheno |
  sed -e 's/1$/CTRL/g' \
  -e 's/2$/SCZ/g' \
  -e 's/-9$/BD/g' |
  join - <(cut -d" " -f2,5 $pfx.fam | sort) |
  cut -d" " -f2- | sort |
  join - out/$pfx.deleterious |
  sed -e 's/ /\t/g' \
  -e 's/1$/MALE/g' \
  -e 's/2$/FEMALE/g') \
  > $pfx.urv.extra.tsv

# generate table with URV variants, annotations, and phenotypes
tail -n+2 $pfx.urv.extra.tsv | cut -f1 |
  annot2table.py \
  -jar ~/res/SnpSift.jar \
  -vcf $pfx.annot.vcf.gz \
  -fld CPG \
  -ann EFFECT GENE FEATUREID HGVS_C HGVS_P \
  -x $pfx.urv.extra.tsv \
  -c PHENOTYPE GENDER \
  -i /dev/stdin |
  sed 's/&.*_\(variant\|deletion\|insertion\|gained\|lost\)\t/\t/g' \
  > $pfx.urv.tsv

# compute association between genotypes and batches
for grp in $(cut -f2 $batch | sort | uniq); do
  awk -F"\t" -v grp="$grp" '$2==grp {print "0",$1,"2"} $2!=grp {print "0",$1,"1"}' $batch |
    plink --bfile $pfx \
    --pheno /dev/stdin \
    --covar kgp/$pfx.lite.all.pca \
    --covar-name $(echo PC{1..20}|tr ' ' ',') \
    --logistic sex \
    --out out/$pfx.$grp
  gzip -f out/$pfx.$grp.assoc.logistic
  gzip -cd out/$pfx.$grp.assoc.logistic.gz |
    awk '$5=="ADD" && $9<5e-8 {print $2}' |
    sort > out/$pfx.$grp.fail
done

# generate list of variants having batch effects and list of variants to exclude
for grp in $(cut -f2 $batch | sort | uniq); do cat out/$pfx.$grp.fail; done | sort | uniq > out/$pfx.fail
cat out/$pfx.{nopass,nofrq,hwe.ex,fail} | sort | uniq > out/$pfx.ex

# compute association with schizophrenia phenotype using logistic regression
awk -F"\t" '{print "0",$1,$2}' $pheno |
  plink --bfile $pfx \
  --remove out/$pfx.rm \
  --pheno /dev/stdin \
  --covar kgp/$pfx.lite.all.pca \
  --covar-name $(echo PC{1..5}|tr ' ' ',') \
  --logistic sex \
  --out out/$pfx
gzip -f out/$pfx.assoc.logistic

# generate Manhattan plot for association with schizophrenia phenotype
gzip -cd out/$pfx.assoc.logistic.gz |
  awk 'NR==FNR {x[$1]++} NR>FNR && (FNR==1 || $5=="ADD" && !($2 in x) && $9!="NA") {print $1"\t"$3"\t"$9}' out/$pfx.ex - \
  > out/$pfx.assoc.logistic.tsv
manplot.R out/$pfx.assoc.logistic.tsv $pfx.assoc.logistic.png

# generate gene sets SetID
(sort ~/res/genes/ens.txt |
  join -t $'\t' out/$pfx.genes.SetID - |
  cut -f2 | sort | uniq | awk '{print "all\t"$1}'
sort ~/res/genes/pLI09.txt |
  join -t $'\t' out/$pfx.genes.SetID - |
  cut -f2 | sort | uniq | awk '{print "lof\t"$1}'
join -a1 -a2 <(sort ~/res/genes/rbfox2.txt) <(sort ~/res/genes/fmrp.txt) |
  join -a1 -a2 - <(sort ~/res/genes/celf4.txt) |
  join -a1 -a2 - <(sort ~/res/genes/synaptome.txt) |
  join -t $'\t' out/$pfx.genes.SetID - |
  cut -f2 | sort | uniq | awk '{print "synaptic\t"$1}') > out/$pfx.sets.SetID

# perform SKAT analyses
mkdir -p burden
awk -F"\t" '{print "0",$1,$2}' $pheno > out/$pfx.phe
default="--bfile $pfx --min-mac 1 --exclude out/$pfx.ex --remove out/$pfx.rm --pheno out/$pfx.phe --covar out/$pfx.skat.cov --weights-beta 1 1"
for class in disruptive damaging deleterious missense; do
  cut -d" " -f1 $pfx.urv.rlist | sort | join out/$pfx.$class - > out/$pfx.$class.urv
  cut -d" " -f1 $pfx.urv.rlist | sort | join -a1 -v1 out/$pfx.$class - > out/$pfx.$class.not.urv
  join -a1 -v1 out/$pfx.$class out/$pfx.nonpsych > out/$pfx.$class.not.nonpsych

  for maxac in 1 5 10; do
    skat.py $default \
      --max-mac $maxac \
      --extract out/$pfx.$class \
      --setid out/$pfx.genes.SetID \
      --out burden/$pfx.genes.mac$maxac.$class
    skat.py $default \
      --max-mac $maxac \
      --extract out/$pfx.$class.not.urv \
      --setid out/$pfx.sets.SetID \
      --out burden/$pfx.sets.mac$maxac.$class
  done

  for maxaf in 001 005; do
    skat.py $default \
      --max-maf .$maxaf \
      --extract out/$pfx.$class \
      --setid out/$pfx.genes.SetID \
      --out burden/$pfx.genes.maf$maxaf.$class
    skat.py $default \
      --max-maf .$maxaf \
      --extract out/$pfx.$class.not.urv \
      --setid out/$pfx.sets.SetID \
      --out burden/$pfx.sets.maf$maxaf.$class
  done

  for grp in genes sets; do
    skat.py $default \
      --max-mac 1 \
      --extract out/$pfx.$class.urv \
      --setid out/$pfx.$grp.SetID \
      --out burden/$pfx.$grp.urv.$class
  done
done

# generate QQ plot for gene burden tests
for thr in urv mac1 mac5 mac10 maf001 maf005; do
  cut -d" " -f2 burden/$pfx.genes.$thr.deleterious.out | sed 's/^P\.value$/P/g' |
    qqplot.R - $pfx.genes.$thr.deleterious.png
done

convert -pointsize 20 -draw 'text 100,90 "Ultra-rare variants (URVs)"' swexm.genes.urv.deleterious.png qqplot.urv.png
convert -pointsize 20 -draw 'text 100,90 "Singletons"' swexm.genes.mac1.deleterious.png qqplot.mac1.png
convert -pointsize 20 -draw 'text 100,90 "Minor allele count ≤5"' swexm.genes.mac5.deleterious.png qqplot.mac5.png
convert -pointsize 20 -draw 'text 100,90 "Minor allele count ≤10"' swexm.genes.mac10.deleterious.png qqplot.mac10.png
convert -pointsize 20 -draw 'text 100,90 "Minor allele frequency <0.1%"' swexm.genes.maf001.deleterious.png qqplot.maf001.png
convert -pointsize 20 -draw 'text 100,90 "Minor allele frequency <0.5%"' swexm.genes.maf005.deleterious.png qqplot.maf005.png

convert qqplot.{urv,mac5,maf001}.png -append qqplot.left.png
convert qqplot.{mac1,mac10,maf005}.png -append qqplot.right.png
convert qqplot.{left,right}.png +append qqplot.png

# run R script for the analyses on the Sweden exome cohort
genovese2016.R

###########################################################################
## DE NOVO ANALYSES                                                      ##
###########################################################################

# generate VCF file from Menachem's table
wget https://fromem03.u.hpc.mssm.edu/denovos/annotated.LIT_de_novos.txt
(echo "##fileformat=VCFv4.1";
 echo "##INFO=<ID=DISEASE,Number=1,Type=String,Description=\"Disease type\">";
 echo "##INFO=<ID=STUDY,Number=1,Type=String,Description=\"Study type\">";
 echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
 gcol LOCUS REF ALT DISEASE STUDY < annotated.LIT_de_novos.txt | tail -n+2 |
   sed 's/^chr//g' | tr ':' '\t' | awk -v OFS="\t" '$3!~"N" && $4!~"N" {split($2,a,".")
   print $1,a[1],".",$3,$4,".",".","DISEASE="$5";STUDY="$6}' |
   sed 's/^X/23/g' | sort -k1,1n -k2,2n | sed 's/^23/X/g') |
   bgzip > annotated.LIT_de_novos.vcf.gz && tabix -f annotated.LIT_de_novos.vcf.gz

# annotate and generate table
vcf2annot.sh annotated.LIT_de_novos.vcf.gz annotated.LIT_de_novos ~/res/human_g1k_v37.fasta ~/res
annot2set.sh annotated.LIT_de_novos.annot.vcf.gz out/denovo ~/res
fld="DISEASE STUDY KGPhase3 RS EXACAC NONPSYCHAC"
annot2table.py -jar ~/res/SnpSift.jar -vcf annotated.LIT_de_novos.annot.vcf.gz -fld $fld -d out/denovo.damaging |
  bgzip > annotated.LIT_de_novos.annot.tsv.gz

join -a1 -a2 out/denovo.synonymous out/denovo.misbenign > out/denovo.nondeleterious
declare -A disease=( ["SSC"]="SIB" ["Fromer"]="SCZ" )
declare -A size=( ["SSC"]="1911" ["Fromer"]="617" )
for class in deleterious nondeleterious; do
  for study in Fromer SSC; do
      n=$(zcat annotated.LIT_de_novos.annot.tsv.gz |
	awk -F"\t" -v disease=${disease[$study]} -v study=$study '$5==disease && $6==study && $9=="."' |
	cut -f1-4 | tr '\t' ':' | sort |
        join - out/denovo.$class | wc -l)
    rate=$(bc -l <<< $n/${size[$study]})
    lowci=$(bc -l <<< \($n-1.96*sqrt\($n\)\)/${size[$study]})
    highci=$(bc -l <<< \($n+1.96*sqrt\($n\)\)/${size[$study]})
    echo "${disease[$study]} $study - $class - n: $n - rate: $rate - CI: $lowci - $highci"
  done
done
