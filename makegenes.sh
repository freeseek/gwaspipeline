#!/bin/bash
###
#  makegenes.sh - generate genesets analyzed in Genovese et al. 2016
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

# preparation steps
mkdir -p ~/res/genes/

# download USCS gene tables
wget -NP ~/res/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/{knownGene,ensGene,ensemblToGeneName}.txt.gz

# download dictionary for mouse genes
wget -O- ftp://ftp.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt |
  tail -n+2 | awk -F"\t" '{print $3"\t"$1}' |
  sort -t $'\t' -k1,1 > ~/res/genes/symbol2mgi.txt

# download dictionary for mouse genes
wget -O- ftp://ftp.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt |
  tail -n+2 | awk -F"\t" '$11!="null" {print $3"\t"$11}' |
  sort -t $'\t' -k1,1 > ~/res/genes/symbol2ensmusg.txt

# download relationship for mouse to human genes
wget -O- 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"?><!DOCTYPE Query><Query formatter = "TSV"><Dataset name = "hsapiens_gene_ensembl" ><Filter name = "with_homolog_mmus" excluded = "0"/><Attribute name = "ensembl_gene_id" /><Attribute name = "mmusculus_homolog_ensembl_gene" /></Dataset></Query>' |
  awk -F"\t" '{print $2"\t"$1}' |
  sort -t $'\t' -k1,1 > ~/res/genes/ensmusg2ensg.txt

# download relationship for mouse to human genes
wget -O- "http://www.genenames.org/cgi-bin/download?col=gd_mgd_id&col=md_ensembl_id&status=Approved&format=text&submit=submit" |
  tail -n+2 | awk -F"\t" '$1!="" {split($1,a,", "); for (i=1; i<=length(a); i++) print a[i]"\t"$2}' |
  sort -t $'\t' -k1,1 > ~/res/genes/mgi2ensg.txt

# download dictionary for human Ensembl genes
gzip -cd ~/res/ensemblToGeneName.txt.gz | cut -f2 | sort | uniq > ~/res/genes/ens.txt
gzip -cd ~/res/ensGene.txt.gz | cut -f2,13 | sort -t$'\t' -k1,1 > ~/res/genes/enst2ensg.txt
gzip -cd ~/res/ensemblToGeneName.txt.gz | sort -t$'\t' -k1,1 > ~/res/genes/enst2symbol.txt
join -t $'\t' ~/res/genes/enst2ensg.txt ~/res/genes/enst2symbol.txt | cut -f2- |
  sort -t $'\t' -k1,1 | uniq > ~/res/genes/ensg2symbol.txt

# download relationship for withdrawn human genes
wget -O- "http://www.genenames.org/cgi-bin/download?col=gd_prev_sym&col=gd_aliases&col=md_ensembl_id&status=Approved&format=text&submit=submit" |
  tail -n+2 | awk -F"\t" '$1!="" && $3!="" {split($1,a,", "); for (i=1; i<=length(a); i++) print a[i]"\t"$3}
  $2!="" && $3!="" {split($2,a,", "); for (i=1; i<=length(a); i++) print a[i]"\t"$3}' |
  sort -t $'\t' -k2,2 |
  join -1 2 -t $'\t' - ~/res/genes/ensg2symbol.txt | cut -f2- |
  sort -t $'\t' -k1,1 |
  join -t $'\t' -a1 -v1 - ~/res/genes/ens.txt > ~/res/genes/homo2symbol.txt

# download relationship for withdrawn mouse genes
wget -O- "ftp://ftp.informatics.jax.org/downloads/reports/MRK_List1.rpt" |
  awk -F"\t" '$9~"^withdrawn, = " {print $7"\t"substr($9,14)}' |
  sort -t $'\t' -k1,1 > ~/res/genes/mus2symbol.txt

# generate bed file with human Ensembl genes
gzip -cd ~/res/ensGene.txt.gz | awk -F"\t" '$7!=$8' | cut -f3,7,8,13 | sed 's/^chr//' | grep -v _ |
  sort -t $'\t' -k4,4 | join -t $'\t' -1 4 - ~/res/genes/ensg2symbol.txt | cut -f2- |
  bedtools sort | uniq > ~/res/ensGene.bed

# generate bed file with human Ensembl coding exons
gzip -cd ~/res/ensGene.txt.gz | awk '{split($10, a, ","); split($11, b, ",");
  for (i=1; i<length(a); i++) {from=$7>a[i]?$7:a[i]; to=$8<b[i]?$8:b[i];
  if (to>from) print $2"\t"$3"\t"from"\t"to}}' |
  sort -t $'\t' -k1,1 |
  join -t $'\t' - ~/res/genes/enst2symbol.txt |
  cut -f2- | bedtools sort | bedtools merge -c 4 -o collapse |
  awk '{split($4, a, ","); for (i=1; i<=length(a); i++) x[a[i]]++; for(i in x) print $1"\t"$2-2"\t"$3+2"\t"i; delete x}' |
  sed -e 's/chrM/MT/g' -e 's/^chr//g' | grep -v _ > ~/res/ensGene.coding.bed

# gene expression table (Fagerberg et al. 2014)
wget -P /tmp http://www.mcponline.org/content/suppl/2013/12/05/M113.035600.DC1/mcp.M113.035600-2.xlsx
xlsx2csv -d tab /tmp/mcp.M113.035600-2.xlsx | (read -r; printf "%s\n" "$REPLY"; sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/ensg2symbol.txt - | cut -f2- | sort -t $'\t' -k1,1 | uniq) \
  > ~/res/genes/fagerberg.tsv

# brain cell types expression table (Cahoy et al. 2008)
wget -P /tmp http://www.jneurosci.org/content/suppl/2008/01/03/28.1.264.DC1/335_Cahoy_S_Table_S3b_dChip_Data_Table_2007-09-11.xls
localc --nologo --convert-to xlsx --outdir /tmp /tmp/335_Cahoy_S_Table_S3b_dChip_Data_Table_2007-09-11.xls
xlsx2csv -d tab /tmp/335_Cahoy_S_Table_S3b_dChip_Data_Table_2007-09-11.xlsx | tail -n+2 | cut -f2-14,16 | tr -d '\r' |
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]=$2} NR>FNR && $14 in x {$14=x[$14]} NR>FNR {gene=$14; NF=13; print gene,$0}' ~/res/genes/mus2symbol.txt - |
  (read -r; printf "%s\n" "$REPLY"; sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/symbol2ensmusg.txt - | cut -f2- | sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/ensmusg2ensg.txt - | cut -f2- | sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/ensg2symbol.txt - | cut -f2- | sort -t $'\t' -k1,1) > ~/res/genes/cahoy.tsv

# neuron cell types expression table (Mo et al. 2015)
wget -P /tmp --user-agent="" http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627315004250/1-s2.0-S0896627315004250-mmc3.xlsx/272195/html/S0896627315004250/b58398a09cc84990c71014436d21a2e6/mmc3.xlsx
xlsx2csv -s 5 -d tab /tmp/mmc3.xlsx |
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]=$2} NR>FNR && $1 in x {$1=x[$1]} NR>FNR' ~/res/genes/mus2symbol.txt - |
  (read -r; printf "%s\n" "$REPLY"; sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/symbol2ensmusg.txt - | cut -f2- | sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/ensmusg2ensg.txt - | cut -f2- | sort -t $'\t' -k1,1 |
  join -t $'\t' ~/res/genes/ensg2symbol.txt - | cut -f2- | sort -t $'\t' -k1,1) > ~/res/genes/mo.tsv

# missense constrained gene list (Samocha et al. 2014)
wget -P /tmp http://www.nature.com/ng/journal/v46/n9/extref/ng.3050-S3.xls
localc --nologo --convert-to xlsx --outdir /tmp /tmp/ng.3050-S3.xls
xlsx2csv -s 3 -d tab /tmp/ng.3050-S3.xlsx | tail -n+2 | cut -f2 | sort |
  join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt |
  cut -f2 | sort | uniq > ~/res/genes/constrained.txt

# haplo-insufficient gene list (Lek et al. 2015)
wget -P /tmp ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
tail -n+2 /tmp/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | awk -F"\t" '$20>.9 {print $2}' |
  sort | join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt |
  cut -f2 | sort | uniq > ~/res/genes/pLI09.txt

# FMRP gene list (Darnell et al. 2011)
wget -O /tmp/mmc2.xls --user-agent="" http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0092867411006556/1-s2.0-S0092867411006556-mmc2.xls/272196/html/S0092867411006556/d8254825323f8fa37904a434e0350c2c/mmc2.xls
localc --nologo --convert-to xlsx --outdir /tmp /tmp/mmc2.xls
xlsx2csv -d tab -s 3 /tmp/mmc2.xlsx |
  tail -n+4 | awk -F"\t" '$29<.1 {print $3}' | sort |
  join -t $'\t' -a1 - ~/res/genes/mus2symbol.txt | cut -f2 | sort |  
  join -t $'\t' - ~/res/genes/symbol2ensmusg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensmusg2ensg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensg2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/fmrp.txt

# CELF4 gene list (Wagnon et al. 2012)
wget -O /tmp/journal.pgen.1003067.s004.XLSX "http://journals.plos.org/plosgenetics/article/asset?unique&id=info:doi/10.1371/journal.pgen.1003067.s004"
xlsx2csv -d tab -s 1 /tmp/journal.pgen.1003067.s004.XLSX |
  tail -n+2 | awk -F"\t" '$2>.2 {print $1}' | sort |
  join -t $'\t' - ~/res/genes/ensmusg2ensg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensg2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/celf4.txt

# RBFOX gene lists (Weyn et al. 2014)
wget -O /tmp/mmc2.xlsx --user-agent="" http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S2211124714000849/1-s2.0-S2211124714000849-mmc2.xlsx/280959/html/S2211124714000849/ce11d03f4c399c44bd5647aa242e6dc3/mmc2.xlsx
xlsx2csv -d tab -s 1 /tmp/mmc2.xlsx |
  tail -n+5 | awk -F"\t" '$13>=4 {print $10}' | sort |
  join -t $'\t' -a1 - ~/res/genes/mus2symbol.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/symbol2ensmusg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensmusg2ensg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensg2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/rbfox2.txt
xlsx2csv -d tab -s 1 /tmp/mmc2.xlsx |
  tail -n+5 | awk -F"\t" '$12+$14>=12 {print $10}' | sort |
  join -t $'\t' -a1 - ~/res/genes/mus2symbol.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/symbol2ensmusg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensmusg2ensg.txt | cut -f2 | sort |
  join -t $'\t' - ~/res/genes/ensg2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/rbfox13.txt

# SynaptomeDB gene list (Pirooznia et al. 2012)
for cat in pre vesicle preactivezone post; do
  wget -O /tmp/cat_$cat.xls http://metamoodics.org/SynaptomeDB/syncats/cat_${cat}_export.php?type=excel
  localc --nologo --convert-to xlsx --outdir /tmp /tmp/cat_$cat.xls
done
(for cat in pre vesicle preactivezone post; do xlsx2csv -d tab /tmp/cat_$cat.xlsx | tail -n+2; done) |
  cut -f2 | sort | join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt |
  cut -f2 | sort | uniq > ~/res/genes/synaptome.txt

# PSD and PSD-95 gene lists (Bayes et al. 2011)
for i in 49 69; do
  if [ $i == 49 ]; then name="psd95";
  elif [ $i == 69 ]; then name="psd"; fi
  wget -O- http://www.genes2cognition.org/db/GeneList/L000000$i |
    grep "<td><em>.*</em></td>" |
    sed 's/^ *<td><em>\(.*\)<\/em><\/td> *$/\1/g' |
    sort | join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt |
    cut -f2 | sort | uniq > ~/res/genes/$name.txt
done

# NMDAR and ARC gene list (Kirov et al. 2012)
wget -P /tmp http://www.nature.com/mp/journal/v17/n2/extref/mp2011154x1.pdf
pdftk /tmp/mp2011154x1.pdf cat 39 40 output /dev/stdout | pdftotext - - |
  grep ^[A-Z][A-Z0-9]*$ | sort |
  join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt |
  cut -f2 | sort | uniq > ~/res/genes/nmdarc.txt

# MIR-137 gene list (Betel et al. 2010)
wget -P /tmp/ http://cbio.mskcc.org/microrna_data/human_predictions_S_C_aug2010.txt.gz
zcat /tmp/human_predictions_S_C_aug2010.txt.gz |
  awk -F"\t" '$2=="hsa-miR-137" {print $4}' |
  sort | uniq | join -a1 -t $'\t' - ~/res/genes/homo2symbol.txt |
  cut -f2 | sort | uniq > ~/res/genes/mir137.txt

# CHD8 list (Cotney et al. 2015)
wget -P /tmp/ http://www.nature.com/article-assets/npg/ncomms/2015/150310/ncomms7404/extref/ncomms7404-s2.xlsx
xlsx2csv -d tab -s1 /tmp/ncomms7404-s2.xlsx | tr -d '\r' | tail -n+2 | cut -f5 |
  sort | join -a1 -t $'\t' - ~/res/genes/homo2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/chd8.hNSC.txt
xlsx2csv -d tab -s2 /tmp/ncomms7404-s2.xlsx | tr -d '\r' | tail -n+2 | cut -f5 |
  sort | join -a1 -t $'\t' - ~/res/genes/homo2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/chd8.hNSC_specific.txt
xlsx2csv -d tab -s3 /tmp/ncomms7404-s2.xlsx | tr -d '\r' | tail -n+2 | cut -f5 |
  sort | join -a1 -t $'\t' - ~/res/genes/homo2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/chd8.human_brain.txt
xlsx2csv -d tab -s4 /tmp/ncomms7404-s2.xlsx | tr -d '\r' | tail -n+2 | cut -f5 |
  sort | join -a1 -t $'\t' - ~/res/genes/homo2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/chd8.hNSC+human_brain.txt
xlsx2csv -d tab -s5 /tmp/ncomms7404-s2.xlsx | tr -d '\r' | tail -n+2 | cut -f5 |
  sort | join -a1 -t $'\t' - ~/res/genes/homo2symbol.txt | cut -f2 | sort | uniq > ~/res/genes/chd8.hNSC+human+mouse.txt

# PGC GWAS gene list (Ripke et al. 2014)
wget -P /tmp https://www.med.unc.edu/pgc/files/resultfiles/scz2.regions.zip
unzip -p /tmp/scz2.regions.zip scz2.anneal.108.txt |
  tail -n+2 | cut -f2,5,6 | sed 's/^chr//g' | awk '{print $1"\t"$2"\t"$3}' |
  bedtools intersect -wb -a ~/res/ensGene.bed -b - | cut -f4- |
  sort -t $'\t' -k2,2n -k3,3n | uniq |
  awk '{x[$2":"$3":"$4]++; y[$2":"$3":"$4]=y[$2":"$3":"$4]"\n"$1}
  END {for (i in x) if (x[i]<=4) print y[i]}' | grep -v ^$ |
  sort | uniq > ~/res/genes/gwas.txt

# X-linked intellectual disability gene lists
for chr in {1..22} X Y; do
  for i in {0..7}; do
    wget -O- --user-agent="M" "http://omim.org/geneMap/$chr?start=$((200*i+1))&limit=$((200*(i+1)))&format=tsv" |
      tail -n+5 | head -n-12
  done
done > /tmp/omim.txt
grep -i "mental retardation" /tmp/omim.txt | grep -v ^X | cut -f3 | cut -d, -f1 | sort |
  join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt | cut -f2 | sort |
  join - ~/res/genes/ens.txt | uniq > ~/res/genes/alid.txt
grep -i "mental retardation" /tmp/omim.txt | grep ^X | cut -f3 | cut -d, -f1 | sort |
  join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt | cut -f2 | sort |
  join - ~/res/genes/ens.txt | uniq > ~/res/genes/xlid.omim.txt
wget -O- http://www.ggc.org/images/Expanded_NGS_XLID_Panel_Gene_List.pdf |
  pdftotext - - | grep [A-Z][A-Z0-9] | grep -v ^Expanded\\\|^$ > ~/res/genes/xlid.gcc.txt
wget -O- http://dnatesting.uchicago.edu/tests/x-linked-non-specific-intellectual-disability-sequencing-panel |
  grep taxonomy | cut -d\> -f3 | cut -d\< -f1 | grep -iv "disability\|^$" |
  sed -e 's/NLGN4Z/NLGN4X/g' -e 's/PHF6PHF8/PHF6\nPHF8/g' > ~/res/genes/xlid.chicago.txt
cat ~/res/genes/xlid.{omim,gcc,chicago}.txt | sort | uniq > ~/res/genes/xlid.txt
wget -P /tmp http://europepmc.org/articles/PMC4053723/bin/gb-2013-14-11-r122-S7.xlsx
xlsx2csv -d tab /tmp/gb-2013-14-11-r122-S7.xlsx |
  awk -F"\t" '$5~"escape" {print $1}' | sort |
  join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt | cut -f2 |
  sort | uniq > ~/res/genes/x.escape.txt

# Developmental disorder gene list (McRae et al. 2016)
wget -P /tmp/ http://www.biorxiv.org/highwire/filestream/13630/field_highwire_adjunct_files/3/049056-4.xlsx
xlsx2csv -d tab /tmp/049056-4.xlsx | tail -n+2 | cut -f1 | sort > ~/res/genes/dd.txt

# SNP and indels denovo gene lists (Fromer et al. 2014)
wget -P /tmp https://fromem03.u.hpc.mssm.edu/denovos/annotated.LIT_de_novos.txt
for disease in AUT CHD EPI ID SCZ; do
  tail -n+2 /tmp/annotated.LIT_de_novos.txt | cut -f7,10,12 |
    awk '$2~"codon" || $2=="frameshift" || $2=="missense" || $2=="nonsense" || $2=="readthrough" || $2~"splice" || $2=="start-lost"' |
    grep ^$disease | cut -f3 | tr ',' '\n' | sort |
    join -t $'\t' -a1 - ~/res/genes/homo2symbol.txt |
    cut -f2 | sort | uniq > ~/res/genes/denovo.${disease,,}.txt
done

# CNV denovo gene lists (Kirov et al. 2012)
wget -P /tmp http://www.nature.com/neuro/journal/v19/n11/extref/nn.4402-S5.xlsx
for class in loss gain; do
  for disease in ASD BD SCZ; do
    xlsx2csv -d tab /tmp/nn.4402-S5.xlsx |
      awk -F"\t" -v class=$class -v disease=$disease '$2==disease && $3==class {print $4"\t"$5"\t"$6}' |
      bedtools window -a - -b ~/res/ensGene.bed | cut -f7 | sort | uniq > ~/res/genes/denovo.$class.${disease,,}.txt
  done
done
