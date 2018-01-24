#!/bin/bash
###
#  install.sh - install basic program and resources for GWAS analyses
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

# install basic tools (Debian specific)
sudo apt-get install bedtools gzip libreoffice-calc pdftk poppler-utils python3 python3-pandas r-base-core r-cran-ggplot2 r-cran-mnormt r-cran-xml samtools tabix unzip wget xlsx2csv

# preparation steps
mkdir -p ~/bin/ ~/res/ && cd /tmp/

# install latest version of bcftools
git clone --branch=develop git://github.com/samtools/bcftools.git
git clone --branch=develop git://github.com/samtools/htslib.git
cd htslib && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/bcftools ~/bin/

# install latest development version of plink
wget https://www.cog-genomics.org/static/bin/plink/plink_linux_x86_64_dev.zip
unzip -od ~/bin/ plink_linux_x86_64_dev.zip plink

# install latest version of GCTA
wget http://cnsgenomics.com/software/gcta/gcta_1.25.2.zip
unzip -od ~/bin/ gcta_1.25.2.zip gcta64

# install latest version of Eagle
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.3.4.tar.gz
tar xzvf Eagle_v2.3.4.tar.gz -C ~/bin/ Eagle_v2.3.4/eagle --strip-components=1
tar xzvf Eagle_v2.3.4.tar.gz -C ~/res/ Eagle_v2.3.4/tables/genetic_map_hg19_withX.txt.gz --strip-components=2

# install latest version of hapi-ur (and fix non-descending order issue)
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/hapi-ur/hapi-ur-1.01.tgz
tar xzvf hapi-ur-1.01.tgz -C ~/bin/ hapi-ur-1.01/hapi-ur hapi-ur-1.01/insert-map.pl --strip-components=1
sed -i -e 's/<= $last_pos/< $last_pos/g' -e 's/ascending/non-descending/g' ~/bin/insert-map.pl

# install latest version of shapeit
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
tar xzvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz -C ~/bin/ bin/shapeit --strip-components=1

# install latest version of impute2
wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
tar xzvf impute_v2.3.2_x86_64_static.tgz -C ~/bin/ impute_v2.3.2_x86_64_static/impute2 --strip-components=1

# install genetic maps for chromsome X
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
tar xzvf 1000GP_Phase3_chrX.tgz -C ~/res/ genetic_map_chrX_{nonPAR,PAR1,PAR2}_combined_b37.txt

# install latest version of data.table
wget https://cran.mtu.edu/src/contrib/data.table_1.9.6.tar.gz
R CMD INSTALL data.table_1.9.6.tar.gz

# install latest version of tidyr
wget https://cran.r-project.org/src/contrib/tidyr_0.4.0.tar.gz
R CMD INSTALL tidyr_0.4.0.tar.gz

# install latest version of psych
wget https://cran.r-project.org/src/contrib/psych_1.5.8.tar.gz
R CMD INSTALL psych_1.5.8.tar.gz

# install latest version of broom
wget https://cran.r-project.org/src/contrib/broom_0.4.0.tar.gz
R CMD INSTALL broom_0.4.0.tar.gz

# install latest version of SKAT
wget https://cran.r-project.org/src/contrib/SKAT_1.1.2.tar.gz
R CMD INSTALL SKAT_1.1.2.tar.gz

# install latest version of snpEff/SnpSift
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip -ojd ~/res/ snpEff_latest_core.zip snpEff/{snpEff.{config,jar},SnpSift.jar}

# install latest version of Beagle
wget https://faculty.washington.edu/browning/beagle/beagle.27Jul16.86a.jar
ln -s beagle.27Jul16.86a.jar ~/res/beagle.jar

# download human genome reference (GRCh37)
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz |
  gzip -d > ~/res/human_g1k_v37.fasta
samtools faidx ~/res/human_g1k_v37.fasta

# download UCSC gene tables
wget -NP ~/res/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/{knownGene,ensGene,ensemblToGeneName}.txt.gz

# download and process dbSNP data
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/00-All{,_papu}.vcf.{gz,gz.md5,gz.tbi}
(tabix -h 00-All.vcf.gz {1..22}
tabix 00-All_papu.vcf.gz X:60000-2699520 | awk -F"\t" -v OFS="\t" '{split($8, x, ";"); $3="rs"substr(x[1], 4); print}'
tabix 00-All.vcf.gz X:2699520-154931044
tabix 00-All_papu.vcf.gz X:154931044-155260560 | awk -F"\t" -v OFS="\t" '{split($8, x, ";"); $3="rs"substr(x[1], 4); print}'
tabix 00-All.vcf.gz Y MT) |
  bgzip > ~/res/dbsnp.vcf.gz &&
  tabix -f ~/res/dbsnp.vcf.gz

# download list of 1000 Genomes project phase 3 list of variants
wget -NP ~/res/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz{,tbi}

# download HRC data (http://www.haplotype-reference-consortium.org/site)

# download gnomAD data (see http://gnomad.broadinstitute.org/downloads)

# download and normalize ExAC data version 0.3
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/{ExAC.r0.3.sites.vep,subsets/ExAC.r0.3.nonpsych.sites}.vcf.gz{,.tbi}
for sfx in sites.vep nonpsych.sites; do
  bcftools annotate -Ou -x ID,INFO/AC_Het,INFO/Het_AFR,INFO/Het_AMR,INFO/Het_EAS,INFO/Het_FIN,INFO/Het_NFE,INFO/Het_OTH,INFO/Het_SAS ExAC.r0.3.$sfx.vcf.gz |
    bcftools norm -Ou -m -any |
    bcftools norm -f ~/res/human_g1k_v37.fasta |
    grep -v ^##bcftools | bgzip > ~/res/ExAC.r0.3.$sfx.norm.vcf.gz
  tabix -f ~/res/ExAC.r0.3.$sfx.norm.vcf.gz
done

# download dbNSFP data version 2.9
wget -O ~/res/dbNSFP.txt.gz https://www.googledrive.com/host/0B7Ms5xMSFMYlSTY5dDJjcHVRZ3M
wget -O ~/res/dbNSFP.txt.gz.tbi https://www.googledrive.com/host/0B7Ms5xMSFMYlOTV5RllpRjNHU2s

# download 1000 Genomes project phase 1 genotype data
pfx="http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase1_vcf/chr"
sfx=".1kg.ref.phase1_release_v3.20101123.vcf.gz"
wget -P ~/res/chrs/ $pfx{{1..22},X}$sfx
for chr in {1..22}; do vcf2plink.py --vcf ~/res/chrs/chr$chr$sfx --ref ~/res/human_g1k_v37.fasta --out ~/res/chrs/kgp.chr$chr; done
vcf2plink.py --vcf ~/res/chrs/$pfx{X}$sfx --ref $ref --out ~/res/chrs/kgp.chrX --impute-sex .9 .9
/bin/cp ~/res/chrs/kgp.chrX.fam ~/res/kgp.fam
cat ~/res/chrs/kgp.chr{{1..22},X}.bim > ~/res/kgp.bim
(echo -en "\x6C\x1B\x01"; for chr in {1..22} X; do tail -c +4 ~/res/chrs/kgp.chr$chr.bed; done) > ~/res/kgp.bed
wget -P ~/res/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
(echo -e "IID\tPOP"; cut -f2,7 ~/res/20130606_g1k.ped | tail -n+2 | sort) > ~/res/kgp.pop
wget -P ~/res/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110426_omni_phased_vcfs/Omni25_genotypes_1362_samples_v2.b37.vcf.gz{,.tbi}
bcftools query -f "%CHROM:%POS:%REF:%ALT\n" ~/res/Omni25_genotypes_1362_samples_v2.b37.vcf.gz | sort > ~/res/kgp.omni25

# build a VCF file with all CpG variants for GRCh37
makecpg.sh
