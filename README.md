gwaspipeline
============

This repositor contains a set of bash and python scripts used for GWAS analyses mostly based on plink.

If you are looking for scripts related to the Sweden exome dataset, you are mostly interested in the following master scripts:

```
install.sh
makegenes.sh
genovese2016.sh
genovese2016.R
```

Current version was updated on October 3rd 2016.

csv2xlsx.py
===========

Basic script to convert text tables into XLSX tables (supports multiple sheets)

genome2ped.py
=============

Identifies trios from a plink genome file and optionally generates PI_HAT x Z0 plot

install.sh
==========

Installs several resources useful for genetic analyses in a ~/res/ directory in your system

markerqc.sh
===========

Perform basic quality control for markers in a plink dataset

plink2lite.sh
=============

Extracts a core set of highly confident non-rare 1000 Genomes project variants to analyze a given dataset befor marker QC

plink2score.py
==============

Compute several polygenic scores from either binary or dosage plink files and optionally adjusts them for covariates (e.g. principal components)

makegenes.sh
============

This scripts installs and generates several gene tables and lists useful for enrichment analyses

sampleqc.sh
===========

Performs basic quality control for samples in a plink dataset

skat.py
=======

Wrapper for performing burden analysis with SKAT (see http://www.hsph.harvard.edu/skat/) with several plink-like options

vcf2annot.sh
============

This script annotates a given VCF using the snpEff/SnpSift software and the dbNSFP database

vcf2plink.py
============

This scripts converts VCF to plink following best practices (see http://bit.ly/1sgOcuP)