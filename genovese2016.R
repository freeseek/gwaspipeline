#!/usr/bin/env Rscript
###
#  genovese2016.R - script for the analyses on the Sweden exome cohort
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

###########################################################################
## LOAD NECESSARY LIBRARIES                                              ##
###########################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

datapath <- paste0(path.expand("~"), .Platform$file.sep, 'swe', .Platform$file.sep, 'data', .Platform$file.sep)
figpath <- paste0(path.expand("~"), .Platform$file.sep, 'swe', .Platform$file.sep, 'figs', .Platform$file.sep)
respath <- paste0(path.expand("~"), .Platform$file.sep, 'res', .Platform$file.sep)
na.zero <- function(x) {replace(x, is.na(x), 0)}
figs <- list() # list of figures
npcs <- 20

###########################################################################
## LOAD DATA ABOUT ULTRA-RARE VARIANTS                                   ##
###########################################################################

setid <- data.frame(setNames(fread(paste0(datapath, 'swe.genes.SetID'), header = FALSE), c('gene', 'name')))
var <- data.frame(setNames(fread(paste0(datapath, 'swe.urv.tsv'), header = FALSE), c('name', 'iid')))

# create individuals table, add information about phenotpye, and information about PCs
fam <- setNames(read.table(paste0(datapath, 'swe.urv.fam'), header = FALSE)[,c('V2','V5')], c('iid', 'sex'))
fam <- merge(fam, setNames(read.table(paste0(datapath, 'pheno.tsv'), header = FALSE), c('iid', 'scz')), all.x = TRUE)
fam$scz <- replace(fam$scz, fam$scz==-9, NA)
fam <- merge(fam, setNames(read.table(paste0(datapath, 'swe.lite.kgp.pc.tsv'), header = TRUE)[,c('IID',paste0('PC',1:20))], c('iid',paste0('pc',1:20))), all.x = TRUE)

# extract kit used
kit <- read.table(paste0(datapath, 'kit.tsv'), header = FALSE)
fam$kit <- is.element(fam$iid, kit$V1[kit$V2==0])
rm(kit)

# extract birth year
birth <- setNames(read.table(paste0(datapath, 'birth.tsv'), header = TRUE), c('iid', 'birth'))
fam <- merge(fam, birth, all.x = TRUE)
rm(birth)

# extract individuals with mismatching sex
sexrm <- read.table(paste0(datapath, 'swe.sex.rm'), header = FALSE)$V2

###########################################################################
## LOAD DATA ABOUT ANNOTATIONS                                           ##
###########################################################################

for (annot in c('all', 'known', 'inframe', 'disruptive', 'missense', 'misdamaging', 'damaging', 'misbenign', 'synonymous', 'deleterious')) {
  var[,annot] <- is.element(var$name, fread(paste0(datapath, 'swe.', annot), sep = '\t', header = FALSE)$V1)
  fam[,paste0('count_', tolower(annot))] <- unname(na.zero(table(var$iid[var[,annot]])[as.character(fam$iid)]))
}
var[,'noncoding'] <- !(var[,'synonymous'] | var[,'misbenign'] | var[,'damaging'] | var[,'deleterious'])
fam[,'count_noncoding'] <- unname(na.zero(table(var$iid[var[,'noncoding']])[as.character(fam$iid)]))

idx <- grepl('[ACGT][ACGT]', var$name)
fam[,'count_indels'] <- na.zero(table(var$iid[idx])[as.character(fam$iid)])
fam[,'count_snps'] <- na.zero(table(var$iid[!idx])[as.character(fam$iid)])
bad <- fam$count_snps > 400  | fam$count_indels > 80 | is.element(fam[,'iid'], c('NA12878', 'NA12891' ,'NA12892'))
keep <- !bad & fam[,'count_all'] <= 100
keepnobp <- keep & !is.na(fam$scz) & !fam$iid %in% sexrm

# extract previous study information
prev <- setNames(read.table(paste0(datapath, 'prev.tsv'), header = TRUE), c('iid', 'prev'))
prev <- is.element(fam$iid, prev$iid[prev$prev==1])

###########################################################################
## ENRICHMENT FUNCTION GIVEN A LIST OF GENESETS                          ##
###########################################################################

get_tidy <- function(model, term, exponentiate = FALSE) {
  estimate <- unname(coef(model)[term])
  std.error <- unname(sqrt(diag(vcov(model)))[term])
  statistic <- unname(coef(model)[term] / sqrt(diag(vcov(model)))[term])
  p.value <- min(pnorm(statistic), pnorm(statistic, lower.tail = FALSE)) * 2
  conf.low <- unname(estimate - 1.96 * std.error)
  conf.high <- unname(estimate + 1.96 * std.error)
  if (exponentiate) {
    return(list(term = term, estimate = exp(estimate), std.error = std.error, statistic = statistic, p.value = p.value, conf.low = exp(conf.low), conf.high = exp(conf.high)))
  } else {
    return(list(term = term, estimate = estimate, std.error = std.error, statistic = statistic, p.value = p.value, conf.low = conf.low, conf.high = conf.high))
  }
}

get_enrichment <- function(setid, var, fam, sets, ctrl, funclasses, flag = FALSE, npcs = 20) {
  df <- data.frame(set = character(),
                   genecount = integer(),
                   funclass = character(),
                   type = character(),
                   estimate = double(),
                   conf.low = double(),
                   conf.high = double(),
                   statistic = double(),
                   p.value = double(),
                   p.adjusted = double(),
                   stringsAsFactors = FALSE)
  
  for (funclass in funclasses) {
    idx <- is.element(var$name, setid$name[is.element(setid$gene,ctrl)]) & var[,funclass]
    fam$count_ctrl <- na.zero(table(var$iid[idx])[as.character(fam$iid)])
    ormod <- glm(as.formula(paste0("(scz-1) ~ count_ctrl + count_all + sex + birth + kit + ", paste0("pc", 1:npcs, collapse = " + "))), family=binomial, fam)

    for (set in names(sets)) {
      geneset <- intersect(sets[[set]], ctrl)
      genecount = length(geneset)
      idx <- is.element(var$name, setid$name[is.element(setid$gene,geneset)]) & var[,funclass]
      if (all(!idx)) next
      fam$count <- na.zero(table(var$iid[idx])[as.character(fam$iid)])
      fit <- glm(as.formula(paste0("(scz-1) ~ count + count_all + sex + birth + kit + ", paste0("pc", 1:npcs, collapse = " + "))), family=binomial, fam)
      stats <- get_tidy(fit, 'count', exponentiate = TRUE)
      if (all(fam$count == fam$count_ctrl)) {
        p <- NA
      } else {
        fit <- glm(as.formula(paste0("(scz-1) ~ count + count_ctrl + count_all + sex + birth + kit + ", paste0("pc", 1:npcs, collapse = " + "))), family=binomial, fam)
        cmp <- anova(ormod, fit, test='Chisq')
        p <- cmp$Pr[2]
      }
      df <- rbind(df, data.frame(set=set, genecount=genecount, funclass=funclass, type='or', estimate=stats[['estimate']], conf.low=stats[['conf.low']], conf.high=stats[['conf.high']], statistic=stats[['statistic']], p.value=stats[['p.value']], p.adjusted=p, row.names=NULL))
      if (flag) {
        fit <- glm(as.formula(paste0("count ~ scz + count_all + sex + birth + kit + ", paste0("pc", 1:npcs, collapse = " + "))), family=gaussian, fam)
        stats <- get_tidy(fit, 'scz')
        if (all(fam$count == fam$count_ctrl)) {
          p <- NA
        } else {
          excmod <- glm(as.formula(paste0("count ~ count_ctrl + count_all + sex + birth + kit + ", paste0("pc", 1:npcs, collapse = " + "))), family=gaussian, fam)
          fit <- glm(as.formula(paste0("count ~ scz + count_ctrl + count_all + sex + birth + kit + ", paste0("pc", 1:npcs, collapse = " + "))), family=gaussian, fam)
          cmp <- anova(excmod, fit, test='Chisq')
          p <- cmp$Pr[2]
        }
        df <- rbind(df, data.frame(set=set, genecount=genecount, funclass=funclass, type='exc', estimate=stats[['estimate']], conf.low=stats[['conf.low']], conf.high=stats[['conf.high']], statistic=stats[['statistic']], p.value=stats[['p.value']], p.adjusted=p, row.names=NULL))
      }
    }
  }

  df$tick <- factor(paste0(df$set, '\n(', prettyNum(df$genecount, big.mark=","), ' genes)'), levels = unique(paste0(df$set, '\n(', prettyNum(df$genecount, big.mark=","), ' genes)')))
  return(df)
}

###########################################################################
## LOAD GENE SETS                                                        ##
###########################################################################

ensgenes <- unique(fread(paste0(respath, 'ensGene.bed'), sep="\t", header = FALSE)$V4)

genesets <- list()
for (set in c('brain', 'neurons', 'xlid', 'gwas', 'mir137', 'psd95', 'psd', 'nmdarc', 'fmrp', 'celf4', 'rbfox2', 'rbfox13', 'synaptome', 'constrained', 'pLI09', 'dd')) {
  genesets[[set]] <- fread(paste0(respath, 'genes', .Platform$file.sep, set, '.txt'), header = FALSE)$V1
}

for (set in c(paste0('denovo.', c('aut', 'chd', 'epi', 'id', 'scz', paste0('loss.', c('asd', 'bd', 'scz')), paste0('gain.', c('asd', 'bd', 'scz')))))) {
  genesets[[set]] <- fread(paste0(respath, 'genes', .Platform$file.sep, set, '.txt'), header = FALSE)$V1
}

###########################################################################
## ANALYSES NECESSARY FOR THE PAPER TEXT                                 ##
###########################################################################

print(paste('Total number of ultra-rare SNPs in wave 6 problematic sample:', fam$count_snps[fam$count_snps > 1000]))
print(paste('Total number of ultra-rare indels in wave 6 problematic sample:', fam$count_indels[fam$count_snps > 1000]))
print(paste('Total number of ultra-rare SNPs in wave 1 problematic samples:', paste(fam$count_snps[fam$count_indels > 80], collapse = ' ')))
print(paste('Total number of ultra-rare indels in wave 1 problematic samples:', paste(fam$count_indels[fam$count_indels > 80], collapse = ' ')))
print(paste('Total number of ultra-rare SNPs in wave 11/12 problematic samples:', paste(fam$count_snps[fam$count_snps > 400 & fam$count_snps < 600], collapse = ' ')))
print(paste('Total number of ultra-rare indels in wave 11/12 problematic samples:', paste(fam$count_indels[fam$count_snps > 400 & fam$count_snps < 600], collapse = ' ')))
print(paste('Range of URVs across non-problematic samples:', min(fam$count_all[!bad]), max(fam$count_all[!bad]), collapse = ' '))
print(paste('Range of ultra-rare SNPs across non-problematic samples:', min(fam$count_snps[!bad]), max(fam$count_snps[!bad]), collapse = ' '))
print(paste('Range of ultra-rare indels across non-problematic samples:', min(fam$count_indels[!bad]), max(fam$count_indels[!bad]), collapse = ' '))

print(paste('Total number of SCZ samples used for URVs:', sum(fam$scz==2, na.rm=TRUE)))
print(paste('Total number of male SCZ samples used for URVs:', sum(fam$scz==2 & fam$sex==1, na.rm=TRUE)))
print(paste('Total number of female SCZ samples used for URVs:', sum(fam$scz==2 & fam$sex==2, na.rm=TRUE)))

print(paste('Total number of CTRL samples used for URVs:', sum(fam$scz==1, na.rm=TRUE)))
print(paste('Total number of male CTRL samples used for URVs:', sum(fam$scz==1 & fam$sex==1, na.rm=TRUE)))
print(paste('Total number of female CTRL samples used for URVs:', sum(fam$scz==1 & fam$sex==2, na.rm=TRUE)))

print(paste('Total number of BD samples used for URVs:', sum(!bad & is.na(fam$scz))))
print(paste('Total number of male BD samples used for URVs:', sum(!bad & is.na(fam$scz) & fam$sex==1)))
print(paste('Total number of female BD samples used for URVs:', sum(!bad & is.na(fam$scz) & fam$sex==2)))

print(paste('Total number of samples used for enrichment analyses:', sum(keepnobp)))
print(paste('Total number of CTRL samples used for enrichment analyses:', sum(keepnobp & fam$scz==1, na.rm=TRUE)))
print(paste('Total number of SCZ samples used for enrichment analyses:', sum(keepnobp & fam$scz==2, na.rm=TRUE)))
print(paste('Total number of samples removed due to sex mismatch or Klinefelter syndrome:', sum(is.element(fam$iid,sexrm) & !bad)))
print(paste('Total number of samples removed due to excess URVs on top of problematic samples:', sum(!keep & !is.element(fam$iid,sexrm) & !bad)))

print(paste('Mean number of dURVS in controls:', mean(fam$count_deleterious[keepnobp & fam$scz==1], na.rm = TRUE)))
print(paste('Mean number of dURVS in cases:', mean(fam$count_deleterious[keepnobp & fam$scz==2], na.rm = TRUE)))
print(paste('Median number of dURVS in controls:', median(fam[keepnobp & fam$scz==1, 'count_deleterious'], na.rm = TRUE)))
print(paste('Median number of dURVS in cases:', median(fam[keepnobp & fam$scz==2, 'count_deleterious'], na.rm = TRUE)))

print(paste('Total number of URVs:', sum(var$iid %in% fam$iid[keepnobp])))
print(paste('Total number of coding URVs:', sum((var$synonymous | var$misbenign | var$damaging | var$disruptive) & var$iid %in% fam$iid[keepnobp])))
print(paste('Total number of synonymous URVs:', sum(var$synonymous & var$iid %in% fam$iid[keepnobp])))
print(paste('Total number of missense non-damaging URVs:', sum(var$misbenign & var$iid %in% fam$iid[keepnobp])))
print(paste('Total number of damaging URVs:', sum(var$damaging & var$iid %in% fam$iid[keepnobp])))
print(paste('Total number of disruptive URVs:', sum(var$disruptive & var$iid %in% fam$iid[keepnobp])))
print(paste('Total number of dURVs:', sum(var$deleterious & var$iid %in% fam$iid[keepnobp])))

fit <- lm(as.formula(paste0('count_synonymous ~ scz + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
print(paste('Exome-wide excess of synonymous URVs:', signif(x$estimate,2), signif(x$conf.low,2), signif(x$conf.high,2), signif(x$p.value,2)))

fit <- lm(as.formula(paste0('count_misbenign ~ scz + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
print(paste('Exome-wide excess of missense non-damaging URVs:', signif(x$estimate,2), signif(x$conf.low,2), signif(x$conf.high,2), signif(x$p.value,2)))

fit <- lm(as.formula(paste0('count_damaging ~ scz + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
print(paste('Exome-wide excess of damaging URVs:', signif(x$estimate,2), signif(x$conf.low,2), signif(x$conf.high,2), signif(x$p.value,2)))

fit <- lm(as.formula(paste0('count_disruptive ~ scz + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
print(paste('Exome-wide excess of disruptive URVs:', signif(x$estimate,2), signif(x$conf.low,2), signif(x$conf.high,2), signif(x$p.value,2)))

fit <- lm(as.formula(paste0('count_deleterious ~ scz + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
print(paste('Exome-wide excess of dURVs:', signif(x$estimate,2), signif(x$conf.low,2), signif(x$conf.high,2), signif(x$p.value,2)))

fit <- glm(as.formula(paste0('(scz-1) ~ count_deleterious + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), family = binomial(link = "logit"), fam[keepnobp,])
x <- get_tidy(fit, 'count_deleterious', exponentiate = TRUE)
print(paste('Odds ratios of dURVs:', signif(x$estimate,3), signif(x$conf.low,3), signif(x$conf.high,3), signif(x$p.value,2)))

for (class in c('all', 'deleterious')) {
  if (class=='all') { opt = '' } else { opt = 'count_all + ' }
  m2 <- lm(as.formula(paste0('count_', class, ' ~ ', opt, 'kit + birth + sex + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
  print(paste('Variance on the number of URVs/dURVs explained by all covariates:', 100*round(summary(m2)$r.squared, 4), '%'))

  if (class=='deleterious') {
    m1 <- lm(as.formula(paste0('count_', class, ' ~ kit + birth + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
    print(paste('Variance uniquely explained by all URV count:', 100*round(summary(m2)$r.squared - summary(m1)$r.squared, 4), '%'))
  }

  m1 <- lm(as.formula(paste0('count_', class, ' ~ ', opt, 'kit + birth + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
  print(paste('Variance uniquely explained by sex:', 100*round(summary(m2)$r.squared - summary(m1)$r.squared, 4), '%'))
  
  m1 <- lm(as.formula(paste0('count_', class, ' ~ ', opt, 'kit + sex + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
  print(paste('Variance uniquely explained by birth year:', 100*round(summary(m2)$r.squared - summary(m1)$r.squared, 4), '%'))
  
  m1 <- lm(as.formula(paste0('count_', class, ' ~ ', opt, 'birth + sex + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
  print(paste('Variance uniquely explained by enrichment kit:', 100*round(summary(m2)$r.squared - summary(m1)$r.squared, 4), '%'))
  
  for (i in 1:20) {
    m1 <- lm(as.formula(paste0('count_', class, ' ~ ', opt, 'scz + kit + birth + sex + ', paste0('pc', (1:npcs)[1:npcs!=i], collapse = ' + '))), fam[keepnobp,])
    print(paste('Variance uniquely explained by', paste0('pc',i), ':', 100*round(summary(m2)$r.squared - summary(m1)$r.squared, 4), '%'))
  }
}

NagelkerkeR2 <- function (rr) 
{
  n <- nrow(rr$model)
  R2 <- (1 - exp((rr$dev - rr$null)/n))/(1 - exp(-rr$null/n))
  RVAL <- list(N = n, R2 = R2)
  return(RVAL)
}

geneset <- union(union(union(genesets[['rbfox2']], genesets[['fmrp']]), genesets[['celf4']]), genesets[['synaptome']])
geneset_idx <- is.element(var$name, setid$name[is.element(setid$gene, geneset)])
for (annot in c('disruptive', 'damaging', 'deleterious')) {
  fam[,paste0('count_set_', tolower(annot))] <- na.zero(table(var$iid[var[,annot] & geneset_idx])[as.character(fam$iid)])
}
remove(geneset, geneset_idx)

m1 <- glm(as.formula(paste0('(scz-1) ~ count_all + birth + sex + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), family = binomial, fam[keepnobp,])
m2 <- glm(as.formula(paste0('(scz-1) ~ count_set_deleterious + count_all + birth + sex + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), family = binomial, fam[keepnobp,])
r2 <- NagelkerkeR2(m2)$R2 - NagelkerkeR2(m1)$R2
x <- get_tidy(m2, 'count_set_deleterious')
print(paste0('Variance (Nagelkerke coefficient of determination) explained by dURVs in potentially synaptic genes: ', signif(r2,2)*100, '% p=', signif(x$p.value,2)))

# estimate for synaptic genes
fit <- lm(as.formula(paste0('count_set_deleterious ~ scz + count_all + birth + sex + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
m1 <- x$estimate; s1 <- ( x$estimate - x$conf.low ) / 1.96
fit <- lm(as.formula(paste0('(count_deleterious - count_set_deleterious) ~ scz + count_all + birth + sex + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
m2 <- x$estimate; s2 <- ( x$estimate - x$conf.low ) / 1.96
n <- 10000000
x <- 1 / (1 + (m2 + s2 * rnorm(n)) / (m1 + s1 * rnorm(n)))
print(paste('Percentage of enrichment explained by synaptic genes:', signif(100*quantile(x, .025),4), '% -', signif(100*quantile(x, .975),4), '%'))

# estimate for synaptic complex genes
geneset <- union(genesets[['psd95']], genesets[['nmdarc']])
geneset_idx <- is.element(var$name, setid$name[is.element(setid$gene, geneset)])
fam[,'count_set_deleterious'] <- na.zero(table(var$iid[var[,annot] & geneset_idx])[as.character(fam$iid)])
fit <- lm(as.formula(paste0('count_set_deleterious ~ scz + count_all + birth + sex + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
m1 <- x$estimate; s1 <- ( x$estimate - x$conf.low ) / 1.96
fit <- lm(as.formula(paste0('(count_deleterious - count_set_deleterious) ~ scz + count_all + birth + sex + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
x <- get_tidy(fit, 'scz')
m2 <- x$estimate; s2 <- ( x$estimate - x$conf.low ) / 1.96
n <- 10000000
x <- 1 / (1 + (m2 + s2 * rnorm(n)) / (m1 + s1 * rnorm(n)))
print(paste('Percentage of enrichment explained by synaptic complex (PSD-95, NMDAR, ARC) genes:', signif(100*quantile(x, .025),2), '% -', signif(100*quantile(x, .975),2), '%'))

###########################################################################
## SEX FIGURE                                                            ##
###########################################################################

df <- read.table(paste0(datapath, 'swe.sexcheck'), header = TRUE)
sex <- setNames(read.table(paste0(datapath, 'sex.tsv'), header = FALSE), c('IID', 'REPSEX'))
df <- merge(df, sex)
df[df$F<.2 & df$YCOUNT>100, 'SNPSEX'] <- 3
figs[['sex']] <- ggplot(df[df$SNPSEX==3 | df$SNPSEX==df$REPSEX,], aes(x=F, y=YCOUNT, fill=as.factor(SNPSEX))) +
  geom_point(alpha=1/2, color = 'black', shape = 21) +
  geom_point(data = df[df$SNPSEX==1 & df$REPSEX==2,], color = 'yellow', shape = 4, show.legend = FALSE) +
  geom_point(data = df[df$SNPSEX==2 & df$REPSEX==1,], color = 'yellow', shape = 4, show.legend = FALSE) +
  scale_x_continuous('Chromosome X inbreeding coefficient (F)') +
  scale_y_continuous('Chromosome Y non-missing genotypes (YCOUNT)') +
  scale_fill_manual('', labels = c('1' = 'Male (XY)', '2' = 'Female (XX)', '3' = 'Klinefelter (XXY)'), values = c('1' = 'blue', '2' = 'red', '3' = 'green')) +
  theme_classic() + theme(panel.grid = element_blank(), legend.justification = c(0, 1), legend.position = c(0, 1.025), legend.key = element_blank(), legend.background = element_blank())
pdf(paste0(figpath, 'sex.pdf')); print(figs[['sex']]); dev.off()

###########################################################################
## GENOME FIGURE                                                         ##
###########################################################################

df <- read.table(paste0(datapath, 'swe.lite.rels.genome'), header = TRUE)
df$TYPE <- factor(NA, levels = c('Parent-child', 'Siblings', 'Duplicates'))
df$TYPE[df$Z1 >= .75] = 'Parent-child'
df$TYPE[df$PI_HAT<.75 & df$Z1 < .75] = 'Siblings'
df$TYPE[df$PI_HAT>=.75] = 'Duplicates'
tr1 <- data.frame(x = c(1,1,.5), y = c(0,1,1))
tr2 <- data.frame(x = c(.35,.35,.5), y = c(.7,1,1))
figs[['genome']] <- ggplot(df, aes(x=PI_HAT, y=Z1, fill=TYPE)) +
  geom_polygon(data = tr1, aes(x=x, y=y), fill='gray', alpha=1/4) + 
  geom_polygon(data = tr2, aes(x=x, y=y), fill='gray', alpha=1/4) +
  geom_point(alpha=1/2, color = 'black', shape = 21) +
  scale_x_continuous('Proportion IBD (PI_HAT)', limits = c(.35,1)) +
  scale_y_continuous('Proportion IBD1 (Z1)', limits = c(0,1)) +
  scale_fill_manual('', values = c('Parent-child' = 'blue', 'Siblings' = 'red', 'Duplicates' = 'green')) +
  theme_classic() + theme(panel.grid = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1.025), legend.key = element_blank(), legend.background = element_blank())
pdf(paste0(figpath, 'genome.pdf')); print(figs[['genome']]); dev.off()

###########################################################################
## FIGURE 1C AND FIGURE S2                                               ##
###########################################################################

lbl <- new.env()
lbl[['missense']] <- 'Missense'
lbl[['disruptive']] <- 'Disruptive'
lbl[['misbenign']] <- 'Missense\nnon-damaging'
lbl[['misdamaging']] <- 'Missense damaging'
lbl[['damaging']] <- 'Damaging'
lbl[['synonymous']] <- 'Synonymous'
lbl[['noncoding']] <- 'Non-coding'
lbl[['inframe']] <- 'In-frame indel'
lbl[['frameshift']] <- 'Frameshift'
lbl[['stop_gained']] <- 'Stop-gained'
lbl[['splice_donor']] <- 'Splice-donor'
lbl[['splice_acceptor']] <- 'Splice-acceptor'
lbl[['protein_protein_contact']] <- 'Protein-protein-contact'
lbl[['mutationassessor']] <- 'Mutation Assessor'
lbl[['mutationtaster']] <- 'Mutation Taster'
lbl[['fathmm']] <- 'FATHMM'
lbl[['provean']] <- 'PROVEAN'
lbl[['lrt']] <- 'LRT'
lbl[['sift']] <- 'SIFT'
lbl[['polyphen2_hvar']] <- 'PolyPhen2 HVAR'
lbl[['polyphen2_hdiv']] <- 'PolyPhen2 HDIV'

# add additional annotations
for (annot in c('inframe', 'frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'protein_protein_contact', 'mutationassessor', 'mutationtaster', 'fathmm', 'provean', 'lrt', 'sift', 'polyphen2_hvar', 'polyphen2_hdiv')) {
  var[,annot] <- is.element(var$name, read.table(paste0(datapath, 'swe.', annot), header = FALSE)$V1) & var[, 'known']
  if (annot %in% c('mutationassessor', 'mutationtaster', 'fathmm', 'provean', 'lrt', 'sift', 'polyphen2_hvar', 'polyphen2_hdiv')) { var[,annot] = var[,annot] & var[,'missense'] }
  fam[,paste0('count_', tolower(annot))] <- na.zero(table(var$iid[var[,annot]])[as.character(fam$iid)])
}

fam[,'count_misdamaging'] <- na.zero(table(var$iid[var[,'misdamaging']])[as.character(fam$iid)])

for (i in 1:3) {
  if (i==1) {lst <- c('disruptive', 'damaging', 'misbenign', 'synonymous'); opt = ' + count_all'; figname1 <- '1c'; figname2 <- '1d'}
  if (i==2) {lst <- c('frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'protein_protein_contact', 'inframe', 'misdamaging', 'mutationassessor', 'mutationtaster', 'fathmm', 'provean', 'lrt', 'sift', 'polyphen2_hvar', 'polyphen2_hdiv'); opt = ' + count_all'; figname1 <- 'preda'; figname2 <- 'predb'}
  if (i==3) {lst <- c('disruptive', 'damaging', 'misbenign', 'synonymous', 'noncoding'); opt = ''; figname1 <- 'noadjusta'; figname2 <- 'noadjustb'}

  df <- data.frame(set = character(), estimate = double(), conf.low = double(), conf.high = double(), statistic = double(), p.value = double())
  for (annot in lst) {
    fit <- lm(as.formula(paste0('count_', annot, ' ~ scz', opt, ' + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
    stats <- get_tidy(fit, 'scz')
    df <- rbind(df, data.frame(set=paste0(lbl[[annot]], '\n(', prettyNum(sum(fam[keepnobp,paste0('count_', annot)]), big.mark=","), ')'), estimate=stats['estimate'], conf.low=stats['conf.low'], conf.high=stats['conf.high'], statistic=stats['statistic'], p.value=stats['p.value'], row.names=NULL))
  }

  figs[[figname1]] <- ggplot(df, aes(x=estimate, y=set, xmin=conf.low, xmax=conf.high)) +
    geom_vline(xintercept = 0, color = 'gray50', linetype = 'dotted') +
    geom_point(color='black', fill='transparent', size=2) +
    geom_errorbarh(color = 'black', height = 0) +
    theme_classic() +
    scale_x_continuous('Excess URVs per case') +
    scale_y_discrete('') +
    geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.value, digits=2))), hjust=.5, vjust=-.5, size=4, parse=TRUE) +
    theme(panel.grid = element_blank(), plot.margin = unit(c(0.4,0,0.4,0), 'lines'))

  df <- data.frame(set = character(), estimate = double(), conf.low = double(), conf.high = double(), statistic = double(), p.value = double())
  for (annot in lst) {
    fit <- glm(as.formula(paste0('(scz-1) ~ count_', annot, opt, ' + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), family=binomial, fam[keepnobp,])
    stats <- get_tidy(fit, paste0('count_', annot), exponentiate = TRUE)
    df <- rbind(df, data.frame(set=paste0(lbl[[annot]], '\n(', prettyNum(sum(fam[keepnobp,paste0('count_', annot)]), big.mark=","), ')'), estimate=stats['estimate'], conf.low=stats['conf.low'], conf.high=stats['conf.high'], statistic=stats['statistic'], p.value=stats['p.value'], row.names=NULL))
  }

  figs[[figname2]] <- ggplot(df, aes(x=estimate, y=set, xmin=conf.low, xmax=conf.high)) +
    geom_vline(xintercept = 1, color = 'gray50', linetype = 'dotted') +
    geom_point(color='black', fill='transparent', size=2) +
    geom_errorbarh(color = 'black', height = 0) +
    theme_classic() + 
    scale_x_continuous('Odds ratios') +
    scale_y_discrete('') +
    geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.value, digits=2))), hjust=.5, vjust=-.5, size=4, parse=TRUE) +
    theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0.4,0.5,0.4,0), 'lines'))
}

figs[['predb']] = figs[['predb']] + scale_x_continuous('Odds ratios', breaks = c(1.0, 1.05, 1.1, 1.15, 1.2))

figs[['1c']] = figs[['1c']] + scale_x_continuous('Excess URVs per case', breaks = c(-.2, -.1, 0, .1, .2, .3), limits = c(-.23, .22))
figs[['1d']] = figs[['1d']] + scale_x_continuous('Odds ratios', breaks = c(1.0, 1.05, 1.1, 1.15), limits = c(.97, 1.11))

###########################################################################
## FIGURE ABOUT ENRICHMENT OF VARIANT CLASSES                            ##
###########################################################################

pdf(paste0(figpath, 'excess.pdf'), width=7, height=7)
grid.arrange(figs[['noadjusta']], figs[['noadjustb']], figs[['preda']], figs[['predb']], nrow=2, ncol=2, widths = c(3,2), heights = c(1,2))
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.98, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('b', x = unit(0.02, 'npc'), y = unit(0.64, 'npc'), gp=gpar(fontface='bold', fontsize=16))
dev.off()

###########################################################################
## EXPLORE THE CORRECT DEFINITION OF MISSENSE DAMAGING VARIANT           ##
###########################################################################

predictor <- c('mutationassessor', 'mutationtaster', 'fathmm', 'provean', 'lrt', 'sift', 'polyphen2_hvar', 'polyphen2_hdiv')
df <- data.frame(p=rep(0, 256), count=rep(0, 256), pred=rep('Predictors\nexcluding FATHMM', 256), stringsAsFactors = FALSE)
for (i in 0:255) {
  idx <- var[,'missense'] & var[, 'known']
  for (j in 0:7) {
    if (floor(i/2^j)%%2==1) {
      idx <- idx & var[,predictor[j+1]]
    }
  }
  fam[, 'count'] <- na.zero(table(var$iid[idx])[as.character(fam$iid)])
  fit <- lm(as.formula(paste0('count ~ scz + count_all + sex + birth + kit + ', paste0('pc', 1:npcs, collapse = ' + '))), fam[keepnobp,])
  stats <- get_tidy(fit, 'scz')
  df[i+1, 'p'] <- stats$p.value
  df[i+1, 'count'] <- sum(fam[keepnobp, 'count'])
}

df[, 'pred'] <- 'Predictors\nexcluding FATHMM'
for (i in 0:255) {if (floor(i/2^2)%%2==1) df[i+1, 'pred'] <- 'Predictors\nincluding FATHMM'}
df[1, 'pred'] <- 'Basic missense\npredictor'
df[249, 'pred'] <- 'Purcell et al. 2014\nlike predictor'
df[252, 'pred'] <- 'Predictor including\nall but FATHMM'
df[256, 'pred'] <- 'Predictor including\nall algorithms'
df[1+2^(0:7), 'pred'] <- 'Single algorithm\npredictors'
df[, 'pred'] <- factor(df[, 'pred'], levels = c('Basic missense\npredictor', 'Single algorithm\npredictors', 'Predictors\nexcluding FATHMM', 'Predictors\nincluding FATHMM', 'Purcell et al. 2014\nlike predictor', 'Predictor including\nall algorithms', 'Predictor including\nall but FATHMM'))

df2 <- df[df$pred=='Single algorithm\npredictors', ]
df2$label <- c('Mutation Assessor', 'Mutation Taster', 'FATHMM', 'PROVEAN', 'LRT', 'SIFT', 'PolyPhen2 HVAR', 'PolyPhen2 HDIV')

figs[['predcomp']] <- ggplot(df, aes(x=-log10(p), y=count/1000, fill=pred)) +
  geom_point(pch=21) +
  scale_x_continuous(expression(-log[10](italic(p)))) +
  scale_y_continuous('Missense URVs predicted damaging (k)') +
  theme_classic() +
  geom_text(data = df2, aes(x=-log10(p), y=count/1000, label=label), hjust=.05, vjust=-.5, size=3, angle = 37) +
  scale_fill_manual(name = '', values = c('Basic missense\npredictor' = 'light blue', 'Predictor including\nall but FATHMM' = 'red', 'Predictor including\nall algorithms' = 'darkred', 'Single algorithm\npredictors' = 'blue', 'Predictors\nincluding FATHMM' = 'dark gray', 'Predictors\nexcluding FATHMM' = 'transparent', 'Purcell et al. 2014\nlike predictor' = 'green')) +
  theme(legend.text=element_text(size = 8), legend.key = element_blank())

pdf(paste0(figpath, 'predcomp.pdf'), width=7, height=3.5)
print(figs[['predcomp']])
dev.off()

###########################################################################
## LOADING FREQUENCY DATA                                                ##
###########################################################################

dt <- fread(paste0('zcat ', datapath, 'swe.skat.frq.counts.gz'))
for (class in c('nonpsych', 'pass', 'known', 'synonymous', 'missense', 'disruptive', 'damaging', 'misbenign', 'polyphen2_hdiv', 'polyphen2_hvar', 'sift', 'lrt', 'mutationassessor', 'mutationtaster', 'fathmm', 'provean')) {
  idx <- fread(paste0(datapath, 'swe.', class), sep = '\t', header=FALSE)$V1
  dt[, class] <- is.element(dt$SNP, idx)
}
dt[, 'psych'] <- !dt[, 'nonpsych', with=FALSE]
m <- 7

###########################################################################
## FIGURE 1                                                              ##
###########################################################################

df <- data.frame(freq = integer(), exac = character(), value = double(), stringsAsFactors = FALSE)
for (i in 1:m) {
  df <- rbind(df, data.frame(freq = i, exac = 'Not in ExAC', value = sum(dt$C1==i & dt$pass & dt$known & dt$psych), row.names=NULL))
  df <- rbind(df, data.frame(freq = i, exac = 'In ExAC', value = sum(dt$C1==i & dt$pass & dt$known & dt$nonpsych), row.names=NULL))
}
figs[['1a']] <- ggplot(df, aes(x = freq, y = value/1e3, fill = reorder(exac, -as.numeric(exac)))) +
  geom_bar(stat = 'identity', color = 'black') +
  theme_classic() +
  scale_x_continuous('Minor allele count', breaks = 1:10) + scale_y_continuous('Count (k)') +
  scale_fill_manual('', values = c('Not in ExAC' = 'transparent', 'In ExAC' = 'dark gray')) +
  guides(fill = guide_legend(override.aes = list(colour = NULL), label.position = 'left')) +
  coord_cartesian(ylim = c(0, max(df[df$exac == 'In ExAC', 'value'] + df[df$exac == 'Not in ExAC', 'value'])/1e3+10)) +
  theme(panel.grid = element_blank(), legend.text = element_text(size = 8), legend.key.height = unit(1, 'line'), legend.justification = c(1, 1), legend.position = c(1.09, 1.27), legend.key = element_blank(), legend.background = element_blank())

df <- data.frame(freq = integer(), annot = character(), value = double(), stringsAsFactors = FALSE)
for (i in 1:m) {
  df <- rbind(df, data.frame(freq = i, annot = 'Synonymous', value = sum(dt$C1 == i & dt$pass & dt$synonymous), row.names=NULL))
  df <- rbind(df, data.frame(freq = i, annot = 'Missense\nnon-damaging', value = sum(dt$C1 == i & dt$pass & dt$misbenign), row.names=NULL))
  df <- rbind(df, data.frame(freq = i, annot = 'Damaging', value = sum(dt$C1 == i & dt$pass & dt$damaging), row.names=NULL))
  df <- rbind(df, data.frame(freq = i, annot = 'Disruptive', value = sum(dt$C1 == i & dt$pass & dt$disruptive), row.names=NULL))
}
figs[['1b']] <- ggplot(df, aes(x = freq, y = value/1e3, fill = reorder(annot,-as.numeric(annot)))) +
  geom_bar(stat='identity', color='black') + theme_classic() +
  scale_x_continuous('Minor allele count', breaks = 1:10) +
  scale_y_continuous('Count (k)') +
  scale_fill_manual('', values = c('Synonymous' = 'transparent', 'Missense\nnon-damaging' = 'light gray', 'Damaging' = 'dark grey', 'Disruptive' = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NULL), label.position = 'left')) +
  coord_cartesian(ylim = c(0, sum(df$value[df$freq==1])/1e3+10)) +
  theme(panel.grid = element_blank(), legend.text = element_text(size = 8), legend.key.height = unit(1, 'line'), legend.justification = c(1, 1), legend.position = c(1.09, 1.27), legend.key = element_blank(), legend.background = element_blank())

pdf(paste0(figpath, 'fig1.pdf'), width=7, height=3.5)
grid.arrange(arrangeGrob(figs[['1a']], figs[['1b']], nrow=2, ncol=1), figs[['1c']], figs[['1d']], nrow=1, ncol=3, widths = c(5,6,4))
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.96, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('b', x = unit(0.02, 'npc'), y = unit(0.46, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('c', x = unit(0.38, 'npc'), y = unit(0.96, 'npc'), gp=gpar(fontface='bold', fontsize=16))
dev.off()

###########################################################################
## PREDICTION FIGURE                                                     ##
###########################################################################

dt[,'misdamaging'] <- dt[,'missense',with=FALSE] & dt[,'damaging',with=FALSE]

df <- data.frame(annot = character(), freq = integer(), exac = character(), value = double(), lo = double(), up = double(), stringsAsFactors = FALSE)
for (class in c('polyphen2_hdiv', 'polyphen2_hvar', 'sift', 'lrt', 'mutationassessor', 'mutationtaster', 'fathmm', 'provean', 'misdamaging')) {
  for (i in 1:m) {
    n <- sum(dt$C1==i & dt$pass & dt$missense & dt$psych)
    value <- sum(dt$C1==i & dt$pass & dt$missense & dt$psych & dt[,class,with = FALSE]) / n
    lo <- value * (1 - 2 / sqrt(n))
    up <- value * (1 + 2 / sqrt(n))
    df <- rbind(df, data.frame(annot = class, freq = i, exac = 'No ExAC', value = value, lo = lo, up = up, row.names=NULL))
    n <- sum(dt$C1==i & dt$pass & dt$missense & dt$nonpsych)
    value <- sum(dt$C1==i & dt$pass & dt$missense & dt$nonpsych & dt[,class,with = FALSE]) / n
    lo <- value * (1 - 2 / sqrt(n))
    up <- value * (1 + 2 / sqrt(n))
    df <- rbind(df, data.frame(annot = class, freq = i, exac = 'In ExAC', value = value, lo = lo, up = up, row.names=NULL))
  }
}

levels(df$annot) <- c('PolyPhen2 HDIV', 'PolyPhen2 HVAR', 'SIFT', 'LRT', 'Mutation Assessor', 'Mutation Taster', 'FATHMM', 'PROVEAN', 'Missense Damaging')
figs[['pred']] <- ggplot(df, aes(x=freq, y=100*value, ymin=100*lo, ymax=100*up, color=exac)) +
  geom_line() +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ annot, scale = 'free') +
  theme_classic() +
  scale_x_continuous('Minor allele count', breaks = 1:10) +
  scale_y_continuous('Percentage (%)') +
  scale_color_manual(guide = FALSE, values = c('No ExAC'='gray','In ExAC'='black')) +
  theme(strip.background = element_blank())

pdf(paste0(figpath, 'pred.pdf')); print(figs[['pred']]); dev.off()

###########################################################################
## FAGERBERG FIGURE                                                      ##
###########################################################################

sets <- list()
sets[['All']] <- ensgenes
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'fagerberg.tsv'), sep="\t", header = TRUE)
medians <- apply(data.matrix(dfexpr[,2:28, with=FALSE]), 1, median)
for (set in names(dfexpr)[2:28]) {
  sets[[paste0(set ,' specific')]] <- dfexpr[as.logical(dfexpr[,set,with=FALSE] > medians * 5), 'Ensembl gene id',with=FALSE][['Ensembl gene id']]
}
rm(dfexpr)
write.table(sets[['brain specific']], paste0(respath, 'genes', .Platform$file.sep, 'brain.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)

df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
df2$tick <- factor(paste0(df2$set, ' (', prettyNum(df2$genecount, big.mark=','), ' genes)'), levels = unique(paste0(df2$set, ' (', prettyNum(df2$genecount, big.mark=','), ' genes)')))
thr_all <- df2$estimate[df2$set=='All' & df2$type=='or' & df2$funclass=='deleterious']

figs[['fagerberg']] <- ggplot(df2, aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5), color='gray50') +
  geom_vline(xintercept = c(1.0), color='gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(thr_all), color='gray50', linetype = 'dotdash') +
  geom_point(color='black', fill='transparent', size=2) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Odds ratios', limits = c(.85,1.25)) +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.3, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), text = element_text(size=10))

pdf(paste0(figpath, 'fagerberg.pdf')); print(figs['fagerberg']); dev.off()

###########################################################################
## CAHOY AND MO FIGURE                                                   ##
###########################################################################

sets <- list()
sets[['All']] <- ensgenes
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'cahoy.tsv'), sep="\t", header = TRUE)
slct <- c(2,3,4,5,6,7,9,11,12,13)
medians <- apply(data.matrix(dfexpr[, slct, with=FALSE]), 1, median)
for (set in names(dfexpr)[slct]) {
  sets[[paste0(set, ' specific')]] <- dfexpr[as.logical(dfexpr[,set,with=FALSE] > medians + .5), 'Gene Symbol', with=FALSE][['Gene Symbol']]
}
rm(dfexpr)
write.table(sets[['Neurons P7n specific']], paste0(respath, 'genes', .Platform$file.sep, 'neurons.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)

df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
thr_all <- df2$estimate[df2$set=='All' & df2$type=='or' & df2$funclass=='deleterious']

figs[['cahoy']] <- ggplot(df2, aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5), color = 'gray50') +
  geom_vline(xintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(thr_all), color = 'gray50', linetype = 'dotdash') +
  geom_point(color='black', fill='transparent', size=2) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Odds ratios') +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.3, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), text = element_text(size=10))

sets <- list()
sets[['All']] <- ensgenes
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'mo.tsv'), sep="\t", header = TRUE)
dfexpr <- data.frame(geneID = dfexpr$geneID, Exc = rowSums(dfexpr[, 2:3, with = FALSE]) / 2, PV = rowSums(dfexpr[, 4:5, with = FALSE]) / 2, VIP = rowSums(dfexpr[, 6:7, with = FALSE]) / 2)
mins <- apply(data.matrix(dfexpr[, 2:4]), 1, min)
sets[['Excitatory\nneuron specific']] <- dfexpr$geneID[dfexpr$Exc > mins * 5]
sets[['Inhibitory PV\nneuron specific']] <- dfexpr$geneID[dfexpr$PV > mins * 5]
sets[['Inhibitory VIP\nneuron specific']] <- dfexpr$geneID[dfexpr$VIP > mins * 5]
sets[['Expressed in\nexcitatory neurons']] <- dfexpr$geneID[dfexpr$Exc > 50]
sets[['Expressed in\ninhibitory PV neurons']] <- dfexpr$geneID[dfexpr$PV > 50]
sets[['Expressed in\ninhibitory VIP neurons']] <- dfexpr$geneID[dfexpr$VIP > 50]
rm(dfexpr)

df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'), flag = TRUE)

figs[['moleft']] <- ggplot(df2[df2$type == 'exc',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5, 4.5), color='gray50') +
  geom_vline(xintercept = 0, color='gray50', linetype = 'dotted') +
  geom_point(size = 1.5) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Excess dURVs per case') +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.1, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0,0.4,0), 'lines'))

figs[['moright']] <- ggplot(df2[df2$type == 'or',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5, 4.5), color='gray50') +
  geom_vline(xintercept = c(1.0), color='gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(df2$estimate[df2$set=='All' & df2$type=='or' & df2$funclass=='deleterious']), color='gray50', linetype = 'dotdash') +
  geom_point(size = 1.5) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Odds ratios') +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.1, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0.5,0.4,0), 'lines'))

pdf(paste0(figpath, 'cahoymo.pdf'))
grid.arrange(figs[['cahoy']], arrangeGrob(figs[['moleft']], figs[['moright']], nrow = 1, ncol = 2, widths = c(4,3)), nrow=2, ncol=1)
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.98, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('b', x = unit(0.02, 'npc'), y = unit(0.48, 'npc'), gp=gpar(fontface='bold', fontsize=16))
dev.off()

###########################################################################
## XLID FIGURE                                                           ##
###########################################################################

ens <- fread(paste0(respath, 'ensGene.bed'), sep="\t", header = FALSE)
sets <- list()
sets[['All']] <- ensgenes
sets[['X linked ID OMIM']] <- as.character(read.table(paste0(respath, 'genes', .Platform$file.sep, 'xlid.omim.txt'), header = FALSE)$V1)
sets[['X linked ID GCC']] <- as.character(read.table(paste0(respath, 'genes', .Platform$file.sep, 'xlid.gcc.txt'), header = FALSE)$V1)
sets[['X linked ID Chicago']] <- as.character(read.table(paste0(respath, 'genes', .Platform$file.sep, 'xlid.chicago.txt'), header = FALSE)$V1)
sets[['X linked ID']] <- as.character(read.table(paste0(respath, 'genes', .Platform$file.sep, 'xlid.txt'), header = FALSE)$V1)
escape <- as.character(read.table(paste0(respath, 'genes', .Platform$file.sep, 'x.escape.txt'), header = FALSE)$V1)
sets[['X linked ID and escape']] <- intersect(sets[['X linked ID']], escape)
sets[['X linked ID not escape']] <- setdiff(sets[['X linked ID']], escape)
sets[['Autosomal linked ID OMIM']] <- as.character(read.table(paste0(respath, 'genes', .Platform$file.sep, 'alid.txt'), header = FALSE)$V1)
sets[['Developmental disorder']] <- genesets[['dd']]

df2 <- get_enrichment(setid, var, fam[keepnobp & fam$sex==1,], sets, ensgenes, funclasses = c('deleterious'))
df2$group <- 'Males'
df3 <- get_enrichment(setid, var, fam[keepnobp & fam$sex==2,], sets, ensgenes, funclasses = c('deleterious'))
df3$group <- 'Females'
df4 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
df4$group <- 'All'
df5 <- rbind(df2, df3, df4)
df5$group <- as.factor(df5$group)

figs[['xlid']] <- ggplot(df5, aes(x = tick, y = estimate, ymin = conf.low, ymax = conf.high, group = group, fill = reorder(group, -as.numeric(group)))) +
  geom_vline(xintercept = c(1.5, 5.5, 7.5), color = 'gray50') +
  geom_hline(yintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_hline(yintercept = c(df5$estimate[df2$set == 'All' & df5$group == 'All']), color = 'gray50', linetype = 'dotdash') +
  geom_errorbar(width = 0, position = position_dodge(width = .8)) +
  geom_point(color = 'black', shape = 21, size = 2, position = position_dodge(width = .8)) +
  theme_classic() + coord_flip() +
  scale_x_discrete('') +
  scale_y_log10('Odds ratios', breaks = c(.4, .6, .8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0)) +
  scale_fill_manual('', values = c('All' = 'red','Males' = 'gray', 'Females' = 'black')) +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.3, size=3, parse=TRUE, position = position_dodge(width = .8)) +
  theme(panel.grid = element_blank(), legend.justification=c(1,1), legend.position=c(1.0,1.05), legend.key = element_blank(), legend.background = element_blank(), text = element_text(size = 10))

pdf(paste0(figpath, 'xlid.pdf')); print(figs['xlid']); dev.off()

###########################################################################
## FIGURE 2                                                              ##
###########################################################################

sets <- list()
sets[['All']] <- ensgenes
sets[['LoF-intolerant']] <- genesets[['pLI09']]
sets[['Missense constrained']] <- genesets[['constrained']]
sets[['microRNA-137']] <- genesets[['mir137']]
sets[['PSD-95']] <- genesets[['psd95']]
sets[['NMDAR or ARC']] <- genesets[['nmdarc']]
sets[['Schizophrenia GWAS']] <- genesets[['gwas']]
sets[['X Linked ID']] <- genesets[['xlid']]
sets[['Developmental disorder']] <- genesets[['dd']]
df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'), flag = TRUE)

figs[['2a']] <- ggplot(df2[df2$type == 'exc',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5, 3.5), color = 'gray50') +
  geom_vline(xintercept = 0, color = 'gray50', linetype = 'dotted') +
  geom_point(size = 1.5) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Excess dURVs per case') +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.1, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0,0.4,0), 'lines'))

figs[['2b']] <- ggplot(df2[df2$type == 'or',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5, 3.5), color = 'gray50') +
  geom_vline(xintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(df2$estimate[df2$set=='All' & df2$type=='or' & df2$funclass=='deleterious']), color = 'gray50', linetype = 'dotdash') +
  geom_point(size = 1.5) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Odds ratios') +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.1, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0.5,0.4,0), 'lines'))

pdf(paste0(figpath, 'fig2.pdf'), width = 7, height = 3)
grid.arrange(figs[['2a']], figs[['2b']], nrow = 1, ncol = 2, widths = c(4,3))
dev.off()

###########################################################################
## FALSE POSITIVE                                                        ##
###########################################################################

sets <- list()
sets[['All']] <- ensgenes
sets[['LoF-intolerant']] <- genesets[['pLI09']]
sets[['non LoF-intolerant']] <- setdiff(ensgenes, genesets[['pLI09']])
sets[['Missense constrained']] <- genesets[['constrained']]
sets[['non missense constrained']] <- setdiff(ensgenes, genesets[['constrained']])

df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious', 'damaging', 'disruptive'), flag = TRUE)

figs[['falsepositivea']] <- ggplot(df2[df2$type == 'exc',], aes(x=tick, y=estimate, ymin=conf.low, ymax=conf.high, group = funclass, fill = reorder(funclass, -as.numeric(funclass)))) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dotted') +
  geom_errorbar(width = 0, position = position_dodge(width = .8)) +
  geom_point(color = 'black', shape = 21, size = 2, position = position_dodge(width = .8)) +
  theme_classic() + coord_flip() +
  scale_x_discrete('') +
  scale_y_continuous('Excess dURVs per case') +
  scale_fill_manual('', values = c('deleterious' = 'red', 'damaging' = 'gray', 'disruptive' = 'black'), labels = c('deleterious' = 'disruptive\nand damaging')) + 
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.value, digits=2))), hjust=-.15, vjust=-.2, size=2.5, parse=TRUE, position = position_dodge(width = .8)) +
  theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size = 10), plot.margin = unit(c(0.4,0,0.4,0), 'lines'))

figs[['falsepositiveb']] <- ggplot(df2[df2$type=='or',], aes(x = tick, y = estimate, ymin = conf.low, ymax = conf.high, group = funclass, fill = reorder(funclass, -as.numeric(funclass)))) +
  geom_hline(yintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(df2$estimate[df2$set=='All' & df2$type=='or' & df2$funclass=='deleterious']), color = 'gray50', linetype = 'dotdash') +
  geom_errorbar(width = 0, position = position_dodge(width = .8)) +
  geom_point(color = 'black', shape = 21, size = 2, position = position_dodge(width = .8)) +
  theme_classic() + coord_flip() +
  scale_x_discrete('') +
  scale_y_log10('Odds ratios', breaks = c(.6, .8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0)) +
  scale_fill_manual('', values = c('deleterious' = 'red', 'damaging' = 'gray', 'disruptive' = 'black'), labels = c('deleterious' = 'disruptive\nand damaging')) + 
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.value, digits=2))), hjust=-.15, vjust=-.2, size=2.5, parse=TRUE, position = position_dodge(width = .8)) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.justification=c(1,0), legend.position=c(1.05,-0.05), legend.key = element_blank(), legend.background = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0.5,0.4,0), 'lines'))

pdf(paste0(figpath, 'falsepositive.pdf'), width = 7, height = 3)
grid.arrange(figs[['falsepositivea']], figs[['falsepositiveb']], nrow = 1, ncol = 2, widths = c(4,3))
dev.off()

###########################################################################
## FIGURE 3                                                              ##
###########################################################################

# first panel
sets <- list()
sets[['All']] <- ensgenes
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'fagerberg.tsv'), sep="\t", header = TRUE)
medians <- apply(data.matrix(dfexpr[,2:28, with=FALSE]), 1, median)
sets[['Kidney specific']] <- dfexpr[as.logical(dfexpr[, 'kidney', with=FALSE] > medians * 5), 'Ensembl gene id', with=FALSE][['Ensembl gene id']]
sets[['Heart specific']] <- dfexpr[as.logical(dfexpr[, 'heart', with=FALSE] > medians * 5), 'Ensembl gene id', with=FALSE][['Ensembl gene id']]
sets[['Colon specific']] <- dfexpr[as.logical(dfexpr[, 'colon', with=FALSE] > medians * 5), 'Ensembl gene id', with=FALSE][['Ensembl gene id']]
sets[['Brain specific']] <- dfexpr[as.logical(dfexpr[, 'brain', with=FALSE] > medians * 5), 'Ensembl gene id', with=FALSE][['Ensembl gene id']]
sets[['Expressed in kidney']] <- dfexpr[as.logical(dfexpr[, 'kidney', with = FALSE] > 5), 'Ensembl gene id', with = FALSE][['Ensembl gene id']]
sets[['Expressed in heart']] <- dfexpr[as.logical(dfexpr[, 'heart', with = FALSE] > 5), 'Ensembl gene id', with = FALSE][['Ensembl gene id']]
sets[['Expressed in colon']] <- dfexpr[as.logical(dfexpr[, 'colon', with = FALSE] > 5), 'Ensembl gene id', with = FALSE][['Ensembl gene id']]
sets[['Expressed in brain']] <- dfexpr[as.logical(dfexpr[, 'brain', with = FALSE] > 5), 'Ensembl gene id', with = FALSE][['Ensembl gene id']]
df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'), flag = TRUE)
df2$group <- 'Tissue type'

# second panel
sets <- list()
sets[['All']] <- ensgenes
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'cahoy.tsv'), sep="\t", header = TRUE)
slct <- c(2,3,4,5,6,7,9,11,12,13)
medians <- apply(data.matrix(dfexpr[, slct, with=FALSE]), 1, median)
sets[['Astrocyte specific']] <- dfexpr[as.logical(dfexpr[, 'Astrocytes  P7-P8', with=FALSE] > medians + .5), 'Gene Symbol', with=FALSE][['Gene Symbol']]
sets[['Oligodendrocyte specific']] <- dfexpr[as.logical(dfexpr[, 'Oligos', with=FALSE] > medians + .5), 'Gene Symbol', with=FALSE][['Gene Symbol']]
sets[['Neuron specific']] <- dfexpr[as.logical(dfexpr[, 'Neurons P7n', with=FALSE] > medians + .5), 'Gene Symbol', with=FALSE][['Gene Symbol']]
sets[['Expressed in astrocytes']] <- dfexpr[as.logical(dfexpr[, 'Astrocytes  P7-P8', with = FALSE] > 9), 'Gene Symbol', with = FALSE][['Gene Symbol']]
sets[['Expressed in oligodendrocytes']] <- dfexpr[as.logical(dfexpr[, 'Oligos', with = FALSE] > 9), 'Gene Symbol', with = FALSE][['Gene Symbol']]
sets[['Expressed in neurons']] <- dfexpr[as.logical(dfexpr[, 'Neurons P7n', with = FALSE] > 9), 'Gene Symbol', with = FALSE][['Gene Symbol']]
rm(dfexpr)
df3 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'), flag = TRUE)
df3$group <- 'Cell type'

# third panel
sets <- list()
sets[['All']] <- ensgenes
sets[['Potentially synaptic']] <- union(union(union(genesets[['rbfox2']], genesets[['fmrp']]), genesets[['celf4']]), genesets[['synaptome']])
sets[['RBFOX2 Weyn']] <- genesets[['rbfox2']]
sets[['RBFOX1/3 Weyn']] <- genesets[['rbfox13']]
sets[['FMRP Darnell']] <- genesets[['fmrp']]
sets[['CELF4 Wagnon']] <- genesets[['celf4']]
sets[['SynaptomeDB']] <- genesets[['synaptome']]
df4 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'), flag = TRUE)
df4$group <- 'Synapse-related'

df5 <- rbind(df2, df3, df4)
df5$group <- factor(df5$group, levels = c('Tissue type', 'Cell type', 'Synapse-related'))
dummy <- data.frame(group = c('Tissue type', 'Cell type', 'Synapse-related'), Z = c(5.5, 4.5, 2.5))

figs[['3a']] <- ggplot(df5[df5$type == 'exc',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5), color = 'gray50') +
  geom_hline(data = dummy, aes(yintercept = Z), color = 'gray50') +
  geom_vline(xintercept = 0, color = 'gray50', linetype = 'dotted') +
  geom_point(size = 1.5) +
  geom_errorbarh(height = 0) +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.1, size=3, parse=TRUE) +
  theme_classic() +
  scale_x_continuous('Excess dURVs per case') +
  scale_y_discrete('') +
  theme(panel.grid = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0,0.4,0), 'lines'), strip.text.y = element_blank()) +
  facet_grid(group ~ ., scale="free_y", space = "free_y")

figs[['3b']] <- ggplot(df5[df5$type == 'or',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_hline(yintercept = c(1.5), color = 'gray50') +
  geom_hline(data = dummy, aes(yintercept = Z), color = 'gray50') +
  geom_vline(xintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(df5$estimate[df5$set=='All' & df5$type=='or' & df5$funclass=='deleterious']), color = 'gray50', linetype = 'dotdash') +
  geom_point(size = 1.5) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Odds ratios', limits = c(.92,1.35)) +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.1, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0.4,0.5,0.4,0), 'lines'), strip.background = element_blank()) +
  facet_grid(group ~ ., scale="free_y", space = "free_y")

pdf(paste0(figpath, 'fig3.pdf'))
grid.arrange(figs[['3a']], figs[['3b']], nrow = 1, ncol = 2, widths = c(3,2))
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.98, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('b', x = unit(0.02, 'npc'), y = unit(0.61, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('c', x = unit(0.02, 'npc'), y = unit(0.33, 'npc'), gp=gpar(fontface='bold', fontsize=16))
dev.off()

###########################################################################
## FIGURE OR                                                             ##
###########################################################################

sets <- list()
sets[['All']] <- ensgenes
sets[['LoF-intolerant']] <- genesets[['pLI09']]
sets[['Missense constrained']] <- genesets[['constrained']]
sets[['mir-137']] <- genesets[['mir137']]
sets[['PSD-95']] <- genesets[['psd95']]
sets[['NMDAR or ARC']] <- genesets[['nmdarc']]
sets[['Schizophrenia GWAS']] <- genesets[['gwas']]
sets[['X Linked ID']] <- genesets[['xlid']]
sets[['Developmental disorder']] <- genesets[['dd']]
sets[['Brain specific']] <- genesets[['brain']]
sets[['Neuron specific']] <- genesets[['neurons']]
sets[['RBFOX2 Weyn']] <- genesets[['rbfox2']]
sets[['RBFOX1/3 Weyn']] <- genesets[['rbfox13']]
sets[['FMRP Darnell']] <- genesets[['fmrp']]
sets[['CELF4 Wagnon']] <- genesets[['celf4']]
sets[['SynaptomeDB']] <- genesets[['synaptome']]

df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious', 'damaging', 'disruptive'))

figs[['or']] <- ggplot(df2[df2$type=='or',], aes(x = tick, y = estimate, ymin = conf.low, ymax = conf.high, group = funclass, fill = reorder(funclass, -as.numeric(funclass)))) +
  geom_vline(xintercept = c(1.5, 9.5, 11.5), color = 'gray50') +
  geom_hline(yintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_hline(yintercept = c(df2$estimate[df2$type=='or' & df2$funclass=='deleterious' & df2$set == 'All']), color = 'gray50', linetype = 'dotdash') +
  geom_errorbar(width = 0, position = position_dodge(width = .8)) +
  geom_point(color = 'black', shape = 21, size = 2, position = position_dodge(width = .8)) +
  theme_classic() + coord_flip() +
  scale_x_discrete('') +
  scale_y_log10('Odds ratios', breaks = c(.6, .8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0)) +
  scale_fill_manual('', values = c('deleterious' = 'red', 'damaging' = 'gray', 'disruptive' = 'black'), labels = c('deleterious' = 'disruptive\nand damaging')) + 
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.15, vjust=-.2, size=2.5, parse=TRUE, position = position_dodge(width = .8)) +
  theme(panel.grid = element_blank(), legend.justification=c(1,1), legend.position=c(1.0,1.05), legend.key = element_blank(), legend.background = element_blank(), text = element_text(size = 10))

df2 <- get_enrichment(setid, var, fam[keepnobp & prev,], sets, ensgenes, funclasses = c('deleterious'))
df2$group <- 'published exomes'
df3 <- get_enrichment(setid, var, fam[keepnobp & !prev,], sets, ensgenes, funclasses = c('deleterious'))
df3$group <- 'new exomes'
df4 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
df4$group <- 'all exomes'
df5 <- rbind(df2, df3, df4)
df5$group <- as.factor(df5$group)

figs[['prev']] <- ggplot(df5, aes(x = tick, y = estimate, ymin = conf.low, ymax = conf.high, group = group, fill = reorder(group, -as.numeric(group)))) +
  geom_vline(xintercept = c(1.5, 9.5, 11.5), color = 'gray50') +
  geom_hline(yintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_hline(yintercept = c(df5$estimate[df5$group=='all exomes' & df5$set == 'All']), color = 'gray50', linetype = 'dotdash') +
  geom_errorbar(width = 0, position = position_dodge(width = .8)) +
  geom_point(color = 'black', shape = 21, size = 2, position = position_dodge(width = .8)) +
  theme_classic() + coord_flip() +
  scale_x_discrete('') +
  scale_y_log10('Odds ratios', breaks = c(.6, .8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0)) +
  scale_fill_manual('', values = c('all exomes' = 'red', 'new exomes' = 'gray', 'published exomes' = 'black')) + 
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.15, vjust=-.2, size=2.5, parse=TRUE, position = position_dodge(width = .8)) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0.4,0.5,0.4,0), 'lines'),
        legend.justification=c(1,1), legend.position=c(1.0,1.05), legend.key = element_blank(), legend.background = element_blank(), text = element_text(size = 10))

pdf(paste0(figpath, 'orprev.pdf'))
grid.arrange(figs[['or']], figs[['prev']], nrow=1, ncol=2, widths = c(3,2))
dev.off()

###########################################################################
## FIGURE 4                                                              ##
###########################################################################

sets <- list()
synsets <- list()
notsynsets <- list()
sets[['All']] <- ensgenes
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'fagerberg.tsv'), sep="\t", header = TRUE)
sets[['Expressed in brain']] <- dfexpr[as.logical(dfexpr[, 'brain', with = FALSE] > 5), 'Ensembl gene id', with = FALSE][['Ensembl gene id']]
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'cahoy.tsv'), sep="\t", header = TRUE)
sets[['Expressed in neurons']] <- dfexpr[as.logical(dfexpr[, 'Neurons P7n', with = FALSE] > 9), 'Gene Symbol', with = FALSE][['Gene Symbol']]
dfexpr <- fread(paste0(respath, 'genes', .Platform$file.sep, 'mo.tsv'), sep="\t", header = TRUE)
sets[['Expressed in\nexcitatory neurons']] <- dfexpr$geneID[rowSums(dfexpr[, 2:3, with = FALSE]) / 2 > 50]
sets[['Expressed in\ninhibitory neurons']] <- dfexpr$geneID[rowSums(dfexpr[, 4:7, with = FALSE]) / 4 > 50]

synapse <- union(union(union(genesets[['celf4']], genesets[['fmrp']]), genesets[['rbfox2']]), genesets[['synaptome']])
for (grp in names(sets)) {
  synsets[[grp]] <- intersect(sets[[grp]], synapse)
  notsynsets[[grp]] <- setdiff(sets[[grp]], synapse)
}

df2 <- get_enrichment(setid, var, fam[keepnobp,], synsets, ensgenes, funclasses = c('deleterious'))
df2$group <- 'Potentially synaptic'
df3 <- get_enrichment(setid, var, fam[keepnobp,], notsynsets, ensgenes, funclasses = c('deleterious'))
df3$group <- 'Not potentially synaptic'
df4 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
df4$group <- 'All'
df4$tick <- factor(paste0(df4$set, '\n(', prettyNum(df3$genecount, big.mark=","), ' + ', prettyNum(df2$genecount, big.mark=","), ' genes)'),
                   levels = unique(paste0(df4$set, '\n(', prettyNum(df3$genecount, big.mark=","), ' + ', prettyNum(df2$genecount, big.mark=","), ' genes)')))

df5 <- rbind(df2, df3, df4)
df5$group <- factor(df5$group)
df5[df5$group=='Not potentially synaptic', c('p.adjusted')] <- NA

figs[['4']] <- ggplot(df5[df5$type=='or' & df5$funclass=='deleterious' & df5$group!='All',], aes(x=set, y=estimate, ymin=conf.low, ymax=conf.high, group = group, fill = group)) +
  geom_vline(xintercept = 1.5, color = 'gray50') +
  geom_hline(yintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_hline(yintercept = c(df5$estimate[df5$set=='All' & df5$type=='or' & df5$funclass=='deleterious' & df5$group == 'All']), color = 'gray50', linetype = 'dotdash') +
  geom_errorbar(width = 0, position = position_dodge(width = -.8)) +
  geom_point(color='black', shape = 21, size = 3, position = position_dodge(width = -.8)) +
  theme_classic() + coord_flip(ylim = c(.9, 1.3)) +
  scale_x_discrete('', label=df5$tick[df5$type=='or' & df5$group == 'All']) +
  scale_y_continuous('Odds ratios', breaks = c(0.9, 1.0, 1.1, 1.2, 1.3)) +
  scale_fill_manual('', values = c('All' = 'red', 'Potentially synaptic' = 'black', 'Not potentially synaptic' = 'gray')) +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.2, size=3, parse=TRUE, position = position_dodge(width = -.8)) +
  theme(panel.grid = element_blank(), legend.position='top', legend.key = element_blank(), legend.background = element_blank(), text = element_text(size = 10))

pdf(paste0(figpath, 'fig4.pdf'), width = 4.5, height = 4); print(figs['4']); dev.off()

###########################################################################
## FIGURE 5                                                             ##
###########################################################################

sets <- list()
sets[['All']] <- ensgenes
sets[['Autism']] <- genesets[['denovo.loss.asd']]
sets[['Bipolar disorder']] <- genesets[['denovo.loss.bd']]
sets[['Schizophrenia']] <- genesets[['denovo.loss.scz']]
tmp <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
tmp$Class <- 'Deletions'
sets[['Autism']] <- genesets[['denovo.gain.asd']]
sets[['Bipolar disorder']] <- genesets[['denovo.gain.bd']]
sets[['Schizophrenia']] <- genesets[['denovo.gain.scz']]
df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, ensgenes, funclasses = c('deleterious'))
df2$Class <- 'Duplications'
df2 <- rbind(tmp, df2)
df2$Class <- factor(df2$Class, levels = c('Duplications', 'Deletions'))

figs[['5a']] <- ggplot(df2[df2$set != 'All',], aes(x = estimate, y = tick, xmin = conf.low, xmax = conf.high, fill = Class)) +
  geom_hline(yintercept = 3.5, color = 'gray50') +
  geom_vline(xintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(df2[df2$set=='All', 'estimate']), color = 'gray50', linetype = 'dotdash') +
  geom_errorbarh(height = 0) +
  geom_point(shape = 21, size = 2) +
  theme_classic() +
  scale_x_continuous('Odds ratios', limits = c(.8,1.8)) +
  scale_y_discrete('') +
  scale_fill_manual('', values = c('Deletions' = 'black', 'Duplications' = 'gray')) +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.3, size=3, parse=TRUE, color = 'black') +
  theme(panel.grid = element_blank(), legend.justification=c(1,1), legend.position=c(1.05,1.1), legend.key = element_blank(), legend.background = element_blank(), text = element_text(size=10)) +
  ggtitle('Genes hit by de-novo CNVs')

sets <- list()
sets[['All']] <- genesets[['pLI09']]
sets[['Autism']] <- intersect(genesets[['denovo.aut']], sets[['All']])
sets[['Epilepsy']] <- intersect(genesets[['denovo.epi']], sets[['All']])
sets[['Congenital\nheart disease']] <- intersect(genesets[['denovo.chd']], sets[['All']])
sets[['Intellectual\ndisability']] <- intersect(genesets[['denovo.id']], sets[['All']])
sets[['Schizophrenia']] <- intersect(genesets[['denovo.scz']], sets[['All']])

df2 <- get_enrichment(setid, var, fam[keepnobp,], sets, sets[['All']], funclasses = c('deleterious'))

figs[['5b']] <- ggplot(df2[df2$set != 'All',], aes(x=estimate, y=tick, xmin=conf.low, xmax=conf.high)) +
  geom_vline(xintercept = c(1.0), color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = c(df2[df2$set=='All', 'estimate']), color = 'gray50', linetype = 'dotdash') +
  geom_point(color='black', fill='transparent', size=2) +
  geom_errorbarh(height = 0) +
  theme_classic() +
  scale_x_continuous('Odds ratios', limits = c(.8,1.6)) +
  scale_y_discrete('') +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p.adjusted, digits=2))), hjust=-.1, vjust=-.3, size=3, parse=TRUE) +
  theme(panel.grid = element_blank(), text = element_text(size=10)) +
  ggtitle('LoF-intolerant genes hit de-novo')

pdf(paste0(figpath, 'fig5.pdf'), width = 7, height = 3)
grid.arrange(figs[['5a']], figs[['5b']], nrow = 1, ncol = 2)
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.96, 'npc'), gp = gpar(fontface = 'bold', fontsize = 16))
grid.text('b', x = unit(0.52, 'npc'), y = unit(0.96, 'npc'), gp = gpar(fontface = 'bold', fontsize = 16))
dev.off()

###########################################################################
## ULTRA-RARE VARIANTS FIGURE                                            ##
###########################################################################

figs[['urva']] <- ggplot(fam[!bad & fam$count_all<=100,], aes(x=count_snps, y=count_indels)) +
  geom_jitter(shape = 'x', alpha = 1/3) +
  geom_jitter(data=fam[!bad & fam$count_all>100,], color='red', shape='x', size=4) +
  scale_x_continuous('Ultra-rare SNPs') +
  scale_y_continuous('Ultra-rare indels') + theme_classic() +
  theme(panel.grid = element_blank())

figs[['urvb']] <- ggplot(fam[!bad & !is.na(fam$birth),], aes(x=as.factor(5*round(birth/5)), fill=as.factor(sex))) +
  geom_bar(position = 'dodge', color = 'black') +
  scale_x_discrete('Birth year') +
  scale_y_continuous('Individuals') +
  theme_classic() +
  coord_cartesian(ylim=c(47,1000)) + 
  scale_fill_manual(name = '', values = c('1' = 'light gray', '2' = 'transparent'), labels = c(paste0('Male\n(n=', sum(fam$sex==1 & !bad & !is.na(fam$birth)), ')'), paste0('Female\n(n=', sum(fam$sex==2 & !bad & !is.na(fam$birth)),')'))) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size = 8), legend.justification=c(1, 1), legend.position=c(1, 1.2), legend.key = element_blank(), legend.background = element_blank())

df <- read.table(paste0(datapath, 'swe.lite.kgp.pc.tsv'), header = TRUE)
figs[['urvc']] <- ggplot(df[df$POP!='SET',], aes(x=PC1, y=PC2, color=POP, shape=POP)) +
  geom_point(alpha = 1/2) +
  scale_colour_discrete('') +
  scale_shape_manual('', values = 0:(length(levels(df[,'POP']))-1)%%25+1) +
  theme_classic() +
  theme(panel.grid = element_blank(), legend.justification=c(1,1), legend.position=c(1.05,1.1), legend.key = element_blank(),legend.background = element_blank()) +
  geom_point(data=df[df$POP=='SET' & is.element(df$IID,fam$iid[!bad & fam$count_all<=100]),], color = 'black', shape = 'x', size=4, alpha = 1/2) +
  geom_point(data=df[df$POP=='SET' & is.element(df$IID,fam$iid[!bad & fam$count_all>100]),], color = 'red', shape = 'x', size=4, alpha = 1) +
  scale_x_continuous('1st PC') +
  scale_y_continuous('2nd PC') +
  coord_cartesian(xlim=c(-1.25,10.75))

fam$phe <- factor(fam$scz, levels=c(1,2,3), labels=c('Control', 'Schizophrenia', 'Bipolar disorder'))
fam$phe[is.na(fam$phe)] <- 'Bipolar disorder'
figs[['urvd']] <- ggplot(fam[!bad & fam$count_all<=100,], aes(x=pc3, y=pc5, color=phe, shape=phe)) +
  geom_point(alpha = 1/4) +
  theme_classic() +
  scale_shape_manual('', values = 0:2) + scale_color_manual('', values = c('black', 'green', 'blue')) +
  theme(panel.grid = element_blank(), legend.justification=c(0,1), legend.position=c(0,1.1), legend.key = element_blank(),legend.background = element_blank()) +
  geom_point(data=fam[!bad & fam$count_all>100,], color='red', shape='x', size=4, alpha = 1) +
  scale_x_continuous('3rd PC (Finnish cline)') +
  scale_y_continuous('5th PC (Northern-Southern Swedish cline)')

y <- 'count_all'
fit <- lm(as.formula(paste0(y, ' ~ pc3 + pc5')), fam[!bad,])
figs[['urve']] <- ggplot(fam[!bad,], aes_string(x='pc3', y=y, color='phe', shape='phe')) +
  geom_hline(yintercept = 100, color='gray') +
  geom_point(alpha=1/2) +
  theme_classic() +
  coord_cartesian(ylim=c(0,max(fam[!bad,y])), xlim=c(-.8,.4)) +
  scale_x_continuous('3rd PC (Finnish cline)') +
  scale_y_continuous('Ultra-rare variants (URVs)') +
  scale_shape_manual('', values = 0:2) +
  scale_color_manual('', values = c('black', 'green', 'blue')) +
  theme(panel.grid = element_blank(), legend.justification=c(0,1), legend.position=c(0,1.1), legend.key = element_blank(),legend.background = element_blank()) +
  annotate("segment", y = 0, yend = 100, x = -fit$coefficients[1]/fit$coefficients[2], xend = (100-fit$coefficients[1])/fit$coefficients[2], colour = "red")

figs[['urvf']] <- ggplot(fam[!bad,], aes_string(x='pc5', y=y, color='phe', shape='phe')) +
  geom_hline(yintercept = 100, color='gray') +
  geom_point(alpha=1/3) +
  theme_classic() +
  coord_cartesian(ylim=c(0,max(fam[!bad,y])), xlim=c(-.25,.45)) +
  scale_x_continuous('5th PC (Northern-Southern Swedish cline)') +
  scale_y_continuous('Ultra-rare variants (URVs)') +
  scale_shape_manual('', values = 0:2) +
  scale_color_manual('', values = c('black', 'green', 'blue')) +
  theme(panel.grid = element_blank(), legend.justification=c(1,1), legend.position=c(1,1.1), legend.key = element_blank(),legend.background = element_blank()) +
  annotate("segment", y = 0, yend = 100, x = -fit$coefficients[1]/fit$coefficients[3], xend = (100-fit$coefficients[1])/fit$coefficients[3], colour = "red")

pdf(paste0(figpath, 'urva.pdf'), width = 4, height = 2); print(figs[['urva']]); dev.off()
pdf(paste0(figpath, 'urvb.pdf'), width = 4, height = 2); print(figs[['urvb']]); dev.off()
pdf(paste0(figpath, 'urvc.pdf'), width = 4, height = 4); print(figs[['urvc']]); dev.off()
pdf(paste0(figpath, 'urvd.pdf'), width = 4, height = 4); print(figs[['urvd']]); dev.off()
pdf(paste0(figpath, 'urve.pdf'), width = 4, height = 4); print(figs[['urve']]); dev.off()
pdf(paste0(figpath, 'urvf.pdf'), width = 4, height = 4); print(figs[['urvf']]); dev.off()
pdf(paste0(figpath, 'urvg.pdf'), width = 5, height = 5); qqnorm(fam$count_all[!bad], ylab='Ultra-rare variants (URVs) count', main=''); dev.off()

pdf(paste0(figpath, 'urv.pdf'), width = 8, height = 10)
grid.arrange(figs[['urva']], figs[['urvb']], figs[['urvc']], figs[['urvd']], figs[['urve']], figs[['urvf']], heights=c(1/5, 2/5, 2/5), nrow=3, ncol=2)
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.99, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('b', x = unit(0.52, 'npc'), y = unit(0.99, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('c', x = unit(0.02, 'npc'), y = unit(0.79, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('d', x = unit(0.52, 'npc'), y = unit(0.79, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('e', x = unit(0.02, 'npc'), y = unit(0.39, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('f', x = unit(0.52, 'npc'), y = unit(0.39, 'npc'), gp=gpar(fontface='bold', fontsize=16))
dev.off()

###########################################################################
## ULTRA-RARE VARIANTS CORRELATIONS FIGURE                               ##
###########################################################################

type <- paste0('count_', c('disruptive', 'damaging', 'misbenign', 'synonymous'))
axislabel <- c('Disruptive URVs', 'Damaging URVs', 'Missense non-damaging URVs', 'Synonymous URVs')
for (i in 2:4) {
  for (j in 1:(i-1)) {
    r2 <- paste0('r^2 == ', signif(cor.test(fam[!bad,type[i]], fam[!bad,type[j]])$estimate^2, 2))
    figs[[paste0('urvcorr', rawToChar(as.raw(96+(i-1)*(i-2)/2+j)))]] <- ggplot(fam[!bad & fam$count_all<=100,], aes_string(x = type[j], y = type[i])) +
      geom_jitter(shape = 'x', alpha = 1/3) +
      geom_jitter(data=fam[!bad & fam$count_all>100,], color='red', shape='x', size=4) +
      scale_x_continuous(axislabel[j]) +
      scale_y_continuous(axislabel[i]) +
      annotate('text', label = r2, x = 0, y = max(fam[!bad,type[i]]), hjust=0, vjust=1, size = 4, parse=TRUE) +
      theme_classic() +
      theme(panel.grid = element_blank())
  }
}
pdf(paste0(figpath, 'urvcorr.pdf'), width = 8, height = 10)
grid.arrange(figs[['urvcorra']], figs[['urvcorrb']], figs[['urvcorrc']], figs[['urvcorrd']], figs[['urvcorre']], figs[['urvcorrf']], nrow=3, ncol=2)
grid.text('a', x = unit(0.02, 'npc'), y = unit(0.99, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('b', x = unit(0.52, 'npc'), y = unit(0.99, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('c', x = unit(0.02, 'npc'), y = unit(0.66, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('d', x = unit(0.52, 'npc'), y = unit(0.66, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('e', x = unit(0.02, 'npc'), y = unit(0.33, 'npc'), gp=gpar(fontface='bold', fontsize=16))
grid.text('f', x = unit(0.52, 'npc'), y = unit(0.33, 'npc'), gp=gpar(fontface='bold', fontsize=16))
dev.off()

###########################################################################
## WAVE ANALYSES                                                         ##
###########################################################################

df <- setNames(read.table(paste0(datapath, 'swe.wave'), header = FALSE), c('iid', 'wave'))
fam <- merge(fam, df, all.x = TRUE)
lst <- paste0('count_', c('all', 'synonymous', 'misbenign', 'damaging', 'disruptive', 'deleterious'))
df <- fam[,c('wave', 'kit', 'scz', lst)]
idx <- keepnobp & (df$wave<=7 | df$wave>=11)
df2 <- reshape(df[idx,], varying = lst, v.names = "count", timevar = "class", times = lst, direction = 'long')
rownames(df2) <- NULL
df2$class <- as.factor(df2$class)
levels(df2$class) <- c('All', 'Damaging', 'Disruptive\nand damaging', 'Disruptive', 'Missense\nnon-damaging', 'Synonymous')
df2$class <- factor(df2$class, levels = c('All', 'Synonymous', 'Missense\nnon-damaging', 'Damaging', 'Disruptive', 'Disruptive\nand damaging'))

lbls = as.character(1:12)
for (i in 1:12) {
  lbls[i] = paste0(lbls[i], '\n(', sum(keepnobp & df$wave==i & df$scz==1), ',', sum(keepnobp & df$wave==i & df$scz==2), ')')
}

figs[['wave']] <- ggplot(df2, aes(x=as.factor(wave), y=count, fill=as.factor(scz), color=as.factor(kit))) +
  stat_summary(fun.data=mean_cl_normal,geom="errorbar", size=.3, width=.2, position = position_dodge(width = 1)) +
  stat_summary(fun.y=mean, geom="point", size=2, pch=21, position=position_dodge(width=1)) +
  theme_classic() + scale_x_discrete('Wave', labels=lbls[c(1:7,11,12)]) + scale_y_continuous('URVs') +
  scale_color_manual(guide=FALSE, values = c('FALSE' = 'black', 'TRUE' = 'red')) +
  scale_fill_manual('', values = c('1' = 'white', '2' = 'dark gray'), labels = c('1' = 'Controls', '2' = 'Cases')) +
  facet_grid('class ~ .', scales = 'free') +
  theme(legend.position='top', legend.key = element_blank(), legend.background = element_blank(), text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size = 0.2), strip.background = element_blank())

pdf(paste0(figpath, 'wave.pdf'))
print(figs[['wave']])
dev.off()
