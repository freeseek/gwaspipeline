#!/usr/bin/env Rscript
###
#  manplot.R - generate Manhattan plot from plink assoc file
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

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop('About:   Generate Manhattan plot from plink assoc file (Oct 3rd 2016)\nError: Usage:   manplot.R <plink assoc> <output pdf/png> [threshold]')
}

if (length(args) > 2) {
  threshold <- as.numeric(args[3])
} else {
  threshold <- 5e-8
}

output <- args[2]
if (regexpr('.pdf$', output)>0) {
  pdf(output)
} else if (regexpr('.png$', output)>0) {
  png(output, width = 960, height = 960, type = 'cairo')
} else {
  stop('Output file must be a pdf or a png file')
}

if (args[1]=='-') {
  df <- read.table(file('stdin'), header=TRUE)
} else {
  df <- read.table(args[1], header=TRUE)
}
if(!all(names(df) %in% c('CHR', 'BP', 'P'))) { stop('Missing headers in the plink assoc file') }

# this coordinates are for GRCh37/hg19 human genome reference
chrlen <- c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067,
            159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
            115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
            59128983, 63026520, 48129895, 51305566, 155270560, 59373566)

ticks <- c(0, cumsum(chrlen))[1:23] + chrlen[1:23]/2

df$CHR[df$CHR==25] <- 23
df <- df[df$CHR %in% 1:23 & !is.na(df$P) & df$P>0 & df$P<=1, ]
df$POS <- c(0, cumsum(chrlen))[df$CHR] + df$BP

p <- ggplot(df, aes(x = POS, y = -log10(P), color = as.factor(CHR %% 2))) + geom_point() +
     scale_x_continuous(name = 'Chromosome', breaks = ticks, labels = c(1:22, 'X'), limits = c(0, sum(chrlen[1:23]))) +
     scale_y_continuous(expression(-log[10](italic(p))), breaks = 0:8, limits = c(0, .5 + -log10(min(df$P)))) +
     scale_colour_manual(guide = FALSE, values = c('gray10', 'gray60')) +
     geom_hline(yintercept = -log10(threshold), color = 'red') +
     theme_bw()

print(p)
