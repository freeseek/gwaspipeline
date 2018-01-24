#!/usr/bin/env Rscript
###
#  qqplot.R - generate QQ plot from plink assoc file
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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop('About:   Generate QQ plot from plink assoc file (Oct 3rd 2016)\nError: Usage:   qqplot.R <plink assoc> <output pdf/png>')
}

output <- args[2]
if (regexpr('.pdf$', output)>0) {
  pdf(output)
} else if (regexpr('.png$', output)>0) {
  png(output, width = 480, height = 480, type = 'cairo')
} else {
  stop('Output file must be a pdf or a png file')
}

if (args[1]=='-') {
  p <- read.table(file('stdin'), header=TRUE)$P
} else {
  p <- read.table(args[1], header=TRUE)$P
}

N <- length(p) # number of p-values

# compute the null distribution 
MAX <- -log10(min(min(p), 1/N))

# compute the confidence intervals
if (N>100) {
  ind <- c(1:100, seq(101, N, (N-101)/round(sqrt(N))))
} else {
  ind <- c(1:N)
}

c95 <- rep(0, length(ind))
c05 <- rep(0, length(ind))

# the jth order statistic from a 
# uniform(0,1) sample 
# has a beta(j,n-j+1) distribution 
# (Casella & Berger, 2002, 
# 2nd edition, pg 230, Duxbury)

for(i in 1:length(ind)) {
  c95[i] <- qbeta(0.95, ind[i], N-ind[i]+1)
  c05[i] <- qbeta(0.05, ind[i], N-ind[i]+1)
}

# plot the two confidence lines
par(mar = c(5,4,4,2) + 0.1 + c(0, 2, 0, 0))
plot(-log10(ind/N), -log10(c95), ylim=c(0,MAX), xlim=c(0,MAX), type='l', axes=FALSE, xlab='', ylab='')
par(new=T)
plot(-log10(ind/N), -log10(c05), ylim=c(0,MAX), xlim=c(0,MAX), type='l', axes=FALSE, xlab='', ylab='')

# add the diagonal
abline(0, 1, col='red')
par(new=T)

p <- sort(p)
idx1 <- rev(!duplicated(p))
p <- sort(p, decreasing=T)
idx2 <- !duplicated(p)

x <- -log10(1-(which(idx1) + which(idx2) - 2)/2/N)
y <- -log10(p[idx2])
plot(x, y, ylim = c(0, MAX), xlim = c(0, MAX), xlab = expression('Expected -log'[10]*'(p)'), ylab = expression('Observed -log'[10]*'(p)'), cex.axis = 1.5, cex.lab = 1.5)
