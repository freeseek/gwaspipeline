#!/bin/bash
###
#  plink2lite.sh - generate a lite version of a plink dataset
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

# generate a lite version of a plink dataset
# requires plink prefix of the input dataset
# requires plink prefix of the output dataset

set -e -o pipefail

if [ $# -lt 3 ]; then
  echo "About:   generate a lite version of a plink dataset (Oct 3rd 2016)"
  echo "Usage:   plink2lite.sh <input prefix> <output prefix> <respath> [--opt opt]"
  exit 1
fi

in="$1"
out="$2"
res="$3"

if [ $# -ge 5 ]; then
  opt="$4 $5"
fi

# path where KGP dataset is kept
kgp="$res/kgp"

# makes sure the ouput directory exists
dir=$(dirname $out)
mkdir -p $dir

# if list of markers with at least 1% frequency in KGP set does not exist, create it
if [ ! -f $kgp.maf.001 ]; then
  if [ ! -f $kgp.frq.counts.gz ]; then
    plink --bfile $kgp --keep-allele-order --freq counts --out $kgp && gzip $kgp.frq.counts
  fi
  gzip -cd $kgp.frq.counts.gz | tail -n+2 | awk '$5/($5+$6)>.01 && $6/($5+$6)>.01 {print $2}' | sort > $kgp.maf.001
fi

# extract KGP markers from dataset
plink --bfile $in $opt \
  --keep-allele-order \
  --extract $kgp.maf.001 \
  --make-bed \
  --out $out

# cleanup
if [ -f $out.nosex ]; then rm $out.nosex; fi
if [ -f $out.hh ]; then rm $out.hh; fi
if [ -f $out.log ]; then rm $out.log; fi
