#!/usr/bin/env ruby

#
# configuration for flaveria pipeline
#
# checks the install location of:
#   fastqc
#   spades
#   trimmomatic
#   khmer
#

require 'which'
require 'trollop'

include Which

fastqc = which "fastqc"
fastqc = which "spades.py"
fastqc = which "trimmomatic"
fastqc = which "normalize-by-median.py"

