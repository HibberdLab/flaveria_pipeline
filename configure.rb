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
spades = which "spades.py"
trimmomatic = which "trimmomatic"
khmer = which "normalize-by-median.py"

memory = 90
threads = 24

