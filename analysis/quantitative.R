#!/usr/bin/Rscript
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(pipeR)
library(ggplot2)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))
#########1#########2#########3#########4#########5#########6#########7#########
.argv = commandArgs(trailingOnly=TRUE)
stopifnot(length(.argv) > 0)
#########1#########2#########3#########4#########5#########6#########7#########

as.bits = function(num, digit=8) {
    vapply(num, function(x) {
        paste(as.integer(intToBits(x)[digit:1]), collapse='')
    }, '')
}
#as.bits(c(15, 42))

as.bins = function(num, digit=8) {
    ldply(num, function(x) {
        t(as.integer(intToBits(x)[digit:1]))
    })
}
#as.bins(c(15, 42))

if (FALSE) {
    .v = sample(255, 1e3, replace=TRUE)
    .u = as.bits(.v) %>>% (? head(.))
    .d = as.bins(.v) %>>% (? head(.))
}

#########1#########2#########3#########4#########5#########6#########7#########

.indir = '.'
.raw = read.csv(file.path(.indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()

.final = .raw %>>% filter(time==max(time))

as.bins(.final$toepad_L)

