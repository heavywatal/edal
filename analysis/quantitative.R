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

expand_table = function(tbl, by='n_TODO') {
    tbl[rep(as.integer(row.names(tbl)), tbl$n),]
}
#expand_table(data.frame(a=1:3, n=3:1))

is.heterozygote = function(tbl, trait) {
    as.bins(tbl[[paste0(trait, '_L')]]) != as.bins(tbl[[paste0(trait, '_R')]])
}
#.heterozygote = is.heterozygote(.final, 'neutral')
#.heterozygote %>>% colMeans()

heterozygosity = function(tbl, trait) {
    .is_h = is.heterozygote(tbl, trait)
    c(mean=mean(.is_h), sd=sd(.is_h))
}

if (FALSE) {
    .theta = 10 ^ seq(-2, 2, length=40)
    data_frame(theta=.theta, H_exp=.theta / (.theta + 1)) %>>%
        ggplot(aes(theta, H_exp))+
        scale_x_log10()+
        geom_line()
}

#########1#########2#########3#########4#########5#########6#########7#########

.indir = '.'
.traits = c('toepad', 'limb', 'height_pref', 'diameter_pref', 'male', 'female', 'choosiness', 'neutral')
names(.traits) = .traits

parse_hetero = function(indir) {
    .conf = read.conf(file.path(indir, 'program_options.conf'))
    .conf = .conf %>>% select(label, carrying_capacity, mu_locus,
                toepad_select, pref_compe, morph_compe, mating_sigma)
    .raw = read.csv(file.path(indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()
    .final = .raw %>>% filter(time==max(time)) %>>% expand_table(n)
    if (FALSE) {
        .raw %>>% group_by(time, row, col) %>>%
            tally(n) %>%
            each(head, tail)(.)
    }
    .conf = mutate(.conf, N=sum(.final$n))
    .h = ldply(.traits, heterozygosity, tbl=.final, .id='trait')
    bind_cols(.conf[rep(1,length(.traits)),], .h)
}

.dirs = list.dirs('~/working/anolis20141210', full.names=TRUE, recursive=FALSE)
names(.dirs) = .dirs
.out = ldply(.dirs, parse_hetero, .parallel=TRUE)

.out %>>%
    group_by(toepad_select, morph_compe, mating_sigma, trait) %>>%
    tally()

.out %>>%
    mutate(theta=4 * N * mu_locus,
           expected=theta / (theta + 1),
           diff=mean - expected) %>>%
    filter(! trait %in% c('limb', 'diameter_pref')) %>>%
    gather(parameter, x, toepad_select, morph_compe, mating_sigma) %>>%
    select(-.id) %>>%
    ggplot(aes(x, diff), colour=N)+
    geom_point()+
    facet_grid(trait ~ parameter, scales='free_x')
