#!/usr/bin/Rscript
library(stringr)
library(tidyverse)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))
.argv = commandArgs(trailingOnly=TRUE)
stopifnot(length(.argv) > 0)
#########1#########2#########3#########4#########5#########6#########7#########
if (FALSE) {
# for each locus/site => far from infinite-allele

as.bits = function(num, digit=8) {
    vapply(num, function(x) {
        paste(as.integer(intToBits(x)[digit:1]), collapse='')
    }, '')
}
#as.bits(c(15, 42))

as.bins = function(num, digit=8) {
    purrr::map_dfr(num, ~{
        t(as.integer(intToBits(.x)[digit:1]))
    })
}
#as.bins(c(15, 42))

if (FALSE) {
    .v = sample(255, 1e3, replace=TRUE)
    .u = as.bits(.v) %>% {print(head(.)); .}
    .d = as.bins(.v) %>% {print(head(.)); .}
}

.is.heterozygote = function(.tbl, trait) {
    as.bins(.tbl[[paste0(trait, '_L')]]) != as.bins(.tbl[[paste0(trait, '_R')]])
}

.heterozygosity = function(.tbl, trait) {
    .is_h = .is.heterozygote(.tbl, trait)
    c(mean=mean(.is_h), sd=sd(.is_h))
}

}
#########1#########2#########3#########4#########5#########6#########7#########
# for each trait (8 loci) => near infinite-allele

is.heterozygote = function(.tbl, trait) {
    .tbl[[paste0(trait, '_L')]] != .tbl[[paste0(trait, '_R')]]
}

heterozygosity = function(.tbl, trait) {
    .is_h = is.heterozygote(.tbl, trait)
    c(mean=mean(.is_h), sd=sd(.is_h))
}

if (FALSE) {
    .theta = 10 ^ seq(-2, 2, length=40)
    data_frame(theta=.theta, H_exp=.theta / (.theta + 1)) %>%
        ggplot(aes(theta, H_exp))+
        scale_x_log10()+
        geom_line()
}

#########1#########2#########3#########4#########5#########6#########7#########

.indir = '.'
.traits = c('toepad', 'limb', 'height_pref', 'diameter_pref', 'male', 'female', 'choosiness', 'neutral') %>%
    {setNames(., .)} %>%
    print()

expand_table = function(.tbl, by='n_TODO') {
    .tbl[rep(as.integer(row.names(.tbl)), .tbl$n),]
}
#expand_table(data.frame(a=1:3, n=3:1))

parse_hetero = function(indir) {
    .conf = wtl::read_boost_ini(file.path(indir, 'program_options.conf'))
    .conf = .conf %>% select(label, carrying_capacity, mu_locus,
                toepad_select, pref_compe, morph_compe, mating_sigma)
    .raw = read_csv(file.path(indir, 'evolution.csv.gz'))
    .final = .raw %>% filter(time==max(time)) %>% expand_table(n)
    if (FALSE) {
        .raw %>% group_by(time, row, col) %>%
            tally(n) %>%
            each(head, tail)(.)
    }
    .conf = mutate(.conf, N=sum(.final$n))
    .h = purrr::map_dfr(.traits, heterozygosity, .tbl=.final, .id='trait')
    bind_cols(.conf[rep(1,length(.traits)),], .h)
}

#.top = '~/working/anolis20141210'
.top = '~/working/anolis20150131'
.dirs = list.dirs(.top, full.names=TRUE, recursive=FALSE) %>% {setNames(., .)}
.out = purrr::map_dfr(.dirs, parse_hetero, .parallel=TRUE)

.p = .out %>%
    mutate(mu_locus = as.numeric(mu_locus),
        theta=4 * N * mu_locus * 8,
        var=sd^2,
        expected=theta / (theta + 1),  # infinite-allele model
#        expected=1 - (1 / (8 * N * mu_locus + 1) ^ 0.5),  # stepwise-mutation model
        diff=mean - expected) %>%
    filter(! trait %in% c('limb', 'diameter_pref')) %>%
    filter(mu_locus > 5e-5) %>%
    select(-.id, -sd) %>%
    gather(yparam, yval, mean, var) %>% (? summary(.)) %>%
    mutate(expected=ifelse(yparam=='mean', expected, NA)) %>%
    ggplot(aes(theta, yval), alpha=0.6)+
    geom_point(aes(colour=as.factor(toepad_select)))+
    geom_line(aes(y=expected))+
    facet_grid(yparam  ~ trait)+
    labs(x='theta = 4Nµ', y='Heterozygosity')+
    scale_colour_discrete(name="selection (s0)")+
    scale_x_log10()+
    theme_bw()
.p
ggsave('heterozygosity.pdf', .p, width=6, height=3, scale=2)
