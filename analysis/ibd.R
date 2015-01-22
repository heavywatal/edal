#!/usr/bin/Rscript
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(pipeR)
library(ggplot2)
library(gridExtra)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))
library(diveRsity)
#########1#########2#########3#########4#########5#########6#########7#########
.argv = commandArgs(trailingOnly=TRUE)
stopifnot(length(.argv) > 0)
#########1#########2#########3#########4#########5#########6#########7#########

#.indir = 'm0.008_s2.0_C1.0_20141204_140953_23036@node46'
#.indir = file.path('~/working/anolis20141203', .indir)
#.raw = read.csv(file.path(.indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()

parse_final = function(tbl) {tbl %>>%
    filter(time==max(time)) %>>%
    group_by(patch=sprintf('%d-%d', row, col)) %>>%
    select(matches('_L$|_R$')) %>>%
    mutate(toepad = sprintf('%03d%03d', toepad_L, toepad_R),
           height_pref = sprintf('%03d%03d', height_pref_L, height_pref_R),
           male = sprintf('%03d%03d', male_L, male_R),
           female = sprintf('%03d%03d', female_L, female_R),
           choosiness = sprintf('%03d%03d', choosiness_L, choosiness_R),
           neutral = sprintf('%03d%03d', neutral_L, neutral_R)) %>>%
    select(-matches('_L$|_R$')) %>>%
    ungroup()# %>>% (? .)
}

#data(Test_data, package = "diveRsity")
#head(Test_data, 41)

as.genepop = function(tbl) {
    .names = colnames(tbl)
    .ncol = ncol(tbl)
    .header = c('BLANK', .names[2:.ncol], rep(NA, .ncol * (.ncol - 1)))
    dim(.header) = c(.ncol, .ncol)
    colnames(.header) = .names
    .boundary = data.frame(t(c('POP', rep(NA, .ncol - 1))), stringsAsFactors=FALSE)
    colnames(.boundary) = .names
    tbl = tbl %>>%
        group_by(patch) %>>%
        do(rbind(.boundary, .))
    rbind(.header, tbl)
}
#.gp = as.genepop(.final) %>>% (? head(., 10))
#.dc = diffCalc(.gp, NULL, fst=TRUE, pairwise=TRUE)
#.dc$pw_locus$Fst

diff_calc = function(gp) {
    dc = gp %>>%
        diffCalc(NULL, fst=TRUE, pairwise=TRUE)
    dc$pw_locus$Fst %>>%
    name_rows() %>>%
    gather(pair, gen_dist, -.rownames) %>>%
    rename(trait = .rownames) %>>%
    mutate(trait = factor(trait, levels=rownames(dc$pw_locus$Fst))) %>>%
    extract(pair, c('row.i', 'col.i', 'row.j', 'col.j'), '(\\d+)-(\\d+).+?(\\d+)-(\\d+)', convert=TRUE) %>>%
    mutate(geo_dist = row.j - row.i + col.j - col.i) %>>%
    select(-matches('^row|^col')) %>>%
    arrange(trait, geo_dist)
}

pipeline = function(infile) {infile %>>%
    read.csv(stringsAsFactors=FALSE) %>>%
    parse_final() %>>%
    as.genepop() %>>%
    diff_calc()
}

.indirs = list.dirs('~/working/anolis20141203', recursive=FALSE)
.infiles = file.path(.indirs, 'program_options.conf')
names(.infiles) = file.path(.indirs, 'evolution.csv.gz')
.conf = ldply(.infiles, wtl::read.conf, .id='path') %>>%
    select(path, toepad_select, morph_compe, migration_rate) %>>%
    mutate(path=as.character(path))

ibd = .conf %>>%
    mdply(function(path, ..., .parallel=TRUE) {
        message(path)
        pipeline(path)
    })
save(ibd, file='~/working/anolis20141203/ibd.rda')

.gp = ibd %>>%
    filter(!migration_rate %in% c(0.001, 0.003, 0.004, 0.008)) %>>%
    mutate(migration_rate=as.factor(migration_rate)) %>>%
    ggplot(aes(geo_dist, gen_dist, group=migration_rate, colour=migration_rate))+
    geom_point(alpha=0.4)+
    geom_smooth(method=glm, se=FALSE, alpha=0.8)+
    facet_grid(toepad_select + morph_compe ~ trait, labeller=label_both)+
    theme_bw()
ggsave('isolation_by_distance.png', .gp, width=6, height=4, scale=1.5)

quit()
################################################################################

library(adegenet)

.genind = df2genind(.final %>>% select(-patch), pop=.final$patch)
pairwise.fst(.genind)

.traits = colnames(.final)[-1]
names(.traits) = .traits
llply(.traits, function(trait){
    .final %>>%
        select_(trait) %>>%
        df2genind(pop=.final$patch) %>>%
        pairwise.fst()
})
