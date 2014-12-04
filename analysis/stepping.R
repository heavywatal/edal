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
#########1#########2#########3#########4#########5#########6#########7#########
.argv = commandArgs(trailingOnly=TRUE)
stopifnot(length(.argv) > 0)
#########1#########2#########3#########4#########5#########6#########7#########

.flags = list(
    h="help",
    v="verbose",
    T="time",
    a="beta_param",
    K="carrying_capacity",
    b="birth_rate",
    p="height_pref",
    P="diameter_pref",
    s="toepad_select",
    S="limb_select",
    c="height_compe",
    C="diameter_compe",
    f="mating_sigma",
    u="mu_locus",
    U="mu_neutral",
    m="migration_rate"
)

.traits = list(
    toepad_P=expression(paste('Toepad size ', italic(x[0]))),
    limb_P=expression(paste('Limb length ', italic(x[1]))),
    height_pref_P=expression(paste('Height pref. ', italic(y[0]))),
    diameter_pref_P=expression(paste('Diameter pref. ', italic(y[1]))),
    male_P=expression(paste('Male trait ', italic(m))),
    female_P=expression(paste('Female trait ', italic(f))),
    choosiness_P=expression(paste('Choosiness ', italic(c))),
    neutral_P='Neutral'
)

.tr = c(.flags, .traits)
.axes = matrix(names(.traits), ncol=2, byrow=TRUE)
.heat_colours = c('#333333', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')

.plot = function(x, data_) {
    .tidy = data_ %>>%
        group_by_('row', 'col', 'time', x) %>>%
        tally(wt=n)# %>>% (? .)
    .p = ggplot(.tidy, aes_string(x=x, y='time'))
    .p = .p + geom_tile(aes(fill=n))
    .p = .p + scale_fill_gradientn(colours=.heat_colours)
    .p = .p + coord_cartesian(c(0, 1), range(.tidy$time))
    .p = .p + xlab(.tr[[x]])
    .p = .p + facet_grid(. ~ col)
    .p = .p + theme(legend.position='none')
    .p = .p + theme(panel.background=element_rect(fill='#000000'))
    .p = .p + theme(panel.grid=element_blank())
    .p = .p + theme(axis.text=element_blank())
    .p
}
#.plot('toepad_P')

main = function(.indir, .force=FALSE) {
    if (!file.info(.indir)$isdir) {return}
    .label = basename(.indir)
    .outfile = paste0(.label, '.png')
    if (.force || !file.exists(.outfile)) {
        message(.outfile)
        .raw = read.csv(file.path(.indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()
        .ys = c('toepad_P', 'height_pref_P', 'male_P', 'female_P', 'choosiness_P', 'neutral_P')
        .pl = llply(.ys, .plot, data_=.raw)
        .grob = do.call(arrangeGrob, c(.pl, list(ncol=1, main=.indir)))
        #print(.grob)
        ggsave(.outfile, .grob, width=4, height=4, scale=3)
    }
}

l_ply(.argv, main, .parallel=TRUE)
quit()
#########1#########2#########3#########4#########5#########6#########7#########

library(diveRsity)

.indir = 'm0.008_s2.0_C1.0_20141204_140953_23036@node46'
.indir = file.path('~/working/anolis20141203', .indir)
.raw = read.csv(file.path(.indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()

.final = .raw %>>%
    filter(time==max(time)) %>>%
    group_by(patch=sprintf('%d-%d', row, col)) %>>%
    select(matches('_L$|_R$')) %>>%
    mutate(toepad = sprintf('%03d%03d', toepad_L, toepad_R)) %>>%
    mutate(height_pref = sprintf('%03d%03d', height_pref_L, height_pref_R)) %>>%
    mutate(male = sprintf('%03d%03d', male_L, male_R)) %>>%
    mutate(female = sprintf('%03d%03d', female_L, female_R)) %>>%
    mutate(choosiness = sprintf('%03d%03d', choosiness_L, choosiness_R)) %>>%
    mutate(neutral = sprintf('%03d%03d', neutral_L, neutral_R)) %>>%
    select(-matches('_L$|_R$')) %>>% (? .)

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
.gp = as.genepop(.final) %>>% (? head(., 10))

.dp = fastDivPart(.gp, NULL, pairwise=TRUE)
.dp$standard
.dp$estimate
.dp$meanPairwise
.dp$pairwise
.dp$pw_locus


as.bits = function(num, digit=8) {
    vapply(num, function(x) {
        paste(as.integer(intToBits(x)[digit:1]), collapse='')
    }, '')
}

as.bits(c(15, 42))

as.bins = function(num, digit=8) {
    ldply(num, function(x) {
        t(as.integer(intToBits(x)[digit:1]))
    })
}

as.bins(c(15, 42))

.v = sample(255, 1e3, replace=TRUE)

.u = as.bits(.v) %>>% (? head(.))
.d = as.bins(.v) %>>% (? head(.))
