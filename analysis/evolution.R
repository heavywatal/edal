#!/usr/bin/Rscript
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(pipeR)
library(ggplot2)
library(gridExtra)
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

.plot = function(x, y, data) {
    .hist_list = data %>>%
        group_by_(x, y) %>>%
        tally(wt=n)
    .p = ggplot(.hist_list, aes_string(x=x, y=y))
    .p = .p + geom_point(aes(colour=n), shape=15, size=8)
    .p = .p + scale_colour_gradientn(colours=.heat_colours)
    .p = .p + xlim(0, 1) + ylim(0, 1)
    .p = .p + xlab(.tr[[x]]) + ylab(.tr[[y]])
#    .p = .p + theme(legend.position='none')
    .p = .p + theme(panel.background=element_rect(fill='#000000'))
    .p = .p + theme(panel.grid=element_blank())
    .p = .p + theme(axis.text=element_blank())
    .p
}
#.plot('toepad_P', 'limb_P', .raw %>>% filter(time==1000))

.plot_time = function(x) {
    .tidy = .raw %>>%
        group_by_('time', x) %>>%
        tally(wt=n)# %>>% (? .)
    .p = ggplot(.tidy, aes_string(x=x, y='time'))
    .p = .p + geom_tile(aes(fill=n))
    .p = .p + scale_fill_gradientn(colours=.heat_colours)
#    .p = .p + xlim(0, 1)
    .p = .p + coord_cartesian(c(0, 1), range(.tidy$time))
    .p = .p + xlab(.tr[[x]])
    .p = .p + theme(legend.position='none')
    .p = .p + theme(panel.background=element_rect(fill='#000000'))
    .p = .p + theme(panel.grid=element_blank())
    .p = .p + theme(axis.text=element_blank())
    .p
}
#.plot_time('toepad_P')

#.indir = '.'
#.indir = 's1_p1_c10_20141110_033830_26387@node47'
.force = FALSE
for (.indir in .argv) {
    if (!file.info(.indir)$isdir) {next}
    .label = basename(.indir)
    .outfile = paste0('trajectory_', .label, '.pdf')
    if (.force || !file.exists(.outfile)) {
        message(.outfile)
        .raw = read.csv(file.path(.indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()
        .ys = c('toepad_P', 'height_pref_P', 'male_P', 'female_P', 'choosiness_P', 'neutral_P')
        .pl = llply(.ys, .plot_time)
        .grob = do.call(arrangeGrob, c(.pl, list(nrow=1, main=.indir)))
        #print(.grob)
        ggsave(.outfile, .grob, width=4, height=1, scale=3)
    }
    next
    .outfile = file.path(.indir, 'evolution.pdf')
    if (!file.exists(.outfile)) {
        message(.outfile)
        .raw = read.csv(file.path(.indir, 'evolution.csv.gz'), stringsAsFactors=FALSE) %>>% tbl_df()
        .done = .raw %>>%
        #    filter(time < 500) %>>%
            group_by(time) %>>%
            do(gpl={
                .pl = plyr::mlply(.axes, .plot, data=.)
                do.call(gridExtra::arrangeGrob, c(.pl,
                    list(nrow=1, main=paste0('T = ', .$time[1]))))
               })

        .mgrob = do.call(marrangeGrob, c(.done$gpl, list(nrow=1, ncol=1, top=NULL)))
        ggsave(.outfile, .mgrob, width=4, height=1, scale=6)
    }
}
