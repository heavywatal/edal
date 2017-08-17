#!/usr/bin/Rscript
library(stringr)
library(tidyverse)
library(gridExtra)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))
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
    .tidy = data_ %>%
        group_by_('row', 'col', 'time', x) %>%
        tally(wt=n)# %>% print()
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
        .raw = read_csv(file.path(.indir, 'evolution.csv.gz'))
        .ys = c('toepad_P', 'height_pref_P', 'male_P', 'female_P', 'choosiness_P', 'neutral_P')
        .pl = purrr::map(.ys, .plot, data_=.raw)
        .grob = do.call(arrangeGrob, c(.pl, list(ncol=1, main=.indir)))
        #print(.grob)
        ggsave(.outfile, .grob, width=3, height=4, scale=3)
    }
}

purrr::walk(.argv, main)
