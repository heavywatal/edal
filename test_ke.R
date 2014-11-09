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
print(.argv)
.indir = .argv[1]
#.indir = '.'
#########1#########2#########3#########4#########5#########6#########7#########

theme_tile = theme(panel.background=element_blank())+
    theme(panel.grid=element_blank())+
    theme(axis.ticks=element_blank())+
    theme(axis.text=element_blank())

.heat_colours = c('#000000', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')
scale_fill_heat = scale_fill_gradientn(colours=.heat_colours)

.tr = c(
    height='u',
    diameter='v',
    toepad='x0',
    limb='x1',
    height_pref='y0',
    diameter_pref='y1'
)

label_tr = function(variable, value) {
    paste0(.tr[variable], ": ", value)
}

#########1#########2#########3#########4#########5#########6#########7#########
## Ke: effective carrying capacity

.conf = wtl::read.conf(file.path(.indir, 'program_options.conf'))
.conf = .conf %>>% select(ends_with('pref'), ends_with('compe'), ends_with('select'), mating_sigma)

.plot_ke_morph_vs_pref = function(.indir, .z='Ke_v3u') {
    .raw = read.csv(file.path(.indir, 'possible_phenotypes.csv.gz'))
    .tidy = .raw %>>%
        select(toepad, limb, ends_with('pref'), matches(.z)) %>>%
        mutate_each(funs(f = . / 16), toepad, limb, height_pref, diameter_pref)

    .p = .tidy %>>% ggplot(aes(x=toepad, y=limb))
    .p = .p + geom_tile(aes_string(fill=.z))
    .p = .p + facet_grid(diameter_pref ~ height_pref, as.table=FALSE, labeller=label_tr)
    .p = .p + labs(title='facet by preference')
    .p = .p + theme_tile + scale_fill_heat + theme(strip.text=element_text(size=rel(0.6)))
    .pl = list(.p)

    .p = .tidy %>>% ggplot(aes(x=height_pref, y=diameter_pref))
    .p = .p + geom_tile(aes_string(fill=.z))
    .p = .p + facet_grid(limb ~ toepad, as.table=FALSE, labeller=label_tr)
    .p = .p + labs(title='facet by morphology')
    .p = .p + theme_tile + scale_fill_heat + theme(strip.text=element_text(size=rel(0.6)))
    .pl = c(.pl, list(.p))

    .pl
}
.pl_p = .plot_ke_morph_vs_pref(.indir)
.grob = do.call(gridExtra::arrangeGrob, c(.pl_p, list(nrow=1, main='Possible phenotypes')))
#print(.grob)

#ggsave(file.path(.indir, 'possible_phenotypes.png'), .grob, width=2, height=1, scale=7)
ggsave(file.path(.indir, 'possible_phenotypes.pdf'), .grob, width=2, height=1, scale=7)

#########1#########2#########3#########4#########5#########6#########7#########
## facet geographic

.plot_geographic = function(.indir) {
    .raw = read.csv(file.path(.indir, 'possible_geographic.csv.gz'))
    .tidy = .raw %>>% tbl_df() %>>%
        filter(toepad==8, limb==8) %>>%
        select(-toepad, -limb) %>>%
        mutate(sojourn=resource * Xi_quad) %>>%
        mutate(sojourn_N=resource * Xi_quad / DI) %>>%
        mutate(Kuv = fitness * sojourn) %>>%
        mutate_each(funs(f = . / 16), height_pref, diameter_pref, height, diameter)

    # see 'sojourn_N' is strange
    .pl = llply(c('resource', 'Xi_quad', 'fitness', 'DI', 'Dxi', 'Kuv'), function(.z) {
        .p = .tidy %>>% ggplot(aes(x=height, y=diameter))
        .p = .p + geom_tile(aes_string(fill=.z))
        .p = .p + facet_grid(diameter_pref ~ height_pref, as.table=FALSE, labeller=label_tr)
        .p = .p + theme_tile + scale_fill_heat + theme(strip.text=element_text(size=rel(0.6)))
        .p
    })
    .pl
}
.pl_g = .plot_geographic(.indir)
.grob = do.call(gridExtra::arrangeGrob, c(.pl_g, list(nrow=2, main='toepad=8, limb=8, facet by preference')))
#print(.grob)

#ggsave(file.path(.indir, 'possible_geographic.png'), .grob, width=2, height=1, scale=10)
ggsave(file.path(.indir, 'possible_geographic.pdf'), .grob, width=2, height=1, scale=10)

#do.call(gridExtra::arrangeGrob, c(.pl_g, .pl_p, list(nrow=2, main='main')))
