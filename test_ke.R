library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(pipeR)
library(ggplot2)
setwd('ignore')
#########1#########2#########3#########4#########5#########6#########7#########

is.odd = function(x) x %% 2 != 0
is.odd(-3:3)

.theme_tile =
    theme(panel.grid=element_blank())+
    theme(panel.background=element_blank())+
    theme(axis.ticks=element_blank())+
    theme(axis.text=element_blank())

heat_colours = c('#000000', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')
scale_fill_heat = scale_fill_gradientn(colours=heat_colours)

#########1#########2#########3#########4#########5#########6#########7#########
## Ke: effective number of competitoes

.working_dir = '~/working/anolis20140826'
setwd(.working_dir)
list.files(.working_dir)

.indir = list.files(.working_dir)[1]
.raw = read.csv(file.path(.indir, 'possible_ke.csv.gz'))
.raw = read.csv('possible_ke.csv.gz') %>>% (? nrow(.)) %>>% (? sample_n(., 20))

.tidy = .raw %>>%
    gather(key, Ke, ends_with('Ke')) %>>% (? sample_n(., 20)) %>>%
    mutate(normalized=ifelse(key=='Ke', 1, 0))

.p = .tidy %>>%
    filter(is.odd(toepad), is.odd(limb)) %>>%
    mutate_each(funs(f = . / 16), toepad, limb, height_pref, diameter_pref) %>>%
    ggplot(aes(x=height_pref, y=diameter_pref))
.p = .p + geom_tile(aes(fill=Ke))
.p = .p + facet_grid(limb ~ normalized + toepad, as.table=FALSE, labeller=label_both)
.p = .p + .theme_tile + scale_fill_heat
.p = .p + theme(strip.text=element_text(size=6))
.p

.p = .tidy %>>%
    filter(is.odd(height_pref), is.odd(diameter_pref)) %>>%
    mutate_each(funs(f = . / 16), toepad, limb, height_pref, diameter_pref) %>>%
    ggplot(aes(x=toepad, y=limb))
.p = .p + geom_tile(aes(fill=Ke))
.p = .p + facet_grid(diameter_pref ~ normalized + height_pref, as.table=FALSE, labeller=label_both)
.p = .p + .theme_tile + scale_fill_heat
.p = .p + theme(strip.text=element_text(size=6))
.p

#########1#########2#########3#########4#########5#########6#########7#########
## sojourn time

.raw = read.csv('sojourn_time.csv.gz') %>>% (? sample_n(., 30))

.p = .raw %>>%
    filter(!is.na(height_pref %>>% match(seq(1, 15, 2) / 16)),
        !is.na(diameter_pref %>>% match(seq(1, 15, 2) / 16))) %>>%
    ggplot(aes(x=height, y=diameter))
.p = .p + geom_tile(aes(fill=time))
.p = .p + facet_grid(diameter_pref ~ normalized + height_pref, as.table=FALSE)
.p = .p + scale_fill_gradientn(colours=c('#000000', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000'))
.p = .p + .theme_tile
.p

.p = .raw %>>%
    filter(!is.na(height %>>% match(seq(1, 15, 2) / 16)),
        !is.na(diameter %>>% match(seq(1, 15, 2) / 16))) %>>%
    ggplot(aes(x=height_pref, y=diameter_pref))
.p = .p + geom_tile(aes(fill=time))
.p = .p + facet_grid(diameter ~ normalized + height, as.table=FALSE)
.p = .p + scale_fill_gradientn(colours=c('#000000', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000'))
.p = .p + .theme_tile
.p


#########1#########2#########3#########4#########5#########6#########7#########
#########1#########2#########3#########4#########5#########6#########7#########
## deprecated

.odd = .raw %>>% filter(is.odd(toepad), is.odd(limb), is.odd(height_pref), is.odd(diameter_pref))
.odd %>>% (? nrow(.)) %>>% sample_n(20)

.p = ggplot(.odd, aes(x=height_pref, y=diameter_pref))
.p = .p + geom_tile(aes(fill=Ke))
.p = .p + facet_grid(limb ~ toepad, as.table=FALSE)
.p = .p + .theme_tile + scale_fill_heat
.p

.p = ggplot(.odd, aes(x=toepad, y=limb))
.p = .p + geom_tile(aes(fill=Ke))
.p = .p + facet_grid(diameter_pref ~ height_pref, as.table=FALSE)
.p = .p + .theme_tile + scale_fill_heat
.p

.p = .raw %>>%
    filter(toepad==7, limb==7) %>>%
    ggplot(aes(height_pref, Ke, group=diameter_pref, colour=diameter_pref))
.p = .p + geom_point()
.p = .p + geom_line()
.p

.p = .raw %>>%
    filter(height_pref==7, diameter_pref==1) %>>%
    ggplot(aes(toepad, Ke, group=limb, colour=limb))
.p = .p + geom_point()
.p = .p + geom_line()
.p

#########1#########2#########3#########4#########5#########6#########7#########

.filename = 'ke.csv'
.tbl = read.csv(.filename)

.p = ggplot(dplyr::filter(.tbl, height_pref==median(height_pref), diameter_pref==median(height_pref)),
            aes(x=toepad, y=limb))
.p = .p + geom_tile(aes(fill=Ke))
.p

.p = ggplot(dplyr::filter(.tbl, toepad==median(toepad), limb==median(limb)),
            aes(x=height_pref, y=diameter_pref))
.p = .p + geom_tile(aes(fill=Ke))
.p

.p = ggplot(dplyr::filter(.tbl, toepad==median(toepad), limb==1),
            aes(x=height_pref, y=diameter_pref))
.p = .p + geom_tile(aes(fill=Ke))
.p

.p = ggplot(.tbl, aes(x=height_pref, y=diameter_pref))
.p = .p + geom_tile(aes(fill=Ke))
.p = .p + facet_grid(limb ~ toepad, as.table=FALSE)
.p = .p + theme(panel.grid=element_blank(), panel.background=element_blank())
.p = .p + theme(axis.ticks=element_blank())
.p = .p + theme(axis.text=element_blank())
.p
ggsave('ke_group_by_morph.pdf', .p, width=16, height=16)

.p = ggplot(.tbl, aes(x=toepad, y=limb))
.p = .p + geom_tile(aes(fill=Ke))
.p = .p + facet_grid(diameter_pref ~ height_pref, as.table=FALSE)
.p = .p + theme(panel.grid=element_blank(), panel.background=element_blank())
.p = .p + theme(axis.ticks=element_blank())
.p = .p + theme(axis.text=element_blank())
.p
ggsave('ke_group_by_preference.pdf', .p, width=16, height=16)

plyr::mlply(unique(.tbl[, c('height_pref', 'diameter_pref')]),
    function(height_pref, diameter_pref){
    .sub = subset(.tbl, height_pref==height_pref && diameter_pref==diameter_pref)
    dplyr::filter(.tbl)
})

#########1#########2#########3#########4#########5#########6#########7#########

.hexize = function(n, maximum=16) {
    as.hexmode(as.integer(n / maximum * 255))
}

.colourize = function(a, b) {
    sprintf('%02x00%02x', .hexize(a), .hexize(b))
}

.peaks = .tbl %.%
    dplyr::group_by(toepad, limb) %.%
    dplyr::filter(Ke==max(Ke)) %.%
    dplyr::mutate(colour=.colourize(height_pref, diameter_pref))

.p = ggplot(.peaks, aes(x=toepad, y=limb))
.p = .p + geom_tile(aes(fill=colour))
.p = .p + theme(panel.grid=element_blank(), panel.background=element_blank())
.p = .p + theme(axis.ticks=element_blank())
.p = .p + theme(axis.text=element_blank())
.p

#########1#########2#########3#########4#########5#########6#########7#########
## too slow

.plot = function(.tbl) {
    .p = ggplot(.tbl, aes(x=height_pref, y=diameter_pref))
    .p = .p + geom_tile(aes(fill=Ke))
    .p = .p + ggtitle(paste0(.tbl$toepad[1], ', ', .tbl$limb[1]))
    .p = .p + theme(panel.grid=element_blank(), panel.background=element_blank())
    .p = .p + theme(axis.ticks=element_blank())
    .p = .p + theme(axis.text=element_blank())
    .p = .p + theme(axis.title=element_blank())
    .p = .p + theme(legend.position='none')
    .p
}

.gpl = plyr::mlply(unique(.tbl[, c('height_pref', 'diameter_pref')]),
    function(height_pref, diameter_pref){
    .plot(subset(.tbl, height_pref==height_pref && diameter_pref==diameter_pref))
})

.gtree = wtl::grid_grob(.gpl, 16, 16)

#########1#########2#########3#########4#########5#########6#########7#########
## group_by() does not consider order

.gpl = .tbl %.%
    dplyr::group_by(toepad, limb) %.%
    dplyr::do(.plot)

