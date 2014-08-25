library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(pipeR)
library(ggplot2)
setwd('ignore')

#########1#########2#########3#########4#########5#########6#########7#########
## Ke: effective number of competitoes

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

