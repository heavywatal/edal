library(plyr)
library(reshape2)
library(ggplot2)

#########1#########2#########3#########4#########5#########6#########7#########

.filename = 'population.csv'
.raw = read.csv(.filename)
head(.raw)

.hist_matrix = reshape2::dcast(.raw, toepad_P ~ limb_P, length)
.hist_list = reshape2::melt(.hist_matrix, value.name='frequency',
                            id.vars='toepad_P', variable.name='limb_P')
.hist_list$limb_P = as.numeric(as.character(.hist_list$limb_P))

.p = ggplot(.hist_list, aes(x=toepad_P, y=limb_P))
.p = .p + geom_tile(aes(fill=frequency))
.p = .p + xlim(0, 1) + ylim(0, 1)
.p

.axes = matrix(c('toepad_P', 'limb_P',
                 'height_pref_P', 'diameter_pref_P',
                 'male_P', 'female_P',
                 'choosiness_P', 'neutral_P'),
               ncol=2, byrow=TRUE)

.plot = function(x, y) {
    .hist_matrix = reshape2::dcast(.raw, get(y) ~ get(x), length)
    .hist_matrix = plyr::rename(.hist_matrix, c('get(y)'=y))
    .hist_list = reshape2::melt(.hist_matrix, value.name='frequency',
                                id.vars=y, variable.name=x)
    .hist_list[, x] = as.numeric(as.character(.hist_list[, x]))
    .p = ggplot(.hist_list, aes_string(x=x, y=y))
    .p = .p + geom_tile(aes(fill=frequency))
    .p = .p + scale_fill_gradient(low='black')
    .p = .p + xlim(0, 1) + ylim(0, 1)
    .p = .p + theme(legend.position='none')
    .p = .p + theme(panel.background=element_rect(fill='black'))
    .p = .p + theme(panel.grid=element_blank())
    .p
}

.gpl = plyr::mlply(.axes, .plot)
.grob = wtl::grid_grob(.gpl, 2, 2)
.grob
ggsave2('phenotype_spectra.png', plot=.grob)
