library(plyr)
library(reshape2)
library(ggplot2)

#########1#########2#########3#########4#########5#########6#########7#########

.filename = 'population.csv'
.raw = read.csv(.filename)

.hist_matrix = reshape2::dcast(.raw, toepad_P ~ limb_P, length)
.hist_list = (.hist_matrix)

.p = ggplot(.hist_list, aes(x=toepad_P, y=limb_P))
.p = .p + geom_tile(fill=frequency)
.p

