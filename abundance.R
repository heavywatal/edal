library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(pipeR)
library(ggplot2)
setwd('ignore')

#########1#########2#########3#########4#########5#########6#########7#########

tile_abundance = function(filename) {
    system('make && make test')
    tbl = read.csv(filename)
    .p = ggplot(tbl, aes(x=height, y=diameter))
    .p = .p + geom_tile(aes(fill=abundance))
    .p
}
tile_abundance('abundance_v3.csv')

#########1#########2#########3#########4#########5#########6#########7#########

system('make && make test')
.draw = function(label) {
    tbl = read.csv(sprintf('abundance_%s.csv', label))
    .p = ggplot(tbl, aes(x=height, y=diameter))
    .p = .p + geom_tile(aes(fill=abundance))
    .p = .p + labs(title=label)
    .p
}
files = c('beta', 'triangle', 'v3')
.gpl = plyr::llply(files, .draw)
.grobs = wtl::grid_grob(.gpl)
wtl::ggsave2('abundance_v3.png', .grobs, width=7, height=2, scale=2)

#########1#########2#########3#########4#########5#########6#########7#########

theta = function(u, c0=2, c1=1) {
    (c0 - c1 * u)
}
u = seq(0, 1, by=0.01)
plot(u, 1 / theta(u, 5, 4))

plot(u, exp(-theta(u, 10, 1)), ylim=c(0, 1))
plot(u, exp(-theta(u, 10, 9)), ylim=c(0, 1))

plot(u, exp(-1 / theta(u, 2, 1)) / theta(u, 2, 1), ylim=c(0, 1))
plot(u, exp(-1 / theta(u, 10, 9)) / theta(u, 10, 9), ylim=c(0, 1))

