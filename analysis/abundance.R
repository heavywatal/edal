library(stringr)
library(tidyverse)

system('make && ./a.out --test=2')
.draw = function(label) {
    .tbl = read_csv(sprintf('ignore/abundance_%s.csv', label))
    .p = ggplot(.tbl, aes(x=height, y=diameter))
    .p = .p + geom_tile(aes(fill=abundance))
    .p = .p + labs(title=label)
    .p
}
.gpl = purrr::map(c('beta', 'triangle', 'v3'), .draw)
.grobs = wtl::grid_grob(.gpl)
wtl::ggsave2('abundance_v3.png', .grobs, width=7, height=2, scale=2)

.gpl = plyr::llply(c('normal', 'exp', 'old'), .draw)
.grobs = wtl::grid_grob(.gpl)
wtl::ggsave2('abundance_old.png', .grobs, width=7, height=2, scale=2)

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

