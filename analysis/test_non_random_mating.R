#!/usr/bin/Rscript
library(tidyverse)

nrm_basic = function(male, female, choosiness, sigma=0.05) {
    exp(-(2 * choosiness - 1) ^ 2 * (female - male) ^ 2 / (2 * sigma ^ 2))
}

nrm_TPH2011 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        nrm_basic(1 - male, female, choosiness, sigma)
    )
}

nrm_TPG2013 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        2 - nrm_basic(male, female, choosiness, sigma)
    )
}

nrm_Debarre2012 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        1 -(2 * choosiness - 1) ^ 2 * (1 - exp(-(female - male) ^ 2 / (2 * sigma ^ 2))),
        1 -(2 * choosiness - 1) ^ 2 * (    exp(-(female - male) ^ 2 / (2 * sigma ^ 2)))
    )
}

.tbl = expand.grid(
    male=seq(0, 1, length=17),
    female=seq(0, 1, length=5),
    choosiness=seq(0, 1, length=17),
    sigma=c(0.02, 0.03, 0.04)) %>%
    mutate(
        TPH2011=nrm_TPH2011(male, female, choosiness, sigma),
        TPG2013=nrm_TPG2013(male, female, choosiness, sigma),
        Debarre2012=nrm_Debarre2012(male, female, choosiness, sigma)) %>%
    gather(method, mating_prob, -male, -female, -choosiness, -sigma)

p = ggplot(.tbl, aes(male, mating_prob, group=female, colour=as.factor(female)))+
    geom_line()+
    facet_grid(sigma + method ~ choosiness)+
    theme_bw()+
    theme(legend.position='top')+
    theme(axis.text.x=element_blank())+
    labs(y='Mating preference, Psi')

ggsave('test_non_random_mating.png', p, width=12, height=8)
