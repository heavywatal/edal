#!/usr/bin/Rscript
library(pipeR)
library(dplyr)
library(tidyr)
library(ggplot2)

nrm_basic = function(male, female, choosiness, sigma=0.05) {
    exp(-(2 * choosiness - 1) ^ 2 * (female - male) ^ 2 / (2 * sigma ^ 2))
}

nrm_TPH2011 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        nrm_basic(1 - male, female, choosiness, sigma)
    )
}

nrm_xavier2013 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        2 - nrm_basic(male, female, choosiness, sigma)
    )
}

nrm_xavier2013a = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        1 - nrm_basic(male, female, choosiness, sigma)
    )
}

nrm_debarre2012 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        1 -(2 * choosiness - 1) ^ 2 * (1 - exp(-(female - male) ^ 2 / (2 * sigma ^ 2))),
        1 -(2 * choosiness - 1) ^ 2 * (    exp(-(female - male) ^ 2 / (2 * sigma ^ 2)))
    )
}

tbl = expand.grid(
    male=seq(0, 1, length=30),
    female=seq(0, 1, length=5),
    choosiness=seq(0.1, 0.9, length=5)) %>>%
    mutate(
        TPH2011=nrm_TPH2011(male, female, choosiness),
        xavier2013=nrm_xavier2013(male, female, choosiness),
        xavier2013a=nrm_xavier2013a(male, female, choosiness),
        debarre2012=nrm_debarre2012(male, female, choosiness)) %>>%
    gather(method, mating_prob, -male, -female, -choosiness)

p = ggplot(tbl, aes(male, mating_prob, group=female, colour=as.factor(female)))+
    geom_line()+
    facet_grid(method ~ choosiness, labeller=label_both)+
    theme_bw()+
    theme(legend.position='top')+
    labs(y='Mating probability, Psi')

ggsave('test_non_random_mating.png', p)
