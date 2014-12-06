#!/usr/bin/Rscript
library(pipeR)
library(dplyr)
library(tidyr)
library(ggplot2)

nrm_basic = function(male, female, choosiness, sigma=0.05) {
    exp(-(2 * choosiness - 1) ^ 2 * (female - male) ^ 2 / (2 * sigma ^ 2))
}

nrm_anolis = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        nrm_basic(1 - male, female, choosiness, sigma)
    )
}

nrm_xavier1 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        1 - nrm_basic(male, female, choosiness, sigma)
    )
}

nrm_xavier2 = function(male, female, choosiness, sigma=0.05) {
    ifelse(choosiness >= 0.5,
        nrm_basic(male, female, choosiness, sigma),
        2 - nrm_basic(male, female, choosiness, sigma)
    )
}

tbl = expand.grid(
    male=seq(0, 1, length=30),
    female=seq(0, 1, length=5),
    choosiness=seq(0, 1, length=5)) %>>%
    mutate(
        anolis=nrm_anolis(male, female, choosiness),
        xavier1=nrm_xavier1(male, female, choosiness),
        xavier2=nrm_xavier2(male, female, choosiness)) %>>%
    gather(method, mating_prob, anolis, xavier1, xavier2)

p = ggplot(tbl, aes(male, mating_prob, group=female, colour=as.factor(female)))+
    geom_line()+
    facet_grid(method ~ choosiness, labeller=label_both)+
    theme_bw()+
    theme(legend.position='top')+
    labs(y='Mating probability, Psi')

ggsave('test_non_random_mating.png', p)
