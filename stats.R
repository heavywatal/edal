library(plyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)

library(doMC)
doMC::registerDoMC(parallel::detectCores())

#########1#########2#########3#########4#########5#########6#########7#########

.flags = list(
    h="help",
    v="verbose",
    T="time",
    m="migration_rate",
    K="carrying_capacity",
    b="birth_rate",
    p="height_pref",
    P="diameter_pref",
    c="height_compe",
    C="diameter_compe",
    s="toepad_select",
    S="limb_select",
    f="mating_sigma",
    u="mu_locus",
    U="mu_neutral",
    a="beta_param"
)

.traits = list(
    toepad_P=expression(paste('Toepad size ', italic(x[0]))),
    limb_P=expression(paste('Limb length ', italic(x[1]))),
    height_pref_P=expression(paste('Height pref. ', italic(y[0]))),
    diameter_pref_P=expression(paste('Diameter pref. ', italic(y[1]))),
    male_P=expression(paste('Male trait ', italic(m))),
    female_P=expression(paste('Female trait ', italic(f))),
    choosiness_P=expression(paste('Choosiness ', italic(c))),
    neutral_P='Neutral'
)

.tr = c(.flags, .traits)
.axes = matrix(names(.traits), ncol=2, byrow=TRUE)

.plot = function(x, y, data) {
    .hist_matrix = reshape2::dcast(data, get(y) ~ get(x), length)
    .hist_matrix = plyr::rename(.hist_matrix, c('get(y)'=y))
    .hist_list = reshape2::melt(.hist_matrix, value.name='frequency',
                                id.vars=y, variable.name=x)
    .hist_list[, x] = as.numeric(as.character(.hist_list[, x]))
    .p = ggplot(.hist_list, aes_string(x=x, y=y))
    .p = .p + geom_tile(aes(fill=frequency))
    .p = .p + scale_fill_gradient(low='white', high='darkcyan')
    .p = .p + xlim(0, 1) + ylim(0, 1)
    .p = .p + xlab(.tr[[x]]) + ylab(.tr[[y]])
    .p = .p + theme(legend.position='none')
    .p = .p + theme(panel.background=element_rect(fill='white', colour='black'))
    .p = .p + theme(panel.grid=element_blank())
    .p = .p + theme(axis.text=element_blank())
    .p = .p + theme(axis.title=element_text(family='Linux Libertine O'))
    .p
}
.run_grob(.rundirs[1])

.run_grob = function(.rundir) {
    cat(.rundir, '\n')
    .raw = read.csv(file.path(.rundir, 'population.csv.gz'))
    .hist_matrix = reshape2::dcast(.raw, toepad_P ~ limb_P, length)
    .hist_list = reshape2::melt(.hist_matrix, value.name='frequency',
                                id.vars='toepad_P', variable.name='limb_P')
    .hist_list$limb_P = as.numeric(as.character(.hist_list$limb_P))
    .gpl = plyr::mlply(.axes, .plot, data=.raw)
    .grob = wtl::grid_grob(.gpl, 1, 4)
    .grob
}

.read_conf = function(.rundir) {
    .files = list.files(.rundir, full.names=TRUE)
    .conf_file = grep('\\.conf$', .files, value=TRUE)
    .ret = wtl::read.conf(.conf_file, id=.rundir)
    dplyr::select(.ret, -help, -test, -verbose, -seed)
}

#########1#########2#########3#########4#########5#########6#########7#########

.topdir = getwd()
.rundirs = list.files(.topdir, full.names=TRUE)

.conf = plyr::ldply(.rundirs, .read_conf)

.aggr = function(x) {
    .label = subset(.conf, id==x[1])$label
    print(.label)
    print(x)
    .gpl = plyr::llply(x, .run_grob, .parallel=TRUE)
    .grob = wtl::grid_grob(.gpl, length(x), 1)
    .key = .tr[[str_extract(.label, '^\\w')]]
    wtl::ggsave2(sprintf('tile_%s_%s.png', .key, .label), .grob)
    .label
}
dplyr::summarise(dplyr::group_by(.conf, label), g=.aggr(id))

