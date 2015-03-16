#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import re
import itertools
import multiprocessing

import torque

#########1#########2#########3#########4#########5#########6#########7#########


def param20140124():
    params = dict()
    params.update(p=[0.2, 0.4, 0.6, 0.8])
    params.update(P=[0.2, 0.4, 0.6, 0.8])
    params.update(c=[0.2, 0.4, 0.6, 0.8])
    params.update(C=[0.2, 0.4, 0.6, 0.8])
    params.update(s=[0.2, 0.4, 0.6, 0.8])
    params.update(S=[0.2, 0.4, 0.6, 0.8])
    params.update(f=[0.1, 0.2, 0.3, 0.4])
    return [x + [make_label(x)] for x in sequential(params)]


def param20140130():
    params = dict()
    params.update(p=[0.2, 0.4, 0.6, 0.8])
    params.update(c=[0.2, 0.4, 0.6, 0.8])
    params.update(s=[0.1, 0.2, 0.3, 0.4])
    args_list = upperlower(params) + sequential(dict(f=[0.05, 0.1, 0.15, 0.2]))
    return [x + [make_label(x)] for x in args_list]


def simple_trait1d_patch0d():
    const = ['-D1', '--row=1', '--col=1', '-m0', '-K10000']
    params = dict()
    params.update(p=[0.1, 1, 10])
    params.update(s=[0.1, 1, 10])
    params.update(c=[0.1, 1, 10])
    return [const + x + [make_label(x)] for x in product(params)]


def adaptive_dynamics():
    const = ['-D1', '--row=1', '--col=1', '-m0', '-K10000', '-p1e6', '-c1e6']
    params = dict()
    params.update(s=[1.0, 2.0, 3.0])
    params.update(C=[1.0, 2.0, 3.0])
    params.update(f=[0.02, 0.03, 0.04])
    return [const + x + [make_label(x)] for x in product(params)]


def stepping_stone():
    const = ['-D1', '--row=1', '--col=6', '-K2000', '-p1e6', '-c1e6', '-f0.03']
    params = dict()
    params.update(s=[2.0, 3.0])
    params.update(C=[1.0, 2.0])
    params.update(m=[0.001, 0.002, 0.003, 0.004])
    return [const + x + [make_label(x)] for x in product(params)]


def mutation_drift():
    const = ['-D1', '--row=1', '--col=1', '-m0', '-p1e6', '-c1e6', '-C1e6', '-f1e6']
    params = dict()
    params.update(s=[0.5, 2.0, 8.0])
    params.update(u=[5e-5, 1e-4, 2e-4])
    params.update(K=[2000, 4000, 8000])
    return [const + x + [make_label(x)] for x in product(params)]


def full_model():
    const = ['--row=6', '--col=6', '-K1000', '--symmetric']
    params = dict()
    params.update(m=[0.005, 0.01])
    params.update(s=[1.0, 2.0, 3.0])
    params.update(p=[1.0, 2.0, 3.0])
    params.update(c=[1.0, 2.0, 3.0])
    params.update(f=[0.01, 0.03, 0.05])
    return [const + x + [make_label(x)] for x in product(params)]


def upperlower(params):
    ret = []
    for (key, vals) in params.items():
        var_args = []
        for value in vals:
            var_args.append('-{}{}'.format(key, value))
            var_args.append('-{}{}'.format(key.upper(), value))
            ret.append(var_args)
    return ret


def sequential(params):
    ret = []
    for (key, vals) in params.items():
        var_args = []
        for value in vals:
            if len(key) > 1:
                var_args.append('--{}={}'.format(key, value))
            else:
                var_args.append('-{}{}'.format(key, value))
            ret.append(var_args)
    return ret


def product(params):
    ret = []
    for vals in itertools.product(*params.values()):
        var_args = []
        for (key, value) in zip(params.keys(), vals):
            if len(key) > 1:
                var_args.append('--{}={}'.format(key, value))
            else:
                var_args.append('-{}{}'.format(key, value))
        ret.append(var_args)
    return ret


def make_label(var_args):
    label = '_'.join([s.lstrip('-') for s in var_args])
    return '--label=' + re.sub('[^\w\._]+', '_', label)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-q', '--queue',
                        choices=['low', 'batch', 'high'], default='batch')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    parser.add_argument('--mode', type=int, default=0)
    parser.add_argument('outfile', nargs='?', default=sys.stdout)
    args = parser.parse_args()

    program = os.path.join(os.path.dirname(__file__), 'a.out')
    constargs = [program]
    constargs.append('--mode={}'.format(args.mode))
    constargs.append('--top_dir=' + os.getcwd())
    constargs.append('-T10000')
    constargs.append('-I100')

    args_list = full_model()
    commands = [constargs + x for x in args_list] * args.repeat

    qargs = dict()
    qargs['-q'] = args.queue
    rex = re.compile('--label=(\S+)')
    for cmd in commands:
        qargs['-N'] = rex.search(' '.join(cmd)).group(1)
        torque.qsub(cmd, args.dry_run, ppn=1, **qargs)
