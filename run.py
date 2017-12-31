#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import re
import itertools
import wtl.options as wopt
now = wopt.now()


def param20140124():
    params = wopt.OrderedDict()
    params['p'] = [0.2, 0.4, 0.6, 0.8]
    params['P'] = [0.2, 0.4, 0.6, 0.8]
    params['c'] = [0.2, 0.4, 0.6, 0.8]
    params['C'] = [0.2, 0.4, 0.6, 0.8]
    params['s'] = [0.2, 0.4, 0.6, 0.8]
    params['S'] = [0.2, 0.4, 0.6, 0.8]
    params['f'] = [0.1, 0.2, 0.3, 0.4]
    for v in wopt.sequential(params):
        x = wopt.make_args(v)
        yield x + [make_label(x)]


def param20140130():
    params = wopt.OrderedDict()
    params['p'] = [0.2, 0.4, 0.6, 0.8]
    params['c'] = [0.2, 0.4, 0.6, 0.8]
    params['s'] = [0.1, 0.2, 0.3, 0.4]
    singles = wopt.sequential(dict(f=[0.05, 0.1, 0.15, 0.2]))
    args_list = itertools.chain(upperlower(params), singles)
    return [x + [make_label(x)] for x in args_list]


def simple_trait1d_patch0d():
    const = ['-D1', '--row=1', '--col=1', '-m0', '-K10000']
    params = wopt.OrderedDict()
    params['p'] = [0.1, 1, 10]
    params['s'] = [0.1, 1, 10]
    params['c'] = [0.1, 1, 10]
    for v in wopt.product(params):
        x = wopt.make_args(v)
        yield const + x + [make_label(x)]


def adaptive_dynamics():
    const = ['-D1', '--row=1', '--col=1', '-m0', '-K10000', '-p1e6', '-c1e6']
    params = wopt.OrderedDict()
    params['s'] = [1.0, 2.0, 3.0]
    params['C'] = [1.0, 2.0, 3.0]
    params['f'] = [0.02, 0.03, 0.04]
    for v in wopt.product(params):
        x = wopt.make_args(v)
        yield const + x + [make_label(x)]


def stepping_stone():
    const = ['-D1', '--row=1', '--col=6', '-K2000', '-p1e6', '-c1e6', '-f0.03']
    params = wopt.OrderedDict()
    params['s'] = [2.0, 3.0]
    params['C'] = [1.0, 2.0]
    params['m'] = [0.001, 0.002, 0.003, 0.004]
    for v in wopt.product(params):
        x = wopt.make_args(v)
        yield const + x + [make_label(x)]


def mutation_drift():
    const = ['-D1', '--row=1', '--col=1', '-m0',
             '-p1e6', '-c1e6', '-C1e6', '-f1e6']
    params = wopt.OrderedDict()
    params['s'] = [0.5, 2.0, 8.0]
    params['u'] = [5e-5, 1e-4, 2e-4]
    params['K'] = [2000, 4000, 8000]
    for v in wopt.product(params):
        x = wopt.make_args(v)
        yield const + x + [make_label(x)]


def full_model():
    const = ['--row=6', '--col=6', '-K1000', '--symmetric']
    params = wopt.OrderedDict()
    params['m'] = [0.005, 0.01]
    params['s'] = [1.0, 2.0, 3.0]
    params['p'] = [1.0, 2.0, 3.0]
    params['c'] = [1.0, 2.0, 3.0]
    params['f'] = [0.01, 0.03, 0.05]
    for v in wopt.product(params):
        x = wopt.make_args(v)
        yield const + x + [make_label(x)]


def upperlower(params):
    for (key, vals) in params.items():
        for value in vals:
            lopt = '-{}{}'.format(key, value)
            uopt = '-{}{}'.format(key.upper(), value)
            yield [lopt, uopt]


def make_label(args):
    return '--label=' + wopt.join(args)


if __name__ == '__main__':
    parser = wopt.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-B', '--batch', action='store_const', const='batch')
    parser.add_argument('-Q', '--torque', action='store_const', const='torque',
                        dest='batch')
    parser.add_argument('-q', '--queue',
                        choices=['low', 'batch', 'high'], default='batch')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    parser.add_argument('--mode', type=int, default=0)
    (args, rest) = parser.parse_known_args()

    project = os.path.dirname(__file__)
    program = os.path.join(project, 'a.out')
    constargs = [program]
    constargs.append('--mode={}'.format(args.mode))
    constargs.append('--top_dir=' + os.getcwd())
    constargs.append('-T10000')
    constargs.append('-I100')

    args_list = full_model()
    commands = [constargs + x for x in args_list] * args.repeat

    if args.batch == 'torque':
        import torque
        qargs = dict()
        qargs['-q'] = args.queue
        for cmd in commands:
            strcmd = ' '.join(cmd)
            qargs['-N'] = re.search('--label=(\S+)', strcmd).group(1)
            torque.qsub(cmd, args.dry_run, ppn=1, **qargs)
    elif args.batch == 'batch':
        wopt.map_async(commands, args.jobs, args.dry_run)
    else:
        for cmd in commands:
            print(' '.join(cmd))
        print('Use -B or -Q for batch run')
