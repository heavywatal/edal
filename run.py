#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import re

import torque

#########1#########2#########3#########4#########5#########6#########7#########


def param20140124():
    ret = []
    params = dict()
    params.update(p=[0.2, 0.4, 0.6, 0.8])
    params.update(P=[0.2, 0.4, 0.6, 0.8])
    params.update(c=[0.2, 0.4, 0.6, 0.8])
    params.update(C=[0.2, 0.4, 0.6, 0.8])
    params.update(s=[0.2, 0.4, 0.6, 0.8])
    params.update(S=[0.2, 0.4, 0.6, 0.8])
    params.update(f=[0.2, 0.4, 0.6, 0.8])
    for (key, vals) in params.items():
        for value in vals:
            var_args = ['-{}{}'.format(key, value)]
            var_args.append('--label=' + make_label(var_args))
            ret.append(var_args)
    return ret


def make_label(var_args):
    label = '_'.join([s.lstrip('-') for s in var_args])
    return re.sub('[^\w\._]+', '_', label)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('--ppn', type=int, default=4)
    parser.add_argument('-q', '--queue',
                        choices=['low', 'batch', 'high'], default='batch')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    parser.add_argument('outfile', nargs='?', default=sys.stdout)
    args = parser.parse_args()

    program = os.path.join(os.path.dirname(__file__), 'a.out')
    constargs = [program]
    constargs.append('--top_dir=' + os.getcwd())
    constargs.append('--ppn={}'.format(args.ppn))
    constargs.append('-T50000')

    args_list = param20140124()
    commands = [constargs + x for x in args_list] * args.repeat

    qargs = dict()
    qargs['-q'] = args.queue
    rex = re.compile('--label=(\S+)')
    for cmd in commands:
        qargs['-N'] = rex.search(' '.join(cmd)).group(1)
        torque.qsub(cmd, args.dry_run, ppn=args.ppn, **qargs)
