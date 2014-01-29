#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import re
import subprocess


def main(exe):
    exe = os.path.abspath(exe)
    help = subprocess.check_output([exe, '--help'])
    print(help)
    rex = re.compile('-(\w).+?--(\S+)')
    print('list(')
    for mobj in rex.finditer(help):
        if mobj:
            print('    {}="{}",'.format(mobj.group(1), mobj.group(2)))
    print(')')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('exe', nargs='?', default='a.out')
    args = parser.parse_args()

    main(args.exe)
