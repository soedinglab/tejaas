#!/usr/bin/env python

import os


def version():
    ## Get the version from version.py without importing the package
    vfn = os.path.join(os.path.dirname(__file__), '../version.py')
    exec(compile(open(vfn).read(), vfn, 'exec'))
    res = locals()['__version__']
    return res


def method_selector(method):
    onlyjpa = False
    jpa = False
    rrg = False
    if method == 'jpa':
        jpa = True
    elif method == 'rr':
        rrg = True
    ## Legacy flag, remove.
    if method == 'jpa-rr':
        rrg = True
    return jpa, rrg


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def need_new_jpanull_file(method_is_jpa, fpath):
    return method_is_jpa and fpath is not None and not is_non_zero_file(fpath)
