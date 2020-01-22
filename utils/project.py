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
    elif method == 'jpa-rr':
        jpa = True
        rrg = True
    if jpa and not rrg:
        jpa = False
        onlyjpa = True
    return jpa, rrg, onlyjpa
