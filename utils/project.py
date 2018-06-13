#!/usr/bin/env python

def version():
    ## Get the version from version.py without importing the package
    exec(compile(open('version.py').read(), 'version.py', 'exec'))
    res = locals()['__version__']
    return res


def method_selector(method):
    jpa = False
    rrg = False
    if method == 'jpa':
        jpa = True
    elif method == 'rr':
        rrg = True
    elif method == 'jpa-rr':
        jpa = True
        rrg = True
    return jpa, rrg