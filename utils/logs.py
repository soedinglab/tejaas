#!/usr/bin/env python

import logging
import sys
from utils import project

class MyLogger(logging.getLoggerClass()):

    def __init__(self, name, format="%(asctime)s | %(name)s | %(levelname)s | %(message)s", level=logging.DEBUG):

        formatter = logging.Formatter(format)
        # log to stdout
        handler = logging.StreamHandler(sys.stdout)
        # log to file
        # handler = logging.FileHandler('hello.log')
        handler.setFormatter(formatter)

        # Complete logging config.
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)
        self.logger.addHandler(handler)


    def info(self, msg, extra=None):
        self.logger.info(msg, extra=extra)


    def error(self, msg, extra=None):
        self.logger.error(msg, extra=extra)


    def debug(self, msg, extra=None):
        self.logger.debug(msg, extra=extra)


    def warn(self, msg, extra=None):
        self.logger.warn(msg, extra=extra)
