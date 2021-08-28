'''
'''

import unittest

class UnittestTester:

    def __init__(self, test_class, verbosity = 1):
        self.test_class = test_class
        self.loader = unittest.TestLoader()
        if not isinstance(test_class, list):
            self.test_class = [test_class]
        else:
            self.test_class = test_class


    def _test_list(self):
        test_list = [self.loader.loadTestsFromTestCase(tclass) \
                        for tclass in self.test_class]
        return test_list


    def _test_suites(self):
        suites = unittest.TestSuite(self._test_list())
        return suites


    def _result_class(self):
        result_class = UnittestResult
        return result_class


    def execute(self):
        self.runner = unittest.TextTestRunner(resultclass = self._result_class())
        self.runner.run(self._test_suites())


'''
Example of overriding classes in the unittest.
Here, the TextTestResult is redefined to 
skip the progress dots.
(print personal log messages instead)
'''
class UnittestResult(unittest.TextTestResult):
    def addSuccess(self, test):
        super(unittest.TextTestResult, self).addSuccess(test)
        if self.showAll:
            self.stream.writeln("ok")
        elif self.dots:
            # All this for removing a dot from output :)
            #self.stream.write('.')
            self.stream.flush()


'''
Example of injecting the overriden class to the main test program.
This can be used as
    import cpydemo.unittest_tester as m_unittest
    m_unittest.main()
'''
class MTestProgram(unittest.TestProgram):
    def runTests(self):
        self.testRunner = unittest.TextTestRunner(resultclass = UnittestResult)
        super(MainTestProgram, self).runTests()

main = MTestProgram
