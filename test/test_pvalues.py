import unittest
  
class OutputTest(unittest.TestCase):

    # def assertMultiLineEqual(self, first, second, msg=None):
    #     """Assert that two multi-line strings are equal.

    #     If they aren't, show a nice diff.

    #     """
    #     self.assertTrue(isinstance(first, str),
    #             'First argument is not a string')
    #     self.assertTrue(isinstance(second, str),
    #             'Second argument is not a string')

    #     if first != second:
    #         message = ''.join(difflib.ndiff(first.splitlines(True),
    #                                             second.splitlines(True)))
    #         if msg:
    #             message += " : " + msg
    #         self.fail("Multi-line strings are unequal:\n" + message)

    def test_output(self):
        first = open("gold/GOLD_result_rr.txt")
        second = open("data/result_rr.txt")
        self.assertMultiLineEqual(first.read(), second.read())
        first.close()
        second.close()

if __name__ == '__main__':
    unittest.main()
