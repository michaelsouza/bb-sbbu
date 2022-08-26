from tkinter import SE
import unittest

from bb import *

class TestNMR(unittest.TestCase):
    def test_segments(self):
        fnmr = '/home/michael/gitrepos/bb-sbbu/DATA_TEST/testA.nmr'
        nmr = NMR(fnmr)
        S = nmr.segments()        
        Sans = [Segment(4, 5), Segment(6, 10), Segment(11,15), Segment(18, 20)]
        self.assertEqual(len(S), len(Sans))
        for i in range(len(S)):
            self.assertTrue(S[i]==Sans[i])


if __name__ == '__main__':
    unittest.main()