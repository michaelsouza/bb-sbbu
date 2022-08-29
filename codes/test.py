# References:
# 1. https://docs.python.org/3/library/unittest.html

from tkinter import SE
import unittest

from bb import *


class TestNMR(unittest.TestCase):
    def test_segments(self):
        fnmr = 'DATA_TEST/testA.nmr'
        nmr = NMR(fnmr)
        S = nmr.segments
        Sans = [Segment(4, 5), Segment(6, 10),
                Segment(11, 15), Segment(18, 20)]
        self.assertEqual(len(S), len(Sans))
        for i in range(len(S)):
            self.assertTrue(S[i] == Sans[i])


class TestBB(unittest.TestCase):
    def test_solveA(self):
        fnmr = 'DATA_TEST/testA.nmr'
        nmr = NMR(fnmr)
        E = {edge.eid: edge for edge in nmr.pruneEdges}
        S = {s.sid: s for s in nmr.segments}
        bb = BB(fnmr, E, S)
        orderOPT, costUB = bb.solve()
        self.assertEqual(costUB)

    def test_solveB(self):
        fnmr = 'DATA_TEST/testB.nmr'
        nmr = NMR(fnmr)
        E = {edge.eid: edge for edge in nmr.pruneEdges}
        S = {s.sid: s for s in nmr.segments}
        bb = BB(fnmr, E, S)
        orderOPT, costUB = bb.solve()
    
    def test_solveC(self):
        fnmr = 'DATA_TEST/testC.nmr'
        nmr = NMR(fnmr)
        E = {edge.eid: edge for edge in nmr.pruneEdges}
        S = {s.sid: s for s in nmr.segments}
        bb = BB(fnmr, E, S)
        orderOPT, costUB = bb.solve()


class TestBBPerm(unittest.TestCase):
    def test_start_from(self):
        fnmr = 'DATA_TEST/testC.nmr'
        nmr = NMR(fnmr)
        E = {edge.eid: edge for edge in nmr.pruneEdges}
        bb = BBPerm(E.keys())
        # warming up BB
        bb.next()
        for i in range(10):
            if i % 5 == 0:
                bb.prune()
            bb.next()
        # save old state
        old_order = bb.order.copy()
        old_idx = bb.idx
        # advance 10 steps
        ORDER1 = []
        for i in range(10):
            if i % 5 == 0:
                bb.prune()
            bb.next()
            ORDER1.append(bb.order.copy())
        # restart from old state
        bb.start_from(idx=old_idx, order=old_order)
        # advance 10 steps
        ORDER2 = []
        for i in range(10):
            if i % 5 == 0:
                bb.prune()
            bb.next()
            ORDER2.append(bb.order.copy())
        # ORDER1 and ORDER2 must be equal
        self.assertEqual(len(ORDER1), len(ORDER2))
        for i in range(len(ORDER1)):
            order1 = ORDER1[i]
            order2 = ORDER2[i]
            self.assertEqual(len(order1), len(order2))
            for j in range(len(order1)):
                self.assertEqual(order1[j], order2[j])


if __name__ == '__main__':
    unittest.main()
