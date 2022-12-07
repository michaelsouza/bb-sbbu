# References:
# 1. https://docs.python.org/3/library/unittest.html

import unittest
from tkinter import SE
from bb import *


class TestNMR(unittest.TestCase):
    def test_segments(self):
        fnmr = "DATA_TEST/testA.nmr"
        nmr = NMR(fnmr)
        S = nmr.segments
        Sans = [
            NMRSegment(4, 5),
            NMRSegment(6, 10),
            NMRSegment(11, 15),
            NMRSegment(18, 20),
        ]
        self.assertEqual(len(S), len(Sans))
        for i in range(len(S)):
            self.assertTrue(S[i] == Sans[i])


class TestBB(unittest.TestCase):
    def test_solveA(self):
        nmr = NMR("DATA_TEST/testA.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        self.assertEqual(costBB, 168)
        for i, eid in enumerate(orderBB):
            if eid == 1:
                id1 = i
            if eid == 2:
                id2 = i
        self.assertLess(id1, id2)

    def test_solveB(self):
        nmr = NMR("DATA_TEST/testB.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])

    def test_solveC(self):
        nmr = NMR("DATA_TEST/testC.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])

    def test_solveD(self):
        nmr = NMR("DATA_TEST/testD.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])

    def test_solveE(self):
        nmr = NMR("DATA_TEST/testE.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])


class TestBBPerm(unittest.TestCase):
    def test_minGT(self):
        E = {9, 15, 5, 20, 18, 7}
        p = BBPerm(E)
        emin = p.minGT(15)
        self.assertEqual(emin, 18)

    def test_start_from(self):
        nmr = NMR("DATA_TEST/testC.nmr")
        E = nmr.E
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


class TestPriorityTree(unittest.TestCase):
    def test_optimality(self):
        wdir = "DATA_TEST"
        for fn in sorted(os.listdir(wdir)):
            if not fn.endswith(".nmr"):
                continue
            nmr = NMR(os.path.join(wdir, fn))
            p = PriorityTree(nmr)
            order, cost = p.solve()
            orderOPT, costOPT = order_brute(nmr)
            self.assertEqual(costOPT, cost)


class TestGreedy(unittest.TestCase):
    def test_optimality(self):
        wdir = "DATA_TEST"
        for c in ["A", "B", "C", "D", "E", "F"]:
            nmr = NMR(os.path.join(wdir, "test%s.nmr" % c))
            orderBF, costBF = order_brute(nmr)
            orderGD, costGD = order_greedy(nmr)
            self.assertEqual(costBF, costGD)


if __name__ == "__main__":
    unittest.main()
