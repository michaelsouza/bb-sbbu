# References:
# 1. https://docs.python.org/3/library/unittest.html

import os
import pandas as pd
import unittest
# from tkinter import SE
from bb import *


class TestNMR(unittest.TestCase):
    def test_segments(self):
        fnmr = "data/nmr_test/testA.nmr"
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
        nmr = NMR("data/nmr_test/testA.nmr")
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
        nmr = NMR("data/nmr_test/testB.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])

    def test_solveC(self):
        nmr = NMR("data/nmr_test/testC.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])

    def test_solveD(self):
        nmr = NMR("data/nmr_test/testD.nmr")
        bb = BB(nmr)
        orderBB, costBB = bb.solve()
        orderBB, costBF = order_brute(nmr)
        self.assertEqual(costBB, costBF)
        for i in range(len(nmr.E)):
            self.assertEqual(orderBB[i], orderBB[i])

    def test_solveE(self):
        nmr = NMR("data/nmr_test/testE.nmr")
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
        nmr = NMR("data/nmr_test/testC.nmr")
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
        wdir = os.path.join("data", "nmr_test")
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
        wdir = os.path.join("data", "nmr_test")
        for c in ["A", "B", "C", "D", "E", "F"]:
            nmr = NMR(os.path.join(wdir, "test%s.nmr" % c))
            orderBF, costBF = order_brute(nmr)
            orderGD, costGD = order_greedy(nmr)
            self.assertEqual(costBF, costGD)


class TestInstances(unittest.TestCase):
    def test_dmdgp_constraints(self):        
        backbone_folder = os.path.join('data','backbones')
        nmr_folder = os.path.join('data','nmr')

        PDB = ['1n6t']

        # get the nmr files for 1n6t
        nmr_files = [fn for fn in os.listdir(nmr_folder) if fn.split('_')[0] in PDB and fn.endswith('.nmr')]

        # take the first 3
        for nmr_file in nmr_files:
            nmr_file = os.path.join(nmr_folder, nmr_file)
            # read edges
            E = {}
            with open(nmr_file, 'r') as f:
                for edge in f:
                    i, j, dij = edge.split()[:3]
                    E[int(i), int(j)] = float(dij)

            # read backbone
            backbone_file = os.path.basename(nmr_file).split('_dmax')[0] + '.csv'
            backbone_file = os.path.join(backbone_folder, backbone_file)
            backbone = pd.read_csv(backbone_file)

            # read dmax
            dmax = float(os.path.basename(nmr_file).split('dmax_')[1].split('.')[0])
            
            # check all possible distances
            for i in range(len(backbone)):
                x_i = backbone.iloc[i][['x', 'y', 'z']].values
                i += 1
                for j in range(i+1, len(backbone)):
                    x_j = backbone.iloc[j][['x', 'y', 'z']].values
                    dij = np.linalg.norm(x_i - x_j)

                    # the edge indices are 1-based
                    j += 1

                    # check discretization distances
                    if (j - 3) <= i < j:
                        self.assertTrue((i, j) in E)
                        self.assertAlmostEqual(E[i, j], dij)
                        continue
                    
                    # check possible pruning edges
                    if dij >= dmax:
                        self.assertTrue((i, j) not in E)
                    else:
                        self.assertTrue((i, j) in E)
                        self.assertAlmostEqual(E[i, j], dij)
            

if __name__ == "__main__":
    unittest.main()
