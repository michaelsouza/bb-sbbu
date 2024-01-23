# References:
# 1. https://docs.python.org/3/library/heapq.html
# 2. https://docs.python.org/3/library/stdtypes.html#set-types-set-frozenset

import os
import sys
import time
import copy
import pickle
import numpy as np
import networkx as nx
from functools import cmp_to_key
from itertools import permutations
from heapq import heapify, heappop, heappush


class NMRSegment:
    SID = 0  # class static variable

    def __init__(self, i, j) -> None:
        NMRSegment.SID += 1
        self.sid = NMRSegment.SID
        self.i = i
        self.j = j
        # set of the indexes of the prune edges that cover this segment
        self.eid = set()
        self.update_weight()

    def resetSID():
        NMRSegment.SID = 0

    def add_eid(self, eid):
        self.eid.add(eid)

    def update_weight(self):
        self.weight = int(2**(self.j - self.i + 1))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, NMRSegment):
            return self.i == other.i and self.j == other.j
        return False


class NMREdge:
    EID = 0  # class static variable

    def __init__(self, i, j) -> None:
        NMREdge.EID += 1
        self.eid = NMREdge.EID  # edge id
        self.i = i
        self.j = j
        self.sid = set()  # set of the indexes of the segments covered by this edge

    def add_sid(self, sid):
        self.sid.add(sid)

    def check_cover(self, s):
        return self.i + 3 <= s.i and s.j <= self.j

    def resetEID():
        NMREdge.EID = 0


class NMR:
    def __init__(self, fnmr: str) -> None:
        self.fnmr = fnmr
        NMREdge.resetEID()
        self.edges = []
        with open(fnmr, 'r') as fid:
            for row in fid:
                i = int(row.split()[0])
                j = int(row.split()[1])
                self.edges.append(NMREdge(i, j))
        self.nnodes = np.max([edge.j for edge in self.edges])
        self.pruneEdges = [edge for edge in self.edges if edge.j > edge.i + 3]
        self.segments = self._segments()
        self.E, self.S = self._ordering_data()

    def _segments(self):
        NMRSegment.resetSID()
        # I: sorted list of all atoms covered by prune edges
        I = set()
        for edge in self.pruneEdges:
            I.update(np.arange(edge.i + 3, edge.j + 1))
        I = sorted(list(I))

        # E[i]: set of the eid's of the edges covering atom 'i'
        E = {i: set() for i in I}
        for edge in self.pruneEdges:
            for i in range(edge.i + 3, edge.j+1):
                E[i].add(edge.eid)

        # Consecutive atoms covered by the same set of edges belong to the same segment
        S = []  # S: list of all segments
        s = NMRSegment(I[0], I[0])
        for j in I:
            if E[s.i] == E[j]:
                s.j = j
            else:
                # the s.i and s.j were possibly updated inside the loop, so it's necessary
                # to update/correct the s.weight.
                s.update_weight()
                S.append(s)
                s = NMRSegment(j, j)
        s.update_weight()
        S.append(s)

        # O(len(S) * len(self.pruneEdges))
        for s in S:
            for edge in self.pruneEdges:
                if edge.check_cover(s):
                    edge.add_sid(s.sid)
                    s.add_eid(edge.eid)
        return S

    def ordering_graph(self, use_weight=False):
        G = nx.Graph()
        E, S = self.E, self.S
        
        # add segment-vertex
        def s_lbl(sid): return '%d:%d' % (S[sid].i, S[sid].j)
        for sid in S:
            G.add_node(s_lbl(sid), weight=S[sid].weight)

        # add edge-vertex
        for eid in E:
            e_lbl = '%d' % eid
            G.add_node(e_lbl)
            # add hyperedge
            for sid in E[eid].sid:
                G.add_edge(e_lbl, s_lbl(sid))
        return G

    def _ordering_data(self):
        E = {e.eid: e for e in self.pruneEdges}
        S = {s.sid: s for s in self.segments}
        return E, S


def order_cost(order, E, S, costUB=np.inf):
    total_cost = 0  # total cost
    # B[sid] is set to true by the first edge that covers it.
    B = set()
    for eid in order:
        edge_cost = 1
        for sid in E[eid].sid:
            # first edge to cover sid
            if sid not in B and sid in S:
                edge_cost *= S[sid].weight
                B.add(sid)
        # the cost of an edge that covers no segment is zero (not one)
        total_cost += edge_cost if edge_cost > 1 else 0
        if total_cost >= costUB:
            return np.inf
    return total_cost


def order_sbbu(nmr):
    E = nmr.E
    order = list(E)  # list of edges eid
    order = sorted(order, key=lambda eid: (E[eid].j, E[eid].i))
    return order, order_cost(order, E, nmr.S)


def order_brute(nmr: NMR):
    '''Evaluate all posible permutations.'''
    E, S = nmr.E, nmr.S
    orderOPT, costOPT = [], np.inf
    for p in permutations(E):
        c = order_cost(p, E, S)
        if c < costOPT:
            costOPT = c
            orderOPT = list(p).copy()
    return orderOPT, costOPT


def order_greedy(nmr:NMR):
    E, S = copy.deepcopy(nmr.E), copy.deepcopy(nmr.S)
    order = []
    while len(E) > 0:
        minEid, minC = None, np.inf        
        for eid in E:
            e: NMREdge = E[eid]
            c = 1
            for sid in e.sid:
                if sid in S:
                    s: NMRSegment = S[sid]
                    c *= s.weight
            if c < minC:
                minEid, minC = eid, c
        order.append(minEid)
        for sid in E[minEid].sid:
            S.pop(sid, None)
        E.pop(minEid)
    return order, order_cost(order, nmr.E, nmr.S)



def cost_relax(U, S):
    total_cost = 0
    for sid in U:
        total_cost += S[sid].weight
    return total_cost


class BBPerm:
    def __init__(self, keys) -> None:
        self.keys = list(keys)
        # index of the element last element inserted
        self.idx = -1
        # current order
        self.order = np.zeros(len(keys), dtype=int)
        # heap of available items
        self.h = [key for key in keys]
        heapify(self.h)
        self.state = 'n'  # n:normal, p:prune

    def next(self):
        if self.state == 'n':
            if len(self.h) > 0:
                emin = heappop(self.h)
                self.idx += 1
                self.order[self.idx] = emin
                return emin
            self.state = 'p'
            return self.next()
        # self.state == 'p'
        if self.idx == -1:
            return None
        emin = self.minGT(self.order[self.idx])
        heappush(self.h, self.order[self.idx])
        if emin is None:
            # set invalid value
            self.order[self.idx] = -1
            self.idx -= 1
            return self.next()
        self.order[self.idx] = emin
        self.state = 'n'
        return emin

    def minGT(self, elem):
        # returns the smallest element bigger than elem
        if len(self.h) == 0:
            return None
        self.buffer = []
        while len(self.h) > 0:
            item = heappop(self.h)
            if item > elem:
                break
            self.buffer.append(item)
        # restore elements on buffer
        for e in self.buffer:
            heappush(self.h, e)
        return item if item > elem else None

    def prune(self):
        self.state = 'p'

    def start_from(self, idx, order):
        self.idx = idx
        self.order[:(idx+1)] = order[:(idx+1)]
        self.state = 'n'
        S = set(order[:(idx+1)])
        self.h = []
        for e in self.keys:
            if e not in S:
                heappush(self.h, e)


class BB:
    def __init__(self, nmr: NMR) -> None:
        self.nmr = nmr
        self.E, self.S = nmr.E, nmr.S
        self.nedges = len(self.E)
        self.idx = -1
        self.perm = BBPerm(nmr.E)
        self.order = np.zeros(self.nedges, dtype=int)
        self.timeout = False

    def order_rem(self, C, U):
        # Returns the total_cost of the eids removed from self.order
        # Remark: The dictionary C is updated.
        idx = self.idx
        total_cost = 0
        while idx >= self.perm.idx and self.perm.order[idx] != self.order[idx]:
            eid = self.order[idx]
            # set invalid value
            self.order[idx] = -1
            eid_cost = 1
            for sid in self.E[eid].sid:
                C[sid] -= 1
                if C[sid] == 0:
                    eid_cost *= self.S[sid].weight
                    U.add(sid)
            if eid_cost > 1:
                total_cost += eid_cost
            idx -= 1
        return total_cost

    def order_add(self, eid, C, U):
        ''' Add eid edge to the self.order and update C and U.
            C[sid] : number of edges on self.order covering segment sid
            U: set of the uncovered segments.
        '''
        self.order[self.perm.idx] = eid
        eid_cost = 1
        for sid in self.E[eid].sid:
            C[sid] += 1
            # the current eid is the only one covering the sid
            if C[sid] == 1:
                eid_cost *= self.S[sid].weight
                U.remove(sid)
        return eid_cost if eid_cost > 1 else 0

    def load(self):
        fname = self.nmr.fnmr.replace('.nmr', '.pkl')
        print('> unpicliking', fname)
        with open(fname, 'rb') as fid:
            data = pickle.load(fid)
        self.idx = data['idx']
        self.order = data['order']
        self.orderOPT = data['orderOPT']
        self.costUB = data['costUB']
        C, U = data['C'], data['U']
        return C, U

    def dump(self, C, U):
        fname = self.nmr.fnmr.replace('.nmr', '.pkl')
        print('> picliking', fname)
        data = {}
        data['idx'] = self.idx
        data['order'] = self.order
        data['orderOPT'] = self.orderOPT
        data['costUB'] = self.costUB
        data['U'] = U
        data['C'] = C
        with open(fname, 'wb') as fid:
            pickle.dump(data, fid)

    def solve(self, unpickling=False, tmax=60):
        tic = time.time()
        if unpickling:
            self.load()
        else:
            # initial optimal solution
            self.orderOPT, self.costUB = order_sbbu(self.nmr)

        # C[sid] : number of edges already included in the order that cover segment sid
        C = {sid: 0 for sid in self.S}

        # U: set of the uncovered segments
        U = set([sid for sid in C])

        # first cost_relax
        costLB = cost_relax(U, self.S)
        if costLB == self.costUB:
            return self.orderOPT, self.costUB

        partial_cost = 0
        eid = self.perm.next()
        # loop through all permutations
        while eid is not None:
            partial_cost -= self.order_rem(C, U)
            partial_cost += self.order_add(eid, C, U)
            self.idx = self.perm.idx
            # when U is empty, the partial_cost is total.
            costLB = partial_cost + cost_relax(U, self.S)
            toc = time.time() - tic
            if toc > tmax:
                self.timeout = True
                print('> timeoutBB %f seconds' % toc)
                self.dump(C, U)
                break
            if costLB >= self.costUB:
                self.perm.prune()
            elif self.perm.idx == (self.nedges - 1) and costLB < self.costUB:
                self.costUB = costLB
                self.orderOPT[:] = self.order
            eid = self.perm.next()
        return self.orderOPT, self.costUB


def write_log(fid, line):
    print(line)
    fid.write(line + '\n')


class PriorityTree:
    def __init__(self, nmr: NMR) -> None:
        self.nmr = nmr
        self.E, self.S = nmr.E, nmr.S
        # sort edges by the number of segments
        E = sorted(self.E, key=lambda eid: len(self.E[eid].sid), reverse=True)
        # sort segments by the number of edges
        b = {sid:False for sid in self.S}
        self.ordS = []
        for eid in E:
            s = []
            for sid in self.E[eid].sid:
                if b[sid]:
                    continue
                b[sid] = True
                s.append(sid)
            self.ordS += sorted(s, key=lambda sid: len(self.S[sid].eid), reverse=True)
        # Ek: number of uncovered segments
        # When Ek is zero, we can calculate the cost of solving the edge eid
        self.Ek = {eid:len(self.E[eid].sid) for eid in self.E}
        # G: graph of priorities
        self.G = nx.DiGraph()
        self.order, self.cost = order_sbbu(self.nmr)
        self.timeout = False

    def check_path(self, eidA, eidB):
        '''eidA: source
           eidB: target
        '''
        try:
            # A > B 
            pAB = nx.shortest_path(self.G, source=eidA, target=eidB)
        except:
            pAB = None
        return pAB is not None
                

    def available_edges(self, E):
        C = {eid:True for eid in E}
        for i in range(len(E)):
            eidA = E[i]
            if (eidA not in self.G) or (C[eidA] == False):
                continue
            for j in range(i+1,len(E)):
                eidB = E[j]
                if self.check_path(eidA, eidB):
                    C[eidA] = False
                    break
                if self.check_path(eidB, eidA):
                    C[eidB] = False
        return sorted([eid for eid in C if C[eid]])

    def add_precedence(self, eidA, E):
        # eidA < eidB
        P = []
        for eidB in E:
            if (eidB != eidA) and (self.G.has_edge(eidB, eidA) == False):
                self.G.add_edge(eidB, eidA)
                P.append((eidB, eidA))
        return P

    def edge_cost(self, c_eid, eid, costUB):
        cost = 1
        for sid in self.E[eid].sid:
            if c_eid[sid] == eid:
                s: NMRSegment = self.S[sid]
                cost *= s.weight
                if cost >= costUB:
                    cost = costUB
                    break
        return 0 if cost == 1 else cost

    def add_cost(self, sid, c_eid, costUB):
        cost = 0
        for eid in self.S[sid].eid:
            if self.Ek[eid] == 0:
                continue
            self.Ek[eid] -= 1
            if self.Ek[eid] == 0 and cost < costUB:
                cost += self.edge_cost(c_eid, eid, costUB)
        return cost

    def rem_cost(self, i, sid, costADD):
        costREM = costADD[i]
        costADD[i] = 0
        for eid in self.S[sid].eid:
            self.Ek[eid] += 1
        return costREM

    def backtracking(self, level, E, P, c_idx, c_eid, cost, costADD):
        while level >= 0:
            sid = self.ordS[level]
            cost -= self.rem_cost(level, sid, costADD)
            c_eid[sid] = None
            self.G.remove_edges_from(P[level])
            P[level] = []
            if c_idx[level] < (len(E[level]) - 1):
                c_idx[level] += 1
                break
            E[level] = []
            c_idx[level] = 0
            level -= 1
        return level if level >= 0 else None, cost

    def save_order(self, c_eid:dict):
        '''Convert from c_eid (dict) to self.order (list)'''
        self.order = np.unique([c_eid[sid] for sid in c_eid])
        # remove all nodes with zero degree
        d = dict(self.G.degree())
        self.G.remove_nodes_from([eid for eid in d if d[eid] == 0])
        # shortest paths
        P = {eid:None for eid in self.order}
        for eid in self.order:
            P[eid] = nx.shortest_path(self.G, source=eid)
        def cmp(eidA, eidB):
            if eidA in P[eidB]:
                return -1
            if eidB in P[eidA]:
                return +1
            return 0
        self.order = sorted(self.order, key=cmp_to_key(cmp))

    def solve(self,tmax=60):
        # init cost_relax
        costLB = cost_relax(self.S, self.S)
        if costLB == self.cost:
            return self.order, self.cost
        # c: vector of each segment choice
        c_eid = {sid:None for sid in self.S}
        c_idx = np.zeros(len(self.ordS), dtype=int)
        costADD = c_idx.copy()  
        level, cost = 0, 0 # index of the current segment
        # E[i]: edges available at level 'i'
        E = [[] for _ in range(len(c_idx))]
        # P[i]: set of pairs precedence (eidA, eidB) added at level 'i'
        P = [[] for _ in range(len(c_idx))]
        tic = time.time()
        while level is not None:
            toc = time.time() - tic
            if toc > tmax:
                self.timeout = True
                print('> timeoutBB %f seconds' % toc)
                return self.order, self.cost
            sid = self.ordS[level]
            if len(E[level]) == 0:
                E[level] = self.available_edges(list(self.S[sid].eid))
            eid = E[level][c_idx[level]]
            c_eid[sid] = eid
            P[level] = self.add_precedence(eid, E[level])
            costADD[level] = self.add_cost(sid, c_eid, self.cost)
            cost += costADD[level]
            # solution found
            if (cost < self.cost) and (level == (len(self.ordS) - 1)):
                self.cost = cost
                self.save_order(c_eid) 
            # next
            if (cost < self.cost) and (level < (len(self.ordS) - 1)):
                level += 1
            else:
                level, cost = self.backtracking(level, E, P, c_idx, c_eid, cost, costADD)
        return self.order, self.cost


def call_solvers(*argv):
    fnmr = '/home/michael/gitrepos/bb-sbbu/DATA_TEST/testC.nmr'
    tmax = 1
    clean_log = False
    for i, arg in enumerate(argv):
        if arg == '-fnmr':
            fnmr = argv[i+1]
        if arg == '-tmax':
            tmax = float(argv[i+1])
        if arg == '-clean_log':
            clean_log = True

    flog = fnmr.replace('.nmr', '.log')
    # check if already has a log file
    if not clean_log and os.path.exists(flog):
        print('> skip (already solved) %s' % fnmr)
        sys.exit(0)
    # create log file
    fid = open(flog, 'w')
    write_log(fid, '> fnmr ' + fnmr)

    # read instance
    nmr = NMR(fnmr)
    E, S = nmr.E, nmr.S

    write_log(fid, '> tmax (secs) ....... %g' % tmax)
    write_log(fid, '> nnodes ............ %d' % nmr.nnodes)
    write_log(fid, '> lenE .............. %d' % len(E))
    write_log(fid, '> lenS .............. %d' % len(S))

    costRELAX = cost_relax(S, S)
    write_log(fid, '> costRX ............ %d' % costRELAX)

    # call order_greedy
    tic = time.time()
    orderGREEDY, costGREEDY = order_greedy(nmr)
    toc = time.time() - tic
    write_log(fid, '> costGD ............ %d' % costGREEDY)
    write_log(fid, '> timeGD (secs) ..... %g' % toc)

    # call order_sbbu
    tic = time.time()
    orderSBBU, costSBBU = order_sbbu(nmr)
    toc = time.time() - tic
    write_log(fid, '> costSB ............ %d' % costSBBU)
    write_log(fid, '> timeSB (secs) ..... %g' % toc)

    # call priority_tree
    tic = time.time()
    pt = PriorityTree(nmr)
    orderPT, costPT = pt.solve()
    toc = time.time() - tic
    write_log(fid, '> costPT ............ %d' % costPT)
    write_log(fid, '> timePT (secs) ..... %g' % toc)

    # call order_bb
    # tic = time.time()
    # bb = BB(nmr)
    # costBB, costBB = bb.solve(tmax=tmax)
    # toc = time.time() - tic
    # write_log(fid, '> timeoutBB ......... %s' % bb.timeout)
    # write_log(fid, '> costBB ............ %d' % costBB)
    # write_log(fid, '> timeBB (secs) ..... %g' % toc)

    fid.close()


if __name__ == '__main__':
    call_solvers(*sys.argv)
