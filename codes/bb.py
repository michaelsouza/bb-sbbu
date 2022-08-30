# References:
# 1. https://docs.python.org/3/library/heapq.html
# 2. https://docs.python.org/3/library/stdtypes.html#set-types-set-frozenset

import os
import sys
import time
import pickle
import numpy as np
import networkx as nx
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

    @property
    def weight(self):
        return int(2**(self.j - self.i + 1))

    def resetSID():
        NMRSegment.SID = 0

    def add_eid(self, eid):
        self.eid.add(eid)

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
    def __init__(self, fnmr) -> None:
        self.fnmr = fnmr
        NMREdge.resetEID()
        self.edges = []
        with open(fnmr, 'r') as fid:
            for row in fid:
                i = int(row.split()[0])
                j = int(row.split()[1])
                self.edges.append(NMREdge(i, j))
        self.nnodes = np.max([np.max([edge.i, edge.j]) for edge in self.edges])
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

        # E[i]: set of the indexes of the edges covering atom 'i'
        E = {i: set() for i in I}
        for eid, edge in enumerate(self.pruneEdges):
            for i in range(edge.i + 3, edge.j+1):
                E[i].add(eid)

        # Consecutive atoms covered by the same set of edges belong to the same segment
        S = []  # S: list of all segments
        s = NMRSegment(I[0], I[0])
        for j in I:
            if E[s.i] == E[j]:
                s.j = j
            else:
                S.append(s)
                s = NMRSegment(j, j)
        S.append(s)

        # O(len(S) * len(self.pruneEdges))
        for s in S:
            for edge in self.pruneEdges:
                if edge.check_cover(s):
                    edge.add_sid(s.sid)
                    s.add_eid(edge.eid)
        return S

    @property
    def ordering_graph(self):
        G = nx.Graph()
        E, S = self._ordering_data
        s_lbl = lambda sid: '%d:%d' % (S[sid].i, S[sid].j)
        for sid in S:                
            G.add_node(s_lbl(sid), weight=S[sid].weight)
        
        for eid in E:
            e_lbl = '%d' % eid
            G.add_node(e_lbl)
            for sid in E[eid].sid:            
                G.add_edge(e_lbl, s_lbl(sid))
        return G
    
    def _ordering_data(self):
        E = {e.eid:e for e in self.pruneEdges}
        S = {s.sid:s for s in self.segments}
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


def cost_relax(U, S):
    total_cost = 0
    for sid in U:
        total_cost += S[sid].weight
    return total_cost


class BBPerm:
    def __init__(self, keys) -> None:
        self.elems = list(keys)
        # index of the element last element inserted
        self.idx = -1
        # current order
        self.order = np.zeros(len(keys), dtype=int)
        # heap of available items
        self.h = [elem for elem in keys]
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
        if emin is None:
            heappush(self.h, self.order[self.idx])
            # set invalid value
            self.order[self.idx] = -1
            self.idx -= 1
            return self.next()
        heappush(self.h, self.order[self.idx])
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
        for e in self.elems:
            if e not in S:
                heappush(self.h, e)


class BB:
    def __init__(self, nmr: NMR) -> None:
        self.nmr = nmr
        self.E, self.S = nmr.E, nmr.S
        self.nedges = len(self.E)
        self.idx = -1
        self.perm = BBPerm(nmr.E)
        self.order = np.zeros(len(self.E), dtype=int)
        self.timeout = False

    def order_rem(self, C, U):
        # return the total_cost of the removed eids
        # Remark: The dictionary C is updated.
        idx = self.idx
        total_cost = 0
        while idx >= 0 and self.perm.order[idx] != self.order[idx]:
            eid = self.order[idx]
            # set invalid value
            self.order[idx] = -1
            cost_edge = 1
            for sid in self.E[eid].sid:
                C[sid] -= 1
                if C[sid] == 0:
                    cost_edge *= self.S[sid].weight
                    U.add(sid)
            if cost_edge > 1:
                total_cost += cost_edge
            idx -= 1
        return total_cost

    def order_add(self, eid, C, U):
        # add eid to
        # remark: order and C are updated
        self.order[self.perm.idx] = eid
        total_cost = 0
        cost_edge = 1
        for sid in self.E[eid].sid:
            C[sid] += 1
            # the current eid is the only one covering the sid
            if C[sid] == 1:
                cost_edge *= self.S[sid].weight
                U.remove(sid)
        if cost_edge > 1:
            total_cost += cost_edge
        return total_cost

    def load(self):
        fname = self.nmr.fnmr.replace('.nmr','.pkl')
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
        fname = self.nmr.fnmr.replace('.nmr','.pkl')
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
        for i in range(self.idx + 1):
            eid = self.order[i]
            for sid in self.E[eid].sid:
                C[sid] += 1

        # U: set of the uncovered segments
        U = set([sid for sid in C if C[sid] == 0])

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
            # whe U is empty, the partial_cost is total.
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

def call_solvers(*argv):
    fnmr = '/home/michael/gitrepos/bb-sbbu/DATA_TEST/testC.nmr'
    tmax = 1
    for i, arg in enumerate(argv):
        if arg == '-fnmr':
            fnmr = argv[i+1]
        if arg == '-tmax':
            tmax = float(argv[i+1])            
    
    flog = fnmr.replace('.nmr', '.log')    
    # check if already has a log file
    if os.path.exists(flog):
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
    write_log(fid, '> costRELAX ......... %d' % costRELAX)

    # call order_sbbu
    tic = time.time()
    orderSBBU, costSBBU = order_sbbu(nmr)
    toc = time.time() - tic
    write_log(fid, '> costSBBU .......... %d' % costSBBU)
    write_log(fid, '> timeSBBU (secs) ... %g' % toc)

    # call order_bb if needed
    tic = time.time()
    bb = BB(nmr)
    costBB, costBB = bb.solve(tmax=tmax)
    toc = time.time() - tic
    write_log(fid, '> timeoutBB ......... %s' % bb.timeout)
    write_log(fid, '> costBB ............ %d' % costBB)
    write_log(fid, '> timeBB (secs) ..... %g' % toc)

    fid.close()

if __name__ == '__main__':
    call_solvers(sys.argv)
