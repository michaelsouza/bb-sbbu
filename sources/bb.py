# References:
# 1. https://docs.python.org/3/library/heapq.html
# 2. https://docs.python.org/3/library/stdtypes.html#set-types-set-frozenset

import sys
import time
import numpy as np
from heapq import heapify, heappop, heappush


class Segment:
    SID = 0  # class static variable

    def __init__(self, i, j) -> None:
        Segment.SID += 1
        self.sid = Segment.SID
        self.i = i
        self.j = j
        # set of the indexes of the prune edges that cover this segment
        self.eid = set()

    @property
    def weight(self):
        return int(2**(self.j - self.i + 1))

    def add_eid(self, eid):
        self.eid.add(eid)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Segment):
            return self.i == other.i and self.j == other.j
        return False


class EdgeNMR:
    EID = 0  # class static variable

    def __init__(self, i, j) -> None:
        EdgeNMR.EID += 1
        self.eid = EdgeNMR.EID  # edge id
        self.i = i
        self.j = j
        self.sid = set()  # set of the indexes of the segments covered by this edge

    def add_sid(self, sid):
        self.sid.add(sid)

    def check_cover(self, s):
        return self.i + 3 <= s.i and s.j <= self.j


class NMR:
    def __init__(self, fnmr) -> None:
        self.edges = []
        with open(fnmr, 'r') as fid:
            for row in fid:
                i = int(row.split()[0])
                j = int(row.split()[1])
                self.edges.append(EdgeNMR(i, j))
        self.nnodes = np.max([np.max([edge.i, edge.j]) for edge in self.edges])
        self.pruneEdges = [edge for edge in self.edges if edge.j > edge.i + 3]

    def segments(self):
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
        s = Segment(I[0], I[0])
        for j in I:
            if E[s.i] == E[j]:
                s.j = j
            else:
                S.append(s)
                s = Segment(j, j)
        S.append(s)

        # O(len(S) * len(self.pruneEdges))
        for s in S:
            for edge in self.pruneEdges:
                if edge.check_cover(s):
                    edge.add_sid(s.sid)
                    s.add_eid(edge.eid)
        return S


def order_cost(order, E, S):
    total_cost = 0  # total cost
    # B[sid] is set to true by the first edge that covers it.
    B = {sid: False for sid in S}
    for eid in order:
        edge_cost = 1
        for sid in E[eid].sid:
            # first edge to cover sid
            if B[sid] == False:
                edge_cost *= S[sid].weight
                B[sid] = True
        # the cost of an edge that covers no segment is zero (not one)
        total_cost += edge_cost if edge_cost > 1 else 0
    return total_cost


def order_sbbu(E, S):
    order = list(E)  # list of edges eid
    order = sorted(order, key=lambda eid: (E[eid].j, E[eid].i))
    return order, order_cost(order, E, S)


def cost_relax(U, S):
    total_cost = 0
    for sid in U:
        total_cost += S[sid].weight
    return total_cost


class BBPerm:
    def __init__(self, elems) -> None:
        # index of the element last element inserted
        self.idx = -1
        # current order
        self.order = np.zeros(len(elems), dtype=int)
        # heap of available items
        self.h = [elem for elem in elems]
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


def order_rem(bb, idx, order,  E, S, C, U):
    # return the total_cost of the removed eids
    # Remark: The dictionary C is updated.
    total_cost = 0
    while idx >= 0 and bb.order[idx] != order[idx]:
        eid = order[idx]
        # set invalid value
        order[idx] = -1
        cost_edge = 1
        for sid in E[eid].sid:
            C[sid] -= 1
            if C[sid] == 0:
                cost_edge *= S[sid].weight
                U.add(sid)
        if cost_edge > 1:
            total_cost += cost_edge
        idx -= 1
    return total_cost


def order_add(bb, eid, order, E, S, C, U):
    # add eid to
    # remark: order and C are updated
    order[bb.idx] = eid
    total_cost = 0
    cost_edge = 1
    for sid in E[eid].sid:
        C[sid] += 1
        # the current eid is the only one covering the sid
        if C[sid] == 1:
            cost_edge *= S[sid].weight
            U.remove(sid)
    if cost_edge > 1:
        total_cost += cost_edge
    return total_cost


def order_bb(E, S):
    # initial optimal solution
    orderOPT, costUB = order_sbbu(E, S)
    # loop through all permutations
    idx = 0  # index of the permutation component
    # C[sid] : number of edges already included in the order that cover segment sid
    C = {sid: 0 for sid in S}
    # U: set of the uncovered sets
    U = set(S)
    bb = BBPerm(E)
    idx = -1
    nedges = len(E)
    order = np.zeros(len(E), dtype=int)
    partial_cost = 0
    eid = bb.next()
    while eid is not None:
        partial_cost -= order_rem(bb, idx, order, E, S, C, U)
        partial_cost += order_add(bb, eid, order, E, S, C, U)
        idx = bb.idx
        # whe U is empty, the partial_cost is total.
        costLB = partial_cost + cost_relax(U, S)
        if costLB >= costUB:
            bb.prune()
        elif bb.idx == (nedges - 1) and costLB < costUB:
            costUB = costLB
            orderOPT[:] = order
        eid = bb.next()
    return orderOPT, costUB


def write_log(fid, line):
    print(line)
    fid.write(line + '\n')


if __name__ == '__main__':
    fnmr = '/home/michael/gitrepos/bb-sbbu/DATA_TEST/testB.nmr'
    if len(sys.argv) > 1:
        fnmr = sys.argv[1]
    # create log file
    flog = fnmr.replace('.nmr', '.log')
    fid = open(flog, 'w')
    write_log(fid, '> fnmr: ' + fnmr)
    
    # read instance
    nmr = NMR(fnmr)
    E = {edge.eid: edge for edge in nmr.pruneEdges}
    S = {s.sid: s for s in nmr.segments()}

    write_log(fid, '> nnodes ............ %d' % nmr.nnodes)
    write_log(fid, '> lenE .............. %d' % len(E))
    write_log(fid, '> lenS .............. %d' % len(S))

    costRELAX = cost_relax(S, S)
    write_log(fid, '> costRELAX ......... %d' % costRELAX)
    
    # call order_sbbu
    tic = time.time()    
    orderSBBU, costSBBU = order_sbbu(E, S)
    toc = time.time() - tic
    write_log(fid, '> costSBBU .......... %d' % costSBBU)
    write_log(fid, '> timeSBBU (secs) ... %g' % toc)

    # call order_bb if needed
    if costRELAX < costSBBU:
        tic = time.time()
        costBB, costBB = order_bb(E, S)
        toc = time.time() - tic
        write_log(fid, '> costBB ............ %d' % costBB)
        write_log(fid, '> timeBB (secs) ..... %g' % toc)
    fid.close()
