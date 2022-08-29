# References:
# 1. https://docs.python.org/3/library/itertools.html

from itertools import permutations, chain
from bb import *


def order_bb_memo(E, S, costUB=np.inf, memo={}):
    # empty subproblem
    if len(S) == 0:
        return tuple(E), 0
    # remember a solution with more relaxed upper bound
    Skey = frozenset(S)
    if Skey in memo:
        m = memo[Skey]
        if m['costUB'] < costUB:
            return m['orderOPT'], m['costUB']
        return [], costUB
    
    # U: set of the uncovered segments
    U = set(S)
    orderOPT = []
    costRELAX = cost_relax(S, S)
    # k: split index 
    k = 3 if len(S) > 3 else len(S)    
    # loop on all prefixes of length k
    bb = BBPerm(E)
    idx = -1
    order = np.zeros(len(E), dtype=int)
    eid = bb.next()
    while eid is not None:
        costP -= order_rem(bb, idx, order, E, S, C, U)
        costP += order_add(bb, eid, order, E, S, C, U)
        if costP >= costUB:
            continue
        F = tuple()
        if k < len(S):
            # C: set of segments covered by P
            C = list(chain.from_iterable(E[e].sid for e in P))
            C = set(C)
            # U: set of segments not covered by P
            U = {s:S[s] for s in S if s not in C}
            costLB = costP + cost_relax(U, S)
            if costLB >= costUB:
                continue
            Pkey = frozenset(P)
            # create suffix problem        
            F = {e:E[e] for e in E if e not in Pkey}
            F, costF = order_bb_memo(F, U, costUB - costP, memo)
            costP = costP + costF
            if costP >= costUB:
                continue
        orderOPT = P + F
        costUB = costP
        # best solution found
        if costRELAX == costUB:
            break
    # add/update memory
    memo[Skey] = {'orderOPT': orderOPT, 'costUB': costUB}
    return orderOPT, costUB


if __name__ == '__main__':
    fnmr = '/home/michael/gitrepos/bb-sbbu/DATA_TEST/testC.nmr'
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

    tic = time.time()
    costBB, costBB = order_bb_memo(E, S, costSBBU)
    toc = time.time() - tic
    write_log(fid, '> costBB ............ %d' % costBB)
    write_log(fid, '> timeBB (secs) ..... %g' % toc)
    fid.close()
