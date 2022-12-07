import os
import tqdm
import numpy as np
import networkx as nx
from bb import order_brute, order_greedy, NMR

if __name__ == "__main__":
    nnodes, nedges = 20, 5
    nsamples = 1000
    np.random.seed(1)
    MESSAGE = []
    for k in tqdm.tqdm(range(nsamples)):
        lenE, E = 0, {}
        while lenE < nedges:
            i = np.random.randint(1, nnodes - 3)
            j = np.random.randint(i + 4, nnodes + 1)
            if i not in E:
                E[i] = set()
            if j not in E[i]:
                E[i].add(j)
                lenE += 1

        fn = "DATA_TEST/testRAND_%04d.nmr" % k
        with open(fn, "w") as fid:
            for i in sorted(E):
                for j in sorted(E[i]):
                    fid.write("%d %d\n" % (i, j))

        nmr = NMR(fn)
        orderBF, costBF = order_brute(nmr)
        orderGD, costGD = order_greedy(nmr)
        G = nmr.ordering_graph

        if costBF == costGD:
            os.remove(fn)
        elif nx.number_connected_components(G) > 1:
            MESSAGE.append("Remove disconnected(fn: %s)." % fn)
            os.remove(fn)
        else:
            msg = "Found interesting instance (BF: %03d, GD: %03d, fn: %s)."
            MESSAGE.append(msg % (costBF, costGD, fn))

    for message in MESSAGE:
        print(message)
