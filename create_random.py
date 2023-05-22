import os
import tqdm
import numpy as np
import networkx as nx
from codes.bb import order_brute, order_greedy, NMR

def create_edges(nnodes, nedges):
    lenE, E = 0, {}
    while lenE < nedges:
        i = np.random.randint(1, nnodes - 3)
        j = np.random.randint(i + 4, nnodes + 1)
        if i not in E:
            E[i] = set()
        if j not in E[i]:
            E[i].add(j)
            lenE += 1
    return E

if __name__ == "__main__":
    nnodes, nedges = 15, 4
    nsamples = 1000
    np.random.seed(1)
    MESSAGE = []

    # create data/nmr_rand folder
    if not os.path.exists("data/nmr_rand"):
        os.makedirs("data/nmr_rand")
        
    for k in tqdm.tqdm(range(nsamples)):
        E = create_edges(nnodes, nedges)

        fn = f"data/nmr_rand/test{k}_chain_A_dmax_5.nmr"
        with open(fn, "w") as fid:
            for i in sorted(E):
                for j in sorted(E[i]):
                    fid.write("%3d %3d 1 1 X X PRO PRO\n" % (i, j))

        nmr = NMR(fn)
        orderBF, costBF = order_brute(nmr)
        orderGD, costGD = order_greedy(nmr)
        G = nmr.ordering_graph

        if costBF == costGD:
            os.remove(fn)
        elif nx.number_connected_components(G) > 1:
            os.remove(fn)
        else:
            msg = "Found interesting instance (BF: %3d, GD: %3d, fn: %s)."
            MESSAGE.append(msg % (costBF, costGD, fn))

    for message in MESSAGE:
        print(message)
