import os
import tqdm
import numpy as np
from bb import order_brute, order_greedy, NMR

if __name__ == '__main__':
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

        fn = 'DATA_TEST/testRAND_%04d.nmr' % k
        with open(fn, 'w') as fid:
            for i in sorted(E):
                for j in sorted(E[i]):
                    fid.write('%d %d\n' % (i, j))

        nmr = NMR(fn)
        orderBF, costBF = order_brute(nmr)
        orderGD, costGD = order_greedy(nmr)

        if costBF != costGD:
            MESSAGE.append('Found interesting instance (BF: %03d, GD: %03d, fn: %s).' % (
                costBF, costGD, fn))
        else:
            os.remove(fn)
    for message in MESSAGE:
        print(message)
