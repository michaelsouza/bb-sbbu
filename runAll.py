import os
import sys
import multiprocessing as mp
from codes.bb import call_solvers

WDIR = ['DATA_EPSD_00_DMAX_50', 'DATA_EPSD_00_DMAX_60']

# default parameters
tmax = 60

# read parameters
for i, arg in enumerate(sys.argv):
    if arg == '-tmax':
        tmax = int(sys.argv[i+1])

ARGV = []
for wdir in sorted(WDIR):
    for fnmr in sorted(os.listdir(wdir)):
        if not fnmr.endswith('.nmr'):
            continue
        arg = '-tmax %g -fnmr %s' % (tmax, os.path.join(wdir, fnmr))
        ARGV.append(arg.split())

pool = mp.Pool(mp.cpu_count() - 1)
pool.starmap(call_solvers, ARGV)
pool.close()
