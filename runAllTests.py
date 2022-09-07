import os
import sys
import multiprocessing as mp
import tqdm

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
        arg = 'build/bb.bin -tmax %g -fnmr %s -clean_log' % (tmax, os.path.join(wdir, fnmr))
        ARGV.append(arg)
        
pool = mp.Pool(mp.cpu_count() - 1)
for _ in tqdm.tqdm(pool.map(os.system, ARGV), total=len(ARGV)):
    pass
pool.close()
