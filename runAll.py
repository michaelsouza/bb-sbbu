import os
import sys
import tqdm
import multiprocessing as mp

WDIR = ['DATA_EPSD_00_DMAX_50', 'DATA_EPSD_00_DMAX_60']

# default parameters
tmax = 3600

# read parameters
for i, arg in enumerate(sys.argv):
    if arg == '-tmax':
        tmax = int(sys.argv[i+1])

FNMR = []
for wdir in sorted(WDIR):
    for fnmr in sorted(os.listdir(wdir)):
        if not fnmr.endswith('.nmr'):
            continue
        FNMR.append(os.path.join(wdir, fnmr))

# sort by file size
FNMR = sorted(FNMR, key=lambda fnmr: os.stat(fnmr).st_size)

ARGV = []
for fnmr in FNMR:
    arg = './build/bb.bin -tmax %d -fnmr %s -clean_log' % (tmax, fnmr)
    # arg = 'python codes/bb.py -tmax %d -fnmr %s -clean_log' % (tmax, fnmr)
    ARGV.append(arg)

print('Saving args.txt')
with open('args.txt', 'w') as fid:
    for arg in ARGV:
        fid.write('%s\n' % arg)

pool = mp.Pool(mp.cpu_count() - 1)
for _ in tqdm.tqdm(pool.imap_unordered(os.system, ARGV), total=len(ARGV)):
    pass
pool.close()
