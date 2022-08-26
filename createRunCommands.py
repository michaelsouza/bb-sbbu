import os
WDIR = ['DATA_EPSD_00_DMAX_50', 'DATA_EPSD_00_DMAX_60']

with open('callBB.sh', 'w') as fid:
    for wdir in WDIR:
        for fnmr in os.listdir(wdir):
            if not fnmr.endswith('.nmr'):
                continue
            fnmr = os.path.join(wdir, fnmr)
            cmd = "python sources/bb.py '%s'\n" % fnmr        
            fid.write(cmd)