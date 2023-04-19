import os
import sys
import tqdm
import multiprocessing as mp


def remove_logs(WDIR: list):
    # clean log files
    for wdir in WDIR:
        print('Cleaning wdir=%s' % wdir)
        for fn in os.listdir(wdir):
            fn = os.path.join(wdir, fn)
            # remove log files
            if fn.endswith('.log'):
                try:
                    os.remove(fn)
                except OSError as e:
                    print('Could not remove %s: %s' % (fn, e))


def get_nmr_files(WDIR: list):
    # get all nmr files
    FNMR = []
    for wdir in WDIR:
        print('Getting nmr files from wdir=%s' % wdir)
        for fn in os.listdir(wdir):
            fn = os.path.join(wdir, fn)
            # add nmr files
            if fn.endswith('.nmr'):
                FNMR.append(fn)

    # sort by file size
    FNMR = sorted(FNMR, key=lambda fnmr: os.stat(fnmr).st_size)

    return FNMR


def get_command_lines(FNMR: str, tmax: int, clean_log: bool, dump_only: bool, solvers:list=['BB'], verbose:bool=False):
    CMD = [] # list of arguments
    for fnmr in FNMR:
        for solver in solvers:
            cmd = f'./build/bb.bin -tmax {tmax} -fnmr {fnmr} -solver {solver}'
            # add dump flag
            if dump_only:
                cmd += ' -dump'
            # add clean_log flag
            if clean_log:
                cmd += ' -clean_log'
            if verbose:
                cmd += ' -verbose'
            # arg = 'python codes/bb.py -tmax %d -fnmr %s -clean_log' % (tmax, fnmr)
            CMD.append(cmd)

    # write all cmd to file
    with open('run_all.cmd', 'w') as fd:
        fd.write('\n'.join(CMD))

    return CMD


if __name__ == "__main__":
    # set default parameters
    tmax = 3600  # 1 hour
    wdir = 'data/nmr_test' # list of directories to run
    dump = False # if True, call the solvers with -dump 
    clean_log = False # if True, call the solvers with -clean_log
    solvers = ['BF']
    verbose = False
    ncpu = mp.cpu_count() - 1

    # read parameters
    for i, arg in enumerate(sys.argv):
        if arg == '-tmax':
            tmax = int(sys.argv[i+1])
        elif arg == '-wdir':
            wdir = sys.argv[i+1].split(',')
        elif arg == '-solvers':
            solvers = sys.argv[i+1].split(',')
        elif arg == '-dump':
            dump = True
        elif arg == '-clean_log':
            clean_log = True
        elif arg == '-verbose':
            verbose = True
        elif arg == '-ncpu':
            ncpu = int(sys.argv[i+1])
        elif arg == '-help':
            print('Usage: python runAll.py [-tmax <int>] [-wdir <str>]')
            print('   -tmax <int>: maximum time to run each problem')
            print('   -wdir <str>: comma separated list of directories to run')
            print('   -solvers <str>: comma separated list of solvers to run')
            print('   -dump: if True, call the solvers with -dump')
            print('   -clean_log: if True, call the solvers with -clean_log')
            print('   -verbose: if True, print more information')
            print('   -ncpu <int>: number of cores to use')
            print('   -help: print this help message')
            sys.exit(0)
    
    # print parameters
    print('Parameters:')
    print('   tmax ........ %d' % tmax)
    print('   WDIR ........ %s' % wdir)
    print('   solver ...... %s' % solvers)
    print('   dump ........ %s' % dump)
    print('   clean_log ... %s' % clean_log)
    print('   verbose ..... %s' % verbose)
    print('   ncpu ........ %d' % ncpu)
    print('')

    # clean log files
    if clean_log:
        remove_logs(wdir)

    # get all nmr files
    FNMR = get_nmr_files(wdir)
    
    # get command lines
    CMD = get_command_lines(FNMR, tmax, clean_log, dump, solvers, verbose)
    
    # run all command lines in CMD in parallel, but leave one core for the OS
    print('Running %d jobs in parallel' % ncpu)
    with mp.Pool(ncpu) as pool:
        for _ in tqdm.tqdm(pool.imap(os.system, CMD), total=len(CMD)):
            pass
