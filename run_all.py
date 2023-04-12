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


def get_command_lines(FNMR: str, tmax: int):
    CMD = [] # list of arguments
    for fnmr in FNMR:
        cmd = './build/bb.bin -tmax %d -fnmr %s -clean_log' % (tmax, fnmr)
        # arg = 'python codes/bb.py -tmax %d -fnmr %s -clean_log' % (tmax, fnmr)
        CMD.append(cmd)

    # write all cmd to file
    with open('runAll.cmd', 'w') as fd:
        fd.write('\n'.join(CMD))

    return CMD


if __name__ == "__main__":
    # set default parameters
    tmax = 3600  # 1 hour
    WDIR = ['DATA_TEST'] # list of directories to run

    # read parameters
    for i, arg in enumerate(sys.argv):
        if arg == '-tmax':
            tmax = int(sys.argv[i+1])
        elif arg == '-wdir':
            WDIR = sys.argv[i+1].split(',')
        elif arg == '-help':
            print('Usage: python runAll.py [-tmax <int>] [-wdir <str>]')
            print('   -tmax <int>: maximum time to run each problem')
            print('   -wdir <str>: comma separated list of directories to run')
            print('   -help: print this help message')
            sys.exit(0)
    
    # print parameters
    print('Parameters:')
    print('   tmax = %d' % tmax)
    print('   WDIR = %s' % WDIR)
    print('')

    # clean log files
    remove_logs(WDIR)

    # get all nmr files
    FNMR = get_nmr_files(WDIR)
    
    # get command lines
    CMD = get_command_lines(FNMR, tmax)
    
    # run all command lines in CMD in parallel, but leave one core for the OS
    ncpu = mp.cpu_count() - 1
    print('Running %d jobs in parallel' % ncpu)
    with mp.Pool(ncpu) as pool:
        for _ in tqdm.tqdm(pool.imap(os.system, CMD), total=len(CMD)):
            pass
