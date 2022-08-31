# To visualize the results use snakeviz
# snakeviz profiling/bb.stats

import os
from codes.bb import *
import cProfile, pstats


if __name__ == "__main__":
    print('Get git log')
    os.system('git log -1 > profiling/git_commit.txt')
    
    print('Call profiler')
    profiler = cProfile.Profile()
    profiler.enable()
    call_solvers('-tmax', 10, '-fnmr', 'profiling/data/4wua.nmr', '-clean_log')
    profiler.disable()
    
    fn = 'profiling/bb.stats'
    print('Save profile stats on', fn)
    stats = pstats.Stats(profiler)
    stats.dump_stats(fn)
    print('To visualize stats, use snakeviz')
    print('$snakeviz ', fn)
    