import os
import sys
import tqdm
import pandas as pd
import re

def get_logs(WDIR: list):
    # get all log files
    FLOG = []
    for wdir in WDIR:
        print('Getting log files from wdir=%s' % wdir)
        for fn in os.listdir(wdir):
            fn = os.path.join(wdir, fn)
            # add log files
            if fn.endswith('.log'):
                FLOG.append(fn)
    return FLOG    


def read_log(flog):
    # check if file exists
    if not os.path.isfile(flog):
        raise FileNotFoundError('File %s not found' % flog)
    
    d = {}
    # read file
    with open(flog, 'r') as fd:
        for line in fd:
            if line.startswith('Run time'):
                text = line
                break
            # parse text by getting the first and the last words
            text = line.strip().split()
            field = text[1]
            value = text[-1]
            # convert to float if possible
            try:
                value = float(value)
            except ValueError:
                pass
            d[field] = value
    return d


if __name__ == "__main__":
    # set default parameters
    WDIR = ['DATA_TEST'] # list of directories to run

    # read parameters
    for i, arg in enumerate(sys.argv):
        if arg == '-wdir':
            WDIR = sys.argv[i+1].split(',') # list of directories to run
        elif arg == '-help':
            print('Usage: python runAll.py [-wdir <str>]')
            print('   -wdir <str>: comma separated list of directories to run')
            print('   -help: print this help message')

    # print parameters
    print('Parameters:')
    print('   WDIR = %s' % WDIR)
    print('')

    # get all log files in each directory in WDIR
    FLOG = get_logs(WDIR)

    # read all log files
    df = []
    for flog in tqdm.tqdm(FLOG):
        df.append(read_log(flog))
    # convert to dataframe
    df = pd.DataFrame(df)
    # convert to correct data types
    df['verbose'] = df['tmax'].astype(bool)
    df['clean_log'] = df['clean_log'].astype(bool)
    df['nnodes'] = df['nnodes'].astype(int)
    df['lenE'] = df['lenE'].astype(int)
    df['lenS'] = df['lenS'].astype(int)

    # save to csv
    fn = 'results.csv'
    print('Saving results to %s' % fn)
    df.to_csv(fn, index=False)
    