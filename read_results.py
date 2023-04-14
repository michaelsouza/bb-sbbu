import os
import sys
import pandas as pd

def get_logs(wdir: list):
    # get all log files
    FLOG = []
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
    wdir = 'data/nmr_test'

    # read parameters
    for i, arg in enumerate(sys.argv):
        if arg == '-wdir':
            wdir = sys.argv[i+1] # list of directories to run
        elif arg == '-help':
            print('Usage: python runAll.py [-wdir <str>]')
            print('   -wdir <str>: directory to run')
            print('   -help: print this help message')

    # print parameters
    print('Parameters:')
    print('   wdir = %s' % wdir)
    print('')

    # get all log files in each directory in WDIR
    FLOG = get_logs(wdir)

    # read all log files
    df = []
    for flog in FLOG:
        df.append(read_log(flog))
    # convert to dataframe
    df = pd.DataFrame(df)
    # convert to correct data types
    df['dump'] = df['dump'].astype(bool)
    df['verbose'] = df['verbose'].astype(bool)
    df['clean_log'] = df['clean_log'].astype(bool)
    df['tmax'] = df['tmax'].astype(int)
    df['|V|'] = df['nnodes'].astype(int)
    df['|E|'] = df['lenE'].astype(int)
    df['|S|'] = df['lenS'].astype(int)    
    # drop unnecessary columns
    df.drop(['nnodes', 'lenE', 'lenS'], axis=1, inplace=True)
    
    # sort columns
    df = df[sorted(df.columns)]
    
    # save to csv
    fn = os.path.join(wdir, 'results.csv')
    print('Saving results to %s' % fn)
    df.to_csv(fn, index=False)
    