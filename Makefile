runAll: 
	python runAll.py

prof_bb:
	python -m cProfile -o profiling/bb.prf codes/bb.py -tmax 60 -fnmr profiling/data/4wua.nmr -clean_log

clean:
	rm -f DATA_EPSD_00_DMAX_50/*.pkl
	rm -f DATA_EPSD_00_DMAX_50/*.log
	rm -f DATA_EPSD_00_DMAX_60/*.pkl
	rm -f DATA_EPSD_00_DMAX_60/*.log