runAll: 
	python runAllTests.py

prof_bb:
	python -m cProfile -o profiling/bb.prf codes/bb.py -tmax 60 -fnmr profiling/data/4wua.nmr -clean_log

bb.bin: codes/bb.h codes/bb.cpp
	g++ -O3 codes/bb.cpp -o bb.bin

bb.dbg: codes/bb.h codes/bb.cpp
	g++ -g -O0 codes/bb.cpp -o bb.dbg

testAll: codes/bb.h codes/tests.cpp
	g++ -g -O0 codes/tests.cpp -I/home/michael/gitrepos/googletest/googletest/include -o tests.bin
	./tests.bin

clean:
	rm -f DATA_EPSD_00_DMAX_50/*.pkl
	rm -f DATA_EPSD_00_DMAX_50/*.log
	rm -f DATA_EPSD_00_DMAX_60/*.pkl
	rm -f DATA_EPSD_00_DMAX_60/*.log