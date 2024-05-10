CXX := mpic++
CC := mpicc

lab-1: lab-1.cpp
	$(CXX) $^ -o lab-1

lab-2: lab-2.cpp
	clang++ -std=c++20 $^ -o lab-2

.PHONY: run-1
run-1: lab-1
	mpirun -np 3 lab-1
	# python3 plot-solution.py

.PHONY: run-2
run-2: lab-2
	time ./lab-2

.PHONY: clean
clean:
	rm -f lab-1 lab-2 *.o *.csv
