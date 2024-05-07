CXX := mpic++
CC := mpicc

main.: lab-1.o
	$(CXX) $^ -o main

main: main.o functions.o
	$(CXX) $^ -o main

.PHONY: run
run: main.o
	mpirun -np 3 main
	python3 plot-solution.py

.PHONY: clean
clean:
	rm main *.o *.csv
