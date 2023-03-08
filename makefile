compile:
	rm *.dat
	g++ simulator.cpp -o sim
	g++ -fopenmp parallel.cpp -o par
