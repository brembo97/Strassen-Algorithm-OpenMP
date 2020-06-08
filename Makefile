all: strassen

strassen:
	g++ -fopenmp strassen.cpp -o strassen
	g++ -w strassen.cpp -o strassen_serial

clean:
	rm -f strassen
