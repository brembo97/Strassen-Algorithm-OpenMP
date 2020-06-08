/**
 * LEAD University
 * Data Science Program
 * BCD-9218: Parallel and Distributed Computing
 * Instructor Diego Jimenez, Eng. (diego.jimenez@ulead.ac.cr)
 * Input/output operations for matrices.
 */

#include <iostream>
#include <fstream>

using namespace std;

// Reads a matrix from text file
int readArrayFile(int *A, int N, char *fileName){
	ifstream file(fileName);
	if(file.is_open()){

		// reading array values
		for(int i=0; i<N; i++){
			file >> A[i];
		}

		// closing file
		file.close();

	} else {
		cout << "Error opening file: " << fileName << endl;
		return 1;
	}
	return 0;
}

// Prints array to standard output
void printArray(int *A, int N){
	for(int i=0; i<N; i++){
		cout << A[i] << "\t";
	}
	cout << endl;
}

// Reads a matrix from text file
int readMatrixFile(int **matrix, int N, char *fileName){
	ifstream file(fileName);
	if(file.is_open()){

		// reading matrix values
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				file >> matrix[i][j];
			}
		}

		// closing file
		file.close();

	} else {
		cout << "Error opening file: " << fileName << endl;
		return 1;
	}
	return 0;
}

// Prints matrix to standard output
void printMatrix(int **matrix, int N){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}
