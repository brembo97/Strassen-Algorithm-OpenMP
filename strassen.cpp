/**
 * LEAD University
 * Data Science Program
 * BCD-9218: Parallel and Distributed Computing
 * Instructor Diego Jimenez, Eng. (diego.jimenez@ulead.ac.cr)
 * OpenMP parallel Strassen algorithm for matrix multiplication.
 */

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include "timer.h"
#include <vector>
#include "io.h"

void add(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int size)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void sub(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int size)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void multiply(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int size)
{
    int i, j,k;

    for(i=0; i < size; ++i){
        for(j=0; j < size; ++j){
            for(k=0; k < size; ++k){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// TODO: Function implementing Strassen's algorithm
void strassen(int **A, int **B, int **C, int N){

	// TODO: YOUR CODE GOES HERE
	int new_size = N/2;
	vector<int> z(new_size);

    // Declarando submatrices
	vector<vector<int>>
            a11(new_size, z), a12(new_size, z), a21(new_size, z), a22(new_size, z),
            b11(new_size, z), b12(new_size, z), b21(new_size, z), b22(new_size, z),
            c11(new_size, z), c12(new_size, z), c21(new_size, z), c22(new_size, z),
            m1(new_size, z), m2(new_size, z), m3(new_size, z), m4(new_size, z),
            m5(new_size, z), m6(new_size, z), m7(new_size, z),
            aResult(new_size, z), bResult(new_size, z);

	//Particion de matrices
	int i ,j;
	for(i = 0;i < new_size; ++i)
	{
	    for(j = 0; j < new_size;++j){
	        a11[i][j] = A[i][j];
	        a12[i][j] = A[i][j + new_size];
	        a21[i][j] = A[i + new_size][j];
	        a22[i][j] = A[i + new_size][j + new_size];

            b11[i][j] = B[i][j];
            b12[i][j] = B[i][j + new_size];
            b21[i][j] = B[i + new_size][j];
            b22[i][j] = B[i + new_size][j + new_size];
	    }
	}

#pragma omp parallel
    {
        //Calculando matrices intermedias
        add(a11, a22, aResult, new_size);
        add(b11, b22, bResult, new_size);
        multiply(aResult,bResult,m1, new_size);

        add(a21, a22,aResult, new_size);
        multiply(aResult,b11,m2, new_size);

        sub(b12, b22, bResult, new_size);
        multiply(a11, bResult, m3, new_size);

        sub(b21, b11, bResult, new_size);
        multiply(a22, bResult, m4, new_size);

        add(a11, a12, aResult, new_size);
        multiply(aResult, b22, m5, new_size);

        sub(a21, a11, aResult, new_size);
        add(b11, b12, bResult, new_size);
        multiply(aResult, bResult, m6, new_size);

        sub(a12, a22, aResult, new_size);
        add(b21, b22, bResult, new_size);
        multiply(aResult, bResult, m7, new_size);

        //Calculando submatrices de C

        add(m3, m5, c12, new_size);
        add(m2, m4, c21, new_size);

        add(m1, m4, aResult, new_size);
        add(aResult, m7, bResult, new_size);
        sub(bResult, m5, c11, new_size);

        sub(m1, m2, aResult, new_size);
        add(aResult, m3, bResult, new_size);
        add(bResult, m6, c22, new_size);
    }

    //Agrupar las submatrices de C
    for(i = 0; i < new_size; ++i){
        for(j = 0; j < new_size; ++j){
            C[i][j] = c11[i][j];
            C[i][j + new_size] = c12[i][j];
            C[i + new_size][j] = c21[i][j];
            C[i + new_size][j + new_size] = c22[i][j];
        }
    }
}




// Main method
int main(int argc, char* argv[]) {
	int N;
	int **A, **B, **C;
	double elapsedTime;

	// checking parameters
	if (argc != 2 && argc != 4) {
		cout << "Parameters: <N> [<fileA> <fileB>]" << endl;
		return 1;
	}
	N = atoi(argv[1]);

	// allocating matrices
	A = new int*[N];
	B = new int*[N];
	C = new int*[N];
	for (int i=0; i<N; i++){
		A[i] = new int[N];
		B[i] = new int[N];
		C[i] = new int[N];
	}

	// reading files (optional)
	if(argc == 4){
		readMatrixFile(A,N,argv[2]);
		readMatrixFile(B,N,argv[3]);
	}

	if(argc == 2){
	    for(int i = 0; i < N; ++i){
	        for(int j = 0; j < N; ++j){
	            A[i][j] = rand() % 11;
                B[i][j] = rand() % 11;
	        }
	    }
	}

	// starting timer
	timerStart();

	// TODO: YOUR CODE GOES HERE
    strassen(A, B, C, N);

	// testing the results is correct
	if(argc == 4 || argc == 2){
		printMatrix(C,N);
	}

	// stopping timer
	elapsedTime = timerStop();

	cout << "Duration: " << elapsedTime << " seconds" << std::endl;

	// releasing memory
	for (int i=0; i<N; i++) {
		delete [] A[i];
		delete [] B[i];
		delete [] C[i];
	}
	delete [] A;
	delete [] B;
	delete [] C;

	return 0;
}
