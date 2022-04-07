#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetri_(int* dim, double* a, int* lda, int* ipiv, double* work,int *lwork, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );