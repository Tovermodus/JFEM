#include "com_simonschmidt_SparseSolver.h"
#include "sparseSolver.h"


JNIEXPORT void JNICALL Java_com_simonschmidt_SparseSolver_sayHello
  (JNIEnv* env, jobject thisObject) {
    std::cout << "Hello from  hh C++ !!" << std::endl;
}

JNIEXPORT jdoubleArray JNICALL Java_com_simonschmidt_SparseSolver_iSolve
  (JNIEnv * env, jobject thisObject, jintArray xs_, jintArray ys_, jdoubleArray vals_, jint entries, jint rows, jint cols, jdoubleArray vect_){

    jint *xs = env->GetIntArrayElements(xs_, 0);
    jint *ys = env->GetIntArrayElements(ys_, 0);
    double *vals = env->GetDoubleArrayElements(vals_, 0);
    double *vect = env->GetDoubleArrayElements(vect_, 0);
    char trans = 'N';
    int nrhs = 1;
    int dim = rows;
    double* a = (double*)malloc(dim*dim*sizeof(double));
    for(int i = 0; i < dim*dim; i++)
        a[i] = 0;
    int LDA = dim;
    int LDB = dim;
    int info;

    std::vector<double> b;
    for(int i = 0; i < entries; i++)
    {
        a[xs[i]*dim+ys[i]] += vals[i];
    }
    env->ReleaseIntArrayElements(xs_, xs, 0);
    env->ReleaseIntArrayElements(ys_, ys, 0);
    env->ReleaseDoubleArrayElements(vals_, vals, 0);
    for(int i = 0; i < dim; i++)
    {
        b.push_back(vect[i]);
    }
    env->ReleaseDoubleArrayElements(vect_, vect, 0);
    int ipiv[dim];
    dgetrf_(&dim, &dim, a, &LDA, ipiv, &info);

    if(info != 0)
    {
        free(a);
        return nullptr;
       }
    dgetrs_(&trans, &dim, &nrhs, a, &LDA, ipiv, & *b.begin(), &LDB, &info);
    free(a);
    if(info != 0)
    {
        return nullptr;
       }
    jdoubleArray result = env->NewDoubleArray(dim);
    jdouble ret[dim];
    for(int i = 0; i < dim; i++)
    {
        ret[i] = b[i];
    }
    env->SetDoubleArrayRegion(result,0,dim,ret);
    return result;
}
JNIEXPORT jdoubleArray JNICALL Java_com_simonschmidt_SparseSolver_iInverse
  (JNIEnv * env, jobject thisObject, jintArray xs_, jintArray ys_, jdoubleArray vals_, jint entries, jint rows, jint cols){

    jint *xs = env->GetIntArrayElements(xs_, 0);
    jint *ys = env->GetIntArrayElements(ys_, 0);
    double *vals = env->GetDoubleArrayElements(vals_, 0);
    char trans = 'N';
    int nrhs = 1;
    int dim = rows;
    double* a = (double*)malloc(dim*dim*sizeof(double));
    double* work = (double*)malloc(dim*dim*sizeof(double));
    int lwork = dim*dim;
    for(int i = 0; i < dim*dim; i++)
        a[i] = 0;
    int LDA = dim;
    int LDB = dim;
    int info;

    std::vector<double> b;
    for(int i = 0; i < entries; i++)
    {
        a[ys[i]*dim+xs[i]] += vals[i];
    }
    env->ReleaseIntArrayElements(xs_, xs, 0);
    env->ReleaseIntArrayElements(ys_, ys, 0);
    env->ReleaseDoubleArrayElements(vals_, vals, 0);
    int ipiv[dim];
    dgetrf_(&dim, &dim, a, &LDA, ipiv, &info);
    if(info != 0)
    {
        free(a);
        free(work);
        return nullptr;
       }
    dgetri_(&dim, a, &dim, ipiv, work, &lwork, &info);
    free(work);
    if(info != 0)
    {
        free(a);
        return nullptr;
        }
    jdoubleArray result = env->NewDoubleArray(dim*dim);
    jdouble *ret = (jdouble*)malloc(dim*dim*sizeof(double));
    for(int i = 0; i < dim*dim; i++)
    {
        ret[i] = a[i];
    }
    env->SetDoubleArrayRegion(result,0,dim*dim,ret);
        free(a);
    free(ret);
    return result;
}
int main()
{
return 0;
}
