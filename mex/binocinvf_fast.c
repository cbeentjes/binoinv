/*==========================================================
 * binocinvf_fast.c - 
 *
 * Computes the inverse binomial CDF efficiently
 *
 * The calling syntax is:
 *
 *		outMatrix = binocinvf_fast(U,N,p)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"  /*Include MEX library*/
#include "matrix.h"

/* Trick to convert math.h library FP_NAN into MEX compatible NaN */
#include <math.h>
#ifdef FP_NAN
#undef FP_NAN
#define FP_NAN mxGetNaN()
#endif

#include "binoinvf.h"

/* Reset FP_NAN macro */
#include <math.h>

/*Computational routine prototypes */
void binocinvf_fast_U(float *U, float *N, float *p, float *Z, mwSize n);
void binocinvf_fast_N(float *U, float *N, float *p, float *Z, mwSize n);
void binocinvf_fast_p(float *U, float *N, float *p, float *Z, mwSize n);
void binocinvf_fast_UN(float *U, float *N, float *p, float *Z, mwSize n);
void binocinvf_fast_Up(float *U, float *N, float *p, float *Z, mwSize n);
void binocinvf_fast_Np(float *U, float *N, float *p, float *Z, mwSize n);
void binocinvf_fast_full(float *U, float *N, float *p, float *Z, mwSize n);

/* Gateway/Main function MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
                 
{
    float *U_matrix, *N_matrix, *p_matrix; /* MxN input matrix */
    const mwSize *dims_U, *dims_N, *dims_p; /* Dimensions of the inputs */
    const mwSize *dims_out;                 /* Dimensions of the output matrix */
    mwSize n_out, n_dim;                    /* Number of dimensions and elements in output matrix */
    float *out_matrix;                     /* Output matrix */
    
    /* check for proper number of arguments */
    if (nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:nrhs","Three inputs required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:nlhs","One output required.");
    }

    /* make sure the input arguments are type float */
    if (!mxIsSingle(prhs[0])
        || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:notSingle","Input U must be type float.");
    }
    if (!mxIsSingle(prhs[1])
        || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:notSingle","Input N must be type float.");
    }
    if (!mxIsSingle(prhs[2])
        || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:notSingle","Input p must be type float.");
    }

    /* get the number of dimensions of the inputs */
    dims_U= mxGetDimensions(prhs[0]);
    dims_N= mxGetDimensions(prhs[1]);
    dims_p= mxGetDimensions(prhs[2]);

    /* make sure the input arguments are the same size or scalar */
    /* compare U and N */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[1])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[1])) {
            mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_U[d]!=dims_N[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare U and p */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_U[d]!=dims_p[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare N and p */
    if (mxGetNumberOfElements(prhs[1])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[1])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[1]); d++) {
                if (dims_N[d]!=dims_p[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocinvf_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }

    /* create a pointer to the real data in the inputs */
    U_matrix = (float *)mxGetData(prhs[0]);
    N_matrix = (float *)mxGetData(prhs[1]);
    p_matrix = (float *)mxGetData(prhs[2]);

    /* get the dimensions of the output matrix */
    if (mxGetNumberOfElements(prhs[0])!=1) {
        n_dim     = mxGetNumberOfDimensions(prhs[0]);
        dims_out  = mxGetDimensions(prhs[0]);
    } else if (mxGetNumberOfElements(prhs[1])!=1) {
        n_dim     = mxGetNumberOfDimensions(prhs[1]);
        dims_out  = mxGetDimensions(prhs[1]);
    } else {
        n_dim     = mxGetNumberOfDimensions(prhs[2]);
        dims_out  = mxGetDimensions(prhs[2]);
    }
    /* create the output matrix */
    plhs[0] = mxCreateNumericArray(n_dim,dims_out,mxSINGLE_CLASS,mxREAL);
    n_out   = mxGetNumberOfElements(plhs[0]);

    /* get a pointer to the real data in the output matrix */
    out_matrix = (float *)mxGetPr(plhs[0]);

    /* call the computational routine */
    if (mxGetNumberOfElements(prhs[0])==1) {
        if (mxGetNumberOfElements(prhs[1])==1) {
            binocinvf_fast_UN(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        } else if (mxGetNumberOfElements(prhs[2])==1) {
            binocinvf_fast_Up(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        } else {
            binocinvf_fast_U(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[1])==1) {
        if (mxGetNumberOfElements(prhs[2])==1) {
            binocinvf_fast_Np(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        } else {
            binocinvf_fast_N(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[2])==1) {
        binocinvf_fast_p(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
    } else {
        binocinvf_fast_full(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
    }
            
}

/*Computational routines */
/* Scalar U */
void binocinvf_fast_U(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[0],N[i],p[i]);
    }
}

/* Scalar N */
void binocinvf_fast_N(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[i],N[0],p[i]);
    }
}

/* Scalar p */
void binocinvf_fast_p(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[i],N[i],p[0]);
    }
}

/* Scalar U & N */
void binocinvf_fast_UN(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[0],N[0],p[i]);
    }
}

/* Scalar U & p */
void binocinvf_fast_Up(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[0],N[i],p[0]);
    }
}

/* Scalar N & p */
void binocinvf_fast_Np(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[i],N[0],p[0]);
    }
}

/* Matrix N, U & p */
void binocinvf_fast_full(float *U, float *N, float *p, float *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinvf(U[i],N[i],p[i]);
    }
}
