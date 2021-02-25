/*==========================================================
 * binocinv_fast.c -
 *
 * Computes the inverse binomial complementary CDF efficiently
 *
 * The calling syntax is:
 *
 *		outMatrix = binocinv_fast(U,N,p)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"  /*Include MEX library*/

/* Trick to convert math.h library FP_NAN into MEX compatible NaN */
#include <math.h>
#ifdef FP_NAN
#undef FP_NAN
#define FP_NAN mxGetNaN()
#endif

#include "binoinv.h"

/* Reset FP_NAN macro */
#include <math.h>

/*Computational routine prototypes */
void binocinv_fast_U(double *U, double *N, double *p, double *Z, mwSize n);
void binocinv_fast_N(double *U, double *N, double *p, double *Z, mwSize n);
void binocinv_fast_p(double *U, double *N, double *p, double *Z, mwSize n);
void binocinv_fast_UN(double *U, double *N, double *p, double *Z, mwSize n);
void binocinv_fast_Up(double *U, double *N, double *p, double *Z, mwSize n);
void binocinv_fast_Np(double *U, double *N, double *p, double *Z, mwSize n);
void binocinv_fast_full(double *U, double *N, double *p, double *Z, mwSize n);

/* Gateway/Main function MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])

{
    double *U_matrix, *N_matrix, *p_matrix; /* MxN input matrix */
    const mwSize *dims_U, *dims_N, *dims_p; /* Dimensions of the inputs */
    const mwSize *dims_out;                 /* Dimensions of the output matrix */
    mwSize n_out, n_dim;                    /* Number of dimensions and elements in output matrix */
    double *out_matrix;                     /* Output matrix */

    /* check for proper number of arguments */
    if (nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:nrhs","Three inputs required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:nlhs","One output required.");
    }

    /* make sure the input arguments are type double */
    if (!mxIsDouble(prhs[0])
        || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:notDouble","Input P must be type double.");
    }
    if (!mxIsDouble(prhs[1])
        || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:notDouble","Input a must be type double.");
    }
    if (!mxIsDouble(prhs[2])
        || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:notDouble","Input b must be type double.");
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
            mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_U[d]!=dims_N[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare U and p */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_U[d]!=dims_p[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare N and p */
    if (mxGetNumberOfElements(prhs[1])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[1])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[1]); d++) {
                if (dims_N[d]!=dims_p[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocinv_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }

    /* create a pointer to the real data in the inputs */
    U_matrix = mxGetPr(prhs[0]);
    N_matrix = mxGetPr(prhs[1]);
    p_matrix = mxGetPr(prhs[2]);

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
    plhs[0] = mxCreateNumericArray(n_dim,dims_out,mxDOUBLE_CLASS,mxREAL);
    n_out   = mxGetNumberOfElements(plhs[0]);

    /* get a pointer to the real data in the output matrix */
    out_matrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    if (mxGetNumberOfElements(prhs[0])==1) {
        if (mxGetNumberOfElements(prhs[1])==1) {
            binocinv_fast_UN(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        } else if (mxGetNumberOfElements(prhs[2])==1) {
            binocinv_fast_Up(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        } else {
            binocinv_fast_U(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[1])==1) {
        if (mxGetNumberOfElements(prhs[2])==1) {
            binocinv_fast_Np(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        } else {
            binocinv_fast_N(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[2])==1) {
        binocinv_fast_p(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
    } else {
        binocinv_fast_full(U_matrix,N_matrix,p_matrix,out_matrix,n_out);
    }

}

/*Computational routines */
/* Scalar U */
void binocinv_fast_U(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[0],N[i],p[i]);
    }
}

/* Scalar N */
void binocinv_fast_N(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[i],N[0],p[i]);
    }
}

/* Scalar p */
void binocinv_fast_p(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[i],N[i],p[0]);
    }
}

/* Scalar U & N */
void binocinv_fast_UN(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[0],N[0],p[i]);
    }
}

/* Scalar U & p */
void binocinv_fast_Up(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[0],N[i],p[0]);
    }
}

/* Scalar N & p */
void binocinv_fast_Np(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[i],N[0],p[0]);
    }
}

/* Matrix N, U & p */
void binocinv_fast_full(double *U, double *N, double *p, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = binocinv(U[i],N[i],p[i]);
    }
}