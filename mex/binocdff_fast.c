/*==========================================================
 * binocdff_fast.c - 
 *
 * Computes the binomial CDF efficiently in single precision
 *
 * The calling syntax is:
 *
 *		outMatrix = binocdff_fast(k,N,p,flag)
 *      
 *      flag    0 - lower tail CDF
 *              1 - upper tail CDF (complementary CDF)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"  /*Include MEX library*/

/* Trick to convert math.h library FP_NAN macro into MEX compatible NaN */
#include <math.h>
#ifdef FP_NAN
#undef FP_NAN
#define FP_NAN mxGetNaN()
#endif

#include "binocdff.h"

/* Reset FP_NAN macro */
#include <math.h>

/*Computational routine prototypes */
void binocdff_fast_k(float *k, float *N, float *p, int flag, float *Z, mwSize n);
void binocdff_fast_N(float *k, float *N, float *p, int flag, float *Z, mwSize n);
void binocdff_fast_p(float *k, float *N, float *p, int flag, float *Z, mwSize n);
void binocdff_fast_kN(float *k, float *N, float *p, int flag, float *Z, mwSize n);
void binocdff_fast_kp(float *k, float *N, float *p, int flag, float *Z, mwSize n);
void binocdff_fast_Np(float *k, float *N, float *p, int flag, float *Z, mwSize n);
void binocdff_fast_full(float *k, float *N, float *p, int flag, float *Z, mwSize n);

/* Gateway/Main function MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
                 
{
    float *k_matrix, *N_matrix, *p_matrix; /* MxN input matrix */
    int flag;                               /* Flag to select tail of distribution */
    const mwSize *dims_k, *dims_N, *dims_p; /* Dimensions of the inputs */
    const mwSize *dims_out;                 /* Dimensions of the output matrix */
    mwSize n_out, n_dim;                    /* Number of dimensions and elements in output matrix */
    float *out_matrix;                     /* Output matrix */
    
    /* check for proper number of arguments */
    if (nrhs<3) {
        mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:nrhs","Three inputs required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:nlhs","One output required.");
    }

    /* make sure the input arguments are type float */
    if (!mxIsSingle(prhs[0])
        || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:notDouble","Input P must be type float.");
    }
    if (!mxIsSingle(prhs[1])
        || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:notDouble","Input a must be type float.");
    }
    if (!mxIsSingle(prhs[2])
        || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:notDouble","Input b must be type float.");
    }

    /* get the number of dimensions of the inputs */
    dims_k= mxGetDimensions(prhs[0]);
    dims_N= mxGetDimensions(prhs[1]);
    dims_p= mxGetDimensions(prhs[2]);

    /* make sure the input arguments are the same size or scalar */
    /* compare k and N */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[1])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[1])) {
            mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_k[d]!=dims_N[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare k and p */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_k[d]!=dims_p[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare N and p */
    if (mxGetNumberOfElements(prhs[1])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[1])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[1]); d++) {
                if (dims_N[d]!=dims_p[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:binocdff_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    
    if ( nrhs!=4 || mxGetNumberOfElements(prhs[3])!=1)
        flag = 0;
    else {
        flag = (int) mxGetPr(prhs[3])[0];
    }

    /* create a pointer to the real data in the inputs */
    k_matrix = (float *)mxGetPr(prhs[0]);
    N_matrix = (float *)mxGetPr(prhs[1]);
    p_matrix = (float *)mxGetPr(prhs[2]);

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
            binocdff_fast_kN(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
        } else if (mxGetNumberOfElements(prhs[2])==1) {
            binocdff_fast_kp(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
        } else {
            binocdff_fast_k(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[1])==1) {
        if (mxGetNumberOfElements(prhs[2])==1) {
            binocdff_fast_Np(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
        } else {
            binocdff_fast_N(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[2])==1) {
        binocdff_fast_p(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
    } else {
        binocdff_fast_full(k_matrix,N_matrix,p_matrix,flag,out_matrix,n_out);
    }
            
}

/*Computational routines */
/* Scalar k */
void binocdff_fast_k(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;    
    for (i=0; i<n; i++) {
        binom_cdff(k[0],N[i],p[i],1-p[i], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;    
    }
}

/* Scalar N */
void binocdff_fast_N(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;    
    for (i=0; i<n; i++) {
        binom_cdff(k[i],N[0],p[i],1-p[i], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;    
    }
}

/* Scalar p */
void binocdff_fast_p(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;    
    for (i=0; i<n; i++) {
        binom_cdff(k[i],N[i],p[0],1-p[0], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;    
    }
}

/* Scalar k & N */
void binocdff_fast_kN(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;
    for (i=0; i<n; i++) {
        binom_cdff(k[0],N[0],p[i],1-p[i], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;    
    }
}

/* Scalar k & p */
void binocdff_fast_kp(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;    
    for (i=0; i<n; i++) {
        binom_cdff(k[0],N[i],p[0],1-p[0], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;
    }
}

/* Scalar N & p */
void binocdff_fast_Np(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;    
    for (i=0; i<n; i++) {
        binom_cdff(k[i],N[0],p[0],1-p[0], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;    
    }
}

/* Matrix N, k & p */
void binocdff_fast_full(float *k, float *N, float *p, int flag, float *Z, mwSize n)
{
    mwSize i;
    float S0, S1;    
    for (i=0; i<n; i++) {
        binom_cdff(k[i],N[i],p[i],1-p[i], &S0, &S1);
        if (flag)
            Z[i] = S1;                        
        else
            Z[i] = S0;    
    }
}
