/*==========================================================
 * bratio_wrapper.c -
 *
 * Computes the incomplete beta function efficiently via R-math library
 *
 * The calling syntax is:
 *
 *		outMatrix = bratio_wrapper(x,z,w,flag)
 *
 *      flag    1 - lower tail
 *              0 - upper tail
 *
 * This is a MEX-file for MATLAB.
 *
 * Compile as mex bratio_wrapper.c -lRmath and make sure R standalone math
 * library is compiled and found by linker by setting LD_LIBRARY_PATH if
 * necessary
 *
 *========================================================*/

#include "mex.h"  /*Include MEX library*/

#define MATHLIB_STANDALONE 1
#include <Rmath.h>


/*Computational routine prototypes */
void betainc_fast_x(double *x, double *z, double *w, int flag, double *Z, mwSize n);
void betainc_fast_z(double *x, double *z, double *w, int flag, double *Z, mwSize n);
void betainc_fast_w(double *x, double *z, double *w, int flag, double *Z, mwSize n);
void betainc_fast_xz(double *x, double *z, double *w, int flag, double *Z, mwSize n);
void betainc_fast_xw(double *x, double *z, double *w, int flag, double *Z, mwSize n);
void betainc_fast_zw(double *x, double *z, double *w, int flag, double *Z, mwSize n);
void betainc_fast_full(double *x, double *z, double *w, int flag, double *Z, mwSize n);

/* Gateway/Main function MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])

{
    double *x_matrix, *z_matrix, *w_matrix; /* MxN input matrix */
    int flag;                               /* Flag to select tail of distribution */
    const mwSize *dims_x, *dims_z, *dims_w; /* Dimensions of the inputs */
    const mwSize *dims_out;                 /* Dimensions of the output matrix */
    mwSize n_out, n_dim;                    /* Number of dimensions and elements in output matrix */
    double *out_matrix;                     /* Output matrix */

//     printf("%s\n",R_VERSION_STRING);

    /* check for proper number of arguments */
    if (nrhs<3) {
        mexErrMsgIdAndTxt("MyToolbox:betainc_fast:nrhs","Three inputs required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("MyToolbox:betainc_fast:nlhs","One output required.");
    }

    /* make sure the input arguments are type double */
    if (!mxIsDouble(prhs[0])
        || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:betainc_fast:notDouble","Input P must be type double.");
    }
    if (!mxIsDouble(prhs[1])
        || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:betainc_fast:notDouble","Input a must be type double.");
    }
    if (!mxIsDouble(prhs[2])
        || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:betainc_fast:notDouble","Input b must be type double.");
    }

    /* get the number of dimensions of the inputs */
    dims_x= mxGetDimensions(prhs[0]);
    dims_z= mxGetDimensions(prhs[1]);
    dims_w= mxGetDimensions(prhs[2]);

    /* make sure the input arguments are the same size or scalar */
    /* compare k and N */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[1])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[1])) {
            mexErrMsgIdAndTxt("MyToolbox:betainc_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_x[d]!=dims_z[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:betainc_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare k and p */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:betainc_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
                if (dims_x[d]!=dims_w[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:betainc_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }
    /* compare N and p */
    if (mxGetNumberOfElements(prhs[1])!=1
        && mxGetNumberOfElements(prhs[2])!=1) {
        if (mxGetNumberOfDimensions(prhs[1])!=mxGetNumberOfDimensions(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:betainc_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.");
        } else {
            for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[1]); d++) {
                if (dims_z[d]!=dims_w[d]) {
                    mexErrMsgIdAndTxt("MyToolbox:betainc_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
                }
            }
        }
    }

    if ( nrhs!=4 || mxGetNumberOfElements(prhs[3])!=1)
        flag = 1;
    else {
        flag = (int) mxGetPr(prhs[3])[0];
    }

    /* create a pointer to the real data in the inputs */
    x_matrix = mxGetPr(prhs[0]);
    z_matrix = mxGetPr(prhs[1]);
    w_matrix = mxGetPr(prhs[2]);

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
            betainc_fast_xz(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
        } else if (mxGetNumberOfElements(prhs[2])==1) {
            betainc_fast_xw(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
        } else {
            betainc_fast_x(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[1])==1) {
        if (mxGetNumberOfElements(prhs[2])==1) {
            betainc_fast_zw(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
        } else {
            betainc_fast_z(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
        }
    } else if (mxGetNumberOfElements(prhs[2])==1) {
        betainc_fast_w(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
    } else {
        betainc_fast_full(x_matrix,z_matrix,w_matrix,flag,out_matrix,n_out);
    }

}

/*Computational routines */
void betainc_fast_x(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[0],z[i],w[i],flag,0);
    }
}

void betainc_fast_z(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[i],z[0],w[i],flag,0);
    }
}

void betainc_fast_w(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[i],z[i],w[0],flag,0);
    }
}

void betainc_fast_xz(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[0],z[0],w[i],flag,0);
    }
}

void betainc_fast_xw(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[0],z[i],w[0],flag,0);
    }
}

void betainc_fast_zw(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[i],z[0],w[0],flag,0);
    }
}

void betainc_fast_full(double *x, double *z, double *w, int flag, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = pbeta(x[i],z[i],w[i],flag,0);
    }
}
