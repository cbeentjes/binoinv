/*==========================================================
 * erfcxf_fast.c - 
 *
 * Computes the scaled complimentary error function
 *
 * The calling syntax is:
 *
 *		outMatrix = erfcxf_fast(x)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"  /*Include MEX library*/

#include "erfcxf.h"



/*Computational routine prototypes */

/* Gateway/Main function MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
                 
{
    float *X_matrix; /* MxN input matrix */
    const mwSize *dims_X; /* Dimensions of the inputs */
    const mwSize *dims_out;                 /* Dimensions of the output matrix */
    mwSize n_out, n_dim;                    /* Number of dimensions and elements in output matrix */
    float *out_matrix;                     /* Output matrix */
    
    /* check for proper number of arguments */
    if (nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:erfcx_fast:nrhs","One input required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("MyToolbox:erfcx_fast:nlhs","One output required.");
    }

    /* make sure the input arguments are type double */
    if (!mxIsSingle(prhs[0])
        || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:erfcx_fast:notSingle","Input X must be type single.");
    }

    /* get the number of dimensions of the inputs */
    dims_X= mxGetDimensions(prhs[0]);

    /* create a pointer to the real data in the inputs */
    X_matrix = (float *)mxGetData(prhs[0]);

    /* get the dimensions of the output matrix */
    n_dim     = mxGetNumberOfDimensions(prhs[0]);
    dims_out  = mxGetDimensions(prhs[0]);    

    /* create the output matrix */
    plhs[0] = mxCreateNumericArray(n_dim,dims_out,mxSINGLE_CLASS,mxREAL);
    n_out   = mxGetNumberOfElements(plhs[0]);

    /* get a pointer to the real data in the output matrix */
    out_matrix = (float *)mxGetPr(plhs[0]);

    /* call the computational routine */
    mwSize i;
    for (i=0; i<n_out; i++) {
        out_matrix[i] = erfcxf(X_matrix[i]);
    }
            
}