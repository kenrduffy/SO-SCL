
/*
* Install: mex -O polarTrans.c
*/

/* C Functions */
#include <mex.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>

/* if power of 2 */
uint8_t Is2Power(uint64_t nNum){
    return nNum > 0 ? ((nNum & (~nNum + 1)) == nNum ? true : false) : false;
}

/* First Odd Then Even (1 3 5 7 ... 2 4 6 8 ...) */
void OddEven(uint8_t *w, uint64_t n){
    uint8_t *temp = (uint8_t *) calloc (n, sizeof(uint8_t));
    for(size_t i = 0; i < n/2; i++){
        temp[i] = w[2*i];
        temp[n/2 + i] = w[2*i + 1];
    }
    for(size_t i = 0; i < n; i++)
        w[i] = temp[i];
    free(temp);
}

/* Bit-Reversal permutation */
void BitReversal(uint8_t *input, uint64_t n){
    if(n <= 2)
        return;
    OddEven(input, n);
    BitReversal(input, n/2);
    BitReversal(&input[n/2], n/2);
}

/* Polar-Transform */
void polarTrans(uint8_t *u, uint8_t *c, uint64_t n){
    uint32_t m = (uint32_t)log2(n);
    for(size_t i = 0; i < n; i++)
        c[i] = u[i];
    for(size_t lambda = 0; lambda < m; lambda++){
        uint64_t sec = (uint64_t) 1 << lambda;
        uint64_t ele = (uint64_t) 1 << (m-lambda-1);
        for(size_t j = 0; j < sec; j++){
            for(size_t i = 0; i < ele; i++)
                c[j*ele*2+i] = c[j*ele*2+i] ^ c[j*ele*2+i+ele];
        }
    }
}

/* Mexfunction Interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    uint8_t ifBitReversal = 0; /* default Bit-Reversal or not (0 for not) */
    /* Check for proper number of arguments */
    if( nrhs < 1 || nrhs > 2 )
        mexErrMsgTxt("Must have either 1 or 2 input arguments.");
    if( nlhs > 1 )
        mexErrMsgTxt("Too many output arguments.");
    if( mxGetN(prhs[0]) != 1 || !mxIsClass(prhs[0], "double"))
        mexErrMsgTxt("First Input must be a column-vector of type double.");
    if( nrhs == 2 ){
        if( mxGetNumberOfElements(prhs[1]) != 1 || !mxIsClass(prhs[1], "double") )
            mexErrMsgTxt("Second Input must be a scalar of type double.");
        else
            ifBitReversal = (uint8_t)*mxGetPr(prhs[1]);
    }
    /* input */
    uint64_t n = mxGetNumberOfElements(prhs[0]);
    if( !Is2Power(n) )
        mexErrMsgTxt("Input Size must be a power of 2.");
    double *u_p = mxGetPr(prhs[0]);
    /* cast the input into a vector of integers */
    uint8_t *u_int = calloc( n, sizeof( uint8_t ) );
		for (size_t i = 0; i < n; i++)
			u_int[i] = (uint8_t) u_p[i];
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *c_p = mxGetPr(plhs[0]);
    uint8_t *c_int = calloc( n, sizeof( uint8_t ) );
    /* use C functions in mexfunction */
    polarTrans(u_int, c_int, n);
    if( ifBitReversal != 0 )
        BitReversal(c_int, n);
    /* cast to output */
    for (size_t i = 0; i < n; i++) 			
		c_p[i] = (double) c_int[i];
    /* Clean up memory */
    free(u_int);
    free(c_int);
}