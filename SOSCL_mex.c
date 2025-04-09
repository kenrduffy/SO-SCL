/*
* All code is subject to license:
* GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf
*/

/*
* References to be cited when using this code.
* P. Yuan, M. Medard, K. Galligan, K. R. Duffy. "Soft-output (SO) GRAND and 
* Iterative Decoding to Outperform LDPC Codes". IEEE Trans. Wireless Commun., 
* 2025.
* P. Yuan, K. R. Duffy & M. Médard. "Near-optimal generalized decoding of 
* Polar-like codes", IEEE ISIT, 2024.
* P. Yuan, K. R. Duffy & M. Médard. "Soft-output successive cancellation 
* list decoding", IEEE Transactions on Information Theory, 71 (2), 
* 1007–1017, 2025.
*/

/*
* Install: mex -O SOSCL_mex.c
*/

/* C Functions */

#include <mex.h>
#include <math.h>
#include <stdint.h>
#define Inf 0x7fffffff

/* u-Node */
struct u_node {
    uint8_t frz;
    uint64_t degree;
    uint64_t *index;
};

/* if power of 2 */
uint8_t Is2Power(uint64_t nNum){
    return nNum > 0 ? ((nNum & (~nNum + 1)) == nNum ? true : false) : false;
}

/* Sign Function */
int8_t sign(double x){
    if(x >= 0)
        return 1;
    else
        return -1;
}

/* Find Min */
double findMin (double a, double b){
    return !(b < a) ? a : b;
}

/* First Odd Then Even (1 3 5 7 ... 2 4 6 8 ...) */
void OddEven(size_t *w, size_t n){
    size_t *temp = (size_t *) calloc (n, sizeof(size_t));
    for(size_t i = 0; i < n/2; i++){
        temp[i] = w[2*i];
        temp[n/2 + i] = w[2*i + 1];
    }
    for(size_t i = 0; i < n; i++)
        w[i] = temp[i];
    free(temp);
}

/* Bit-Reversal permutation */
void BitReversal(size_t *input, size_t n){
    if(n <= 2)
        return;
    OddEven(input, n);
    BitReversal(input, n/2);
    BitReversal(&input[n/2], n/2);
}

/* Jacobian logarithm (approximated) */
double JacLog(double x){
    if(x > 100.0)
        return x;
    else if(x < -100.0)
        return 0.0;
    else
        return log(1+exp(x));
}

/* log of Sum PM */
double logSumPM(double a, double b){
    return findMin(a, b) - JacLog(-fabs(a-b));
}
    
/* Check-Note Update Function */
double CN(double a, double b){
//     return findMin(fabs(a), fabs(b)) * sign(a) * sign(b); /* min-sum approximation */
    // return 2 * atanh( tanh(a/2) * tanh(b/2) ); /* tanh rule */
    return findMin(fabs(a), fabs(b)) * sign(a) * sign(b) + JacLog(-fabs(a+b)) - JacLog(-fabs(a-b)); /* tanh rule (approximated) */
}

/* Path-Metric Update Function */
double updatePM(double PM_in, double L, uint8_t u){
    return PM_in + JacLog(-(1-2*u)*L);
//     if ( 2 * u == 1 - sign(L) ) /* approximation */
//         return PM_in;
//     else
//         return PM_in + fabs(L);
}

/* Quick Sort */
void QuickSort (double *a, uint64_t n) {
    uint64_t i, j;
    double t, p;
    if (n < 2)
        return;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    QuickSort(a, i);
    QuickSort(a + i, n - i);
}

/* Get the median of the EVEN length array */
double median(double *array, uint64_t n){
    double *x = (double *) calloc (n, sizeof(double));
    for (size_t i = 0; i < n; i++)
        x[i] = array[i];
    QuickSort(x, n);
    double m = ( x[n/2-1] + x[n/2] )/2;
    free(x);
    return m;
}

/* Recursively Update Matrix B (Bits) */
void recursivelyUpdateB(uint32_t lambda, uint64_t phi, uint64_t n, uint32_t m,  uint8_t *C){
    uint64_t psi = (uint64_t)(phi/2); // int psi = floor(phi/2);
    for(size_t beta = 0; beta < (uint64_t)pow(2.0, m-lambda); beta++){
        C[(lambda-1)*n*2+2*beta*2+(psi%2)] = C[lambda*n*2+beta*2+((phi-1)%2)] ^ C[lambda*n*2+beta*2+(phi%2)];
        C[(lambda-1)*n*2+(2*beta+1)*2+(psi%2)] = C[lambda*n*2+beta*2+(phi%2)];
    }
    if (psi%2 == 1)
        recursivelyUpdateB(lambda-1, psi, n, m, C);
}

/* Recursively Update Matrix P (LLRs) */
void recursivelyCalcP(uint32_t lambda, uint64_t phi, uint64_t n, uint32_t m, double *P, uint8_t *C){
    if (lambda == 0)
        return;
    uint64_t psi = (uint64_t)(phi/2); // int psi = floor(phi/2);
    if (phi%2 == 0)
        recursivelyCalcP(lambda-1, psi, n, m, P, C);
    for(size_t beta = 0; beta < (uint64_t)pow(2.0,m-lambda); beta++){
        if (phi%2 == 0)
            P[lambda*n+beta] = CN(P[(lambda-1)*n+2*beta],P[(lambda-1)*n+2*beta+1]);
        else
            P[lambda*n+beta] = (1-2* C[lambda*n*2+beta*2+(phi-1)%2] )*P[(lambda-1)*n+2*beta]+P[(lambda-1)*n+2*beta+1];
    }
}

/* Main Function */
void dec(uint8_t *uhat, uint8_t *chat, double *PM, double *llr, struct u_node *u_nodes, uint64_t n, uint64_t L, double *log_minus_pNL, uint8_t ifBR){
    uint32_t m = (uint32_t)log2(n);
    uint64_t NumOfCan = 1;
    /* For Lazy Copy */
    uint64_t beta;
    uint32_t t;
    size_t *stack = (size_t *) calloc (L, sizeof(size_t)); 
    /* Matrices */
    double *P = (double *) calloc (L*(m+1)*n, sizeof(double)); // P: L*(m+1)*n
    uint8_t *C = (uint8_t *) calloc (L*(m+1)*n*2, sizeof(uint8_t)); //C: L*(m+1)*n*2
    double *llrPrime = (double *) calloc (L, sizeof(double)); // llrPrime: 1*L
    double *PM_temp = (double *) calloc (2*L, sizeof(double)); // PM_temp: 2*L
    size_t *brp = (size_t *)calloc(n, sizeof(size_t));
    /* for pY */
    log_minus_pNL[0] = Inf;
    uint64_t nFutureFrz = 0;
    for(size_t i = 0; i < n; i++){
        if( u_nodes[i].frz != 0 )
            nFutureFrz++;
    }
    /* Initialize */
    for(size_t i = 0; i < L; i++)
        PM[i] = 0.0;
    if( ifBR == 0 ){
        /* if not Bit-reversal version */
        for(size_t i = 0; i < n; i++)
            brp[i] = i;
        BitReversal(brp, n); /* change the permutation */
        for(size_t i = 0; i < n; i++)
            P[i] = llr[brp[i]];
    }
    else{
        /* if Bit-reversal version */
        for(size_t i = 0; i < n; i++)
            P[i] = llr[i];
    }
    /* Main Loop */
    for(size_t phi = 0; phi < n; phi++){
        for(size_t l = 0; l < NumOfCan; l++){
            recursivelyCalcP(m, phi, n, m, &P[l*(m+1)*n], &C[l*(m+1)*n*2]);
            llrPrime[l] = P[l*(m+1)*n+m*n]; // P(l,m,0)
        }
        if( u_nodes[phi].frz == 1 ){
            /* Frozen */
            for(size_t l = 0; l < NumOfCan; l++){
                C[l*(m+1)*n*2 + m*n*2 + (phi%2)] = 0; // C(l,m,0,phi%2)
                uhat[l*n + phi] = 0; // uhat(l,phi)
                PM[l] = updatePM(PM[l], llrPrime[l], 0);
            }
        }
        else if( u_nodes[phi].frz == 2 ){
            /* Dynamical Frozen */
            for(size_t l = 0; l < NumOfCan; l++){
                uint8_t dfrzbit = 0;
                for(size_t i = 0; i < u_nodes[phi].degree; i++)
                    dfrzbit ^= uhat[l*n + u_nodes[phi].index[i]];
                C[l*(m+1)*n*2 + m*n*2 + (phi%2)] = dfrzbit; // C(l,m,0,phi%2)
                uhat[l*n + phi] = dfrzbit; // uhat(l,phi)
                PM[l] = updatePM(PM[l], llrPrime[l], dfrzbit);
            }
        }
        else{
            /* Unfrozen */
            for(size_t l = 0; l < NumOfCan; l++){
                /* calculate all PMs */
                for(uint8_t u_temp = 0; u_temp < 2; u_temp++)
                    PM_temp[u_temp*L+l] = updatePM(PM[l], llrPrime[l], u_temp); // PM_temp(u_temp,l)
            }
            if(NumOfCan < L){
                /* List is not full */
                for(size_t l = 0; l < NumOfCan; l++){
                    PM[NumOfCan + l] = PM[l];
                    /* begin copy P */
                    for(size_t j = 0; j < m+1; j++){
                        beta = j*n + n/(uint64_t)(pow(2.0,j));
                        for(size_t i = j*n; i < beta; i++)
                            P[(NumOfCan + l)*(m+1)*n + i] = P[l*(m+1)*n + i];
                    }
                    /* begin copy C */
                    if(phi <= 1)
                        t = 1;
                    else
                        t = (uint32_t)floor(log2(phi));
                    for(size_t j = m - t; j < m + 1; j++){
                        beta = j*n + n/ (uint64_t)(pow(2.0,j));
                        for(size_t i = j*n; i < beta; i++)
                            C[(NumOfCan + l)*(m+1)*n*2 + 2*i ] = C[l*(m+1)*n*2 + 2*i];
                    }
                    t = t - 1;
                    for(size_t j = m-t; j < m+1; j++){
                        beta = j*n + n/ (uint64_t)(pow(2.0,j));
                        for(size_t i = j*n; i < beta; i++)
                            C[(NumOfCan + l)*(m+1)*n*2 + 2*i + 1] = C[l*(m+1)*n*2 + 2*i + 1];
                    }
                    /* begin copy uhat */
                    for(size_t i = 0; i < phi; i++)
                        uhat[(NumOfCan + l)*n + i] = uhat[l*n + i]; // uhat(NumOfCan+l,i)=uhat(l,i)
                    PM[l] = PM_temp[l]; // PM_temp(0,l)
                    C[l*(m+1)*n*2 + m*n*2 + (phi%2)] = 0; // C(l,m,0,phi%2)
                    uhat[l*n + phi] = 0; // uhat(l,phi);
                    PM[NumOfCan + l] = PM_temp[L + l]; // PM_temp(1,l)
                    C[(NumOfCan + l)*(m+1)*n*2 + m*n*2 + (phi%2)] = 1; // C(NumOfCan+l,m,0,phi%2)
                    uhat[(NumOfCan + l)*n + phi] = 1; // uhat(NumOfCan+l,phi)
                }
                NumOfCan = NumOfCan * 2;
            }
            else{
                /* List is full */
                double tao = median(PM_temp, 2*L);
                uint64_t a = 0, b = 0;
                for(size_t l = 0; l < NumOfCan; l++){
                    if (PM_temp[l] > tao && PM_temp[L+l] > tao)
                        stack[L/2 + b++] = l;
                    else if (PM_temp[l] < tao && PM_temp[L+l] < tao)
                        stack[a++] = l;
                    else if (PM_temp[l] < tao && PM_temp[L+l] > tao){
                        C[l*(m+1)*n*2 + m*n*2 + (phi%2)] = 0; // C(l,m,0,phi%2)
                        uhat[l*n + phi] = 0; // uhat(l,phi)
                        PM[l] = PM_temp[l]; // PM_temp(0,l)
                    }
                    else{
                        C[l*(m+1)*n*2 + m*n*2 + (phi%2)] = 1; // C(l,m,0,phi%2)
                        uhat[l*n + phi] = 1; // uhat(l,phi)
                        PM[l] = PM_temp[L+l]; // PM_temp(0,l)
                    }
                }
                /* for pY: approximate sum of the prob. of valid codewords in the subtree */
                for(size_t l = 0; l < 2*NumOfCan; l++){
                    if (PM_temp[l] > tao)
                        log_minus_pNL[0] = logSumPM( log_minus_pNL[0], PM_temp[l] + (double)nFutureFrz*log(2) );
                }
                //if(a != b)
                    //mexErrMsgTxt(" Stack Error ");                
                if(a != 0){
                    /* lazy copy */
                    for (size_t j = 0; j < a; j++){
                        size_t i_act = stack[j];
                        size_t i_inact = stack[L/2 + j];
                        /* copy PM */
                        PM[i_inact] = PM[i_act];
                        /* begin copy P */
                        for(size_t j = 1; j < m+1; j++){
                            uint64_t beta = j*n + n/ (uint64_t)(pow(2.0,j));
                            for(size_t i = j*n; i < beta; i++)
                                P[i_inact*(m+1)*n + i] = P[i_act*(m+1)*n + i];
                        }
                        /* begin copy C */
                        if(phi <= 1)
                            t = 1;
                        else
                            t = (uint32_t)floor(log2(phi));
                        for(size_t j = m-t; j < m+1; j++){
                            beta = j*n + n/ (uint64_t)(pow(2.0,j));
                            for(size_t i = j*n; i < beta; i++)
                                C[i_inact*(m+1)*n*2 + 2*i ] = C[i_act*(m+1)*n*2 + 2*i];
                        }
                        t = t - 1;
                        for(size_t j = m-t; j < m+1; j++){
                            beta = j*n+n/(uint64_t)(pow(2.0,j));
                            for(size_t i = j*n; i < beta; i++)
                                C[i_inact*(m+1)*n*2 + 2*i + 1] = C[i_act*(m+1)*n*2 + 2*i + 1];
                        }
                        /* begin copy uhat */
                        for(size_t i = 0; i < phi; i++)
                            uhat[i_inact*n + i] = uhat[i_act*n + i]; // uhat(i_inact,i)=uhat(i_act,i)
                        /* continue both path */
                        PM[i_inact] = PM_temp[i_act];
                        C[i_inact*(m+1)*n*2 + m*n*2 + (phi%2)] = 0;
                        uhat[i_inact*n + phi] = 0; // uhat(i_inact,phi)
                        PM[i_act] = PM_temp[L + i_act];
                        C[i_act*(m+1)*n*2 + m*n*2 + (phi%2)] = 1;
                        uhat[i_act*n + phi] = 1; // uhat(i_act,phi)
                    }
                } //end if (!act.empty())
            } // end if List is full
        } //end if info bits
        if(phi%2 == 1){
            for(size_t l = 0; l != NumOfCan; l++)
                recursivelyUpdateB(m, phi, n, m, &C[l*(m+1)*n*2]);
        }
        if( u_nodes[phi].frz != 0 )
            nFutureFrz--;
    } // end for phi
    /* output */
    if( ifBR == 0 ){
        for(size_t l = 0; l < L; l++){
            for(size_t i = 0; i < n; i++)
                chat[l*n + i] = C[l*(m+1)*n*2 + brp[i]*2];
        }
    }
    else{
        for(size_t l = 0; l < L; l++){
            for(size_t i = 0; i < n; i++)
                chat[l*n + i] = C[l*(m+1)*n*2 + i*2];
        }
    }
    /* Release */
    free(brp);
    free(stack);
    free(C);
    free(P);
    free(llrPrime);
    free(PM_temp);
}

/* Mexfunction Interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    uint8_t ifBitReversal = 0; /* default Bit-Reversal or not (0 for not) */
    /* Check for proper number of arguments */
    if( nrhs < 4 || nrhs > 5 )
        mexErrMsgTxt("Must have either 5 or 6 input arguments.");
    if( nlhs > 4 )
        mexErrMsgTxt("Too many output arguments.");
    if( mxGetN(prhs[0]) != 1 || !mxIsClass(prhs[0], "double"))
        mexErrMsgTxt("First Input must be a column-vector of type double.");
    if( mxGetN(prhs[1]) != 1 || !mxIsClass(prhs[1], "logical"))
        mexErrMsgTxt("Second Input must be a column-vector of type logical.");
    if( mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]) )
        mexErrMsgTxt("First Input and Second Input must be same length.");
    if( mxGetN(prhs[2]) != 2 || !mxIsClass(prhs[2], "uint64"))
        mexErrMsgTxt("Third Input must be a (x,2) Matrix of type uint64.");
    if( mxGetNumberOfElements(prhs[3]) != 1 || !mxIsClass(prhs[3], "double"))
        mexErrMsgTxt("Fourth Input must be a double scalar.");
    if( nrhs == 5 ){
        if( mxGetNumberOfElements(prhs[4]) != 1 || !mxIsClass(prhs[4], "double") )
            mexErrMsgTxt("Fifth Input must be a double scalar.");
        else
            ifBitReversal = (uint8_t)*mxGetPr(prhs[4]);
    }
    /* input */
    uint64_t n = mxGetNumberOfElements(prhs[0]);
    if( !Is2Power(n) )
        mexErrMsgTxt("Input Size must be a power of 2.");
    double *llr_p = mxGetPr(prhs[0]);
    mxLogical *frz_p = (mxLogical *)mxGetData(prhs[1]);
    uint64_t *dfrzCons = (uint64_t *)mxGetData(prhs[2]);
    uint64_t weight = mxGetM(prhs[2]);
    uint64_t L = (uint64_t)mxGetScalar(prhs[3]);
    if( !Is2Power(L) )
        mexErrMsgTxt(" L must be a power of 2 ! ");
    /* cast the input into a vector of integers */
    double *llr_f = calloc( n, sizeof( double ) );
    for (size_t i = 0; i < n; i++)
        llr_f[i] = llr_p[i];
    /* Creat struct for FRZ and DFRZ */
    struct u_node *u_nodes = calloc( n, sizeof( struct u_node ) );
    for (size_t i = 0; i < n; i++){
        u_nodes[i].frz = (uint8_t)frz_p[i];
        u_nodes[i].degree = 0;
    }
    for (size_t i = 0; i < weight; i++){
        u_nodes[ dfrzCons[i] - 1 ].frz = 2;
        u_nodes[ dfrzCons[i] - 1 ].degree++;
    }
    for (size_t i = 0; i < n; i++){
        if (u_nodes[i].frz != 2)
            continue;
        u_nodes[i].index = calloc( u_nodes[i].degree, sizeof( uint64_t ) );
    }
    uint64_t count = 0;
    for (size_t i = 0; i < n; i++){
        if (u_nodes[i].frz != 2)
            continue;
        for (size_t j = 0; j < u_nodes[i].degree; j++)
            u_nodes[i].index[j] = dfrzCons[ weight + count++ ] - 1;
    }
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(n, L, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, L, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, L, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *uhat_p          = mxGetPr(plhs[0]);
    double *chat_p          = mxGetPr(plhs[1]);
    double *PM_p            = mxGetPr(plhs[2]);
    double *log_minus_pNL   = mxGetData(plhs[3]);
    uint8_t *uhat_int = calloc( n*L, sizeof( uint8_t ) );
    uint8_t *chat_int = calloc( n*L, sizeof( uint8_t ) );
    double *PM_f = calloc( L, sizeof( double ) );
    /* use C functions in mexfunction */
    dec(uhat_int, chat_int, PM_f, llr_f, u_nodes, n, L, log_minus_pNL, ifBitReversal);
    /* cast to output */
    for (size_t i = 0; i < n*L; i++) 			
		uhat_p[i] = (double) uhat_int[i];
    for (size_t i = 0; i < n*L; i++) 			
		chat_p[i] = (double) chat_int[i];
    for (size_t i = 0; i < L; i++) 			
		PM_p[i] = PM_f[i];
    /* Clean up memory */
    free(llr_f);
    free(uhat_int);
    free(chat_int);
    free(PM_f);
    for (size_t i = 0; i < n; i++){
        if (u_nodes[i].frz != 2)
            continue;
        free(u_nodes[i].index);
    }
    free( u_nodes );
}