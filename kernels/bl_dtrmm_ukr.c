#include <bl_config.h>
#include "bl_dgemm_kernel.h"

//micro-panel a is stored in column major, lda=DGEMM_MR.
#define a(i,j) a[ (j)*DGEMM_MR + (i) ]
//micro-panel b is stored in row major, ldb=DGEMM_NR.
#define b(i,j) b[ (i)*DGEMM_NR + (j) ]
//result      c is stored in column major.
#define c(i,j) c[ (j)*ldc + (i) ]


void bl_dtrmm_ukr( int    k,
                   double *a,
                   double *b,
                   double *c,
                   unsigned long long ldc,
                   aux_t* data,
                   int    xa,
                   int    ya )
{
    int l, j, i;

    for ( l = 0; l < k; ++l )                               // Loop 0.1, column [l] of A, times, row [l] of B
    {                 
        // if column[l] of a is in upper triangular, then it is all-zeros, no need to calculate.
        // also, all columns following column[l], are all-zeros too.
        //
        // lowest element of column[l] of A is: (xa+aux.m-1, ya+l)
        if ( xa + data->m - 1 < ya + l )
            break;

        for ( j = 0; j < DGEMM_NR; ++j )                    // Loop 0.2, walk through row[l] of B
        { 
            for ( i = 0; i < DGEMM_MR; ++i )                // Loop 0.3, walk through column[l] of A
            { 
                c( i, j ) += a( i, l ) * b( l, j );
            }
        }
    }

}