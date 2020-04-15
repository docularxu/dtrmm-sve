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
                   int    offset )
{
    int l, j, i;
    register double a_il_reg;
    register double b_lj0_reg, b_lj1_reg, b_lj2_reg, b_lj3_reg;
    register double c_ij0_reg, c_ij1_reg, c_ij2_reg, c_ij3_reg;
    double *b0_pntr;
    double *a_pntr;
    double *c_pntr;   // TODO: c_pntr can be serialized too.

    b0_pntr = b;

    for ( l = 0; l < k; ++l )                               // Loop 0.1, column [l] of A, times, row [l] of B
    {                 
        // if column[l] of a is in upper triangular, then it is all-zeros, no need to calculate.
        // also, all columns following column[l], are all-zeros too.
        //
        // lowest element of column[l] of A is: (xa+aux.m-1, ya+l)
        if ( offset < l )
            break;

        for ( j = 0; j < DGEMM_NR; j+=4 )                    // Loop 0.2, walk through row[l] of B
        { 
            b_lj0_reg = *b0_pntr ++ ;
            b_lj1_reg = *b0_pntr ++ ;
            b_lj2_reg = *b0_pntr ++ ;
            b_lj3_reg = *b0_pntr ++ ;

            a_pntr = a + DGEMM_MR * l;

            for ( i = 0; i < DGEMM_MR; ++i )                // Loop 0.3, walk through column[l] of A
            {
                a_il_reg = *a_pntr ++;

                c_ij0_reg = (double) 0.0;
                c_ij1_reg = (double) 0.0;
                c_ij2_reg = (double) 0.0;
                c_ij3_reg = (double) 0.0;

                // TODO: b(l, j+0:3) can be loaded into a length-4 vector using one instruction
                //       but for now, let's read them from b_pntr.
                c_ij0_reg += a_il_reg * b_lj0_reg;  // b( l, j + 0 );
                c_ij1_reg += a_il_reg * b_lj1_reg;  // b( l, j + 1 );
                c_ij2_reg += a_il_reg * b_lj2_reg;  // b( l, j + 2 );
                c_ij3_reg += a_il_reg * b_lj3_reg;  // b( l, j + 3 );

                // memory store
                c( i, j + 0 ) += c_ij0_reg;
                c( i, j + 1 ) += c_ij1_reg;
                c( i, j + 2 ) += c_ij2_reg;
                c( i, j + 3 ) += c_ij3_reg;
            }  // end of i
        }  // end of j
    }  // end of l

}
