/*
 * --------------------------------------------------------------------------
 * BLISLAB 
 * --------------------------------------------------------------------------
 * Copyright (C) 2016, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * bl_dtrmm.c
 *
 *
 * Purpose:
 * this is the main file of blislab dtrmm.
 *
 * Todo:
 *
 *
 * Modification:
 *
 * 
 * */

#include <stdio.h>

#include "bl_dgemm_kernel.h"
#include "bl_dtrmm.h"

static inline void packA_zeros_mcxkc_d(
        int    m,
        int    k,
        double *XA,
        int    ldXA,
        int    offseta,
        double *packA
        )
{
    int    i, p;
    double *a_pntr[ DGEMM_MR ];

    for ( i = 0; i < m; i ++ ) {
        a_pntr[ i ] = XA + ( offseta + i );
    }

    for ( i = m; i < DGEMM_MR; i ++ ) {
        a_pntr[ i ] = XA + ( offseta + 0 );
    }

    for ( p = 0; p < k; p ++ ) {
        for ( i = 0; i < DGEMM_MR; i ++ ) {
            *packA = (double) 0.0; 
            packA ++;
            a_pntr[ i ] = a_pntr[ i ] + ldXA;
        }
    }
}

static inline void packA_mcxkc_d(
        int    m,
        int    k,
        double *XA,
        int    ldXA,
        int    offseta,
        double *packA
        )
{
    int    i, p;
    double *a_pntr[ DGEMM_MR ];

    for ( i = 0; i < m; i ++ ) {
        a_pntr[ i ] = XA + ( offseta + i );
    }

    for ( i = m; i < DGEMM_MR; i ++ ) {
        a_pntr[ i ] = XA + ( offseta + 0 );
    }

    for ( p = 0; p < k; p ++ ) {
        for ( i = 0; i < DGEMM_MR; i ++ ) {
            *packA = *a_pntr[ i ];
            packA ++;
            a_pntr[ i ] = a_pntr[ i ] + ldXA;
        }
    }
}

static inline void packA_diag_mcxkc_d(
        int    m,
        int    k,
        double *XA,
        int    ldXA,
        int    offseta,
        int    offseta_y,
        double *packA
        )
{
    int    i, p;
    int    px, py;
    double *a_pntr[ DGEMM_MR ];

    for ( i = 0; i < m; i ++ ) {
        a_pntr[ i ] = XA + ( offseta + i );
    }

    for ( i = m; i < DGEMM_MR; i ++ ) {
        a_pntr[ i ] = XA + ( offseta + 0 );
    }

    for ( p = 0, py = offseta_y; p < k; p ++, py ++ ) {
        for ( i = 0, px = offseta; i < DGEMM_MR; i ++, px ++ ) {
            // TODO: make these loops in reverse order,
            //      so, when at diagonal, we can 'break', instead of 'contrinue'.
            if ( px < py )     // 0 for upper triangular
                *packA = (double) 0.0;
            else
                *packA = *a_pntr[ i ];
            // printf("packing A: (%d, %d) = %f\n", px, py, *packA);
            packA ++;
            a_pntr[ i ] = a_pntr[ i ] + ldXA;
        }
    }
}

static inline void zero_matrix(
        int m, 
        int n,
        double *A,
        int lda
        )
{
    int i, j;
    double *p;

    p = A;
    for ( j = 0; j < n; j ++ ) {
        for ( i = 0; i < m; i ++ )
            *p ++ = (double) 0.0;
        p += ( lda - m );
    }
}
/*
 * --------------------------------------------------------------------------
 */

static inline void packB_kcxnc_d(
        int    n,
        int    k,
        double *XB,
        int    ldXB, // ldXB is the original k
        int    offsetb,
        double *packB
        )
{
    int    j, p; 
    double *b_pntr[ DGEMM_NR ];

    for ( j = 0; j < n; j ++ ) {
        b_pntr[ j ] = XB + ldXB * ( offsetb + j );
    }

    for ( j = n; j < DGEMM_NR; j ++ ) {
        b_pntr[ j ] = XB + ldXB * ( offsetb + 0 );
    }

    for ( p = 0; p < k; p ++ ) {
        for ( j = 0; j < DGEMM_NR; j ++ ) {
            *packB ++ = *b_pntr[ j ] ++;
        }
    }
}

/*
 * --------------------------------------------------------------------------
 */
void bl_macro_kernel(
        int    m,
        int    n,
        int    k,
        double *packA,
        double *packB,
        double *C,
        int    ldc,
        int    xa,    /* X offset of packA in original XA */
        int    ya     /* Y offset of packA in orininal XA */
        )
{
    int    i, j;
    aux_t  aux;

    aux.b_next = packB;

    // printf("in bl_macro_kernel()\n");
    // printf("   xa=%d, ya=%d\n", xa, ya);
    for ( j = 0; j < n; j += DGEMM_NR ) {                        // 2-th loop around micro-kernel
        aux.n  = min( n - j, DGEMM_NR );
        for ( i = 0; i < m; i += DGEMM_MR ) {                    // 1-th loop around micro-kernel
            aux.m = min( m - i, DGEMM_MR );
            // if this A is in upper triangular, then it's all zero, skip.
            // bottom-left element index is: ( xa + i + aux.m - 1, ya )
            // printf("micro kernel for A: (%d, %d), height=%d, width=%d\n",
            //       xa + i, ya, aux.m, k);
            // printf("             packB: column=%d, width=%d\n", j, aux.n);

            // Note: skip at the same mr*kc sub-A where it's also skipped in Loop 3.
            if ( xa + i + aux.m - 1 < ya ) {     
                 // printf("all ZEROs, skip\n");
                 continue;
            }
            if ( i + DGEMM_MR >= m ) {
                aux.b_next += DGEMM_NR * k;
            }

            ( *bl_trmm_micro_kernel ) (
                    k,
                    &packA[ i * k ],
                    &packB[ j * k ],
                    &C[ j * ldc + i ],
                    (unsigned long long) ldc,
                    &aux,
                    xa + i,
                    ya
                    );
        }                                                        // 1-th loop around micro-kernel
    }                                                            // 2-th loop around micro-kernel
}

#if 1
// C must be aligned
void bl_dtrmm(
        int    m,
        int    n,
        double *XA,
        int    lda,
        double *XB,
        int    ldb
        )
{
    int    i, j, p;
    int    ic, ib, jc, jb, pc, pb;
    int    wb;
    int    ir, jr;
    double *packA, *packB_NC, *packB;
    char   *str;
    int    k = m;

    // Early return if possible
    if ( m == 0 || n == 0 ) {
        printf( "bl_dgemm(): early return\n" );
        return;
    }
    // printf("m=%d, n=%d, lda=%d, ldb=%d\n", m, n, lda, ldb);

    // Allocate packing buffers
    packA  = bl_malloc_aligned( DGEMM_KC, ( DGEMM_MC + 1 ) , sizeof(double) );
    // Allocate packB_NC to hold all nc column of B
    packB_NC  = bl_malloc_aligned( k, ( DGEMM_NC + 1 ) , sizeof(double) );

    for ( jc = 0; jc < n; jc += DGEMM_NC ) {                                       // 5-th loop around micro-kernel
        jb = min( n - jc, DGEMM_NC );

        // TODO: packing all k rows of XB in one loop destroys packB's preload in
        //       cache hierarchy of memory system.
        //  Reason I chnage this is because unlike GEMM, TRMM stores the results 
        //    back into XB. So, this part of XB has to be moved then erased to zero.
        //
        // pack jb column of XB
        packB = packB_NC;
        for ( pc = 0; pc < k; pc += DGEMM_KC ) {                                   // pack B, kc row each step
            pb = min( k - pc, DGEMM_KC );
            for ( j = 0; j < jb; j += DGEMM_NR ) {
                #if 1
                packB_kcxnc_d(
                        min( jb - j, DGEMM_NR ),
                        pb,
                        &XB[ pc ],
                        ldb, // should be ldXB instead
                        jc + j,
                        &packB[ j * pb ]
                        );
                #endif
            }
            packB += pb * ( DGEMM_NC + 1 );
        }
        // zero k*jb of XB
        zero_matrix( k, jb, XB + jc * ldb, ldb );

        packB = packB_NC;
        for ( pc = 0; pc < k; pc += DGEMM_KC ) {                                   // 4-th loop around micro-kernel
            pb = min( k - pc, DGEMM_KC );

            for ( ic = 0; ic < m; ic += DGEMM_MC ) {                               // 3-rd loop around micro-kernel
                // printf("in loop 3.\n");
                ib = min( m - ic, DGEMM_MC );


                for ( i = 0; i < ib; i += DGEMM_MR ) {                             // Loop 3.1, to packA 
                    wb = min( ib - i, DGEMM_MR );
                    // left-top element index:    ( ic + i, pc)
                    // height: wb, length: pb
                    //
                    // left-bottom element index: ( ic + i + wb - 1, pc )
                    // top-right element index: ( ic + i, pc + pb -1 )
                    #if 0
                    printf("in loop 3.1 to packA. height=%d, length=%d\n", wb, pb);
                    printf("  top-left element index:    (%d, %d)\n", ic + i, pc);  
                    printf("  top-right element index:   (%d, %d)\n", ic + i, pc + pb -1 );
                    printf("  left-bottom element index: (%d, %d)\n", ic + i + wb - 1, pc );
                    #endif
                    // if 'left-bottom' is in upper triangular, then it's all zeros, no need to pack.
                    if ( ic + i + wb - 1 < pc ) {
                        // printf("all zeros, skip\n"); 
                        #if 0
                        packA_zeros_mcxkc_d(
                                wb,
                                pb,
                                &XA[ pc * lda ],
                                lda,
                                ic + i,
                                &packA[ 0 * DGEMM_MC * pb + i * pb ]
                                );
                        #endif
                        continue;
                    }
                    // else if 'top-right' is in lower triangular or diagonal, then it should be packed the same way as gemm.
                    else if ( ic + i >= pc + pb - 1 ) {
                        // printf("call packA_mcxkc_d()\n");
                        packA_mcxkc_d(
                                wb,
                                pb,
                                &XA[ pc * lda ],
                                lda,
                                ic + i,
                                &packA[ 0 * DGEMM_MC * pb + i * pb ]
                                );
                    }
                    // else, it crosses the diagonal, and should be handled carefully as triangular way. upper triangular is zero.
                    else {
                        // printf("call packA_diag_mcxkc_d()\n");
                        #if 1
                        packA_diag_mcxkc_d(
                                wb,
                                pb,
                                &XA[ pc * lda ],
                                lda,
                                ic + i,
                                pc,
                                &packA[ 0 * DGEMM_MC * pb + i * pb ]
                                );
                        #endif
                    }
                }

                // printf("call bl_macro_kernel(), top-left index: (%d, %d)\n", ic, pc);
                // printf("    m=%d, n=%d, k=%d\n", ib, jb, pb);
                #if 1
                bl_macro_kernel(
                        ib,
                        jb,
                        pb,
                        packA,
                        packB,
                        &XB[ jc * ldb + ic ], 
                        ldb,
                        ic,
                        pc
                        );
                #endif
            }                                                                     // End 3.rd loop around micro-kernel
            packB += pb * ( DGEMM_NC + 1 );
        }                                                                         // End 4.th loop around micro-kernel
    }                                                                             // End 5.th loop arouncolumn
    free( packA );
    free( packB_NC );
}

#else 

void bl_dtrmm(
        int    m,
        int    n,
        double *XA,
        int    lda,
        double *XB,
        int    ldb
        )
{
    // Local variables.
    int    i, j, p;
    double alpha = 1.0;
    double *XC;

    // Sanity check for early return.
    if ( m == 0 || n == 0 ) return;

    // Reference TRMM implementation.

    XC = bl_malloc_aligned( m, n, sizeof(double) );
    // zero out XC
    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            XC[ j * m + i ] =  (double)( 0.0 );
        }
    }

    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            for ( p = 0; p <= i; p ++ ) {
                // for lower triangular matrix XA, upper triangle is 0
                XC[ j * m + i ] += XA[ p * lda + i ] * XB[ j * ldb + p ];
            }
        }
    }

    // copy result to XB
    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            XB[ j * ldb + i ] = XC[ j * m + i ];
        }
    }
    free( XC );
}
#endif
