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
 * bl_dtrmm_ref.c
 *
 *
 * Purpose:
 *
 * Todo:
 *
 *
 * Modification:
 *
 *
 * */

#ifdef USE_BLAS
/*
 * dtrmm prototype
 *
 */
#include <cblas.h>
extern void cblas_dtrmm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
#undef DGEMM_NR /* the same name is used in cblas.h and bl_dtrmm.h but with different meanings */
#undef DGEMM_NC /* the same name is used in cblas.h and bl_dtrmm.h but with different meanings */
#endif

#include <bl_dtrmm.h>

void bl_dtrmm_ref(
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

#ifdef USE_BLAS
    cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower,
            CblasNoTrans, CblasNonUnit,
            m, n, alpha, XA, lda, XB, ldb);

#else
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
#endif

}

