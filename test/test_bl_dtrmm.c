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
 * test_bl_dtrmm.c
 *
 *
 * Purpose:
 * test driver for BLISLAB dtrmm routine and reference dtrmm routine.
 *
 * Todo:
 *
 *
 * Modification:
 *
 *
 * */


#include "bl_dtrmm.h"

#define TOLERANCE 1E-11
// copy from A to AT
void copy_matrix(
        int    m,
        int    n,
        int    lda,
        double *A,
        double *AT
        )
{
    int i, j;

    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            AT[ j * lda + i ] = A[ j * lda + i ];
        }
    }
}

void computeError(
        int    ldc,
        int    ldc_ref,
        int    m,
        int    n,
        double *C,
        double *C_ref
        )
{
    int    i, j;
    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            if ( fabs( C( i, j ) - C_ref( i, j ) ) > TOLERANCE ) {
                printf( "C[ %d ][ %d ] != C_ref, %E, %E\n", i, j, C( i, j ), C_ref( i, j ) );
                break;
            }
        }
    }

}

void test_bl_dtrmm(
        const int m,
        const int n
        )
{
    int    i, j, p;
    double *A, *B, *B_MY, *B_REF;
    double tmp, error, flops;
    double ref_beg, ref_time, bl_dtrmm_beg, bl_dtrmm_time;
    int    nrepeats;
    int    lda, ldb, ldc, ldc_ref;
    double ref_rectime, bl_dtrmm_rectime;
    int    k = m;

    lda = m;
    // TODO: because of bl_micro_kernel() write beyond B's normal m*n, I have to increase ldb and column as well.
    ldb = k + 8;  // 8 refer to DGEMM_MR is 8

    A    = (double*)malloc( sizeof(double) * lda * k );
    // TODO: because of bl_micro_kernel() write beyond B's normal m*n, I have to increase ldb and column as well.
    B    = (double*)malloc( sizeof(double) * ldb * ( n + 8) );     // 4 refer to DGEMM_NR is 4
    B_MY    = (double*)malloc( sizeof(double) * ldb * ( n + 8) );
    B_REF   = (double*)malloc( sizeof(double) * ldb * ( n + 8) );

    if( A == NULL || B == NULL || B_MY == NULL || B_REF == NULL ) {
        printf("Allocation of memory failed. Abort.\n");
        return;
    }
    nrepeats = 1;

    srand48 (time(NULL));

    // Randonly generate points in [ 0, 1 ].
    for ( p = 0; p < k; p ++ ) {
        for ( i = 0; i < lda; i ++ ) {
            // A( i, p ) = (double)( 1.0 );
            A( i, p ) = (double)( drand48() );
        }
    }
    for ( j = 0; j < n; j ++ ) {
        for ( p = 0; p < ldb; p ++ ) {
            // B( p, j ) = (double)( 1.0 );
            B( p, j ) = (double)( drand48() );
        }
    }

    for ( j = 0; j < n; j ++ ) {
        for ( i = 0; i < ldb; i ++ ) {
            B_MY( i, j ) = (double)( 0.0 );
            B_REF( i, j ) = (double)( 0.0 );
        }
    }

    for ( i = 0; i < nrepeats; i ++ ) {
        // copy B to B_MY
        copy_matrix( k, n, ldb, B, B_MY );
        bl_dtrmm_beg = bl_clock();
        {
            #if 1 
            bl_dtrmm(
                    m,
                    n,
                    A,
                    lda,
                    B_MY,
                    ldb
                    );
            #endif
        }
        bl_dtrmm_time = bl_clock() - bl_dtrmm_beg;

        if ( i == 0 ) {
            bl_dtrmm_rectime = bl_dtrmm_time;
        } else {
            bl_dtrmm_rectime = bl_dtrmm_time < bl_dtrmm_rectime ? bl_dtrmm_time : bl_dtrmm_rectime;
        }
    }

    #if 1
    for ( i = 0; i < nrepeats; i ++ ) {
        // copy B to B_REF
        copy_matrix( k, n, ldb, B, B_REF );
        ref_beg = bl_clock();
        {
            bl_dtrmm_ref(
                    m,
                    n,
                    A,
                    lda,
                    B_REF,
                    ldb
                    );
        }
        ref_time = bl_clock() - ref_beg;

        if ( i == 0 ) {
            ref_rectime = ref_time;
        } else {
            ref_rectime = ref_time < ref_rectime ? ref_time : ref_rectime;
        }
    }

    computeError(
            ldb,
            ldb,
            m,
            n,
            B_MY,
            B_REF
            );
    #endif

    // Compute overall floating point operations.
    flops = ( m * n / ( 1000.0 * 1000.0 * 1000.0 ) ) * ( 2 * k );

    printf( "%5d\t %5d\t %5d\t %5.2lf\t %5.2lf\n", 
            m, n, k, flops / bl_dtrmm_rectime, flops / ref_rectime );

    #if 1
    free( A     );
    free( B     );
    free( B_MY  );
    #endif
    #if 1
    free( B_REF );
    #endif
}

int main( int argc, char *argv[] )
{
    int    m, n;

    if ( argc != 3 ) {
        printf( "Error: require 2 arguments, but only %d provided.\n", argc - 1 );
        exit( 0 );
    }

    sscanf( argv[ 1 ], "%d", &m );
    sscanf( argv[ 2 ], "%d", &n );

    test_bl_dtrmm( m, n );

    return 0;
}

