/* TODO: add copyright */


/* C function prototype */
/* void bl_dtrmm_asm_sve_8x8( int     k,
                              double  *a,
                              double  *b,
                              double  *c,
                              unsigned long ldc,
                              aux_t   *data,
                              int     offset
                            )
 */

#include "bl_config.h"

#define FUNC bl_dtrmm_asm_sve_8x8

#define origK       x0
#define origPA      x1
#define origPB      x2
#define origPC      x3
#define LDC         x4
#define auxData     x5
#define offset      x6
#define feedback    x7

#define counterL    x8
#define a_pntr      x9
#define c_pntr      x10
#define b_pntr      x11

#define tmp_counter        x15
#define tmp2        x16
// 18 must save
// 19 must save
// 20 must save
// 21 must save
// 22 must save
// 23 must save
// 24 must save
// 25 must save
// 26 must save
// 27 must save
// 28 must save
// 29 frame
// 30 link
// 31 sp

//v00 ~ v07: b ( l,0:7 ), dup
//v08 must save
//v09 must save
//v10 must save
//v11 must save
//v12 must save
//v13 must save
//v14 must save, a ( 0:3,l )
//v15 must save, a ( 4:7,l )
//v16 ~ v31: c
//v16 must save
//v17 must save

/******************************************************************************
* Macro definitions
******************************************************************************/
.macro SAVE_REGS
	add	sp, sp, #-(11 * 16)
	// stp	d8, d9, [sp, #(0 * 16)]
	// stp	d10, d11, [sp, #(1 * 16)]
	// stp	d12, d13, [sp, #(2 * 16)]
	stp	d14, d15, [sp, #(3 * 16)]
	stp	d16, d17, [sp, #(4 * 16)]
	// stp	x18, x19, [sp, #(5 * 16)]
	// stp	x20, x21, [sp, #(6 * 16)]
	// stp	x22, x23, [sp, #(7 * 16)]
	// stp	x24, x25, [sp, #(8 * 16)]
	// stp	x26, x27, [sp, #(9 * 16)]
	str	x28, [sp, #(10 * 16)]
.endm

.macro RESTORE_REGS
	mov	x0, #0				// set return value
	// ldp	d8, d9, [sp, #(0 * 16)]
	// ldp	d10, d11, [sp, #(1 * 16)]
	// ldp	d12, d13, [sp, #(2 * 16)]
	ldp	d14, d15, [sp, #(3 * 16)]
	ldp	d16, d17, [sp, #(4 * 16)]
	// ldp	x18, x19, [sp, #(5 * 16)]
	// ldp	x20, x21, [sp, #(6 * 16)]
	// ldp	x22, x23, [sp, #(7 * 16)]
	// ldp	x24, x25, [sp, #(8 * 16)]
	// ldp	x26, x27, [sp, #(9 * 16)]
	ldr	x28, [sp, #(10 * 16)]
	add	sp, sp, #(11*16)
.endm
/******************************************************************************
* End of macro definitions
******************************************************************************/

    .arch   armv8.2-a+sve
    .text
    .align  4
    .global FUNC
    .type   FUNC, %function

FUNC:
    SAVE_REGS
#if 1

    ptrue   p0.d, all
    /* copy from c (0:7,0) to z16, z17 */
    mov     x10, origPC    /* x10 as base pointer of c */
    ldr     z16, [x10]
    ldr     z17, [x10, #1, MUL VL]
    /* copy from c (0:7,1) to z18, z19 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z18, [x10]
    ldr     z19, [x10, #1, MUL VL]
    /* copy from c (0:7,2) to z20, z21 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z20, [x10]
    ldr     z21, [x10, #1, MUL VL]
    /* copy from c (0:7,3) to z22, z23 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z22, [x10]
    ldr     z23, [x10, #1, MUL VL]
    /* copy from c (0:7,4) to z24, z25 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z24, [x10]
    ldr     z25, [x10, #1, MUL VL]
    /* copy from c (0:7,5) to z26, z27 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z26, [x10]
    ldr     z27, [x10, #1, MUL VL]
    /* copy from c (0:7,6) to z28, z29 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z28, [x10]
    ldr     z29, [x10, #1, MUL VL]
    /* copy from c (0:7,7) to z30, z31 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    ldr     z30, [x10]
    ldr     z31, [x10, #1, MUL VL]


    /* initialize a_pntr, b_pntr */
    mov     a_pntr, origPA
    mov     b_pntr, origPB

    /* for counterL=k:0 */
    mov     counterL, origK

.Ldtrmm_kernel_loop_BEGIN:
    subs    counterL, counterL, #1
    blt     .Ldtrmm_kernel_loop_END
    /* crossed diagonal? */
    cmp     offset, #0
    blt     .Ldtrmm_kernel_loop_END

    /* load a ( 0:7, l ) into z14, z15 */
    ldr     z14, [a_pntr]
    ldr     z15, [a_pntr, #1, MUL VL]
    add     a_pntr, a_pntr, #64
    /* load b ( l,0:3 ) into 1st double of z0, z1, z2, z3 respectively */

/* TODO: for b  matrix, we can */
/*        use  LD1RQD, FMLA indexed which has no predicate. */
/* TODO: to replace ld4 and dup, fmla. */

    ld4     { v0.d, v1.d, v2.d, v3.d }[0], [b_pntr]
    add     b_pntr, b_pntr, #32
    /* load b ( l,4:7 ) into 1st double of z4, z5, z6 and z7 respectively */
    ld4     { v4.d, v5.d, v6.d, v7.d }[0], [b_pntr]
    add     b_pntr, b_pntr, #32

    /* dup b ( l,0 ) into all 4 elem. of z0 */
    /* dup b ( l,1 ) into all 4 elem. of z1 */
    /* dup b ( l,2 ) into all 4 elem. of z2 */
    /* dup b ( l,3 ) into all 4 elem. of z3 */
    dup     z0.d, z0.d[0]
    dup     z1.d, z1.d[0]
    dup     z2.d, z2.d[0]
    dup     z3.d, z3.d[0]
    /* dup b ( l,4 ) into all 4 elem. of z4 */
    /* dup b ( l,5 ) into all 4 elem. of z5 */
    /* dup b ( l,6 ) into all 4 elem. of z6 */
    /* dup b ( l,7 ) into all 4 elem. of z7 */
    dup     z4.d, z4.d[0]
    dup     z5.d, z5.d[0]
    dup     z6.d, z6.d[0]
    dup     z7.d, z7.d[0]

    /* calculate c */
    /* c( 0:3,0 ) += a ( 0:3,l ) * b( l,0 )
     *  VFMLA z16 += z14 * z0
     */
    fmla    z16.d, z14.d, z0.d[0]
    /* c( 4:7,0 ) += a ( 4:7,l ) * b( l,0 )
     *  VFMLA z17 += z15 * z0
     */
    fmla    z17.d, z15.d, z0.d[0]

    /* c( 0:7,1 ) += a ( 0:7,l ) * b( l,1 ) */
    fmla    z18.d, z14.d, z1.d[0]
    fmla    z19.d, z15.d, z1.d[0]

    /* c( 0:7,2 ) += a ( 0:7,l ) * b( l,2 ) */
    fmla    z20.d, z14.d, z2.d[0]
    fmla    z21.d, z15.d, z2.d[0]

    /* c( 0:7,3 ) += a ( 0:7,l ) * b( l,3 ) */
    fmla    z22.d, z14.d, z3.d[0]
    fmla    z23.d, z15.d, z3.d[0]

    /* c( 0:7,4 ) += a ( 0:7,l ) * b( l,4 ) */
    fmla    z24.d, z14.d, z4.d[0]
    fmla    z25.d, z15.d, z4.d[0]

    /* c( 0:7,5 ) += a ( 0:7,l ) * b( l,5 ) */
    fmla    z26.d, z14.d, z5.d[0]
    fmla    z27.d, z15.d, z5.d[0]

    /* c( 0:7,6 ) += a ( 0:7,l ) * b( l,6 ) */
    fmla    z28.d, z14.d, z6.d[0]
    fmla    z29.d, z15.d, z6.d[0]

    /* c( 0:7,7 ) += a ( 0:7,l ) * b( l,7 ) */
    fmla    z30.d, z14.d, z7.d[0]
    fmla    z31.d, z15.d, z7.d[0]

    /* decrease offset for next counterL */
    sub    offset, offset, #1
    b       .Ldtrmm_kernel_loop_BEGIN

.Ldtrmm_kernel_loop_END:

/*********************************************
 * store C
 *********************************************/

    /* store z16 ~ z31 to c (0:7, 0:7) in memory */
       /*   assumption: p0.d is all true
        *   ptrue   p0.d, all
        */
    mov     x10, origPC    /* x10 as base pointer of c */
    str     z16, [x10]
    str     z17, [x10, #1, MUL VL]
    /* column c (0:7,1) from z18, z19 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z18, [x10]
    str     z19, [x10, #1, MUL VL]
    /* column c (0:7,2) from z20, z21 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z20, [x10]
    str     z21, [x10, #1, MUL VL]
    /* column c (0:7,3) from z22, z23 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z22, [x10]
    str     z23, [x10, #1, MUL VL]
    /* column c (0:7,4) from z24, z25 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z24, [x10]
    str     z25, [x10, #1, MUL VL]
    /* column c (0:7,5) from z26, z27 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z26, [x10]
    str     z27, [x10, #1, MUL VL]
    /* column c (0:7,6) from z28, z29 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z28, [x10]
    str     z29, [x10, #1, MUL VL]
    /* column c (0:7,7) from z30, z31 */
    add     x10, x10, LDC, LSL #3   /* x10 points to next col. */
    str     z30, [x10]
    str     z31, [x10, #1, MUL VL]

    RESTORE_REGS

#else  // debugging code

/*************************************************************
 * save registers
**************************************************************/
	add	sp, sp, #-(11 * 16)
	stp	d8, d9, [sp, #(0 * 16)]
	stp	d10, d11, [sp, #(1 * 16)]
	stp	d12, d13, [sp, #(2 * 16)]
	stp	d14, d15, [sp, #(3 * 16)]
	stp	d16, d17, [sp, #(4 * 16)]
	stp	x18, x19, [sp, #(5 * 16)]
	stp	x20, x21, [sp, #(6 * 16)]
	stp	x22, x23, [sp, #(7 * 16)]
	stp	x24, x25, [sp, #(8 * 16)]
	stp	x26, x27, [sp, #(9 * 16)]
	str	x28, [sp, #(10 * 16)]

/*************************************************************
 * code body
**************************************************************/


    mov     tmp_counter, xzr
/* code samples for testing */
/* loading C to z16 ~ z31 */
    pfalse  p1.b
    ptrue   p0.d, all

.Ldtrmm_loading_c_BEGIN:
    pnext   p1.d, p0, p1.d
    b.none  .Ldtrmm_loading_c_END
    // do something
    //
    add     tmp_counter, tmp_counter, #1    // total l loops
    b       .Ldtrmm_loading_c_BEGIN

.Ldtrmm_loading_c_END:
    str     tmp_counter, [feedback]    // returns total l loops

/*************************************************************
 * restore registers
**************************************************************/
	mov	x0, #0				// set return value
	ldp	d8, d9, [sp, #(0 * 16)]
	ldp	d10, d11, [sp, #(1 * 16)]
	ldp	d12, d13, [sp, #(2 * 16)]
	ldp	d14, d15, [sp, #(3 * 16)]
	ldp	d16, d17, [sp, #(4 * 16)]
	ldp	x18, x19, [sp, #(5 * 16)]
	ldp	x20, x21, [sp, #(6 * 16)]
	ldp	x22, x23, [sp, #(7 * 16)]
	ldp	x24, x25, [sp, #(8 * 16)]
	ldp	x26, x27, [sp, #(9 * 16)]
	ldr	x28, [sp, #(10 * 16)]
	add	sp, sp, #(11*16)

#endif

    ret

