/* TODO: add copyright */


/* C function prototype */
/* void bl_dtrmm_asm_sve_8x8( int     k,
                              double  *a,
                              double  *b,
                              double  *c,
                              unsigned long ldc,
                              aux_t   *data,
                              int     xa,
                              int     ya )
 */

#include "bl_config.h"

#define FUNC bl_dtrmm_asm_sve_8x8

#define origK       x0
#define origPA      x1
#define origPB      x2
#define origPC      x3
#define LDC         x4
#define auxData     x5


/******************************************************************************
* Macro definitions
******************************************************************************/

.macro INIT8x4
    fmov    d16, xzr
    fmov    d17, xzr
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
    INIT8x4
    ret

