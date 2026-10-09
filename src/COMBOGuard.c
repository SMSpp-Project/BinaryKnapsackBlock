/*--------------------------------------------------------------------------*/
/*-------------------------- File COMBOGuard.c -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The only translation unit that sees COMBO. combo.c is the copy made at
 * configure time from the one in COMBO_ROOT (see CMakeLists.txt), with the
 * product type widened to 128-bit integers: as distributed it is a double,
 * which loses the sign of the determinants of profits and weights past
 * 2^53. combo.c ends a run with exit() on an internal error (state space
 * or memory exhausted, inconsistent data): here exit() jumps back to
 * combo_guarded(), which reports the error to the caller, and the message
 * COMBO prints before it is suppressed. The memory COMBO holds at that
 * point is not freed.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

static jmp_buf combo_guard_env;

static void combo_guard_exit( int code )
{
 (void) code;
 longjmp( combo_guard_env , 1 );
 }

#define exit( code ) combo_guard_exit( code )
#define printf( ... ) ( (void) 0 )
#define vprintf( ... ) ( (void) 0 )

#include "combo.c"

#undef exit
#undef printf
#undef vprintf

/*--------------------------------------------------------------------------*/

int combo_guarded( item * f , item * l , stype c , stype lb , stype * z )
{
 if( setjmp( combo_guard_env ) )
  return( -1 );
 *z = combo( f , l , c , lb , 0 , 1 , 0 );
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File COMBOGuard.c ---------------------------*/
/*--------------------------------------------------------------------------*/
