/*--------------------------------------------------------------------------*/
/*-------------------------- File COMBOBridge.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The bridge between COMBOBinaryKnapsackSolver and the COMBO solver for the
 * 0-1 knapsack. COMBO is a C source file (combo.c) that ends a run with
 * exit() on an internal error: it is compiled alone in COMBOGuard.c, which
 * turns that exit() into an error code, and is reached only through the
 * function declared here.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#ifndef __COMBOBridge
 #define __COMBOBridge

#include <vector>

namespace SMSpp_di_unipi_it {
namespace combo_bridge {

/// solve a 0-1 knapsack with COMBO
/** Solves max sum_i p_i x_i s.t. sum_i w_i x_i <= C, x binary, for the @p n
 * items with the given integer profits and weights (all positive, each
 * weight <= C), writing the optimal x in @p x and returning the optimal
 * value. With @p lb >= 0 only the solutions worth at least @p lb are sought:
 * if there is none, the value returned is below @p lb and @p x is
 * meaningless. Throws std::runtime_error if COMBO stops on an internal error
 * (e.g., its state space is exhausted). */

long long solve_01( int n , const long long * p , const long long * w ,
                    long long C , std::vector< char > & x ,
                    long long lb = -1 );

}  // end( namespace combo_bridge )
}  // end( namespace SMSpp_di_unipi_it )

#endif  /* COMBOBridge.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File COMBOBridge.h -------------------------*/
/*--------------------------------------------------------------------------*/
