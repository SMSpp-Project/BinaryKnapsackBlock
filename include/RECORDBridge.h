/*--------------------------------------------------------------------------*/
/*------------------------- File RECORDBridge.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The bridge between RECORDBinaryKnapsackSolver and the RECORD solver for the
 * 0-1 knapsack (gitlab.com/renanfernandofranco/record). RECORD comes as a
 * single source file that puts its types, `using namespace std` and short
 * macros in the global scope, so it is compiled alone in RECORDBridge.cpp,
 * which includes no SMS++ header, and is reached only through the function
 * declared here.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#ifndef __RECORDBridge
 #define __RECORDBridge

#include <vector>

namespace SMSpp_di_unipi_it {
namespace record_bridge {

/// solve a 0-1 knapsack with RECORD
/** Solves max sum_i p_i x_i s.t. sum_i w_i x_i <= C, x binary, for the @p n
 * items with the given integer profits and weights (all positive, each
 * weight <= C), writing the optimal x in @p x and returning the optimal
 * value. With @p lb >= 0 only the solutions worth at least @p lb are sought
 * (RECORD prunes against it from the start): if there is none, the value
 * returned is below @p lb and @p x is a solution of that value, not an
 * optimal one. Throws std::runtime_error if RECORD stops on its time limit.
 */

long long solve_01( int n , const long long * p , const long long * w ,
                    long long C , std::vector< char > & x ,
                    long long lb = -1 );

}  // end( namespace record_bridge )
}  // end( namespace SMSpp_di_unipi_it )

#endif  /* RECORDBridge.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File RECORDBridge.h -------------------------*/
/*--------------------------------------------------------------------------*/
