/*--------------------------------------------------------------------------*/
/*------------------ File RECORDBinaryKnapsackSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the RECORDBinaryKnapsackSolver class.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>

#include "RECORDBinaryKnapsackSolver.h"

#include "RECORDBridge.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- FACTORY REGISTRATION -------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( RECORDBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

double RECORDBinaryKnapsackSolver::solve_integer_core(
                                                    std::vector< char > & in )
{
 // RECORD works on integers: a fractional profit leaves the core to the
 // core enumeration (warm started as well, see intReopt)
 const std::size_t n = v_w.size();
 std::vector< long long > p( n ) , w( n );
 for( std::size_t i = 0 ; i < n ; ++i ) {
  if( v_p[ i ] != std::floor( v_p[ i ] ) )
   return( CoreDPBinaryKnapsackSolver::solve_integer_core( in ) );
  p[ i ] = static_cast< long long >( v_p[ i ] );
  w[ i ] = v_w[ i ];
  }

 // warm start (intReopt): the previous solution, repaired to the current
 // data, is the incumbent RECORD has to beat, which it prunes against from
 // the start; if it finds nothing better, that solution is the optimum
 if( f_reopt && f_prev_valid ) {
  std::vector< char > win;
  const double lb = warm_incumbent( win );
  if( lb > - Inf< double >() ) {
   const long long zr = record_bridge::solve_01(
                     int( n ) , p.data() , w.data() , f_C , in ,
                     static_cast< long long >( lb ) + 1 );
   if( double( zr ) > lb )
    return( double( zr ) );
   in.swap( win );
   return( lb );
   }
  }

 return( double( record_bridge::solve_01( int( n ) , p.data() , w.data() ,
                                          f_C , in ) ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------- End File RECORDBinaryKnapsackSolver.cpp -----------------*/
/*--------------------------------------------------------------------------*/
