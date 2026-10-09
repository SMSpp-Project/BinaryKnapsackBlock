/*--------------------------------------------------------------------------*/
/*------------------------ File RECORDBridge.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the bridge to RECORD (see RECORDBridge.h). RECORD.hpp is
 * RECORD's source file with its main() cut away, generated at configure time
 * from the RECORD checkout given by RECORD_ROOT (see CMakeLists.txt): the
 * main() is where the source as distributed ignores the multiplicity it
 * reads and solves a bounded knapsack, so the 0-1 problem is set up here.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#include "RECORD.hpp"

#include <stdexcept>

#include "RECORDBridge.h"

/*--------------------------------------------------------------------------*/

long long SMSpp_di_unipi_it::record_bridge::solve_01( int n ,
                                                     const long long * p ,
                                                     const long long * w ,
                                                     long long C ,
                                                     std::vector< char > & x ,
                                                     long long lb )
{
 x.assign( n , 0 );
 long long sw = 0 , sp = 0;
 for( int i = 0 ; i < n ; ++i )
  { sw += w[ i ]; sp += p[ i ]; }
 if( sw <= C ) {               // everything fits (RECORD expects it handled)
  x.assign( n , 1 );
  return( sp );
  }

 clock_gettime( CLOCK_MONOTONIC , &t0 );      // RECORD's time-limit clock
 ::Solver s( n );
 for( int i = 0 ; i < n ; ++i )
  s.add_item( p[ i ] , w[ i ] , i , 1 );      // multiplicity 1: 0-1 items
 if( lb >= 0 )                 // as RECORD does for its own sub-problems
  s.lb = lb;
 const long long z = s.Solve( C );
 if( z < 0 )
  throw( std::runtime_error(
   "record_bridge::solve_01: RECORD stopped on its time limit" ) );

 for( const auto & it : s.items )              // back in insertion order
  x[ it.id ] = char( it.val );
 return( z );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File RECORDBridge.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
