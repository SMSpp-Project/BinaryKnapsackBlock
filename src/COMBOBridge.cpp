/*--------------------------------------------------------------------------*/
/*------------------------- File COMBOBridge.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the bridge to COMBO (see COMBOBridge.h). COMBO reorders
 * the items it is given and keeps no index of them: the solution is mapped
 * back by profit and weight, which identical items share, so that any
 * assignment of the chosen copies among them is an optimal solution as
 * well.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <numeric>
#include <stdexcept>

extern "C" {
#include "combo.h"

int combo_guarded( item * f , item * l , stype c , stype lb , stype * z );
}

#include "COMBOBridge.h"

/*--------------------------------------------------------------------------*/

long long SMSpp_di_unipi_it::combo_bridge::solve_01( int n ,
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
 if( sw <= C ) {               // everything fits (COMBO expects it handled)
  x.assign( n , 1 );
  return( sp );
  }

 std::vector< item > it( n );
 for( int i = 0 ; i < n ; ++i )
  it[ i ] = { itype( p[ i ] ) , itype( w[ i ] ) , 0 };

 // COMBO only updates the solution when it beats its lower bound, which is
 // a value already attained: one less than the value sought
 const stype zlb = lb > 0 ? stype( lb - 1 ) : 0;
 stype z;
 if( combo_guarded( it.data() , it.data() + n - 1 , stype( C ) , zlb , &z ) )
  throw( std::runtime_error(
   "combo_bridge::solve_01: COMBO stopped on an internal error" ) );
 if( z <= zlb )                // nothing beats lb (or the optimum is 0)
  return( z );

 // back to the input order by (profit, weight)
 auto key = []( long long a , long long b , long long c , long long d ) {
  return( ( a < c ) || ( ( a == c ) && ( b < d ) ) );
  };
 std::vector< int > ord( n );
 std::iota( ord.begin() , ord.end() , 0 );
 std::sort( ord.begin() , ord.end() , [ & ]( int a , int b ) {
  return( key( p[ a ] , w[ a ] , p[ b ] , w[ b ] ) );
  } );
 std::sort( it.begin() , it.end() , [ & ]( const item & a , const item & b ) {
  return( key( a.p , a.w , b.p , b.w ) );
  } );
 for( int k = 0 ; k < n ; ) {  // per group, the chosen copies first
  int e = k;
  int ch = 0;
  while( ( e < n ) && ( it[ e ].p == it[ k ].p ) &&
         ( it[ e ].w == it[ k ].w ) )
   ch += it[ e++ ].x ? 1 : 0;
  for( int j = k ; j < k + ch ; ++j )
   x[ ord[ j ] ] = 1;
  k = e;
  }

 long long xp = 0 , xw = 0;
 for( int i = 0 ; i < n ; ++i )
  if( x[ i ] )
   { xp += p[ i ]; xw += w[ i ]; }
 if( ( xp != z ) || ( xw > C ) )
  throw( std::runtime_error(
   "combo_bridge::solve_01: the solution of COMBO is inconsistent" ) );
 return( z );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File COMBOBridge.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
