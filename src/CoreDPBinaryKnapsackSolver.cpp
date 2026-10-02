/*--------------------------------------------------------------------------*/
/*------------------- File CoreDPBinaryKnapsackSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class CoreDPBinaryKnapsackSolver.
 *
 * Break-item core enumeration of the 0-1 knapsack with dominance (Pareto
 * frontier) and an LP upper bound (per-state pruning + global early-exit), i.e.
 * the MINKNAP scheme: only the deviations from the break solution are
 * enumerated, keeping a small state set in memory. On top of it: a strong
 * initial heuristic pool seeding the incumbent, the Martello-Toth U2 ceiling,
 * the surrogate / cardinality relaxation - both as a fractional bound and
 * solved exactly by recursion, so that its optimum certifies and possibly
 * raises the incumbent - weight- and profit-gcd divisibility reductions,
 * Dembo-Hammer per-item fixing, fixing-by-dominance, the guarded DP extension,
 * item aggregation with multiplicity reduction (binary-decomposed copy
 * batches with early termination), and the PH / TPH / SSPH / SPH / GCH primal
 * heuristics running during the enumeration. Correctness is validated against
 * :MILPSolver over the full feature matrix and against the Pisinger reference
 * optima (smallcoeff / largecoeff / hardinstances).
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <type_traits>

#include "CoreDPBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef CORE_STATS
#define CORE_STATS 0             // stderr counters per solve (dev only)
#endif

#if CORE_STATS
static std::size_t core_rec_states = 0;  // states of the recursive calls
static std::size_t cs_merge , cs_nogrow , cs_cand , cs_b2 , cs_b3 , cs_b4;
#define CSTAT( x ) ( x )
#else
#define CSTAT( x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

// length of the longest prefix of v, in the order `before`, whose sum stays
// <= limit: quickselect halving (expected linear) instead of a full sort,
// only the last few elements are sorted; v is permuted
template< class T , class Before >
int prefix_fit( std::vector< T > & v , Before before , T limit )
{
 int lo = 0 , hi = int( v.size() ) , k = 0;
 while( hi - lo > 16 ) {
  const int mid = lo + ( hi - lo ) / 2;
  std::nth_element( v.begin() + lo , v.begin() + mid , v.begin() + hi ,
                    before );
  T s = 0;
  for( int i = lo ; i < mid ; ++i )
   s += v[ i ];
  if( s <= limit ) { limit -= s; k += mid - lo; lo = mid; }
  else hi = mid;
  }
 std::sort( v.begin() + lo , v.begin() + hi , before );
 while( ( lo < hi ) && ( v[ lo ] <= limit ) )
  { limit -= v[ lo++ ]; ++k; }
 return( k );
 }

// a value computed in floating point which is mathematically an integer can
// land just below or just above it, and rounding it would then lose a whole
// unit: bounds are rounded past a relative tolerance (well above the
// machine precision, well below one unit for values up to 1e12), in the
// direction that can only loosen them (a looser bound is still valid, a
// tighter one not)
inline double floor_safe( double x )
{
 return( std::floor( x + 1e-12 * std::max( 1.0 , std::abs( x ) ) ) );
 }

inline double ceil_safe( double x )
{
 return( std::ceil( x - 1e-12 * std::max( 1.0 , std::abs( x ) ) ) );
 }

// an allocator that default-initialises: on a vector of a trivial type,
// resize() then leaves the new elements as they are instead of zeroing them
// (for buffers that are written before being read)
template< class T >
struct NoInit : std::allocator< T > {
 template< class U > struct rebind { using other = NoInit< U >; };
 NoInit( void ) = default;
 template< class U > NoInit( const NoInit< U > & ) {}
 template< class U , class... A >
 void construct( U * q , A &&... a ) {
  if constexpr( sizeof...( A ) == 0 )
   ::new( static_cast< void * >( q ) ) U;
  else
   ::new( static_cast< void * >( q ) ) U( std::forward< A >( a )... );
  }
 };

}  // end( unnamed namespace )

/*--------------------------------------------------------------------------*/
/*--------------------------- FACTORY REGISTRATION -------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( CoreDPBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int CoreDPBinaryKnapsackSolver::compute( bool changedvars )
{
 lock();

 if( update_instance() )       // process all the pending modifications
  extract_instance();          // anything changed: rebuild the cores

 if( f_C < 0 ) {                 // capacity exhausted by the fixed-to-1 items
  f_obj = - Inf< double >();
  unlock();
  return( kInfeasible );
  }

 try {
  enumerate_states();
  }
 catch( ... ) {                // an external solver gave up: not left locked
  unlock();
  throw;
  }

 v_prev_x = f_x;               // the seed of the next (warm started) solve
 f_prev_valid = true;

 unlock();
 return( kOK );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

void CoreDPBinaryKnapsackSolver::extract_instance( void )
{
 normalize_instance();         // base: raw mirror -> normalized double cores

 v_w.clear();  v_p.clear();  v_orig.clear();  v_comp.clear();
 v_cw.clear(); v_cp.clear(); v_corig.clear(); v_ccomp.clear();

 if( ! f_Block ) {              // detaching the Solver
  f_C = 0;
  f_obj = - Inf< double >();
  return;
  }

 f_C = long( std::floor( f_Cd ) );

 // the core enumeration requires integer weights
 auto to_long = []( double w ) -> long {
  const long l = long( w >= 0 ? w + 0.5 : w - 0.5 );   // nearest integer
  if( std::abs( double( l ) - w ) > WeightIntegrality )
   throw( std::invalid_argument(
    "CoreDPBinaryKnapsackSolver::extract_instance: weights must be integers"
    ) );
  return( l );
  };

 // integer items that cannot fit the capacity are pre-fixed: a normal one to
 // 0, a complement one to 1 (its toggle could never be undone)
 v_w.reserve( c_w.size() );  v_p.reserve( c_w.size() );
 v_orig.reserve( c_w.size() );  v_comp.reserve( c_w.size() );
 for( std::size_t k = 0 ; k < c_w.size() ; ++k ) {
  const long w = to_long( c_w[ k ] );
  if( w > f_C ) {
   f_x[ c_orig[ k ] ] = c_comp[ k ] ? 1 : 0;
   continue;
   }
  v_w.push_back( w ); v_p.push_back( c_p[ k ] );
  v_orig.push_back( c_orig[ k ] ); v_comp.push_back( c_comp[ k ] );
  }

 for( std::size_t k = 0 ; k < cc_w.size() ; ++k ) {
  v_cw.push_back( to_long( cc_w[ k ] ) ); v_cp.push_back( cc_p[ k ] );
  v_corig.push_back( cc_orig[ k ] ); v_ccomp.push_back( cc_comp[ k ] );
  }
 }

/*--------------------------------------------------------------------------*/

void CoreDPBinaryKnapsackSolver::enumerate_states( void )
{
 // add the core optimum to the base profit and write the core variables; the
 // pure-binary case uses the fast MINKNAP break-item enumeration, the case
 // with continuous variables the Pareto-frontier + fractional-fill scheme
 f_obj = f_base;
 std::vector< char > in;
 if( v_cw.empty() ) {
  double z = - Inf< double >();
  f_outcome = ( f_reopt && f_prev_valid ) ? 1 : 0;
  if( ( f_reopt >= 2 ) && f_prev_valid )
   z = last_still_optimal( in );
  if( z == - Inf< double >() )
   z = solve_integer_core( in );
  f_obj += z;
  if( f_reopt >= 2 ) {         // the core and solution the next solve checks
   f_last_C = f_C;
   v_last_w = v_w;  v_last_p = v_p;
   v_last_orig = v_orig;  v_last_comp = v_comp;
   v_last_in = in;
   }
  }
 else {
  std::vector< double > cx;
  f_outcome = 0;                // not warm started (see intReopt)
  f_obj += solve_with_continuous( in , cx );
  for( std::size_t j = 0 ; j < v_cw.size() ; ++j )
   f_x[ v_corig[ j ] ] = v_ccomp[ j ] ? ( 1.0 - cx[ j ] ) : cx[ j ];
  }

 // complemented items contribute x_orig = 1 - ( chosen in the core )
 for( std::size_t i = 0 ; i < v_w.size() ; ++i )
  f_x[ v_orig[ i ] ] = v_comp[ i ] ? double( 1 - in[ i ] ) : double( in[ i ] );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::solve_integer_core(
                                                    std::vector< char > & in )
{
 // the extracted integer core, solved by the full break-item core enumeration;
 // when warm started, the repaired previous solution floors the incumbent and
 // is itself the answer if nothing strictly beats it
 if( ! ( f_reopt && f_prev_valid ) )
  return( core_enumerate( v_w , v_p , f_C , in , - Inf< double >() , false ) );

 std::vector< char > win;
 const double lb = warm_incumbent( win );
 const double z = core_enumerate( v_w , v_p , f_C , in , lb , false );
 if( ( z == - Inf< double >() ) && ( lb > - Inf< double >() ) ) {
  in.swap( win );
  return( lb );
  }
 return( z );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::last_still_optimal(
                                            std::vector< char > & in ) const
{
 // the last solution y is optimal for the last core; with the same items, a
 // capacity no larger and no weight of a taken item decreased, y is feasible
 // if it still fits, and every new feasible z that takes no untaken item
 // whose weight decreased is feasible for the last core as well. For such a
 // z, p'z - p'y = ( pz - py ) + sum_k ( p'_k - p_k )( z_k - y_k ), where the
 // first term is <= 0 and so is every term of the sum but those of the items
 // with y_k = 0 and p'_k > p_k, or y_k = 1 and p'_k < p_k, with z_k != y_k.
 // Hence, a z beating y flips some suspect item: one of those, or an untaken
 // item whose weight decreased. With none, y is optimal; with intReopt 3, it
 // is also when, for each suspect item k, a Lagrangian bound of the new core
 // with x_k = 1 - y_k cannot beat p'y, and with intReopt 4 when either that
 // bound or the Martello-Toth one cannot
 const std::size_t m = v_w.size();
 if( ( f_C > f_last_C ) || ( m != v_last_w.size() ) ||
     ( v_last_in.size() != m ) ) {
  f_lambda = -1;
  return( - Inf< double >() );
  }
 if( f_C != f_last_C )         // the break item moves with the capacity
  f_lambda = -1;

 double z = 0;
 long wy = 0;
 bool wchg = false;
 std::vector< std::size_t > sus;
 for( std::size_t k = 0 ; k < m ; ++k ) {
  if( ( v_orig[ k ] != v_last_orig[ k ] ) ||
      ( v_comp[ k ] != v_last_comp[ k ] ) ) {
   f_lambda = -1;
   return( - Inf< double >() );
   }
  if( v_w[ k ] != v_last_w[ k ] )
   wchg = true;
  if( v_last_in[ k ] ) {
   if( v_w[ k ] < v_last_w[ k ] ) {
    f_lambda = -1;
    return( - Inf< double >() );
    }
   if( v_p[ k ] < v_last_p[ k ] )
    sus.push_back( k );
   z += v_p[ k ];
   wy += v_w[ k ];
   }
  else
   if( ( v_p[ k ] > v_last_p[ k ] ) || ( v_w[ k ] < v_last_w[ k ] ) )
    sus.push_back( k );
  }
 if( wchg )                     // the break item moves with the weights
  f_lambda = -1;
 if( wy > f_C )                 // y does not fit the new data
  return( - Inf< double >() );

 bool by_mt = false;           // some suspect item cleared by Martello-Toth
 if( ! sus.empty() ) {
  if( f_reopt < 3 )
   return( - Inf< double >() );
  if( f_cert_wait > 0 ) {        // backing off after consecutive failures
   --f_cert_wait;
   return( - Inf< double >() );
   }
#if CORE_STATS
  fprintf( stderr , "CERT try sus=%zu m=%zu\n" , sus.size() , m );
#endif

  // the multiplier lambda of the Lagrangian relaxation of the capacity, any
  // lambda >= 0 giving the valid bound U = lambda C + sum_j max( 0 , r_j ),
  // r_j = p_j - lambda w_j, and fixing x_k to v costing max( 0 , r_k ) -
  // v r_k more. The best lambda is the efficiency of the break item, found
  // by quickselect (as in core_enumerate()) on the first check and kept
  // while the capacity and the items stay: a stale one is still valid
  if( f_lambda < 0 ) {
   auto more_eff = [ this ]( std::size_t a , std::size_t b ) {
    return( v_p[ a ] * double( v_w[ b ] ) > v_p[ b ] * double( v_w[ a ] ) );
    };
   std::vector< std::size_t > ord( m );
   std::iota( ord.begin() , ord.end() , 0 );
   std::size_t lo = 0 , hi = m;
   long c = f_C;
   while( hi - lo > 16 ) {
    const std::size_t mid = lo + ( hi - lo ) / 2;
    std::nth_element( ord.begin() + lo , ord.begin() + mid ,
                      ord.begin() + hi , more_eff );
    long sw = 0;
    for( std::size_t t = lo ; t < mid ; ++t )
     sw += v_w[ ord[ t ] ];
    if( sw <= c ) { c -= sw; lo = mid; }
    else hi = mid;
    }
   std::sort( ord.begin() + lo , ord.begin() + hi , more_eff );
   while( ( lo < hi ) && ( v_w[ ord[ lo ] ] <= c ) )
    c -= v_w[ ord[ lo++ ] ];
   f_lambda = lo < m ? v_p[ ord[ lo ] ] / double( v_w[ ord[ lo ] ] ) : 0;
   }

  bool intp = true;
  double u = f_lambda * double( f_C );
  for( std::size_t k = 0 ; k < m ; ++k ) {
   if( v_p[ k ] < 0 )            // (the bound needs p >= 0)
    return( - Inf< double >() );
   if( intp && ( v_p[ k ] != double( long( v_p[ k ] ) ) ) )
    intp = false;
   u += std::max( 0.0 , v_p[ k ] - f_lambda * double( v_w[ k ] ) );
   }

  const double tol = 1e-12 * std::max( 1.0 , std::abs( z ) );
  auto beats = [ & ]( double b ) {
   return( intp ? ( floor_safe( b ) > z ) : ( b > z + tol ) );
   };

  // with intReopt 4, a suspect item k that the Lagrangian bound does not
  // clear gets the Martello-Toth bound with x_k = 1 - y_k: the larger of
  // the two continuous bounds with the critical item also fixed to 0 and
  // to 1, which, unlike any continuous bound, sees that an item (say, one
  // as large as the capacity) is either taken whole or not at all. The
  // items are sorted by efficiency on the first such item only
  std::vector< std::size_t > ord;
  auto lp = [ & ]( std::size_t a , int va , std::size_t b , int vb ,
                   std::size_t & crit ) {
   double cap = double( f_C ) - va * double( v_w[ a ] );
   double val = va * v_p[ a ];
   if( b < m ) {
    cap -= vb * double( v_w[ b ] );
    val += vb * v_p[ b ];
    }
   crit = m;
   if( cap < 0 )
    return( - Inf< double >() );
   for( auto j : ord ) {
    if( ( j == a ) || ( j == b ) )
     continue;
    if( double( v_w[ j ] ) <= cap ) {
     cap -= double( v_w[ j ] );
     val += v_p[ j ];
     }
    else {
     val += v_p[ j ] * cap / double( v_w[ j ] );
     crit = j;
     break;
     }
    }
   return( val );
   };
  auto mt_bound = [ & ]( std::size_t k ) {
   if( ord.empty() ) {
    ord.resize( m );
    std::iota( ord.begin() , ord.end() , 0 );
    std::sort( ord.begin() , ord.end() ,
               [ this ]( std::size_t a , std::size_t b ) {
                return( v_p[ a ] * double( v_w[ b ] ) >
                        v_p[ b ] * double( v_w[ a ] ) );
                } );
    }
   const int v = v_last_in[ k ] ? 0 : 1;
   std::size_t c , c2;
   const double l = lp( k , v , m , 0 , c );
   if( c == m )                // the continuous solution is integer
    return( l );
   return( std::max( lp( k , v , c , 0 , c2 ) , lp( k , v , c , 1 , c2 ) ) );
   };

  for( auto k : sus ) {
   const double r = v_p[ k ] - f_lambda * double( v_w[ k ] );
   const double uk = u - ( v_last_in[ k ] ? std::max( 0.0 , r )
                                          : std::max( 0.0 , r ) - r );
   if( beats( uk ) && ( f_reopt >= 4 ) && ( ! beats( mt_bound( k ) ) ) ) {
    by_mt = true;
    continue;
    }
   if( beats( uk ) ) {
#if CORE_STATS
    fprintf( stderr , "CERT fail sus=%zu m=%zu\n" , sus.size() , m );
#endif
    // from the 4th failure in a row, the next 1, 3, 7, ..., 127 checks are
    // skipped: where it never succeeds, its O(m) cost goes away
    if( ++f_cert_fail >= 4 )
     f_cert_wait = ( 1 << std::min( f_cert_fail - 3 , 7 ) ) - 1;
    return( - Inf< double >() );
    }
   }
  }

#if CORE_STATS
 fprintf( stderr , "CERT ok sus=%zu m=%zu\n" , sus.size() , m );
#endif
 f_cert_fail = 0;
 f_outcome = sus.empty() ? 2 : ( by_mt ? 4 : 3 );
 in = v_last_in;
 return( z );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::warm_incumbent( std::vector< char > & in )
 const
{
 // the items outside the core are decided by the fixings and the sign-based
 // pre-fixing of the current data, so only the core image of the previous
 // solution matters: y = x, or y = 1 - x for a complemented item, rounded
 // down (a continuous variable may have become integer)
 const std::size_t m = v_w.size();
 in.assign( m , 0 );
 if( v_prev_x.size() != f_N )
  return( - Inf< double >() );

 long wsum = 0;
 double psum = 0;
 for( std::size_t k = 0 ; k < m ; ++k ) {
  const double x = v_prev_x[ v_orig[ k ] ];
  const double y = v_comp[ k ] ? 1.0 - x : x;
  if( y > 1.0 - 1e-9 ) {
   in[ k ] = 1;
   wsum += v_w[ k ];
   psum += v_p[ k ];
   }
  }

 // repair: drop the least efficient taken items until the capacity holds
 if( wsum > f_C ) {
  std::vector< std::size_t > tk;
  for( std::size_t k = 0 ; k < m ; ++k )
   if( in[ k ] )
    tk.push_back( k );
  std::sort( tk.begin() , tk.end() ,
             [ this ]( std::size_t a , std::size_t b ) {
              return( v_p[ a ] * double( v_w[ b ] ) <
                      v_p[ b ] * double( v_w[ a ] ) );
              } );
  for( auto k : tk ) {
   if( wsum <= f_C )
    break;
   in[ k ] = 0;
   wsum -= v_w[ k ];
   psum -= v_p[ k ];
   }
  }

 return( psum );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::core_enumerate(
                       const std::vector< long > & iw ,
                       const std::vector< double > & ip , long C ,
                       std::vector< char > & in , double lb , bool relx )
{
 const int m = int( iw.size() );
 if( ( ! f_lazy_core ) || ( m < 64 ) || relx )   // (not in the surrogate)
  return( core_enumerate_full( iw , ip , C , in , lb , relx ) );

 // the break item by quickselect on the efficiencies (Balas-Zemel): after the
 // halving, key[ 0 , b ) are the b most efficient items, which all fit, and
 // key[ b ] the most efficient of the others (the break item)
 struct Key { double e; long w; double p; int i; };
 std::vector< Key > key( m );
 double emin = Inf< double >() , emax = - Inf< double >();
 for( int i = 0 ; i < m ; ++i ) {
  key[ i ] = Key{ ip[ i ] / double( iw[ i ] ) , iw[ i ] , ip[ i ] , i };
  emin = std::min( emin , key[ i ].e );
  emax = std::max( emax , key[ i ].e );
  }
 if( emin == emax )            // all equally efficient: nothing can be fixed
  return( core_enumerate_full( iw , ip , C , in , lb , relx ) );
 const auto before = []( const Key & a , const Key & c ) {
  return( a.e > c.e );
  };
 long rc = C;
 int lo = 0 , hi = m;
 while( hi - lo > 16 ) {
  const int mid = lo + ( hi - lo ) / 2;
  std::nth_element( key.begin() + lo , key.begin() + mid , key.begin() + hi ,
                    before );
  long wh = 0;
  for( int i = lo ; i < mid ; ++i )
   wh += key[ i ].w;
  if( wh <= rc ) { rc -= wh; lo = mid; }
  else hi = mid;
  }
 std::sort( key.begin() + lo , key.begin() + hi , before );
 while( ( lo < hi ) && ( key[ lo ].w <= rc ) )
  rc -= key[ lo++ ].w;
 const int b = lo;

 in.assign( m , 0 );
 long wsumb = 0;
 double psumb = 0;
 for( int i = 0 ; i < b ; ++i )
  { wsumb += key[ i ].w; psumb += key[ i ].p; }
 if( b >= m ) {                // everything fits
  for( int i = 0 ; i < m ; ++i )
   in[ i ] = 1;
  return( psumb > lb ? psumb : - Inf< double >() );
  }

 // the break efficiency, the Dantzig bound and the pruning margin of an
 // improving solution (integer profits: it must gain at least 1), before
 // the selections below move the break item from its place
 bool intp = true;
 for( int i = 0 ; ( i < m ) && intp ; ++i )
  if( key[ i ].p != std::floor( key[ i ].p ) )
   intp = false;
 const double eb = key[ b ].e;
 const double base = psumb + double( C - wsumb ) * eb;
 const double ztol = 1e-12 * std::max( 1.0 , std::abs( base ) );
 const double margin = intp ? ( 1 - ztol ) : ztol;

 // stage 1, a primal heuristic: the small core of the (up to) 32 least
 // efficient items of the break solution and the 32 most efficient of the
 // others, key[ b - kl , b + kr ), solved exactly with every other item at
 // its break value. Its optimum is feasible for the whole instance and
 // usually close to the optimum: it floors both the fixing and the
 // enumeration of stage 2
 const int kl = std::min( 32 , b ) , kr = std::min( 32 , m - b );
 if( kl < b )
  std::nth_element( key.begin() , key.begin() + ( b - kl ) ,
                    key.begin() + b , before );
 if( kr < m - b )
  std::nth_element( key.begin() + b , key.begin() + ( b + kr ) , key.end() ,
                    before );
 const int s0 = b - kl , s1 = b + kr;
 long offw = 0;
 double offp = 0;
 for( int i = 0 ; i < s0 ; ++i )
  { offw += key[ i ].w; offp += key[ i ].p; }
 std::vector< long > cw( s1 - s0 );
 std::vector< double > cp( s1 - s0 );
 for( int i = s0 ; i < s1 ; ++i )
  { cw[ i - s0 ] = key[ i ].w; cp[ i - s0 ] = key[ i ].p; }
 // (on a budget of 4 m states: on the hardest instances even this small
 // core is hard, and then its best solution found is taken)
 std::vector< char > cin;
 const long cap0 = f_state_cap;
 f_state_cap = 4 * long( m );
 const double z1 = offp + core_enumerate_full( cw , cp , C - offw , cin ,
                                               - Inf< double >() , relx );
 f_state_cap = cap0;
 // the stage-1 solution, on the input indices: the answer if nothing better
 auto stage1 = [ & ]( void ) -> double {
  if( z1 <= lb )
   return( - Inf< double >() );
  in.assign( m , 0 );
  for( int i = 0 ; i < s0 ; ++i )
   in[ key[ i ].i ] = 1;
  for( int i = s0 ; i < s1 ; ++i )
   in[ key[ i ].i ] = cin[ i - s0 ];
  return( z1 );
  };
 const double z0 = std::max( z1 , lb );

 // stage 2, the fixed core: the stage-1 solution is feasible for the
 // problem restricted to any core containing the small one, so the core
 // optimum is at least z0; an item whose Dembo-Hammer bound (the LP bound
 // with the item flipped, linearised at the break efficiency) stays below z0
 // plus the improving margin can then be fixed at its break value, and the
 // core is made of all the others and of the small core. A core larger than
 // half of the instance saves nothing: the whole of it is solved instead
 // (and as soon as this is clear), still floored at z0
 std::vector< int > core;
 offw = 0;
 offp = 0;
 for( int i = 0 ; i < m ; ++i ) {
  const double U = ( i < b ) ? base - key[ i ].p + double( key[ i ].w ) * eb
                             : base + key[ i ].p - double( key[ i ].w ) * eb;
  if( ( U >= z0 + margin ) || ( ( i >= s0 ) && ( i < s1 ) ) ) {
   core.push_back( i );
   if( 2 * core.size() > std::size_t( m ) ) {
    std::vector< char > fin;
    const double zf = core_enumerate_full( iw , ip , C , fin , z0 , relx );
    if( zf == - Inf< double >() )
     return( stage1() );
    in.swap( fin );
    return( zf );
    }
   }
  else
   if( i < b )
    { offw += key[ i ].w; offp += key[ i ].p; }
  }

 const int mc = int( core.size() );
 cw.resize( mc );
 cp.resize( mc );
 for( int k = 0 ; k < mc ; ++k )
  { cw[ k ] = key[ core[ k ] ].w; cp[ k ] = key[ core[ k ] ].p; }
 std::vector< char > fin;
 const double zc = core_enumerate_full( cw , cp , C - offw , fin , z0 - offp ,
                                        relx );
 if( zc == - Inf< double >() )
  return( stage1() );
 in.assign( m , 0 );
 for( int i = 0 ; i < b ; ++i )         // fixed break items, the core ones
  in[ key[ i ].i ] = 1;                 // overwritten below
 for( int k = 0 ; k < mc ; ++k )
  in[ key[ core[ k ] ].i ] = fin[ k ];
 return( offp + zc );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::core_enumerate_full(
                       const std::vector< long > & iw ,
                       const std::vector< double > & ip , long C ,
                       std::vector< char > & in , double lb , bool relx )
{
 // integer state profits when all the profits are integers whose sum a
 // double holds exactly
 double ps = 0;
 bool intp = true;
 for( std::size_t i = 0 ; ( i < ip.size() ) && intp ; ++i ) {
  intp = ( ip[ i ] == std::floor( ip[ i ] ) );
  ps += std::abs( ip[ i ] );
  }
 if( intp && ( ps < 9e15 ) )
  return( core_enumerate_engine< long >( iw , ip , C , in , lb , relx ) );
 return( core_enumerate_engine< double >( iw , ip , C , in , lb , relx ) );
 }

/*--------------------------------------------------------------------------*/

template< class PT >
double CoreDPBinaryKnapsackSolver::core_enumerate_engine(
                       const std::vector< long > & iw ,
                       const std::vector< double > & ip , long C ,
                       std::vector< char > & in , double lb , bool relx )
{
 // MINKNAP-style enumeration. The items are sorted by efficiency p/w; the
 // "break solution" takes the most efficient ones until the capacity is (about
 // to be) exceeded. We then enumerate only the deviations from that break
 // solution, expanding the core outwards (add the next right item / drop the
 // next left item, alternately), keeping the Pareto frontier of (weight,
 // profit) and pruning every state whose LP upper bound cannot beat the
 // incumbent, which is floored at the caller's @p lb. Identical items are
 // aggregated into one item with a multiplicity, an item that doubles another
 // ( w' , p' ) = ( 2w , 2p ) is folded into it as two further copies, and the
 // copies are expanded in binary decomposition ( 1 , 2 , 4 , ... ) with early
 // termination as soon as one expansion yields no new state. Returns the
 // optimum and fills @p in (in the items' input order), or returns -Inf when
 // no solution strictly beats @p lb; @p relx marks the (surrogate) recursive
 // call, which must not recurse further.

 const int m = int( iw.size() );
#if CORE_STATS
 if( ! relx )
  cs_merge = cs_nogrow = cs_cand = cs_b2 = cs_b3 = cs_b4 = 0;
#endif
 in.assign( m , 0 );
 if( m == 0 )                  // no free item to decide
  return( 0 > lb ? 0.0 : - Inf< double >() );

 // aggregation: identical ( w , p ) items collapse into one registry entry
 // with a multiplicity. Duplicates are found by grouping a ( w , p )-sorted
 // permutation - no hashing and no per-item containers, the group's input
 // indices are the byk[ mstart , mstart + mcount ) span
 // (with small weights, a counting sort on them and a sort by profit inside
 // each weight, which then has few items, in place of the full sort)
 std::vector< int > byk( m );
 const auto byp = [ & ]( int a , int b ) { return( ip[ a ] < ip[ b ] ); };
 const auto [ wlo , whi ] = std::minmax_element( iw.begin() , iw.end() );
 const long wspan = *whi - *wlo + 1;
 if( wspan <= 4 * long( m ) ) {
  std::vector< int > start( wspan + 1 , 0 );
  for( int i = 0 ; i < m ; ++i )
   ++start[ iw[ i ] - *wlo + 1 ];
  for( long v = 0 ; v < wspan ; ++v )
   start[ v + 1 ] += start[ v ];
  for( int i = 0 ; i < m ; ++i )
   byk[ start[ iw[ i ] - *wlo ]++ ] = i;
  // start[ v ] is now the end of the bucket of weight *wlo + v
  for( long v = 0 , b0 = 0 ; v < wspan ; b0 = start[ v++ ] )
   if( start[ v ] - b0 > 1 )
    std::sort( byk.begin() + b0 , byk.begin() + start[ v ] , byp );
  }
 else {
  for( int i = 0 ; i < m ; ++i )
   byk[ i ] = i;
  std::sort( byk.begin() , byk.end() , [ & ]( int a , int b ) {
   return( ( iw[ a ] != iw[ b ] ) ? ( iw[ a ] < iw[ b ] ) : byp( a , b ) );
   } );
  }

 struct Agg {
  long w;
  double p;
  long d;                      ///< multiplicity (own units + 2x each child's)
  int child;                   ///< registry index of the folded double, or -1
  int mstart;                  ///< first own unit in byk
  int mcount;                  ///< number of own units
  };
 std::vector< Agg > reg;      // ( w , p )-ascending by construction
 reg.reserve( m );
 for( int i = 0 ; i < m ; ) {
  int j = i;
  while( ( j < m ) && ( iw[ byk[ j ] ] == iw[ byk[ i ] ] ) &&
         ( ip[ byk[ j ] ] == ip[ byk[ i ] ] ) )
   ++j;
  reg.push_back( Agg{ iw[ byk[ i ] ] , ip[ byk[ i ] ] , long( j - i ) , -1 ,
                      i , j - i } );
  i = j;
  }

 // multiplicity reduction: in decreasing weight order, an item whose exact
 // half ( w/2 , p/2 ) exists (binary search in the sorted registry) is
 // removed and contributes two copies per unit to the half - chains
 // 8 -> 4 -> 2 unwind transitively; the doubling halts once the expanded
 // unit count would blow past 2m (degenerate deep chains)
 {
  auto find_reg = [ & ]( long ww , double pp ) -> int {
   int lo = 0 , hi = int( reg.size() ) - 1;
   while( lo <= hi ) {
    const int mid = ( lo + hi ) / 2;
    if( ( reg[ mid ].w < ww ) ||
        ( ( reg[ mid ].w == ww ) && ( reg[ mid ].p < pp ) ) )
     lo = mid + 1;
    else if( ( reg[ mid ].w == ww ) && ( reg[ mid ].p == pp ) )
     return( mid );
    else
     hi = mid - 1;
    }
   return( -1 );
   };
  long units = m;
  for( int j = int( reg.size() ) - 1 ; j >= 0 ; --j ) {
   if( reg[ j ].w % 2 )
    continue;
   const int h = find_reg( reg[ j ].w / 2 , reg[ j ].p / 2 );
   if( h < 0 )
    continue;
   if( units + reg[ j ].d > 2 * long( m ) )
    break;
   units += reg[ j ].d;        // 2 half-copies replace each unit (net + d)
   reg[ h ].d += 2 * reg[ j ].d;
   reg[ h ].child = j;
   reg[ j ].d = 0;             // folded away: inactive
   }
  }

 // active aggregated items, sorted by non-increasing efficiency p/w
 std::vector< int > regIdx;
 regIdx.reserve( reg.size() );
 for( std::size_t j = 0 ; j < reg.size() ; ++j )
  if( reg[ j ].d > 0 )
   regIdx.push_back( int( j ) );
 std::sort( regIdx.begin() , regIdx.end() , [ & ]( int a , int b ) {
  return( reg[ a ].p * double( reg[ b ].w ) > reg[ b ].p * double( reg[ a ].w ) );
  } );
 const int mA = int( regIdx.size() );
 std::vector< long > w( mA ) , d( mA );
 std::vector< double > p( mA );
 long mEl = 0;
 for( int j = 0 ; j < mA ; ++j ) {
  w[ j ] = reg[ regIdx[ j ] ].w;
  p[ j ] = reg[ regIdx[ j ] ].p;
  d[ j ] = reg[ regIdx[ j ] ].d;
  mEl += d[ j ];
  }

 // expanded per-unit view (copies adjacent, efficiency order preserved): the
 // break / U2 / initial heuristic / surrogate machinery is unit-based; with
 // no multiplicity anywhere (the common case) the unit view coincides with
 // the aggregated one and is aliased rather than copied
 const int mE = int( mEl );
 const bool noagg = ( mE == mA );
 std::vector< long > wEb;
 std::vector< double > pEb;
 std::vector< int > unit_item;
 if( ! noagg ) {
  wEb.resize( mE ); pEb.resize( mE ); unit_item.resize( mE );
  for( int j = 0 , u = 0 ; j < mA ; ++j )
   for( long c = 0 ; c < d[ j ] ; ++c , ++u ) {
    wEb[ u ] = w[ j ];
    pEb[ u ] = p[ j ];
    unit_item[ u ] = j;
    }
  }
 const std::vector< long > & wE = noagg ? w : wEb;
 const std::vector< double > & pE = noagg ? p : pEb;
 // unit -> aggregated item index
 auto uitem = [ & ]( int u ) { return( noagg ? u : unit_item[ u ] ); };

 // integer profits enable the tighter "cannot reach z + 1" pruning margin;
 // more generally every attainable value is a multiple of g_p = gcd of the
 // (integer) profits, so an improving solution must reach z + g_p: profit
 // ceiling instances ( p_i = d ceil( w_i / d ) ) have g_p = d, which both
 // multiplies the pruning margin and snaps the global ceiling down to the
 // attainable grid
 bool intp = true;
 for( int i = 0 ; ( i < mA ) && intp ; ++i )
  if( std::abs( p[ i ] - std::round( p[ i ] ) ) > 1e-9 )
   intp = false;
 long gp = 1;
 if( intp ) {
  gp = std::lround( p[ 0 ] );
  for( int i = 1 ; ( i < mA ) && ( gp > 1 ) ; ++i )
   gp = std::gcd( gp , std::lround( p[ i ] ) );
  if( gp < 1 )
   gp = 1;
  }
 // the margin allows for the rounding errors of the bounds, which are of
 // the order of the machine precision times the value of the solutions: a
 // relative tolerance on an upper bound of the latter, the least of the sum
 // of the profits and of the capacity filled at the highest efficiency
 // (items by non-increasing efficiency: the first one)
 double ptot = 0;
 for( int i = 0 ; i < mA ; ++i )
  ptot += std::abs( p[ i ] ) * double( d[ i ] );
 if( mA > 0 )
  ptot = std::min( ptot ,
                   double( C ) * std::abs( p[ 0 ] ) / double( w[ 0 ] ) );
 const double ztol = 1e-12 * std::max( 1.0 , ptot );
 const double zmargin = intp ? ( double( gp ) - ztol ) : ztol;
 auto snap = [ intp , gp ]( double ub ) -> double {
  return( ( intp && ( gp > 1 ) )
          ? floor_safe( ub / double( gp ) ) * double( gp ) : ub );
  };

 // A3: divisibility-reduced effective capacity (equals C when gcd == 1); all
 // state weights are multiples of the gcd, so it is exact for feasibility and
 // only tightens the bounds
 const long Ceff = divisibility_capacity( w , mA , C );

 // break unit bE: units [ 0 , bE ) form the (feasible) break solution, taking
 // e copies of the (partial) break item bA
 long wsumb = 0;
 double psumb = 0;
 int bE = 0;
 while( ( bE < mE ) && ( wsumb + wE[ bE ] <= Ceff ) ) {
  wsumb += wE[ bE ];
  psumb += pE[ bE ];
  ++bE;
  }
 const int bA = ( bE < mE ) ? uitem( bE ) : mA;
 std::vector< long > cntb( mA , 0 );    // per-item copies in the break solution
 for( int u = 0 ; u < bE ; ++u )
  ++cntb[ uitem( u ) ];

 // completion data: prem[ j ] is the profit of the copies of the items
 // j , j + 1 , ... still out of the break solution, wrem[ j + 1 ] the weight
 // of those of the items 0 , ... , j still in it; a state can gain at most
 // the former from the unprocessed right items, and can shed at most the
 // latter by dropping the unprocessed left ones (fixed items included: the
 // bounds are only looser for it)
 std::vector< double > prem( mA + 1 , 0 );
 std::vector< long > wrem( mA + 1 , 0 );
 for( int j = mA - 1 ; j >= 0 ; --j )
  prem[ j ] = prem[ j + 1 ] + p[ j ] * double( d[ j ] - cntb[ j ] );
 for( int j = 0 ; j < mA ; ++j )
  wrem[ j + 1 ] = wrem[ j ] + w[ j ] * cntb[ j ];

 // A2: global upper bound used as an early-exit ceiling. The Martello-Toth U2
 // bound (break item wholly in or out) is provably <= the plain Dantzig bound,
 // so it is always valid and usually tighter.
 double dantzig =
  snap( upper_bound_u2( wE , pE , mE , bE , wsumb , psumb , Ceff , intp ) );

 // A1: seed the incumbent with a strong initial heuristic (>= break solution),
 // so the Dantzig pruning below starts from a higher z; the caller's lb floors
 // the pruning level, and `beat` tracks whether any solution actually beat it
 // (only then are z / best_in meaningful for the caller). The unit-based
 // heuristic selection is folded into per-item copy counts
 std::vector< long > best_in( mA , 0 );
 std::vector< char > selE( mE , 0 );
 const double z0 = initial_heuristic( wE , pE , mE , bE , wsumb , psumb , Ceff ,
                                      selE );
 for( int u = 0 ; u < mE ; ++u )
  if( selE[ u ] )
   ++best_in[ uitem( u ) ];
 bool beat = ( z0 > lb );
 double z = std::max( z0 , lb );
#if CORE_STATS
 const double cs_z0 = z;
 std::size_t cs_zstep = 0 , cs_zimp = 0;   // last improvement, improvements
#endif


 // the surrogate / cardinality bound is computed lazily, only once the Pareto
 // frontier grows past intSurrTrigger (i.e. the instance is actually hard,
 // COMBO's MINSET gate): this keeps it off the easy instances where it would
 // only add sorting overhead. It uses the (by then larger) incumbent z for
 // the N_min side and for the exact-solve pruning, so one retry is allowed
 // when the incumbent has improved after the first attempt (a stronger z can
 // turn a truncated solve into a certifying one). A solve that ran out of
 // its budget (bit 1 of intSurrAdapt) is retried, without counting as a try,
 // once the enumeration has doubled the states it had generated at the time
 // (surr_redo): each retry has a budget twice as large, so a certifying solve
 // is eventually let through at a geometrically bounded total cost
 int surr_tries = 2;
 double surr_z = - Inf< double >();
 std::size_t surr_redo = std::numeric_limits< std::size_t >::max();

 // B4 (WB item-fixing): the break efficiency is the fill rate of the
 // Dembo-Hammer reduction bound WB( i , 1 )
 const double eb = ( bE < mE ) ? ( pE[ bE ] / double( wE[ bE ] ) ) : 0.0;

 // states are deviations from the break solution: unprocessed left copies are
 // implicitly in, unprocessed right copies implicitly out (so the initial
 // weight is wsumb). `front` is the current Pareto frontier, sorted by
 // strictly ascending weight, and `cur` the buffer the next one is merged
 // into (the two are swapped, never reallocated once large enough);
 // step_item[ k ] / step_mult[ k ] are the item toggled by the k-th stored
 // step and its signed copy count (binary-decomposed batch, < 0 = dropped
 // left copies). For the backtracking each state carries the decisions of
 // the steps since the last checkpoint, one bit each (bit k - C - 1 for the
 // step k > C, C the largest multiple of 64 below k), and the frontier of
 // every checkpoint C > 0 is kept as the step C + 1 read it (chk[ C / 64 -
 // 1 ]): from a state, its mask gives the last decisions, undoing them the
 // weight of its ancestor at the checkpoint, which the weight identifies in
 // that frontier, and so on back to the break state. This replaces a
 // per-state history of every step with one copy of the frontier every 64
 struct WP { long wsum; PT psum; std::uint64_t mask; };
 using WPvec = std::vector< WP , NoInit< WP > >;   // resize() writes nothing
 WPvec front( 1 , WP{ wsumb , PT( psumb ) , 0 } );
 std::vector< PT > pt( mA );   // the item profits, as state profits
 for( int j = 0 ; j < mA ; ++j )
  pt[ j ] = PT( p[ j ] );
 WPvec cur;
 std::vector< WPvec > chk;
 std::vector< int > step_item( 1 , -1 );
 std::vector< long > step_mult( 1 , 0 );
 std::size_t nstates = 1;      // states of all the stored steps

 // record an improving incumbent: rebuild the achieving per-item copy counts
 // by walking the parent chain through the (already final) earlier steps. The
 // state is given as ( pos , item / delta , nx extra unit toggles in xs )
 // where pos indexes the LAST stored step: a merge candidate passes its parent
 // there plus its own batch as ( item , signed copies ), a pairing-heuristic
 // match the state itself plus the paired item(s) in xs, each entry signed:
 // +( j + 1 ) adds one copy of item j, -( j + 1 ) drops one
 auto record_incumbent = [ & ]( double psum , int pos , int item , long delta ,
                                const int * xs , int nx ) {
  z = psum;
  beat = true;
  CSTAT( ( cs_zstep = nstates , ++cs_zimp ) );
  best_in = cntb;
  if( delta )
   best_in[ item ] += delta;
  for( int e = 0 ; e < nx ; ++e ) {
   const int j = std::abs( xs[ e ] ) - 1;
   best_in[ j ] += ( xs[ e ] > 0 ) ? 1 : -1;
   }
  long sw = front[ pos ].wsum;
  std::uint64_t mask = front[ pos ].mask;
  for( std::size_t S = step_item.size() - 1 ; S > 0 ; ) {
   const std::size_t C = ( ( S - 1 ) / 64 ) * 64;   // checkpoint below S
   for( std::size_t k = S ; k > C ; --k )
    if( ( mask >> ( k - C - 1 ) ) & 1u ) {
     best_in[ step_item[ k ] ] += step_mult[ k ];
     sw -= step_mult[ k ] * w[ step_item[ k ] ];
     }
   if( C == 0 )
    break;
   const auto & f = chk[ C / 64 - 1 ];    // the ancestor, by its weight
   const auto at = std::lower_bound( f.begin() , f.end() , sw ,
                                     []( const WP & a , long v ) {
                                      return( a.wsum < v );
                                      } );
   assert( ( at != f.end() ) && ( at->wsum == sw ) );
   mask = at->mask;
   S = C;
   }
  };

 // B2 (guarded DP-extension) per-side state: consecutive expansions that
 // produced no new state, the read-only-check trigger threshold (doubled
 // whenever a check turns out wasteful) and the consecutive read-only skips
 // (an expansion is forced after 40, to let the state pruning operate)
 int nonewR = 0 , nonewL = 0;
 int guardR = 10 , guardL = 10;
 int roSkipR = 0 , roSkipL = 0;

 // B3 (fixing-by-dominance) per-side reference: the last item whose expansion
 // generated no new state (valid when the profit is >= 0); with
 // intDominanceFix = 2 every such item instead marks at once all the still
 // unprocessed items of its side that it dominates
 long   domRw = 0 , domLw = 0;
 double domRp = -1 , domLp = -1;
 std::vector< char > domfix( f_dom_fix > 1 ? mA : 0 , 0 );

 // B5 (pairing heuristics) state: the PH/TPH/SSPH trigger is a threshold
 // on the states generated since the last call, doubled per call (so that
 // they run early, while the incumbent decides how much is pruned); SPH runs
 // after every expansion and uses a deterministic xorshift32 (reproducible
 // solves) to sample one item per block
 std::size_t ph_thresh = 10 * std::size_t( m );
 std::size_t ph_work = 0;
 unsigned rng = 0x9E3779B9u;
 auto rnd = [ & ]( unsigned k ) -> unsigned {
  rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
  return( rng % k );
  };

 // best frontier state under a weight cap: the frontier is sorted by weight
 // with profit increasing along it, so it is the last state below the cap
 auto best_under = [ & ]( const WPvec & f , long wcap ) -> int {
  int lo = 0 , hi = int( f.size() ) - 1 , res = -1;
  while( lo <= hi ) {
   const int mid = ( lo + hi ) / 2;
   if( f[ mid ].wsum <= wcap ) { res = mid; lo = mid + 1; }
   else hi = mid - 1;
   }
  return( res );
  };

 int s = std::min( bA , mA - 1 );  // next left item to (consider to) drop from
 int t = bA;                       // next right item to (consider to) add to
 while( ( s >= 0 ) && ( cntb[ s ] == 0 ) )
  --s;
 while( ( t < mA ) && ( d[ t ] - cntb[ t ] == 0 ) )
  ++t;
 bool goRight = true;
 bool dead = false;            // the frontier emptied: certified, stop all

 // the items fixed at their break value before being reached: fixd[ j ] is
 // set by the reduction below as the incumbent grows, or by dominance;
 // wfixU / pfixU are the weight (left items) and the profit (right items)
 // of those still to be reached, which the completion bounds must not count
 std::vector< char > fixd( mA , 0 );
 long wfixU = 0;
 double pfixU = 0;
 auto fix_item = [ & ]( int j , bool right ) {
  if( fixd[ j ] )
   return;
  fixd[ j ] = 1;
  if( right )
   pfixU += p[ j ] * double( d[ j ] - cntb[ j ] );
  else
   wfixU += w[ j ] * cntb[ j ];
  };

 // reduction: whenever the incumbent has grown, every item still to be
 // reached whose Dembo-Hammer bound (see B4 below) cannot beat it is fixed
 // at once, rather than when its turn comes
 double z_reduced = - Inf< double >();

 // the check of the dominance after an expansion that grew the frontier
 // (see below) is amortised over 10 m generated states
 std::size_t dom_work = 0;
 const std::size_t dom_thresh = 10 * std::size_t( m );
 auto reduce = [ & ]( void ) {
  if( ( ! f_red_fix ) || ( z <= z_reduced ) )
   return;
  z_reduced = z;
  for( int j = s ; j >= 0 ; --j )
   if( cntb[ j ] && ( ! fixd[ j ] ) &&
       ( psumb - p[ j ] + double( Ceff - wsumb + w[ j ] ) * eb <
         z + zmargin ) )
    fix_item( j , false );
  for( int j = t ; j < mA ; ++j )
   if( ( d[ j ] > cntb[ j ] ) && ( ! fixd[ j ] ) &&
       ( psumb + p[ j ] + double( Ceff - wsumb - w[ j ] ) * eb <
         z + zmargin ) )
    fix_item( j , true );
  };

 // the frontier is closed under one more copy of `i` (and stays so under any
 // later expansion), so is it under any item that `i` dominates: mark all
 // those still to be processed on the side of `i`
 auto mark_dominated = [ & ]( int i , bool right ) {
  if( right ) {
   for( int j = t ; j < mA ; ++j )
    if( ( w[ j ] >= w[ i ] ) &&
        ( p[ j ] <= p[ i ] * double( w[ j ] / w[ i ] ) ) )
     { domfix[ j ] = 1; fix_item( j , true ); }
   }
  else
   for( int j = s ; j >= 0 ; --j )
    if( ( w[ j ] <= w[ i ] ) && ( p[ j ] >= p[ i ] ) )
     { domfix[ j ] = 1; fix_item( j , false ); }
  };
 while( ( ( s >= 0 ) || ( t < mA ) ) && ( z < dantzig - 1e-9 ) && ! dead ) {

  // pick the next item to expand, alternating right / left; its copies still
  // on the break side are the batch budget of the inner loop below
  int it;
  bool isRight;
  long avail;
  if( ( t < mA ) && ( goRight || ( s < 0 ) ) ) {
   it = t; isRight = true; avail = d[ it ] - cntb[ it ];
   do ++t; while( ( t < mA ) && ( d[ t ] - cntb[ t ] == 0 ) );
   }
  else {
   it = s; isRight = false; avail = cntb[ it ];
   do --s; while( ( s >= 0 ) && ( cntb[ s ] == 0 ) );
   }
  goRight = ! goRight;

  reduce();
  if( fixd[ it ] ) {           // fixed before its turn: no longer pending
   if( isRight )
    pfixU -= p[ it ] * double( avail );
   else
    wfixU -= w[ it ] * avail;
   CSTAT( ++cs_b4 );
   continue;
   }

  // B4: Dembo-Hammer reduction bound. WB( it , 1 ) bounds the value of any
  // solution deviating on `it` (the LP value function is concave with slope
  // <= eb past the break point); if it cannot beat the incumbent the item is
  // fixed at its break-side value and never expanded
  if( f_red_fix ) {
   const double WB = isRight
    ? psumb + p[ it ] + double( Ceff - wsumb - w[ it ] ) * eb
    : psumb - p[ it ] + double( Ceff - wsumb + w[ it ] ) * eb;
   if( WB < z + zmargin )
    { CSTAT( ++cs_b4 ); continue; }
   }

  // B3: fixing-by-dominance. once expanding item i gave no new state, any
  // heavier right item with p' <= p_i * floor( w' / w_i ) (resp. lighter left
  // item with p' >= p_i) cannot yield an improved solution either - floor(
  // w' / w_i ) copies of i dominate it state-wise - so skip it outright
  if( ( f_dom_fix > 1 ) && domfix[ it ] )
   { CSTAT( ++cs_b3 ); continue; }
  if( f_dom_fix == 1 ) {
   if( isRight ) {
    if( ( domRp >= 0 ) && ( w[ it ] >= domRw ) &&
        ( p[ it ] <= domRp * double( w[ it ] / domRw ) ) )
     { CSTAT( ++cs_b3 ); continue; }
    }
   else
    if( ( domLp >= 0 ) && ( w[ it ] <= domLw ) && ( p[ it ] >= domLp ) )
     { CSTAT( ++cs_b3 ); continue; }
   }

  // LP bound fill rates of the still-unprocessed units. While copies of the
  // current item remain uncommitted they are themselves the most efficient
  // addable (resp. the least efficient droppable) units, so until the last
  // batch the bound must use the current item's own efficiency
  // (the next items not fixed: the fixed ones will never move)
  int tn = t , sn = s;
  while( ( tn < mA ) && ( fixd[ tn ] || ( d[ tn ] == cntb[ tn ] ) ) )
   ++tn;
  while( ( sn >= 0 ) && ( fixd[ sn ] || ( cntb[ sn ] == 0 ) ) )
   --sn;
  const double eff_t  = ( tn < mA ) ? ( p[ tn ] / double( w[ tn ] ) ) : 0.0;
  const double eff_s  = ( sn >= 0 ) ? ( p[ sn ] / double( w[ sn ] ) ) : 0.0;
  const double eff_it = p[ it ] / double( w[ it ] );

  int & nonew  = isRight ? nonewR  : nonewL;
  int & guard  = isRight ? guardR  : guardL;
  int & roSkip = isRight ? roSkipR : roSkipL;

  // the item's copies are expanded in binary-decomposed batches 1, 2, 4, ...
  // plus the residual: any count in [ 0 , avail ] is reachable, with O(log)
  // merges instead of one per copy. As soon as one batch yields no new state
  // none of the following can either (each translate of the frontier is then
  // a translate of an already-merged one), so the item is closed early
  for( long mdone = 0 , mult = 1 ; mdone < avail ; mdone += mult , mult *= 2 ) {
   if( mult > avail - mdone )
    mult = avail - mdone;
   const bool lastB = ( mdone + mult >= avail );

   double eff_nr , eff_nl;
   bool hasL;
   if( isRight ) {
    eff_nr = lastB ? eff_t : eff_it;
    hasL = ( sn >= 0 );
    eff_nl = eff_s;
    }
   else {
    eff_nr = eff_t;
    hasL = ( ! lastB ) || ( sn >= 0 );
    eff_nl = lastB ? eff_s : eff_it;
    }

   const WPvec & prev = front;
   const int np = int( prev.size() );
   const long   dmul = isRight ? mult : - mult;
   const long   dw = dmul * w[ it ];
   const PT dp = PT( dmul ) * pt[ it ];

   // what is left after this batch: the profit the right items can still
   // add, the weight the left ones can still shed
   const long   rest = avail - mdone - mult;
   const double pleft = prem[ t ] - pfixU +
                        ( isRight ? p[ it ] * double( rest ) : 0 );
   const long   wleft = wrem[ s + 1 ] - wfixU +
                        ( isRight ? 0 : w[ it ] * rest );

   // the per-state survival test: not bound-pruned (fractional fill for
   // feasible states, fractional removal for over-capacity ones, each capped
   // by what the unprocessed items can still give: the completion bounds)
   auto survives = [ & ]( long nw , double npp ) -> bool {
    double B;
    if( nw <= Ceff )
     B = std::min( npp + double( Ceff - nw ) * eff_nr , npp + pleft );
    else
     B = ( hasL && ( nw - Ceff <= wleft ) )
         ? npp - double( nw - Ceff ) * eff_nl : - Inf< double >();
    return( B >= z + zmargin );
    };

   // B2: after `guard` consecutive no-new-state expansions on this side, run
   // a read-only check (a scan of the would-be merge that writes nothing): if
   // no toggled state would survive, skip the whole item. A wasteful check (a
   // survivor exists) doubles the threshold; 40 consecutive skips force an
   // expansion so that the state pruning still operates periodically
   if( f_dp_ext && ( ! mdone ) && ( nonew > guard ) && ( roSkip < 40 ) ) {
    bool wouldAdd = false;
    double maxp = - Inf< double >();      // max untoggled profit at <= weight
    for( int ia = 0 , ib = 0 ; ib < np ; ) {
     if( ( ia < np ) && ( prev[ ia ].wsum <= prev[ ib ].wsum + dw ) )
      { maxp = prev[ ia ].psum; ++ia; continue; }
     const long   nw  = prev[ ib ].wsum + dw;
     const double npp = prev[ ib ].psum + dp;
     ++ib;
     if( ( npp > maxp ) && survives( nw , npp ) )
      { wouldAdd = true; break; }
     }
    if( ! wouldAdd ) {
     ++roSkip;
     CSTAT( ++cs_b2 );
     break;
     }
    guard *= 2;
    }
   roSkip = 0;

   // expansion: linear merge of the frontier with its ( +dw , +dp ) translate
   // (both sorted by ascending weight), keeping the merged Pareto frontier;
   // candidates with equal weight are merged higher-profit-first
   // the output is sized for the worst case and filled by index, then cut
   // to the `nc` states kept; the decision of this step goes in the bit
   // K mod 64 of the masks (K the step of prev), which restart at 0 after a
   // checkpoint
   // (cur holds nothing worth keeping: emptied before growing, so that a
   // reallocation copies nothing, and grown geometrically, so that it is
   // rare and the freed blocks are reused rather than handed back)
   const long wdead = Ceff + ( hasL ? wleft : 0 );
   cur.clear();
   if( cur.capacity() < 2 * std::size_t( np ) )
    cur.reserve( 4 * std::size_t( np ) );
   cur.resize( 2 * std::size_t( np ) );
   const std::size_t K = step_item.size() - 1;
   const unsigned bit = unsigned( K % 64 );
   const std::uint64_t keep_mask = ( bit == 0 ) ? 0 : ~std::uint64_t( 0 );
   std::size_t nc = 0;
   // with integer profits the completion bound is folded in the dominance
   // test: a candidate is dropped unless it beats the profit that, plus all
   // the right ones still to add, reaches the improving threshold (taken at
   // the start of the merge, i.e., a weaker but valid one; such a candidate
   // cannot improve the incumbent either), and is left to survives_m() the
   // continuous bound only
   PT frontier = std::numeric_limits< PT >::lowest();
   if constexpr( std::is_integral_v< PT > )
    frontier = PT( ceil_safe( z + zmargin - pleft ) ) - 1;
   // the threshold of the bounds and, with integer profits, the incumbent
   // as a state profit, refreshed when the incumbent improves
   double zthr = z + zmargin;
   PT zst = std::numeric_limits< PT >::lowest();
   if constexpr( std::is_integral_v< PT > )
    if( z > double( std::numeric_limits< PT >::lowest() ) )
     zst = PT( std::floor( z ) );
   const auto survives_m = [ & ]( long nw , double npp ) -> bool {
    if constexpr( std::is_integral_v< PT > ) {
     const double B = ( nw <= Ceff )
      ? npp + double( Ceff - nw ) * eff_nr
      : ( ( hasL && ( nw - Ceff <= wleft ) )
          ? npp - double( nw - Ceff ) * eff_nl : - Inf< double >() );
     return( B >= zthr );
     }
    else
     return( survives( nw , npp ) );
    };
   int newStates = 0;
   bool improved = false;

   int ia = 0 , ib = 0;
   CSTAT( ( ++cs_merge , cs_cand += 2 * std::size_t( np ) ) );
   // one candidate, by ascending weight: false once one is too heavy to be
   // brought back under the capacity by the left items still to drop (then
   // so are all the others)
   const auto take = [ & ]( long nw , PT npp , int par , unsigned took ,
                            std::uint64_t pmask ) -> void {
    if( npp <= frontier )       // Pareto-dominated (weight already >=)
     return;
    frontier = npp;

    // update the incumbent from feasible states, recording the achieving
    // solution now (the state may be bound-pruned right after); with
    // integer profits a feasible improving state beats z by a whole g_p,
    // so its bound passes the threshold and it is enough to look among the
    // surviving states, in integers
    const auto improve = [ & ]( void ) {
     record_incumbent( double( npp ) , par , it , took ? dmul : 0 , nullptr ,
                       0 );
     improved = true;
     zthr = z + zmargin;
     if constexpr( std::is_integral_v< PT > )
      zst = PT( std::floor( z ) );
     };
    if constexpr( ! std::is_integral_v< PT > )
     if( ( nw <= Ceff ) && ( double( npp ) > z ) )
      improve();

    if( survives_m( nw , double( npp ) ) ) {
     if constexpr( std::is_integral_v< PT > )
      if( ( nw <= Ceff ) && ( npp > zst ) )
       improve();
     cur[ nc++ ] = WP{ nw , npp , ( pmask & keep_mask ) |
                                  ( std::uint64_t( took ) << bit ) };
     newStates += int( took );
     }
    };

   // both streams, then the rest of either, with no index check in the loop;
   // the candidates come by ascending weight, so once one is too heavy to be
   // brought back under the capacity by the left items still to drop, so
   // are all the others: the merge stops there
   [ & ]( void ) {
    while( ( ia < np ) && ( ib < np ) ) {
     const WP & a = prev[ ia ];
     const WP & b = prev[ ib ];
     const long wb = b.wsum + dw;
     const PT   pb = b.psum + dp;
     if( ( wb < a.wsum ) || ( ( wb == a.wsum ) && ( pb > a.psum ) ) ) {
      if( wb > wdead )
       return;
      take( wb , pb , ib , 1 , b.mask );
      ++ib;
      }
     else {
      if( a.wsum > wdead )
       return;
      take( a.wsum , a.psum , ia , 0 , a.mask );
      ++ia;
      }
     }
    for( ; ( ia < np ) && ( prev[ ia ].wsum <= wdead ) ; ++ia )
     take( prev[ ia ].wsum , prev[ ia ].psum , ia , 0 , prev[ ia ].mask );
    for( ; ( ib < np ) && ( prev[ ib ].wsum + dw <= wdead ) ; ++ib )
     take( prev[ ib ].wsum + dw , prev[ ib ].psum + dp , ib , 1 ,
           prev[ ib ].mask );
    }();
   cur.resize( nc );

   bool grown;
   if( newStates ) {
    if( ( bit == 0 ) && ( K > 0 ) )      // prev is the checkpoint K
     chk.push_back( front );
    std::swap( front , cur );
    nstates += nc;
    step_item.push_back( it );
    step_mult.push_back( dmul );
    nonew = 0;
    grown = true;

    // dominance after an expansion that grew the frontier: the items still
    // to be reached that `it` dominates are fixed if the frontier is closed
    // under one more copy of `it`, and otherwise (right side) those that,
    // added to the lightest state from which `it` still yields a surviving
    // state, would be too heavy to ever come back under the capacity
    if( ( f_dom_fix > 1 ) && ( ! mdone ) && ( p[ it ] >= 0 ) &&
        ( ( dom_work += front.size() ) >= dom_thresh ) ) {
     dom_work = 0;
     std::vector< int > dom;
     if( isRight ) {
      for( int j = t ; j < mA ; ++j )
       if( ( ! fixd[ j ] ) && ( w[ j ] >= w[ it ] ) &&
           ( p[ j ] <= p[ it ] * double( w[ j ] / w[ it ] ) ) )
        dom.push_back( j );
      }
     else
      for( int j = s ; j >= 0 ; --j )
       if( ( ! fixd[ j ] ) && ( w[ j ] <= w[ it ] ) && ( p[ j ] >= p[ it ] ) )
        dom.push_back( j );
     if( dom.size() > 2 ) {
      // read-only scan of the merge with one more copy of `it`: the weight
      // of the first state that still yields something (-1: none)
      long min_req = -1;
      const int nf = int( front.size() );
      PT maxp = std::numeric_limits< PT >::lowest();
      for( int ia = 0 , ib = 0 ; ib < nf ; ) {
       if( ( ia < nf ) && ( front[ ia ].wsum <= front[ ib ].wsum + dw ) )
        { maxp = front[ ia ].psum; ++ia; continue; }
       const long nw = front[ ib ].wsum + dw;
       const PT npp = front[ ib ].psum + dp;
       if( ( npp > maxp ) && survives( nw , double( npp ) ) )
        { min_req = front[ ib ].wsum; break; }
       ++ib;
       }
      if( min_req < 0 )
       for( int j : dom )
        fix_item( j , isRight );
      else
       if( isRight )
        for( int j : dom )
         if( min_req + w[ j ] - ( w[ j ] / w[ it ] - 1 ) * w[ it ] > wdead )
          fix_item( j , true );
      }
     }
    if( relx ) {               // the surrogate recursion runs on a budget
     f_rec_left -= long( front.size() );
     if( f_rec_left < 0 )      // exhausted: give up, the caller knows
      return( Inf< double >() );
     }
    if( ( f_state_cap >= 0 ) && ( long( nstates ) > f_state_cap ) ) {
     dead = true;              // capped: stop with the best solution so far
     break;
     }
    }
   else {
    // the frontier did not grow: do not store a step (nothing to backtrack
    // through), the merged frontier is the last one pruned under the current
    // bound (a state dropped as dominated by a pruned translate would have
    // been pruned too, being heavier and worth less) and takes its place,
    // the masks of its states untouched - but where they restart (a
    // checkpoint) the last one is compacted instead; on a first-batch
    // failure also remember the item as the dominance reference (B3)
    ++nonew;
    CSTAT( ++cs_nogrow );
    if( bit != 0 ) {
     nstates -= front.size() - nc;
     std::swap( front , cur );
     }
    else {
     std::size_t kept = 0;
     for( std::size_t i = 0 ; i < front.size() ; ++i )
      if( survives( front[ i ].wsum , front[ i ].psum ) )
       front[ kept++ ] = front[ i ];
     nstates -= front.size() - kept;
     front.resize( kept );
     }
    if( ! mdone ) {
     if( isRight ) { domRw = w[ it ]; domRp = p[ it ]; }
     else          { domLw = w[ it ]; domLp = p[ it ]; }
     // the frontier is closed under one more copy of `it` (and stays so
     // under any later expansion), so is it under any item that `it`
     // dominates: mark all those still to be processed on this side
     if( ( f_dom_fix > 1 ) && ( p[ it ] >= 0 ) )
      mark_dominated( it , isRight );
     }
    grown = false;
    }

   if( front.empty() ) {
    dead = true;
    break;
    }

   if( f_heur ) {
    // B5 (GCH): a new incumbent may leave slack capacity; greedily complete it
    // with any still-unselected copies that fit (always feasible)
    if( improved ) {
     long wcur = 0;
     for( int i = 0 ; i < mA ; ++i )
      wcur += best_in[ i ] * w[ i ];
     for( int i = 0 ; ( i < mA ) && ( wcur < Ceff ) ; ++i )
      while( ( best_in[ i ] < d[ i ] ) && ( wcur + w[ i ] <= Ceff ) )
       { ++best_in[ i ]; wcur += w[ i ]; z += p[ i ]; }
     }

    // B5 (SPH): after every expansion, pair a sparse sample of the unprocessed
    // items with the frontier: one random item per block of size gamma on each
    // side, for ~ beta = ceil( |S| / 50 ) pairings in all - much cheaper than
    // the O( |S| ) expansion itself, so it runs unconditionally
    {
     const WPvec & f = front;
     const long beta = long( f.size() / 50 ) + 1;
     const long gamma = std::max( 1L , long( mA ) / beta );
     const long gl = std::min( long( s + 1 ) , gamma );
     const long gr = std::min( long( mA - t ) , gamma );
     if( gl >= 3 )
      for( long st0 = 0 ; st0 <= s ; st0 += gl ) {
       const int i = int( st0 + rnd( unsigned(
                            std::min( gl , long( s + 1 ) - st0 ) ) ) );
       const int j = best_under( f , Ceff + w[ i ] );
       if( ( j >= 0 ) && ( f[ j ].psum - p[ i ] > z ) ) {
        const int xs[ 1 ] = { - ( i + 1 ) };
        record_incumbent( f[ j ].psum - p[ i ] , j , -1 , 0 , xs , 1 );
        }
       }
     if( gr >= 3 )
      for( long st0 = t ; st0 < mA ; st0 += gr ) {
       const int i = int( st0 + rnd( unsigned(
                            std::min( gr , long( mA ) - st0 ) ) ) );
       const int j = best_under( f , Ceff - w[ i ] );
       if( ( j >= 0 ) && ( f[ j ].psum + p[ i ] > z ) ) {
        const int xs[ 1 ] = { i + 1 };
        record_incumbent( f[ j ].psum + p[ i ] , j , -1 , 0 , xs , 1 );
        }
       }
     }

    // B5 (PH / TPH / SSPH): when the states generated since the last call
    // have passed the (doubling) threshold, pair every still-unprocessed item
    // with its best frontier state via binary search; if the number of
    // ( left , right ) pairs is small wrt |S| / log2 |S| also try every pair
    // (TPH); on huge frontiers additionally enumerate all subsets of k
    // random unprocessed items with alpha k 2^k <= |S| , alpha = 20 (SSPH)
    if( ( ph_work += front.size() ) >= ph_thresh ) {
     ph_thresh *= 2;
     ph_work = 0;
     const WPvec & f = front;
     for( int i = t ; i < mA ; ++i ) {      // pair with an unprocessed right
      const int j = best_under( f , Ceff - w[ i ] );
      if( ( j >= 0 ) && ( f[ j ].psum + p[ i ] > z ) ) {
       const int xs[ 1 ] = { i + 1 };
       record_incumbent( f[ j ].psum + p[ i ] , j , -1 , 0 , xs , 1 );
       }
      }
     for( int i = s ; i >= 0 ; --i ) {      // pair with an unprocessed left
      const int j = best_under( f , Ceff + w[ i ] );
      if( ( j >= 0 ) && ( f[ j ].psum - p[ i ] > z ) ) {
       const int xs[ 1 ] = { - ( i + 1 ) };
       record_incumbent( f[ j ].psum - p[ i ] , j , -1 , 0 , xs , 1 );
       }
      }

     // TPH: drop an unprocessed left item AND add an unprocessed right one
     const double l2 = std::log2( double( f.size() ) + 2 );
     if( double( s + 1 ) * double( mA - t ) < double( f.size() ) / l2 )
      for( int il = 0 ; il <= s ; ++il )
       for( int ir = t ; ir < mA ; ++ir ) {
        const int j = best_under( f , Ceff + w[ il ] - w[ ir ] );
        const double v = ( j >= 0 ) ? f[ j ].psum - p[ il ] + p[ ir ]
                                    : - Inf< double >();
        if( v > z ) {
         const int xs[ 2 ] = { - ( il + 1 ) , ir + 1 };
         record_incumbent( v , j , -1 , 0 , xs , 2 );
         }
        }

     // SSPH: largest k with alpha k 2^k <= |S|; k random unprocessed items,
     // all their non-singleton subsets paired with the frontier
     constexpr long alpha = 20;
     int k = 0;
     while( ( k + 1 <= 24 ) &&
            ( alpha * long( k + 1 ) * ( 1L << ( k + 1 ) ) <=
              long( f.size() ) ) )
      ++k;
     const int navail = ( s + 1 ) + ( mA - t );
     if( ( k >= 2 ) && ( navail >= k ) ) {
      std::vector< int > pick( k );
      for( int e = 0 ; e < k ; ++e ) {
       const int r = int( rnd( unsigned( navail ) ) );
       pick[ e ] = ( r <= s ) ? r : ( t + ( r - s - 1 ) );
       }
      // duplicates would unbalance the ( dw , dp ) sums vs the count toggles
      std::sort( pick.begin() , pick.end() );
      pick.erase( std::unique( pick.begin() , pick.end() ) , pick.end() );
      k = int( pick.size() );
      std::vector< int > xs;
      for( unsigned msk = 1 ; ( k >= 2 ) && ( msk < ( 1u << k ) ) ; ++msk ) {
       if( ! ( msk & ( msk - 1 ) ) )        // singletons: PH already did them
        continue;
       long dw2 = 0; double dp2 = 0;
       xs.clear();
       for( int e = 0 ; e < k ; ++e )
        if( msk & ( 1u << e ) ) {
         const int i = pick[ e ];
         if( i <= s )
          { dw2 -= w[ i ]; dp2 -= p[ i ]; xs.push_back( -( i + 1 ) ); }
         else         { dw2 += w[ i ]; dp2 += p[ i ]; xs.push_back( i + 1 ); }
         }
       const int j = best_under( f , Ceff - dw2 );
       if( ( j >= 0 ) && ( f[ j ].psum + dp2 > z ) )
        record_incumbent( f[ j ].psum + dp2 , j , -1 , 0 , xs.data() ,
                          int( xs.size() ) );
       }
      }
     }
    }

   if( ! grown )               // L10 / T11: close the item early
    break;
   }                           // end of the binary-decomposed batch loop

  if( dead )
   break;

  // hard instance detected: tighten the ceiling with the surrogate bound,
  // computed on the per-unit expanded view; with intSurrogate = 2 the exact
  // optimum of the surrogate subproblem both certifies and (when it attains
  // the forced cardinality) raises the incumbent
  if( ( f_surrogate > 0 ) && ( ! relx ) && surr_tries &&
      ( int( front.size() ) > f_surr_trigger ) &&
      ( ( z > surr_z ) || ( nstates >= surr_redo ) ) ) {
   // with bit 1 of intSurrAdapt the recursion of the exact surrogate solve
   // may generate at most twice the states of the enumeration so far: a
   // solve that certifies does it well within that, one that does not would
   // otherwise cost up to orders of magnitude more than the whole solve
   const std::size_t mainst = nstates;
   f_rec_left = ( f_surr_adapt & 2 ) ? 2 * long( mainst )
                                     : std::numeric_limits< long >::max();
#if CORE_STATS
   const double dbef = dantzig;
   core_rec_states = 0;
#endif
   if( f_surrogate > 1 ) {
    const double zb4 = z;
    std::fill( selE.begin() , selE.end() , 0 );
    const double su = surrogate_solve( wE , pE , mE , bE , Ceff , z , selE );
    if( z > zb4 ) {            // improved: fold the unit selection to counts
     beat = true;
     std::fill( best_in.begin() , best_in.end() , 0 );
     for( int u = 0 ; u < mE ; ++u )
      if( selE[ u ] )
       ++best_in[ uitem( u ) ];
     }
    dantzig = std::min( dantzig , snap( su ) );
    }
   else
    dantzig = std::min( dantzig ,
                        snap( surrogate_bound( wE , pE , mE , bE , Ceff ,
                                               z ) ) );
#if CORE_STATS
   fprintf( stderr , "SURRTRY front=%zu step=%zu z=%.0f before=%.0f "
            "bound=%.0f mainst=%zu recst=%zu\n" , front.size() ,
            step_item.size() , z , dbef , dantzig , mainst ,
            core_rec_states );
#endif
   if( ( f_surr_adapt & 2 ) && ( f_rec_left < 0 ) )
    surr_redo = 2 * mainst;    // out of budget: retry at twice the work
   else {
    --surr_tries;
    surr_redo = std::numeric_limits< std::size_t >::max();
    }
   surr_z = z;
   }
  }

#if CORE_STATS
 {
  std::size_t mxf = front.size() , totf = nstates;
  for( const auto & f : chk )
   mxf = std::max( mxf , f.size() );
  if( relx )
   core_rec_states += totf;
  fprintf( stderr ,
           "CORE_STATS m=%d mA=%d mE=%d bE=%d steps=%zu maxfront=%zu "
           "totstates=%zu z=%.0f dantzig=%.0f closed=%s left=%d right=%d\n" ,
           m , mA , mE , bE , step_item.size() , mxf , totf , z , dantzig ,
           ( z >= dantzig - 1e-9 ) ? "bound" : "exhaust" , s + 1 , mA - t );
  if( ! relx )
   fprintf( stderr , "CORE_INC z0=%.0f z=%.0f improvements=%zu "
            "at_states=%zu of %zu\n" , cs_z0 , z , cs_zimp , cs_zstep ,
            nstates );
  if( ! relx )
   fprintf( stderr , "CORE_EXP merges=%zu nogrow=%zu cand=%zu skipB2=%zu "
            "skipB3=%zu skipB4=%zu\n" , cs_merge , cs_nogrow , cs_cand ,
            cs_b2 , cs_b3 , cs_b4 );
  }
#endif

 if( ! beat )                  // nothing strictly above the caller's lb
  return( - Inf< double >() );

 // unfold the per-item copy counts onto the original input units, walking the
 // multiplicity-reduction chain: each level keeps own units first, with the
 // parity adjusted so that the ( 2w , 2p ) child receives an integral count
 for( int j = 0 ; j < mA ; ++j ) {
  int r = regIdx[ j ];
  long c = best_in[ j ];
  while( ( r >= 0 ) && ( c > 0 ) ) {
   const Agg & a = reg[ r ];
   long take = std::min( c , long( a.mcount ) );
   if( ( c - take ) % 2 )
    --take;
   for( long k = 0 ; k < take ; ++k )
    in[ byk[ a.mstart + k ] ] = 1;
   c = ( c - take ) / 2;
   r = a.child;
   }
  }
 return( z );
 }

/*--------------------------------------------------------------------------*/
/*------------------- A3: TRIVIAL DIVISIBILITY BOUND -----------------------*/
/*--------------------------------------------------------------------------*/

long CoreDPBinaryKnapsackSolver::divisibility_capacity(
                       const std::vector< long > & w , int m , long C ) const
{
 // g = gcd of the item weights; every reachable weight is a multiple of g, so
 // the usable capacity is floor( C / g ) * g. g == 1 (the common case) leaves
 // the capacity unchanged. C >= 0 here (compute() returns kInfeasible first).
 if( m == 0 )
  return( C );
 long g = w[ 0 ];
 for( int i = 1 ; ( i < m ) && ( g > 1 ) ; ++i )
  g = std::gcd( g , w[ i ] );
 return( g > 1 ? ( C / g ) * g : C );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ A1: INITIAL HEURISTIC ---------------------------*/
/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::initial_heuristic(
                       const std::vector< long > & w ,
                       const std::vector< double > & p , int m , int b ,
                       long wsumb , double psumb , long C ,
                       std::vector< char > & in )
{
 // all arguments are in efficiency-sorted core order; @p in is the achieving
 // selection. We return the best of a few O(n) feasible extensions of the
 // break solution; the break solution itself is feasible, so the result is
 // always a valid lower bound on the core optimum.

 // H1 - forward greedy fill: keep the break items in, then add the most
 // efficient right items that still fit
 std::vector< char > h( m , 0 );
 for( int i = 0 ; i < b ; ++i )
  h[ i ] = 1;
 long wcur = wsumb;
 double best = psumb;
 for( int t = b ; t < m ; ++t )
  if( wcur + w[ t ] <= C ) { h[ t ] = 1; wcur += w[ t ]; best += p[ t ]; }
 in = h;

 // H2 - swap: drop the least efficient break item ( b - 1 ), then greedily
 // refill from the right items; the freed weight can let a better-profit
 // right item in where H1's plain fill could not fit it
 if( ( b > 0 ) && ( b < m ) ) {
  std::vector< char > g( m , 0 );
  for( int i = 0 ; i < b - 1 ; ++i )
   g[ i ] = 1;
  long wg = wsumb - w[ b - 1 ];
  double pg = psumb - p[ b - 1 ];
  for( int t = b ; t < m ; ++t )
   if( wg + w[ t ] <= C ) { g[ t ] = 1; wg += w[ t ]; pg += p[ t ]; }
  if( pg > best ) { best = pg; in = g; }
  }

 // H3 - drop the least efficient left items until the break item b fits, take
 // it, then greedily fill the rest; keep it if it beats H1 / H2
 if( b < m ) {
  std::vector< char > g( m , 0 );
  for( int i = 0 ; i < b ; ++i )
   g[ i ] = 1;
  long wg = wsumb;
  double pg = psumb;
  int li = b - 1;
  while( ( li >= 0 ) && ( wg + w[ b ] > C ) ) {
   g[ li ] = 0; wg -= w[ li ]; pg -= p[ li ]; --li;
   }
  if( wg + w[ b ] <= C ) {
   g[ b ] = 1; wg += w[ b ]; pg += p[ b ];
   for( int t = b + 1 ; t < m ; ++t )
    if( wg + w[ t ] <= C ) { g[ t ] = 1; wg += w[ t ]; pg += p[ t ]; }
   if( pg > best ) { best = pg; in = g; }
   }
  }

 return( best );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- A2: MARTELLO-TOTH U2 BOUND -------------------------*/
/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::upper_bound_u2(
                       const std::vector< long > & w ,
                       const std::vector< double > & p , int m , int b ,
                       long wsumb , double psumb , long C , bool intp ) const
{
 if( b >= m )                  // the whole core fits: break solution is exact
  return( psumb );

 const long r = C - wsumb;     // residual capacity, 0 <= r < w[ b ]

 // U' : the break item b is left OUT; fill r with the next item b+1 (the most
 // efficient remaining), rounded down for integer profits
 double Uex = psumb;
 if( b + 1 < m ) {
  const double add = double( r ) * ( p[ b + 1 ] / double( w[ b + 1 ] ) );
  Uex += intp ? floor_safe( add ) : add;
  }

 // U'' : the break item b is taken; ( w[ b ] - r ) of weight must be removed
 // from the previous (least efficient included) item b-1, rounded up
 double Uin = psumb + p[ b ];
 if( b >= 1 ) {
  const double rem = double( w[ b ] - r ) *
                     ( p[ b - 1 ] / double( w[ b - 1 ] ) );
  Uin -= intp ? ceil_safe( rem ) : rem;
  }

 // both U' and U'' are <= the plain Dantzig bound, so their max is a valid,
 // tighter ceiling
 return( std::max( Uex , Uin ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------- A2: SURROGATE / CARDINALITY BOUND -----------------------*/
/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::surrogate_card_bound(
                       const std::vector< long > & w ,
                       const std::vector< double > & p , int m ,
                       int card , long C , long slo , long shi ,
                       long * out_sur ) const
{
 // tightest surrogate fractional bound for the cardinality value `card`,
 // minimised over the integer multiplier s in [ slo , shi ]: the surrogate
 // problem is max sum p_i x_i s.t. sum (w_i + s) x_i <= C + s*card, x in [0,1];
 // ANY s in the valid sign range gives a valid upper bound, so the min is too.
 // combo.c surbin's gradient binary search is used to locate the minimiser;
 // with bit 0 of intSurrAdapt it first probes the end of the range next to
 // s = 0, where the surrogate is the plain continuous relaxation: when the
 // gradient points to that end the search is over in one or two steps
 bool probe = ( f_surr_adapt & 1 ) && ( ( slo == 0 ) || ( shi == 0 ) );
 double best = Inf< double >();
 long bestsur = 0;
 std::vector< int > ord( m );
#if CORE_STATS
 const long slo0 = slo;
 int evals = 0;
#endif

 while( slo <= shi ) {
  const long s = probe ? ( ( slo == 0 ) ? 0 : std::max( slo , -1L ) )
                       : slo + ( shi - slo ) / 2;  // floor towards slo
  probe = false;
#if CORE_STATS
  ++evals;
#endif

  long csur = C + s * ( long ) card;
  if( csur < 0 ) csur = 0;

  // items with non-positive surrogate weight ( w_i + s <= 0 ) are always in:
  // they add profit and free capacity; the rest are filled by efficiency
  double psum = 0;
  double cap = double( csur );
  ord.clear();
  for( int i = 0 ; i < m ; ++i ) {
   const long ww = w[ i ] + s;
   if( ww <= 0 ) { psum += p[ i ]; cap -= double( ww ); }
   else ord.push_back( i );
   }
  const auto before = [ & ]( int a , int c ) {
   return( p[ a ] * double( w[ c ] + s ) > p[ c ] * double( w[ a ] + s ) );
   };

  // greedy fractional fill of the surrogate capacity; the critical item is
  // located by quickselect on the efficiencies (Balas-Zemel), each halving
  // taking whole the upper half when it fits, so only a short tail is sorted
  int d = int( m - ord.size() );               // items already taken (ww <= 0)
  double r = cap;
  int bi = -1;
  int lo = 0 , hi = int( ord.size() );
  while( hi - lo > 16 ) {
   const int mid = lo + ( hi - lo ) / 2;
   std::nth_element( ord.begin() + lo , ord.begin() + mid ,
                     ord.begin() + hi , before );
   double wh = 0 , ph = 0;
   for( int t = lo ; t < mid ; ++t )
    { wh += double( w[ ord[ t ] ] + s ); ph += p[ ord[ t ] ]; }
   if( wh <= r ) { psum += ph; r -= wh; d += mid - lo; lo = mid; }
   else hi = mid;
   }
  std::sort( ord.begin() + lo , ord.begin() + hi , before );
  for( int t = lo ; t < hi ; ++t ) {
   const double ww = double( w[ ord[ t ] ] + s );
   if( ww <= r ) { psum += p[ ord[ t ] ]; r -= ww; ++d; }
   else { bi = ord[ t ]; break; }
   }

  double ua, gr;
  if( bi < 0 ) {                                // everything fits: exact
   ua = psum; gr = 0;
   }
  else {
   const double wwb = double( w[ bi ] + s );
   const double pb = p[ bi ];
   ua = psum + r * ( pb / wwb );
   // bound at multiplier s+1 (same break), to get the gradient sign
   const double ub = psum + ( r + double( card - d ) ) * ( pb / ( wwb + 1 ) );
   gr = ub - ua;
   }

  if( ua < best ) { best = ua; bestsur = s; }
  if( gr > 0 ) shi = s - 1; else slo = s + 1;
  }

#if CORE_STATS
 fprintf( stderr , "SURCARD card=%d range=%s evals=%d sur=%ld best=%.0f\n" ,
          card , ( slo0 < 0 ) ? "neg" : "pos" , evals , bestsur , best );
#endif
 if( out_sur )
  *out_sur = bestsur;
 return( best );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::surrogate_bound(
                       const std::vector< long > & w ,
                       const std::vector< double > & p , int m , int b ,
                       long C , double z ) const
{
 // COMBO surrogate (combo.c surrelax): only worth it when the cardinality is
 // essentially forced, i.e. card1 (max items that fit C = N_max) or card2 (min
 // items whose profit beats z = N_min) is within 1 of the break cardinality b.
 if( ( b <= 0 ) || ( b >= m ) )
  return( Inf< double >() );

 long maxw = 0;
 for( int i = 0 ; i < m ; ++i )
  if( w[ i ] > maxw ) maxw = w[ i ];
 const long minsur = - maxw, maxsur = maxw;

 // card1 = N_max: largest k with the k smallest weights fitting C
 std::vector< long > ws( w );
 const int card1 = prefix_fit( ws , std::less< long >() , C );

 // card2 = N_min: smallest k with the k largest profits exceeding z
 std::vector< double > ps( p );
 const int card2 = ( z < 0 ) ? 0 :
  std::min( prefix_fit( ps , std::greater< double >() , z ) + 1 , m );

 // apply the surrogate only for a forced cardinality (combo.c surrelax order);
 // a cardinality dichotomy at b is a complete, always-valid case split
 double u;
 if( card2 == b + 1 ) {
  u = surrogate_card_bound( w , p , m , b + 1 , C , minsur , 0 );
  if( u < z ) u = z;                       // bound for an IMPROVED solution
  }
 else if( card1 == b ) {
  u = surrogate_card_bound( w , p , m , b , C , 0 , maxsur );
  }
 else if( card1 == b + 1 ) {
  u = std::max( surrogate_card_bound( w , p , m , b + 1 , C , minsur , 0 ) ,
                surrogate_card_bound( w , p , m , b , C , 0 , maxsur ) );
  }
 else if( card2 == b ) {
  u = std::max( surrogate_card_bound( w , p , m , b , C , 0 , maxsur ) ,
                surrogate_card_bound( w , p , m , b + 1 , C , minsur , 0 ) );
  }
 else
  return( Inf< double >() );                    // wide window: no tightening

 return( u );
 }

/*--------------------------------------------------------------------------*/
/*------------------ A2: SURROGATE-SOLVE (COMBO solvesur) ------------------*/
/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::surrogate_solve(
                       const std::vector< long > & w ,
                       const std::vector< double > & p , int m , int b ,
                       long C , double & z , std::vector< char > & best_in )
{
 if( ( b <= 0 ) || ( b >= m ) )
  return( Inf< double >() );

 long maxw = 0;
 for( int i = 0 ; i < m ; ++i )
  if( w[ i ] > maxw ) maxw = w[ i ];
 const long minsur = - maxw, maxsur = maxw;

 // card1 = N_max: largest k with the k smallest weights fitting C, so every
 // feasible solution has <= card1 items; card2 = N_min: smallest k with the k
 // largest profits exceeding z, so every IMPROVING solution has >= card2 items
 std::vector< long > ws( w );
 const int card1 = prefix_fit( ws , std::less< long >() , C );
 std::vector< double > ps( p );
 const int card2 = ( z < 0 ) ? 0 :
  std::min( prefix_fit( ps , std::greater< double >() , z ) + 1 , m );

 // upper bound on every solution of one forced-cardinality case: the tightest
 // surrogate fractional bound, improved to the EXACT optimum of the surrogate
 // subproblem when the latter can be computed (the surrogate is a relaxation,
 // so its optimum is itself a valid bound - usually the one that certifies z);
 // an exact optimum attaining the cardinality is also feasible for the
 // original problem, and raises the incumbent when improving
 auto case_bound = [ & ]( int card , long slo , long shi ) -> double {
  long sur = 0;
  const double uf = surrogate_card_bound( w , p , m , card , C , slo , shi ,
                                          & sur );
  if( ( sur == 0 ) || ( uf <= z ) )    // nothing improving in the case anyway
   return( uf );

  // shift the weights by the multiplier; items with non-positive surrogate
  // weight are always taken (profit + freed capacity), the rest go to the
  // RECURSIVE core enumeration (with relx set: no surrogate in the surrogate),
  // floored at the incumbent so it only ever works towards an improvement
  std::vector< long > mw;  std::vector< double > mp;
  std::vector< int > mi;
  std::vector< char > chosen( m , 0 );
  double ps = 0;
  long csur = C + sur * ( long ) card;
  for( int i = 0 ; i < m ; ++i ) {
   const long ww = w[ i ] + sur;
   if( ww <= 0 ) { chosen[ i ] = 1; ps += p[ i ]; csur -= ww; }
   else { mw.push_back( ww ); mp.push_back( p[ i ] ); mi.push_back( i ); }
   }
  if( csur < 0 )                       // surrogate infeasible: nothing beats z
   return( z );

  std::vector< char > sub;
  const double val = core_enumerate( mw , mp , csur , sub , z - ps , true );
  if( val == Inf< double >() )         // out of budget: the fractional bound
   return( uf );
  if( val == - Inf< double >() )       // proved: nothing in the case beats z
   return( z );

  // exact surrogate optimum ( > z ): a valid bound for the case; and when it
  // attains the cardinality and the true capacity it is feasible for the
  // original problem, raising the incumbent
  for( std::size_t e = 0 ; e < mi.size() ; ++e )
   if( sub[ e ] )
    chosen[ mi[ e ] ] = 1;
  long got = 0, wsum = 0;
  double psum = 0;
  for( int i = 0 ; i < m ; ++i )
   if( chosen[ i ] ) { ++got; wsum += w[ i ]; psum += p[ i ]; }
  if( ( got == card ) && ( wsum <= C ) && ( psum > z ) ) {
   z = psum;
   best_in = chosen;
   }
  return( ps + val );
  };

#if CORE_STATS
 fprintf( stderr , "SURR card1=%d card2=%d b=%d m=%d z=%.0f\n" ,
          card1 , card2 , b , m , z );
#endif

 // forced-cardinality dichotomy at b (combo.c surrelax case order): a
 // multiplier >= 0 relaxes all solutions with <= card items, <= 0 those with
 // >= card items; the two-case splits cover every (improving) solution
 if( card2 == b + 1 )                  // improving ==> >= b + 1 items
  return( std::max( case_bound( b + 1 , minsur , 0 ) , z ) );
 if( card1 == b )                      // feasible ==> <= b items
  return( case_bound( b , 0 , maxsur ) );
 if( card1 == b + 1 )                  // <= b or exactly b + 1
  return( std::max( case_bound( b + 1 , minsur , 0 ) ,
                    case_bound( b , 0 , maxsur ) ) );
 if( card2 == b )                      // improving ==> >= b: split at b
  return( std::max( case_bound( b , 0 , maxsur ) ,
                    case_bound( b + 1 , minsur , 0 ) ) );

 return( Inf< double >() );            // wide window: no tightening
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::solve_with_continuous(
                       std::vector< char > & in , std::vector< double > & cx )
{
 // With continuous variables the integer part can no longer be optimised in
 // isolation: a lighter integer solution leaves more room for the (fractional)
 // continuous fill. Since the dominance relation is unchanged, we enumerate the
 // full feasible integer Pareto frontier (dominance only, no Dantzig pruning)
 // and pick the (weight, profit) state maximising profit + g( C - weight ),
 // where g is the fractional knapsack value of the continuous items.

 const int m = int( v_w.size() );
 in.assign( m , 0 );

 // integer core sorted by non-increasing efficiency (keeps the frontier small)
 std::vector< int > ord( m );
 for( int i = 0 ; i < m ; ++i )
  ord[ i ] = i;
 std::sort( ord.begin() , ord.end() , [ & ]( int a , int b ) {
  return( v_p[ a ] * double( v_w[ b ] ) > v_p[ b ] * double( v_w[ a ] ) );
  } );
 std::vector< long > w( m );
 std::vector< double > p( m );
 for( int i = 0 ; i < m ; ++i ) {
  w[ i ] = v_w[ ord[ i ] ];
  p[ i ] = v_p[ ord[ i ] ];
  }

 // forward non-dominated DP: steps[ k ] is the feasible ( weight <= C ) Pareto
 // frontier after the first k items, sorted by ascending weight; each State
 // keeps its parent index in steps[ k - 1 ] for backtracking
 std::vector< std::vector< State > > steps;
 steps.push_back( { State{ 0 , 0.0 , -1 , false } } );
 for( int k = 0 ; k < m ; ++k ) {
  const std::vector< State > & prev = steps.back();
  std::vector< State > cand;
  cand.reserve( 2 * prev.size() );
  for( int j = 0 ; j < int( prev.size() ) ; ++j ) {
   cand.push_back( State{ prev[ j ].wsum , prev[ j ].psum , j , false } );
   if( prev[ j ].wsum + w[ k ] <= f_C )
    cand.push_back( State{ prev[ j ].wsum + w[ k ] ,
                           prev[ j ].psum + p[ k ] , j , true } );
   }
  std::sort( cand.begin() , cand.end() ,
             []( const State & a , const State & b ) {
              return( a.wsum != b.wsum ? a.wsum < b.wsum : a.psum > b.psum );
              } );
  std::vector< State > cur;
  cur.reserve( cand.size() );
  double frontier = - Inf< double >();
  for( const State & st : cand )
   if( st.psum > frontier ) { cur.push_back( st ); frontier = st.psum; }
  steps.push_back( std::move( cur ) );
  }

 // continuous items sorted by non-increasing efficiency (for the fractional
 // fill g( R ) = profit of the best fractional packing of capacity R)
 const int mc = int( v_cw.size() );
 std::vector< int > cord( mc );
 for( int i = 0 ; i < mc ; ++i )
  cord[ i ] = i;
 std::sort( cord.begin() , cord.end() , [ & ]( int a , int b ) {
  return( v_cp[ a ] * double( v_cw[ b ] ) > v_cp[ b ] * double( v_cw[ a ] ) );
  } );
 std::vector< long > cw( mc );
 std::vector< double > cp( mc );
 for( int i = 0 ; i < mc ; ++i ) {
  cw[ i ] = v_cw[ cord[ i ] ];
  cp[ i ] = v_cp[ cord[ i ] ];
  }

 // fractional knapsack value of the continuous items for a (double) capacity:
 // the capacity is kept fractional, matching DPBinaryKnapsackSolver's greedy
 auto cont_g = [ & ]( double R ) -> double {
  double val = 0;
  double rem = R;
  for( int k = 0 ; ( k < mc ) && ( rem > 0 ) ; ++k ) {
   if( double( cw[ k ] ) <= rem ) { val += cp[ k ]; rem -= cw[ k ]; }
   else { val += cp[ k ] * ( rem / double( cw[ k ] ) ); break; }
   }
  return( val );
  };

 // combine: best frontier state by profit + continuous fill of the rest
 const std::vector< State > & frontier = steps.back();
 double best = - Inf< double >();
 int bestpos = 0;
 for( int j = 0 ; j < int( frontier.size() ) ; ++j ) {
  const double val = frontier[ j ].psum + cont_g( f_Cd - frontier[ j ].wsum );
  if( val > best ) { best = val; bestpos = j; }
  }

 // backtrack the integer selection of the best state into in (core order)
 std::vector< char > best_in( m , 0 );
 int pos = bestpos;
 for( int k = m ; k >= 1 ; --k ) {
  const State & st = steps[ k ][ pos ];
  if( st.took )
   best_in[ k - 1 ] = 1;
  pos = st.parent;
  }
 for( int i = 0 ; i < m ; ++i )
  in[ ord[ i ] ] = best_in[ i ];

 // fractional fill of the residual capacity into cx (continuous-core order)
 cx.assign( mc , 0.0 );
 double rem = f_Cd - frontier[ bestpos ].wsum;
 for( int k = 0 ; ( k < mc ) && ( rem > 0 ) ; ++k ) {
  if( double( cw[ k ] ) <= rem ) { cx[ cord[ k ] ] = 1.0; rem -= cw[ k ]; }
  else { cx[ cord[ k ] ] = rem / double( cw[ k ] ); break; }
  }

 return( best );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void CoreDPBinaryKnapsackSolver::get_var_solution( Configuration * solc )
{
 if( pending_changes() )
  throw( std::invalid_argument(
   "CoreDPBinaryKnapsackSolver::get_var_solution: compute() must be called "
   "first" ) );

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 BKB->set_x( f_x.begin() );
 }

/*--------------------------------------------------------------------------*/
/*-------------- End File CoreDPBinaryKnapsackSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
