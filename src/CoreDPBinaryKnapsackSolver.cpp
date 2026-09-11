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
#include <cmath>
#include <functional>
#include <numeric>

#include "CoreDPBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef SURROGATE_BOUND
#define SURROGATE_BOUND 1
#endif
/* A2 surrogate / cardinality bound: tightens the global ceiling on
 * forced-cardinality (correlated) instances. Validated vs :MILPSolver and the
 * Pisinger reference optima. */

#ifndef SURROGATE_SOLVE
#define SURROGATE_SOLVE 1
#endif
/* COMBO surrogate-solve: solve the surrogate subproblem exactly to (possibly)
 * raise the incumbent. Requires SURROGATE_BOUND. Set to 0 to A/B against the
 * bound-only behaviour. */

#ifndef CORE_B2
#define CORE_B2 1                // guarded DP-extension (read-only check)
#endif
#ifndef CORE_B3
#define CORE_B3 1                // fixing-by-dominance
#endif
#ifndef CORE_B4
#define CORE_B4 1                // WB Dembo-Hammer per-item fixing
#endif
#ifndef CORE_B5
#define CORE_B5 1                // primal heuristics (PH/TPH/SSPH/SPH/GCH)
#endif
#ifndef CORE_STATS
#define CORE_STATS 0             // stderr counters per solve (dev only)
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

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

 enumerate_states();

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
  const long l = std::lround( w );
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
 if( v_cw.empty() )
  f_obj += solve_integer_core( in );
 else {
  std::vector< double > cx;
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
 // the extracted integer core, solved by the full break-item core enumeration
 return( core_enumerate( v_w , v_p , f_C , in , - Inf< double >() , false ) );
 }

/*--------------------------------------------------------------------------*/

double CoreDPBinaryKnapsackSolver::core_enumerate(
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
 in.assign( m , 0 );
 if( m == 0 )                  // no free item to decide
  return( 0 > lb ? 0.0 : - Inf< double >() );

 // aggregation: identical ( w , p ) items collapse into one registry entry
 // with a multiplicity. Duplicates are found by grouping a ( w , p )-sorted
 // permutation - no hashing and no per-item containers, the group's input
 // indices are the byk[ mstart , mstart + mcount ) span
 std::vector< int > byk( m );
 for( int i = 0 ; i < m ; ++i )
  byk[ i ] = i;
 std::sort( byk.begin() , byk.end() , [ & ]( int a , int b ) {
  return( ( iw[ a ] != iw[ b ] ) ? ( iw[ a ] < iw[ b ] )
                                 : ( ip[ a ] < ip[ b ] ) );
  } );

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
 const double zmargin = intp ? ( double( gp ) - 1e-9 ) : 1e-9;
 auto snap = [ intp , gp ]( double ub ) -> double {
  return( ( intp && ( gp > 1 ) )
          ? std::floor( ub / double( gp ) ) * double( gp ) : ub );
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

#if SURROGATE_BOUND
 // the surrogate / cardinality bound is computed lazily, only once the Pareto
 // frontier grows past SURR_TRIGGER (i.e. the instance is actually hard,
 // COMBO's MINSET gate): this keeps it off the easy instances where it would
 // only add sorting overhead. It uses the (by then larger) incumbent z for
 // the N_min side and for the exact-solve pruning, so one retry is allowed
 // when the incumbent has improved after the first attempt (a stronger z can
 // turn a truncated solve into a certifying one)
 int surr_tries = 2;
 double surr_z = - Inf< double >();
 constexpr int SURR_TRIGGER = 2000;
#endif

 // B4 (WB item-fixing): the break efficiency is the fill rate of the
 // Dembo-Hammer reduction bound WB( i , 1 )
 const double eb = ( bE < mE ) ? ( pE[ bE ] / double( wE[ bE ] ) ) : 0.0;

 // states are deviations from the break solution: unprocessed left copies are
 // implicitly in, unprocessed right copies implicitly out (so the initial
 // weight is wsumb). steps[ k ] is the Pareto frontier after k expansions,
 // sorted by ascending weight; each State keeps its parent index in
 // steps[ k - 1 ], and step_item[ k ] / step_mult[ k ] the toggled item and
 // its signed copy count (binary-decomposed batch, < 0 = dropped left copies)
 std::vector< std::vector< State > > steps;
 std::vector< int > step_item;
 std::vector< long > step_mult;
 steps.push_back( { State{ wsumb , psumb , -1 , false } } );
 step_item.push_back( -1 );
 step_mult.push_back( 0 );

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
  best_in = cntb;
  if( delta )
   best_in[ item ] += delta;
  for( int e = 0 ; e < nx ; ++e ) {
   const int j = std::abs( xs[ e ] ) - 1;
   best_in[ j ] += ( xs[ e ] > 0 ) ? 1 : -1;
   }
  for( int kk = int( steps.size() ) - 1 ; kk >= 1 ; --kk ) {
   const State & ps = steps[ kk ][ pos ];
   if( ps.took )
    best_in[ step_item[ kk ] ] += step_mult[ kk ];
   pos = ps.parent;
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
 // generated no new state (valid when the profit is >= 0)
 long   domRw = 0 , domLw = 0;
 double domRp = -1 , domLp = -1;

 // B5 (pairing heuristics) state: the PH/TPH/SSPH trigger is a frontier size
 // threshold doubled per call; SPH runs after every expansion and uses a
 // deterministic xorshift32 (reproducible solves) to sample one item per block
 std::size_t ph_thresh = 10 * std::size_t( m );
 unsigned rng = 0x9E3779B9u;
 auto rnd = [ & ]( unsigned k ) -> unsigned {
  rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
  return( rng % k );
  };

 // best frontier state under a weight cap: the frontier is sorted by weight
 // with profit increasing along it, so it is the last state below the cap
 auto best_under = [ & ]( const std::vector< State > & f , long wcap ) -> int {
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

#if CORE_B4
  // B4: Dembo-Hammer reduction bound. WB( it , 1 ) bounds the value of any
  // solution deviating on `it` (the LP value function is concave with slope
  // <= eb past the break point); if it cannot beat the incumbent the item is
  // fixed at its break-side value and never expanded
  {
   const double WB = isRight
    ? psumb + p[ it ] + double( Ceff - wsumb - w[ it ] ) * eb
    : psumb - p[ it ] + double( Ceff - wsumb + w[ it ] ) * eb;
   if( WB < z + zmargin )
    continue;
   }
#endif

#if CORE_B3
  // B3: fixing-by-dominance. once expanding item i gave no new state, any
  // heavier right item with p' <= p_i * floor( w' / w_i ) (resp. lighter left
  // item with p' >= p_i) cannot yield an improved solution either - floor(
  // w' / w_i ) copies of i dominate it state-wise - so skip it outright
  if( isRight ) {
   if( ( domRp >= 0 ) && ( w[ it ] >= domRw ) &&
       ( p[ it ] <= domRp * double( w[ it ] / domRw ) ) )
    continue;
   }
  else
   if( ( domLp >= 0 ) && ( w[ it ] <= domLw ) && ( p[ it ] >= domLp ) )
    continue;
#endif

  // LP bound fill rates of the still-unprocessed units. While copies of the
  // current item remain uncommitted they are themselves the most efficient
  // addable (resp. the least efficient droppable) units, so until the last
  // batch the bound must use the current item's own efficiency
  const double eff_t  = ( t < mA ) ? ( p[ t ] / double( w[ t ] ) ) : 0.0;
  const double eff_s  = ( s >= 0 ) ? ( p[ s ] / double( w[ s ] ) ) : 0.0;
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
    hasL = ( s >= 0 );
    eff_nl = eff_s;
    }
   else {
    eff_nr = eff_t;
    hasL = ( ! lastB ) || ( s >= 0 );
    eff_nl = lastB ? eff_s : eff_it;
    }

   const std::vector< State > & prev = steps.back();
   const int np = int( prev.size() );
   const long   dmul = isRight ? mult : - mult;
   const long   dw = dmul * w[ it ];
   const double dp = double( dmul ) * p[ it ];

   // the per-state survival test: not bound-pruned (fractional fill for
   // feasible states, fractional removal for over-capacity ones)
   auto survives = [ & ]( long nw , double npp ) -> bool {
    const double B = ( nw <= Ceff )
     ? npp + double( Ceff - nw ) * eff_nr
     : ( hasL ? npp - double( nw - Ceff ) * eff_nl : - Inf< double >() );
    return( B >= z + zmargin );
    };

#if CORE_B2
   // B2: after `guard` consecutive no-new-state expansions on this side, run
   // a read-only check (a scan of the would-be merge that writes nothing): if
   // no toggled state would survive, skip the whole item. A wasteful check (a
   // survivor exists) doubles the threshold; 40 consecutive skips force an
   // expansion so that the state pruning still operates periodically
   if( ( ! mdone ) && ( nonew > guard ) && ( roSkip < 40 ) ) {
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
     break;
     }
    guard *= 2;
    }
   roSkip = 0;
#endif

   // expansion: linear merge of the frontier with its ( +dw , +dp ) translate
   // (both sorted by ascending weight), keeping the merged Pareto frontier;
   // candidates with equal weight are merged higher-profit-first
   std::vector< State > cur;
   cur.reserve( 2 * np );
   double frontier = - Inf< double >();
   int newStates = 0;
   bool improved = false;

   int ia = 0 , ib = 0;
   while( ( ia < np ) || ( ib < np ) ) {
    bool tog;
    if( ia >= np )
     tog = true;
    else if( ib >= np )
     tog = false;
    else {
     const long wa = prev[ ia ].wsum , wb = prev[ ib ].wsum + dw;
     tog = ( wb < wa ) ||
           ( ( wb == wa ) && ( prev[ ib ].psum + dp > prev[ ia ].psum ) );
     }
    long nw; double npp; int par; bool took;
    if( tog )
     { nw = prev[ ib ].wsum + dw; npp = prev[ ib ].psum + dp;
       par = ib; took = true; ++ib; }
    else
     { nw = prev[ ia ].wsum; npp = prev[ ia ].psum;
       par = ia; took = false; ++ia; }

    if( npp <= frontier )       // Pareto-dominated (weight already >=)
     continue;
    frontier = npp;

    // update the incumbent from feasible states, recording the achieving
    // solution now (the state may be bound-pruned right after)
    if( ( nw <= Ceff ) && ( npp > z ) ) {
     record_incumbent( npp , par , it , took ? dmul : 0 , nullptr , 0 );
     improved = true;
     }

    if( ! survives( nw , npp ) )
     continue;

    cur.push_back( State{ nw , npp , par , took } );
    if( took )
     ++newStates;
    }

   bool grown;
   if( newStates ) {
    steps.push_back( std::move( cur ) );
    step_item.push_back( it );
    step_mult.push_back( dmul );
    nonew = 0;
    grown = true;
    }
   else {
    // the frontier did not grow: do not store a step (nothing to backtrack
    // through), just compact the last one in place under the current bound
    // (parents of its entries are untouched); on a first-batch failure also
    // remember the item as the dominance reference (B3)
    ++nonew;
    auto & last = steps.back();
    last.erase( std::remove_if( last.begin() , last.end() ,
                                [ & ]( const State & st ) {
                                 return( ! survives( st.wsum , st.psum ) );
                                 } ) , last.end() );
    if( ! mdone ) {
     if( isRight ) { domRw = w[ it ]; domRp = p[ it ]; }
     else          { domLw = w[ it ]; domLp = p[ it ]; }
     }
    grown = false;
    }

   if( steps.back().empty() ) {
    dead = true;
    break;
    }

#if CORE_B5
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
    const std::vector< State > & f = steps.back();
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

   // B5 (PH / TPH / SSPH): when the frontier has grown past the (doubling)
   // threshold, pair every still-unprocessed item with its best frontier
   // state via binary search; if the number of ( left , right ) pairs is
   // small wrt |S| / log2 |S| also try every pair (TPH); on huge frontiers
   // additionally enumerate all subsets of k random unprocessed items with
   // alpha k 2^k <= |S| , alpha = 20 (SSPH)
   if( steps.back().size() >= ph_thresh ) {
    ph_thresh *= 2;
    const std::vector< State > & f = steps.back();
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
#endif

   if( ! grown )               // L10 / T11: close the item early
    break;
   }                           // end of the binary-decomposed batch loop

  if( dead )
   break;

#if SURROGATE_BOUND
  // hard instance detected: tighten the ceiling with the surrogate bound,
  // computed on the per-unit expanded view; with SURROGATE_SOLVE the exact
  // optimum of the surrogate subproblem both certifies and (when it attains
  // the forced cardinality) raises the incumbent
  if( ( ! relx ) && surr_tries &&
      ( int( steps.back().size() ) > SURR_TRIGGER ) && ( z > surr_z ) ) {
#if SURROGATE_SOLVE
   const double zb4 = z;
   std::fill( selE.begin() , selE.end() , 0 );
   const double su = surrogate_solve( wE , pE , mE , bE , Ceff , z , selE );
   if( z > zb4 ) {             // improved: fold the unit selection to counts
    beat = true;
    std::fill( best_in.begin() , best_in.end() , 0 );
    for( int u = 0 ; u < mE ; ++u )
     if( selE[ u ] )
      ++best_in[ uitem( u ) ];
    }
   dantzig = std::min( dantzig , snap( su ) );
#else
   dantzig = std::min( dantzig ,
                       snap( surrogate_bound( wE , pE , mE , bE , Ceff ,
                                              z ) ) );
#endif
   --surr_tries;
   surr_z = z;
   }
#endif
  }

#if CORE_STATS
 {
  std::size_t mxf = 0 , totf = 0;
  for( const auto & st : steps )
   { mxf = std::max( mxf , st.size() ); totf += st.size(); }
  fprintf( stderr ,
           "CORE_STATS m=%d mA=%d mE=%d bE=%d steps=%zu maxfront=%zu "
           "totstates=%zu z=%.0f dantzig=%.0f closed=%s left=%d right=%d\n" ,
           m , mA , mE , bE , steps.size() , mxf , totf , z , dantzig ,
           ( z >= dantzig - 1e-9 ) ? "bound" : "exhaust" , s + 1 , mA - t );
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
  Uex += intp ? std::floor( add ) : add;
  }

 // U'' : the break item b is taken; ( w[ b ] - r ) of weight must be removed
 // from the previous (least efficient included) item b-1, rounded up
 double Uin = psumb + p[ b ];
 if( b >= 1 ) {
  const double rem = double( w[ b ] - r ) *
                     ( p[ b - 1 ] / double( w[ b - 1 ] ) );
  Uin -= intp ? std::ceil( rem ) : rem;
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
 // combo.c surbin's gradient binary search is used to locate the minimiser.
 double best = Inf< double >();
 long bestsur = 0;
 std::vector< int > ord( m );

 while( slo <= shi ) {
  const long s = slo + ( shi - slo ) / 2;     // floor towards slo (>=0 step)

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
  std::sort( ord.begin() , ord.end() , [ & ]( int a , int c ) {
   return( p[ a ] * double( w[ c ] + s ) > p[ c ] * double( w[ a ] + s ) );
   } );

  // greedy fractional fill of the surrogate capacity
  int d = int( m - ord.size() );               // items already taken (ww <= 0)
  double r = cap;
  int bi = -1;
  for( int t = 0 ; t < int( ord.size() ) ; ++t ) {
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
 std::sort( ws.begin() , ws.end() );
 int card1 = 0; long acc = 0;
 while( ( card1 < m ) && ( acc + ws[ card1 ] <= C ) )
  { acc += ws[ card1 ]; ++card1; }

 // card2 = N_min: smallest k with the k largest profits exceeding z
 std::vector< double > ps( p );
 std::sort( ps.begin() , ps.end() , std::greater< double >() );
 int card2 = 0; double pacc = 0;
 while( ( card2 < m ) && ( pacc <= z ) ) { pacc += ps[ card2 ]; ++card2; }

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
 int card1 = 0; long acc = 0;
 { std::vector< long > ws( w ); std::sort( ws.begin() , ws.end() );
   while( ( card1 < m ) && ( acc + ws[ card1 ] <= C ) )
  { acc += ws[ card1 ]; ++card1; } }
 int card2 = 0; double pacc = 0;
 { std::vector< double > ps( p );
   std::sort( ps.begin() , ps.end() , std::greater< double >() );
   while( ( card2 < m ) && ( pacc <= z ) ) { pacc += ps[ card2 ]; ++card2; } }

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
