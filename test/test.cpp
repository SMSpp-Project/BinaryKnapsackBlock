/*--------------------------------------------------------------------------*/
/*---------------------------- File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of BinaryKnapsackBlock and of the Solver of its module, on
 * instances built in memory and needing nothing but the core.
 *
 * The problem encoded by a BinaryKnapsackBlock is
 *
 *     max ( or min ) sum_i P_i x_i   s.t.   sum_i W_i x_i <= C ,
 *                                           x_i in { 0 , 1 } ( or [ 0 , 1 ] )
 *
 * with no restriction on the signs of W, P and C, some of the x possibly
 * fixed. The instances are small enough for the optimum to be found by
 * enumerating the 0-1 items, the continuous ones being filled greedily by
 * efficiency for each assignment of the others (which is the exact optimum
 * of the continuous part), and the same enumeration with every item
 * continuous gives the continuous relaxation.
 *
 * Every exact Solver of the module (DPBinaryKnapsackSolver, also with the
 * reoptimization on, ParallelDPBinaryKnapsackSolver with each of its engines,
 * CoreDPBinaryKnapsackSolver with all its reoptimization on, and
 * RECORDBinaryKnapsackSolver and COMBOBinaryKnapsackSolver the same way when
 * the build has them) has to find that optimum, with a Solution
 * feasible and worth it, or say that there is none; the relaxation Solver
 * (GreedyRelaxationBinaryKnapsackSolver and its incremental variant) have to
 * find the value of the relaxation, a bound on the right side of the optimum
 * and a feasible rounded Solution worth the other bound. The cases are the
 * corners of the data (capacity 0, no item fitting, every item fitting,
 * items of weight 0, negative data, no item at all), the fixings, the two
 * senses, a sequence of changes on one Block with every Solver attached, the
 * same changes made through the abstract representation, the netCDF round
 * trip, the abstract representation itself, the R3 copy, the conditional
 * bounds of the Block and the branching of the relaxation Solver.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "BinaryKnapsackBlock.h"
#include "CoreDPBinaryKnapsackSolver.h"
#include "DPBinaryKnapsackSolver.h"
#include "GreedyRelaxationBinaryKnapsackSolver.h"
#include "IncrementalGreedyRelaxationBinaryKnapsackSolver.h"
#include "ParallelDPBinaryKnapsackSolver.h"

// assert() checks in every build type, the Release one included, as the
// TestAssert.h of the core does: it comes after every header of the library,
// which has to be read as the library was compiled
#undef NDEBUG
#include <cassert>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- DATA -------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

using Index = Block::Index;
using dVec = std::vector< double >;

const double tol = 1e-7;   // relative accuracy asked of the values

int failures = 0;          // number of checks that failed

/// records the outcome of a check, printing it if it failed

bool check( bool ok , const std::string & what )
{
 if( ! ok ) {
  std::cout << "  KO: " << what << std::endl;
  ++failures;
  }
 return( ok );
 }

bool close( double a , double b )
{
 if( std::isinf( a ) || std::isinf( b ) )
  return( a == b );
 return( std::abs( a - b ) <= tol * std::max( 1.0 , std::abs( b ) ) );
 }

std::string str( double v ) { return( std::to_string( v ) ); }

/*--------------------------------------------------------------------------*/
/// the data of one instance, as the BinaryKnapsackBlock holds it

struct Inst {
 double C = 0;
 dVec W , P;
 std::vector< bool > I;              // true for a 0-1 item
 std::vector< unsigned char > fxd;   // 0 free, 1 fixed to 0, 2 fixed to 1
 bool max = true;

 Index n( void ) const { return( W.size() ); }
 };

/// the data of the Block, read back from it

Inst read( const BinaryKnapsackBlock & b )
{
 Inst in;
 in.C = b.get_Capacity();
 in.W = b.get_Weights();
 in.P = b.get_Profits();
 in.I = b.get_Integrality();
 in.fxd = b.get_fxd();
 in.max = ( b.get_objective_sense() == Objective::eMax );
 return( in );
 }

/*--------------------------------------------------------------------------*/
/// the optimum of the instance, - INF (in the max sense) if it is empty
/** The 0-1 free items are enumerated; for each assignment the continuous free
 * items are taken by the greedy fill of the residual capacity, once the
 * sign-trivial ones are settled and the ones with negative weight and profit
 * are complemented, which is the optimum of the continuous part. With
 * \p relax every item is continuous. The value is in the max sense, the
 * profits being changed in sign for a min problem. */

double opt_max( const Inst & in , bool relax = false )
{
 const Index n = in.n();
 std::vector< Index > bin;
 std::vector< Index > cnt;
 double C0 = in.C , P0 = 0;
 auto p = [ & ]( Index i ) { return( in.max ? in.P[ i ] : - in.P[ i ] ); };

 for( Index i = 0 ; i < n ; ++i ) {
  const unsigned char f = i < in.fxd.size() ? in.fxd[ i ] : 0;
  if( f == 1 )
   continue;
  if( f == 2 ) {
   C0 -= in.W[ i ];
   P0 += p( i );
   continue;
   }
  if( in.I[ i ] && ( ! relax ) )
   bin.push_back( i );
  else
   cnt.push_back( i );
  }

 // the continuous part: settle and complement, then the candidates
 double Cc = 0 , Pc = 0;
 std::vector< std::pair< double , double > > cand;   // ( w , p ) , both > 0
 for( auto i : cnt ) {
  const double w = in.W[ i ] , pi = p( i );
  if( ( w <= 0 ) && ( pi >= 0 ) ) { Cc -= w; Pc += pi; continue; }
  if( ( w >= 0 ) && ( pi <= 0 ) ) continue;
  if( w < 0 ) {   // and pi < 0: taken, its complement is the candidate
   Cc -= w; Pc += pi;
   cand.emplace_back( - w , - pi );
   }
  else
   cand.emplace_back( w , pi );
  }
 std::sort( cand.begin() , cand.end() , []( auto & a , auto & b ) {
  return( a.second * b.first > b.second * a.first ); } );

 double best = - INFINITY;
 const Index nb = bin.size();
 for( unsigned long mask = 0 ; mask < ( 1UL << nb ) ; ++mask ) {
  double R = C0 + Cc , V = P0 + Pc;
  for( Index k = 0 ; k < nb ; ++k )
   if( mask & ( 1UL << k ) ) {
    R -= in.W[ bin[ k ] ];
    V += p( bin[ k ] );
    }
  if( R < 0 )
   continue;
  for( auto & [ w , pi ] : cand ) {
   if( R <= 0 )
    break;
   const double t = std::min( 1.0 , R / w );
   V += t * pi;
   R -= t * w;
   }
  best = std::max( best , V );
  }
 return( best );
 }

/// the optimum in the sense of the instance, NaN if it is empty

double optimum( const Inst & in , bool relax = false )
{
 const double v = opt_max( in , relax );
 if( std::isinf( v ) )
  return( NAN );
 return( in.max ? v : - v );
 }

/// the value of a solution, and whether it is feasible up to the tolerance

std::pair< double , bool > evaluate( const Inst & in , const dVec & x )
{
 if( x.size() != in.n() )
  return( std::make_pair( NAN , false ) );
 double v = 0 , w = 0;
 bool ok = true;
 for( Index i = 0 ; i < in.n() ; ++i ) {
  if( ( x[ i ] < - tol ) || ( x[ i ] > 1 + tol ) )
   ok = false;
  if( in.I[ i ] && ( std::abs( x[ i ] - std::round( x[ i ] ) ) > tol ) )
   ok = false;
  const unsigned char f = i < in.fxd.size() ? in.fxd[ i ] : 0;
  if( ( ( f == 1 ) && ( std::abs( x[ i ] ) > tol ) ) ||
      ( ( f == 2 ) && ( std::abs( x[ i ] - 1 ) > tol ) ) )
   ok = false;
  v += in.P[ i ] * x[ i ];
  w += in.W[ i ] * x[ i ];
  }
 if( w > in.C + tol * std::max( 1.0 , std::abs( in.C ) ) )
  ok = false;
 return( std::make_pair( v , ok ) );
 }

/*--------------------------------------------------------------------------*/
/// the Solver of the module, each with its parameters

struct Slv {
 std::string name;
 Solver * s;
 bool exact;
 };

void set_int( Solver * s , const std::string & name , int value )
{
 s->set_par( s->int_par_str2idx( name ) , value );
 }

void set_dbl( Solver * s , const std::string & name , double value )
{
 s->set_par( s->dbl_par_str2idx( name ) , value );
 }

/// every Solver of the module, registered to the Block

std::vector< Slv > attach( BinaryKnapsackBlock * b , bool integer = true )
{
 std::vector< Slv > v;
 auto add = [ & ]( const std::string & name , const std::string & cls ,
                   bool exact ) -> Solver * {
  auto s = Solver::new_Solver( cls );
  assert( s );
  b->register_Solver( s );
  v.push_back( { name , s , exact } );
  return( s );
  };

 // the DP ones need integer weights, the relaxation ones do not
 if( integer ) {
  add( "DP" , "DPBinaryKnapsackSolver" , true );
  set_dbl( add( "DP reopt 0.5" , "DPBinaryKnapsackSolver" , true ) ,
           "dblReopt" , 0.5 );
  set_dbl( add( "DP reopt 1" , "DPBinaryKnapsackSolver" , true ) ,
           "dblReopt" , 1 );
  add( "ParallelDP auto" , "ParallelDPBinaryKnapsackSolver" , true );
  for( int e = 0 ; e <= 3 ; ++e ) {
   auto s = add( "ParallelDP engine " + std::to_string( e ) ,
                 "ParallelDPBinaryKnapsackSolver" , true );
   set_int( s , "intWhichParallel" , e );
   set_int( s , "intMaxThread" , 2 );
   }
  // the core DP with all its reoptimization, and the external solvers the
  // same way, only built if their sources were given
  for( std::string cls : { "CoreDPBinaryKnapsackSolver" ,
                           "RECORDBinaryKnapsackSolver" ,
                           "COMBOBinaryKnapsackSolver" } )
   if( Solver::has_Solver( cls ) )
    set_int( add( cls.substr( 0 , cls.find( "Binary" ) ) , cls , true ) ,
             "intReopt" , 3 );
  }
 add( "Greedy" , "GreedyRelaxationBinaryKnapsackSolver" , false );
 add( "IncrementalGreedy" , "IncrementalGreedyRelaxationBinaryKnapsackSolver" ,
      false );
 return( v );
 }

/// a Block holding the instance

BinaryKnapsackBlock * build( const Inst & in )
{
 auto b = new BinaryKnapsackBlock();
 b->load( in.n() , in.C , in.W , in.P , in.I , in.fxd );
 if( ! in.max )
  b->set_objective_sense( false );
 return( b );
 }

/*--------------------------------------------------------------------------*/
/// solves with every Solver and checks each against the enumeration

void solve_all( BinaryKnapsackBlock * b , std::vector< Slv > & slv ,
                const std::string & label )
{
 const Inst in = read( *b );
 const double opt = optimum( in );
 const double rel = optimum( in , true );
 const bool empty = std::isnan( opt );
 const std::string at = label + ", optimum " + str( opt );

 // the Block has to know whether it is empty, and its conditional bounds
 // have to hold on the two sides of the optimum
 check( b->is_empty() == empty , at + ": is_empty() says " +
        std::to_string( b->is_empty() ) );
 if( ! empty ) {
  const double lb = b->get_valid_lower_bound( true );
  const double ub = b->get_valid_upper_bound( true );
  check( lb <= opt + tol * std::max( 1.0 , std::abs( opt ) ) ,
         at + ": the lower bound of the Block is " + str( lb ) );
  check( ub >= opt - tol * std::max( 1.0 , std::abs( opt ) ) ,
         at + ": the upper bound of the Block is " + str( ub ) );
  }

 for( auto & sl : slv ) {
  const std::string who = at + ", " + sl.name;
  int status;
  try {
   status = sl.s->compute( false );
   }
  catch( std::exception & e ) {
   check( false , who + ": compute() throws " + e.what() );
   continue;
   }

  if( empty ) {
   check( status == Solver::kInfeasible , who + ": the instance is empty, "
          "status " + std::to_string( status ) );
   continue;
   }
  if( ! check( status == Solver::kOK , who + ": status " +
               std::to_string( status ) ) )
   continue;

  const double value = sl.s->get_var_value();
  const double lb = sl.s->get_lb();
  const double ub = sl.s->get_ub();

  if( sl.exact ) {
   check( close( value , opt ) , who + ": value " + str( value ) );
   check( close( lb , opt ) && close( ub , opt ) , who + ": bounds " +
          str( lb ) + " , " + str( ub ) );
   }
  else {
   // the value of the relaxation, and the bounds on the two sides of the
   // optimum, one of them being the value of the rounded solution
   check( close( value , rel ) , who + ": relaxation value " + str( value ) +
          ", not " + str( rel ) );
   const double s = tol * std::max( 1.0 , std::abs( opt ) );
   check( ( lb <= opt + s ) && ( ub >= opt - s ) , who + ": bounds " +
          str( lb ) + " , " + str( ub ) );
   check( in.max ? close( ub , value ) : close( lb , value ) , who +
          ": the bound of the relaxation is not its value" );
   }

  // the Solution: feasible, and worth what the Solver says
  Solution * sol = nullptr;
  if( sl.exact )
   sol = sl.s->get_Solution();
  else
   sol = dynamic_cast< RelaxationSolver * >( sl.s )->get_true_solution();
  auto ks = dynamic_cast< BinaryKnapsackSolution * >( sol );
  // with no item a Solver of the BinaryKnapsackSolver family has no
  // Solution to give, which it says with nullptr
  if( ( ! ks ) && ( in.n() == 0 ) )
   continue;
  if( check( ks != nullptr , who + ": no BinaryKnapsackSolution" ) ) {
   auto [ v , ok ] = evaluate( in , ks->get_x() );
   check( ok , who + ": the Solution is not feasible" );
   const double bound = in.max ? lb : ub;
   check( close( v , bound ) , who + ": the Solution is worth " + str( v ) +
          ", not " + str( bound ) );
   // the check of the Block agrees, on a solution with no fraction in it
   if( std::all_of( in.I.begin() , in.I.end() , []( bool i ) { return i; } ) )
    check( b->is_sol_feasible( ks ) , who + ": is_sol_feasible() says no" );
   }
  delete sol;
  }
 }

/// a random instance: n items, weights in [ wl , wu ], profits in [ pl , pu ]

Inst random_inst( std::mt19937 & rg , Index n , int wl , int wu ,
                  double pl , double pu , double cfrac , double ffrac ,
                  bool max )
{
 std::uniform_int_distribution< int > uw( wl , wu );
 std::uniform_real_distribution< double > up( pl , pu ) , u01( 0 , 1 );
 Inst in;
 in.max = max;
 double tw = 0;
 for( Index i = 0 ; i < n ; ++i ) {
  in.W.push_back( uw( rg ) );
  in.P.push_back( std::round( 2 * up( rg ) ) / 2 );   // halves: ties happen
  in.I.push_back( u01( rg ) >= cfrac );
  const double f = u01( rg );
  in.fxd.push_back( f < ffrac / 2 ? 1 : ( f < ffrac ? 2 : 0 ) );
  tw += std::abs( in.W.back() );
  }
 in.C = std::round( u01( rg ) * tw / 2 );
 return( in );
 }

/// one Block with every Solver attached, which each instance is loaded into

BinaryKnapsackBlock * shared_block = nullptr;
std::vector< Slv > shared_slv;

/// loads the instance in the shared Block and solves it with every Solver
/** The load() reaches the Solver as an NBModification, which they have to
 * take as a new instance; with \p integer false the instance goes in a Block
 * of its own, with the Solver that do not need integer weights only. */

void run( const Inst & in , const std::string & label , bool integer = true )
{
 if( ! integer ) {
  auto b = build( in );
  auto slv = attach( b , false );
  solve_all( b , slv , label );
  b->unregister_Solvers( true );
  delete b;
  return;
  }

 if( ! shared_block ) {
  shared_block = new BinaryKnapsackBlock();
  shared_slv = attach( shared_block );
  }
 shared_block->load( in.n() , in.C , in.W , in.P , in.I , in.fxd );
 shared_block->set_objective_sense( in.max );
 solve_all( shared_block , shared_slv , label );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- CASES ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// random instances of every kind: signs, continuous items, fixings, senses

void test_random( void )
{
 std::mt19937 rg( 17 );
 for( int k = 0 ; k < 120 ; ++k ) {
  const Index n = 1 + k % 12;
  const bool max = k % 2;
  // positive data first, then mixed signs, then with fixings and
  // continuous items
  Inst in;
  if( k < 40 )
   in = random_inst( rg , n , 1 , 20 , 1 , 30 , 0 , 0 , max );
  else if( k < 80 )
   in = random_inst( rg , n , -10 , 20 , -15 , 30 , 0 , 0 , max );
  else
   in = random_inst( rg , n , -10 , 20 , -15 , 30 , 0.3 , 0.3 , max );
  run( in , "random " + std::to_string( k ) );
  }
 }

/// capacity 0: only what weighs nothing or less fits

void test_capacity_zero( void )
{
 Inst in;
 in.C = 0;
 in.W = { 3 , 5 , 2 , 7 };
 in.P = { 4 , 6 , 1 , 9 };
 in.I.assign( 4 , true );
 run( in , "capacity 0, positive weights" );

 in.W = { 3 , -2 , 0 , 7 };
 in.P = { 4 , 1 , 2 , 9 };
 run( in , "capacity 0, a negative and a zero weight" );

 in.max = false;
 in.P = { -4 , 6 , 1 , -9 };
 run( in , "capacity 0, min" );
 }

/// every item weighs more than the capacity, so the optimum is 0

void test_all_too_heavy( void )
{
 Inst in;
 in.C = 4;
 in.W = { 5 , 9 , 6 , 12 , 7 };
 in.P = { 10 , 3 , 8 , 20 , 1 };
 in.I.assign( 5 , true );
 run( in , "every item above the capacity" );

 in.I = { true , false , true , false , true };
 run( in , "every item above the capacity, some continuous" );
 }

/// the capacity takes every item, so the optimum takes all the profitable ones

void test_all_fit( void )
{
 Inst in;
 in.C = 100;
 in.W = { 5 , 9 , 6 , 12 , 7 };
 in.P = { 10 , 3 , 8 , 20 , 1 };
 in.I.assign( 5 , true );
 run( in , "every item fitting" );

 in.C = 39;   // exactly the total weight
 in.P = { 10 , -3 , 8 , 20 , -1 };
 run( in , "every item fitting exactly, some unprofitable" );
 }

/// items of weight 0, taken if and only if profitable

void test_zero_weights( void )
{
 Inst in;
 in.C = 6;
 in.W = { 0 , 4 , 0 , 3 , 0 , 5 };
 in.P = { 3 , 5 , -2 , 4 , 0 , 6 };
 in.I.assign( 6 , true );
 run( in , "zero weights" );

 in.max = false;
 run( in , "zero weights, min" );

 in.I = { false , true , false , true , false , false };
 run( in , "zero weights, some continuous" );
 }

/// negative weights, profits and capacity, empty or not

void test_negative( void )
{
 Inst in;
 in.C = -3;
 in.W = { -4 , 2 , -1 , 5 , -2 };
 in.P = { -6 , 3 , 2 , 7 , -1 };
 in.I.assign( 5 , true );
 run( in , "negative capacity, feasible" );

 in.C = -8;
 run( in , "negative capacity, empty" );

 in.C = 5;
 in.W = { -3 , -5 , 4 , -1 };
 in.P = { -2 , -7 , -1 , -4 };
 in.I.assign( 4 , true );
 run( in , "every profit negative" );

 in.max = false;
 run( in , "every profit negative, min" );
 }

/// no item at all

void test_no_items( void )
{
 Inst in;
 in.C = 5;
 run( in , "no items" );

 in.C = -1;
 run( in , "no items, negative capacity" );
 }

/// items fixed to 0 and to 1, in the load() and by fix_x(), then unfixed

void test_fixed( void )
{
 Inst in;
 in.C = 10;
 in.W = { 6 , 5 , 4 , 3 , -2 , 7 };
 in.P = { 9 , 7 , 5 , 3 , 4 , -1 };
 in.I.assign( 6 , true );
 in.fxd = { 1 , 0 , 2 , 0 , 1 , 2 };
 run( in , "fixed in the load()" );

 in.fxd = { 2 , 2 , 0 , 0 , 0 , 0 };
 run( in , "fixed to 1 beyond the capacity" );

 // the conditional bounds of the Block take the fixed items at their value:
 // a profitable item of negative weight fixed to 0 is not in the lower bound,
 // an unprofitable one fixed to 1 is in the upper bound
 Inst sm;
 sm.C = 3;
 sm.W = { -1 , 2 , 1 };
 sm.P = { 5 , 1 , -2 };
 sm.I.assign( 3 , true );
 sm.fxd = { 1 , 0 , 2 };
 run( sm , "bounds with fixed items" );

 // the same fixings made on a Block with every Solver attached
 in.fxd.clear();
 auto b = build( in );
 auto slv = attach( b );
 solve_all( b , slv , "before the fixings" );
 b->fix_x( false , 0 );
 b->fix_x( true , 5 );
 solve_all( b , slv , "fix_x( i )" );
 std::vector< bool > val = { true , false };
 b->fix_x( val.begin() , Block::Range( 2 , 4 ) );
 solve_all( b , slv , "fix_x( Range )" );
 b->unfix_x( 0 );
 solve_all( b , slv , "unfix_x( i )" );
 // an unordered Subset: nms[ h ] takes value[ h ], whatever the order
 b->fix_x( val.begin() , Block::Subset{ 4 , 1 } );
 check( ( read( *b ).fxd[ 4 ] == 2 ) && ( read( *b ).fxd[ 1 ] == 1 ) ,
        "fix_x( unordered Subset ): item 4 fixed to " +
        std::to_string( read( *b ).fxd[ 4 ] - 1 ) + " and item 1 to " +
        std::to_string( read( *b ).fxd[ 1 ] - 1 ) + ", not to 1 and 0" );
 solve_all( b , slv , "fix_x( Subset )" );
 b->unfix_x( Block::Subset{ 2 , 5 } );
 solve_all( b , slv , "unfix_x( Subset )" );
 b->unfix_x();
 solve_all( b , slv , "unfix_x()" );
 check( read( *b ).fxd == std::vector< unsigned char >( 6 , 0 ) ,
        "unfix_x() leaves something fixed" );
 b->unregister_Solvers( true );
 delete b;
 }

/// both senses, and the change from one to the other

void test_sense( void )
{
 Inst in;
 in.C = 9;
 in.W = { 4 , -3 , 5 , 2 , 6 };
 in.P = { 3 , -2 , -6 , 5 , 4 };
 in.I = { true , true , true , false , true };
 auto b = build( in );
 auto slv = attach( b );
 solve_all( b , slv , "max" );
 b->set_objective_sense( false );
 solve_all( b , slv , "max changed to min" );
 b->set_objective_sense( true );
 solve_all( b , slv , "min changed back to max" );
 b->unregister_Solvers( true );
 delete b;
 }

/// a sequence of changes on one Block, every Solver re-solving after each

void test_reoptimization( void )
{
 std::mt19937 rg( 23 );
 std::uniform_real_distribution< double > u01( 0 , 1 );
 std::uniform_int_distribution< int > uw( -4 , 15 ) , up( -6 , 20 );
 for( int r = 0 ; r < 4 ; ++r ) {
  const Index n = 6 + 2 * r;
  Inst in = random_inst( rg , n , 1 , 15 , 1 , 20 , 0 , 0 , r % 2 );
  auto b = build( in );
  auto slv = attach( b );
  solve_all( b , slv , "run " + std::to_string( r ) + " start" );
  std::uniform_int_distribution< Index > ui( 0 , n - 1 );

  for( int step = 0 ; step < 40 ; ++step ) {
   const std::string at = "run " + std::to_string( r ) + " step " +
                          std::to_string( step );
   const Index i = ui( rg ) , j = ui( rg );
   const Index lo = std::min( i , j ) , hi = std::max( i , j ) + 1;
   dVec vals;
   for( Index k = lo ; k < hi ; ++k )
    vals.push_back( up( rg ) );
   std::string what;
   switch( step % 10 ) {
    case( 0 ): b->chg_profit( up( rg ) , i ); what = "chg_profit"; break;
    case( 1 ): b->chg_weight( uw( rg ) , i ); what = "chg_weight"; break;
    case( 2 ):
     b->chg_capacity( std::round( u01( rg ) * 40 ) - 5 );
     what = "chg_capacity";
     break;
    case( 3 ):
     b->chg_profits( Block::MF_dbl_sp( vals ) , Block::Range( lo , hi ) );
     what = "chg_profits( Range )";
     break;
    case( 4 ): {
     for( auto & v : vals ) v = uw( rg );
     b->chg_weights( Block::MF_dbl_sp( vals ) , Block::Range( lo , hi ) );
     what = "chg_weights( Range )";
     break;
     }
    case( 5 ):
     if( b->is_fixed( i ) ) { b->unfix_x( i ); what = "unfix_x"; }
     else { b->fix_x( u01( rg ) < 0.5 , i ); what = "fix_x"; }
     break;
    case( 6 ): {
     // an unordered Subset: nms[ h ] takes the h-th value
     if( i == j ) { b->chg_profit( up( rg ) , i ); what = "chg_profit"; break; }
     dVec v2 = { double( up( rg ) ) , double( up( rg ) ) };
     const Index hi2 = std::max( i , j ) , lo2 = std::min( i , j );
     b->chg_profits( Block::MF_dbl_sp( v2 ) , Block::Subset{ hi2 , lo2 } );
     check( ( b->get_Profit( hi2 ) == v2[ 0 ] ) &&
            ( b->get_Profit( lo2 ) == v2[ 1 ] ) , at + ": chg_profits( "
            "unordered Subset ) pairs the values with other items" );
     what = "chg_profits( Subset )";
     break;
     }
    case( 7 ): {
     dVec v2 = { double( uw( rg ) ) };
     b->chg_weights( Block::MF_dbl_sp( v2 ) , Block::Subset{ i } );
     what = "chg_weights( Subset )";
     break;
     }
    case( 8 ):
     b->set_objective_sense( b->get_objective_sense() != Objective::eMax );
     what = "set_objective_sense";
     break;
    case( 9 ):
     b->chg_integrality( ! b->get_Integrality( i ) , i );
     what = "chg_integrality";
     break;
    }
   solve_all( b , slv , at + " " + what );
   }
  b->unregister_Solvers( true );
  delete b;
  }
 }

/// the same changes made through the abstract representation

void test_abstract_modifications( void )
{
 Inst in;
 in.C = 12;
 in.W = { 5 , 4 , 6 , 3 , 2 , 7 };
 in.P = { 8 , 6 , 9 , 4 , 3 , 10 };
 in.I.assign( 6 , true );
 auto b = build( in );
 b->generate_abstract_constraints();
 b->generate_objective();
 auto slv = attach( b );
 solve_all( b , slv , "abstract, start" );

 auto obj = b->get_objective< FRealObjective >();
 auto lfo = static_cast< LinearFunction * >( obj->get_function() );
 auto cns = b->get_static_constraint< FRowConstraint >( 0 );
 assert( cns );
 auto lfc = static_cast< LinearFunction * >( cns->get_function() );

 // a C05FunctionModLinRngd on the objective, as FrankWolfeSolver issues
 lfo->modify_coefficient( 2 , -3 );
 check( b->get_Profit( 2 ) == -3 , "the profit does not follow its "
        "coefficient in the objective" );
 solve_all( b , slv , "abstract, a coefficient of the objective" );

 // a C05FunctionModLinSbst on the objective
 lfo->modify_coefficients( dVec{ 1 , 12 } , Block::Subset{ 0 , 4 } );
 solve_all( b , slv , "abstract, a Subset of the objective" );

 // a coefficient and a Range of the constraint
 lfc->modify_coefficient( 1 , 9 );
 check( b->get_Weight( 1 ) == 9 , "the weight does not follow its "
        "coefficient in the constraint" );
 solve_all( b , slv , "abstract, a coefficient of the constraint" );
 lfc->modify_coefficients( dVec{ 1 , -2 } , Block::Range( 3 , 5 ) );
 solve_all( b , slv , "abstract, a Range of the constraint" );

 // the right-hand side of the constraint
 cns->set_rhs( 7 );
 check( b->get_Capacity() == 7 , "the capacity does not follow the RHS" );
 solve_all( b , slv , "abstract, the RHS" );

 // a ColVariable fixed to 1 and one to 0, and one unfixed again
 b->get_Var( 5 )->set_value( 1 );
 b->get_Var( 5 )->is_fixed( true );
 b->get_Var( 0 )->set_value( 0 );
 b->get_Var( 0 )->is_fixed( true );
 check( b->is_fixed( 5 ) && b->is_fixed( 0 ) , "the Block does not see "
        "the fixed ColVariable" );
 solve_all( b , slv , "abstract, ColVariable fixed" );
 b->get_Var( 0 )->is_fixed( false );
 solve_all( b , slv , "abstract, ColVariable unfixed" );

 // the sense of the Objective
 obj->set_sense( Objective::eMin );
 check( b->get_objective_sense() == Objective::eMin , "the Block does not "
        "follow the sense of the Objective" );
 solve_all( b , slv , "abstract, the sense" );

 b->unregister_Solvers( true );
 delete b;
 }

/// the abstract representation says the same as the physical one

void test_abstract_representation( void )
{
 Inst in;
 in.C = 11;
 in.W = { 5 , 4 , -6 , 3 , 2 };
 in.P = { 8 , 6 , -9 , 4 , 3 };
 in.I = { true , true , true , false , true };
 in.max = false;
 for( auto & p : in.P ) p = - p;
 auto b = build( in );
 b->generate_abstract_constraints();
 b->generate_objective();

 auto & x = * b->get_static_variable_v< ColVariable >( 0 );
 check( x.size() == in.n() , "not one ColVariable per item" );
 for( Index i = 0 ; i < in.n() ; ++i )
  check( x[ i ].is_integer() == bool( in.I[ i ] ) , "ColVariable " +
         std::to_string( i ) + " of the wrong integrality" );

 auto cns = b->get_static_constraint< FRowConstraint >( 0 );
 auto lfc = static_cast< LinearFunction * >( cns->get_function() );
 auto obj = b->get_objective< FRealObjective >();
 auto lfo = static_cast< LinearFunction * >( obj->get_function() );
 check( cns->get_rhs() == in.C , "the RHS is not the capacity" );
 check( obj->get_sense() == Objective::eMin , "the Objective is not min" );
 for( Index i = 0 ; i < in.n() ; ++i )
  check( ( lfc->get_coefficient( i ) == in.W[ i ] ) &&
         ( lfo->get_coefficient( i ) == in.P[ i ] ) , "coefficient " +
         std::to_string( i ) + " is not the data" );

 // the solution written in the ColVariable, feasible either way and worth
 // the optimum in the Objective
 auto s = Solver::new_Solver( "DPBinaryKnapsackSolver" );
 b->register_Solver( s );
 check( s->compute() == Solver::kOK , "abstract: the DP does not solve" );
 s->get_var_solution();
 check( b->is_feasible( true ) && b->is_feasible( false ) ,
        "abstract: the solution is not feasible" );
 const double opt = optimum( in );
 check( close( b->get_objective_value() , opt ) , "abstract: the Objective "
        "is worth " + str( b->get_objective_value() ) + ", not " + str( opt ) );

 // the Solution of the Block reads the ColVariable back
 auto sol = b->get_Solution( nullptr , false );
 auto ks = dynamic_cast< BinaryKnapsackSolution * >( sol );
 if( check( ks != nullptr , "abstract: no BinaryKnapsackSolution" ) )
  check( close( evaluate( in , ks->get_x() ).first , opt ) ,
         "abstract: the Solution of the Block is not the optimum" );
 delete sol;

 // an infeasible point is seen as such
 for( Index i = 0 ; i < in.n() ; ++i )
  b->set_x( i , in.W[ i ] > 0 ? 1 : 0 );
 check( ! b->is_feasible( true ) , "abstract: an infeasible point passes" );

 b->unregister_Solvers( true );
 delete b;
 }

/// serialize and deserialize a Block through a netCDF file

void test_netcdf( void )
{
 const auto file = std::filesystem::temp_directory_path() /
  ( "BinaryKnapsackBlock_test_" + std::to_string( std::random_device()() ) +
    ".nc4" );

 auto round_trip = [ & ]( const Inst & in , const std::string & label ) {
  auto b = build( in );
  {
   netCDF::NcFile f( file.string() , netCDF::NcFile::replace );
   f.putAtt( "SMS++_file_type" , netCDF::NcInt() , eBlockFile );
   b->Block::serialize( f , eBlockFile );
   }
  delete b;

  auto nb = dynamic_cast< BinaryKnapsackBlock * >(
                                     Block::deserialize( file.string() ) );
  std::filesystem::remove( file );
  if( ! check( nb != nullptr , label + ": no BinaryKnapsackBlock back" ) )
   return;
  const Inst out = read( *nb );
  check( ( out.C == in.C ) && ( out.W == in.W ) && ( out.P == in.P ) &&
         ( out.I == in.I ) && ( out.max == in.max ) , label +
         ": the data are not the same" );
  auto slv = attach( nb );
  solve_all( nb , slv , label );
  nb->unregister_Solvers( true );
  delete nb;
  };

 Inst in;
 in.C = 13;
 in.W = { 5 , 4 , -6 , 3 , 2 , 0 };
 in.P = { 8.5 , 6 , -9 , 4 , 3 , 1 };
 in.I.assign( 6 , true );
 round_trip( in , "netCDF, max" );

 in.max = false;
 in.P = { -8.5 , 6 , -9 , -4 , 3 , -1 };
 in.I = { true , false , true , true , false , true };
 round_trip( in , "netCDF, min with continuous items" );

 in.C = 0;
 round_trip( in , "netCDF, capacity 0" );
 }

/// the R3 copy is the same problem, and the solutions go back and forth

void test_R3_copy( void )
{
 Inst in;
 in.C = 10;
 in.W = { 5 , 4 , 6 , 3 , -2 };
 in.P = { 8 , 6 , 9 , 4 , -3 };
 in.I = { true , true , false , true , true };
 in.fxd = { 0 , 0 , 0 , 2 , 0 };
 auto b = build( in );
 auto c = dynamic_cast< BinaryKnapsackBlock * >( b->get_R3_Block() );
 if( ! check( c != nullptr , "R3: no BinaryKnapsackBlock copy" ) ) {
  delete b;
  return;
  }
 const Inst out = read( *c );
 check( ( out.C == in.C ) && ( out.W == in.W ) && ( out.P == in.P ) &&
        ( out.I == in.I ) && ( out.fxd == in.fxd ) && ( out.max == in.max ) ,
        "R3: the copy has other data" );

 auto s = Solver::new_Solver( "DPBinaryKnapsackSolver" );
 c->register_Solver( s );
 check( s->compute() == Solver::kOK , "R3: the DP does not solve the copy" );
 check( close( s->get_var_value() , optimum( in ) ) , "R3: the copy has "
        "another optimum" );
 s->get_var_solution();
 b->map_back_solution( c );
 dVec x( in.n() );
 b->get_x( x.begin() );
 check( close( evaluate( in , x ).first , optimum( in ) ) ,
        "R3: map_back_solution() does not bring the optimum back" );

 // and forward again, onto a copy whose free ColVariable are cleared (a
 // fixed one keeps its value [see ColVariable::set_value()])
 for( Index i = 0 ; i < in.n() ; ++i )
  if( ! c->is_fixed( i ) )
   c->set_x( i , 0 );
 b->map_forward_solution( c );
 dVec y( in.n() );
 c->get_x( y.begin() );
 check( x == y , "R3: map_forward_solution() does not bring it forward" );

 c->unregister_Solvers( true );
 delete c;
 delete b;
 }

/// the DP Solver want integer weights, the relaxation ones do not

void test_fractional_weights( void )
{
 Inst in;
 in.C = 7.5;
 in.W = { 2.5 , 3 , 4.25 , 1 };
 in.P = { 5 , 4 , 7 , 1 };
 in.I.assign( 4 , true );
 auto b = build( in );
 for( std::string cls : { "DPBinaryKnapsackSolver" ,
                          "CoreDPBinaryKnapsackSolver" } ) {
  auto s = Solver::new_Solver( cls );
  bool thrown = false;
  try {
   b->register_Solver( s );
   s->compute();
   }
  catch( std::invalid_argument & ) { thrown = true; }
  check( thrown , cls + " takes weights that are not integer" );
  b->unregister_Solver( s , true );
  }
 delete b;

 run( in , "fractional weights" , false );
 }

/// bounds rounded on values that are integers in exact arithmetic
/** The core enumeration rounds its upper bounds to integers: with weights
 * that are multiples of 84 (or of 42) the values computed in floating point
 * once landed just off the integer, the rounding lost a unit and the
 * optimum 5202 came out as 5201; the same weights divided by 84 never did. */

void test_rounding( void )
{
 const std::vector< std::pair< int , int > > wp = {
  { 84 , 110 } , { 84 , 110 } , { 84 , 110 } , { 84 , 110 } , { 84 , 109 } ,
  { 84 , 109 } , { 84 , 109 } , { 336 , 434 } , { 420 , 542 } ,
  { 420 , 542 } , { 252 , 325 } , { 336 , 433 } , { 252 , 324 } ,
  { 168 , 216 } , { 252 , 324 } , { 252 , 324 } , { 84 , 108 } ,
  { 84 , 108 } , { 252 , 324 } , { 420 , 539 } };
 for( int s : { 1 , 2 , 84 } ) {
  Inst in;
  for( const auto & [ w , p ] : wp ) {
   in.W.push_back( w / s );
   in.P.push_back( p );
   in.I.push_back( true );
   in.fxd.push_back( 0 );
   }
  in.C = 4032 / s;
  run( in , "rounding, weights / " + std::to_string( s ) );
  }
 }

/// what each kind of change lets the core DP reuse of the previous solve

void test_reopt_outcome( void )
{
 // two efficient items fill the capacity, the third one is the break item
 Inst in;
 in.C = 2;
 in.W = { 1 , 1 , 2 };
 in.P = { 10 , 10 , 1 };
 in.I.assign( 3 , true );
 in.fxd.assign( 3 , 0 );
 auto b = build( in );
 auto s = dynamic_cast< CoreDPBinaryKnapsackSolver * >(
			   Solver::new_Solver( "CoreDPBinaryKnapsackSolver" ) );
 b->register_Solver( s );
 set_int( s , "intReopt" , 3 );
 check( s->compute() == Solver::kOK , "reopt: the first solve fails" );
 check( s->get_reopt_outcome() == 0 , "reopt: the first solve is warm" );

 // a taken item gains profit: no solve at all, no suspect item
 auto chg = [ & ]( Index i , double p , double z , int outcome ,
		   const std::string & what ) {
  std::vector< double > P = { b->get_Profit( 0 ) , b->get_Profit( 1 ) ,
			      b->get_Profit( 2 ) };
  P[ i ] = p;
  b->chg_profits( Block::MF_dbl_sp( P ) , Block::Range( 0 , 3 ) );
  check( s->compute() == Solver::kOK , "reopt: " + what + " fails" );
  check( std::abs( s->get_var_value() - z ) < 1e-9 ,
	 "reopt: " + what + " gives a wrong optimum" );
  check( s->get_reopt_outcome() == outcome ,
	 "reopt: " + what + " reuses the wrong amount" );
  };
 chg( 0 , 11 , 21 , 2 , "a taken item gaining profit" );

 // the third item gains a little: a suspect, which the Lagrangian bound
 // with the multiplier of the break item clears (it stays at 21)
 chg( 2 , 2 , 21 , 3 , "an untaken item gaining a little" );

 // the third item gains a lot: the bound cannot clear it, and the optimum
 // becomes the third item alone
 chg( 2 , 30 , 30 , 1 , "an untaken item gaining enough to enter" );

 b->unregister_Solvers( true );
 delete b;
 }

/// the two children of branch() cover the relaxation, then undo

void test_branch( void )
{
 Inst in;
 in.C = 10;
 in.W = { 5 , 4 , 6 , 3 , 2 };
 in.P = { 8 , 6 , 9 , 4 , 1 };
 in.I.assign( 5 , true );
 for( std::string cls : { "GreedyRelaxationBinaryKnapsackSolver" ,
                          "IncrementalGreedyRelaxationBinaryKnapsackSolver" } ) {
  auto b = build( in );
  auto s = Solver::new_Solver( cls );
  b->register_Solver( s );
  auto rs = dynamic_cast< RelaxationSolver * >( s );
  assert( rs );
  check( s->compute() == Solver::kOK , cls + ": the root does not solve" );
  const double root = s->get_var_value();

  auto children = rs->branch();
  check( children.size() == 2 , cls + ": not two children" );
  double best = - INFINITY;
  for( Index c = 0 ; c < children.size() ; ++c ) {
   auto chg = dynamic_cast< BinaryKnapsackBlockChange * >( children[ c ] );
   if( ! check( chg && ( chg->type() == BinaryKnapsackBlockChange::eFixX ) &&
                ( chg->num_items() == 1 ) , cls + ": a child is not a fixing" ) )
    continue;
   Inst ch = in;
   ch.fxd.assign( in.n() , 0 );
   ch.fxd[ chg->item( 0 ) ] = chg->data()[ 0 ] ? 2 : 1;

   auto undo = rs->apply( children[ c ] , true );
   check( s->compute() == Solver::kOK , cls + ": a child does not solve" );
   const double v = s->get_var_value();
   check( close( v , optimum( ch , true ) ) , cls + ": child " +
          std::to_string( c ) + " is worth " + str( v ) + ", not " +
          str( optimum( ch , true ) ) );
   check( v <= root + tol , cls + ": a child beats its root" );
   best = std::max( best , optimum( ch ) );
   check( read( *b ).fxd == std::vector< unsigned char >( in.n() , 0 ) ,
          cls + ": the fixing of the child reaches the Block" );

   rs->apply( undo );
   delete undo;
   check( s->compute() == Solver::kOK , cls + ": the root does not solve "
          "again" );
   check( close( s->get_var_value() , root ) , cls + ": the undo does not "
          "bring the root back" );
   }
  for( auto c : children )
   delete c;

  // the optimum is in one of the two children
  check( close( best , optimum( in ) ) , cls + ": the children lose the "
         "optimum" );
  b->unregister_Solvers( true );
  delete b;
  }
 }

}  // namespace

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::vector< std::pair< std::string , std::function< void() > > >
  cases = {
  { "random instances" , test_random } ,
  { "capacity 0" , test_capacity_zero } ,
  { "every item above the capacity" , test_all_too_heavy } ,
  { "every item fitting" , test_all_fit } ,
  { "items of weight 0" , test_zero_weights } ,
  { "negative data" , test_negative } ,
  { "no items" , test_no_items } ,
  { "fixed items" , test_fixed } ,
  { "objective sense" , test_sense } ,
  { "reoptimization" , test_reoptimization } ,
  { "abstract Modification" , test_abstract_modifications } ,
  { "abstract representation" , test_abstract_representation } ,
  { "netCDF round trip" , test_netcdf } ,
  { "R3 copy" , test_R3_copy } ,
  { "fractional weights" , test_fractional_weights } ,
  { "rounding of the bounds" , test_rounding } ,
  { "reuse of the previous solve" , test_reopt_outcome } ,
  { "branching" , test_branch } };

 for( auto & [ name , f ] : cases ) {
  const int before = failures;
  try {
   f();
   }
  catch( std::exception & e ) {
   check( false , name + ": exception " + e.what() );
   }
  std::cout << ( failures == before ? "OK  " : "KO  " ) << name << std::endl;
  }

 if( shared_block ) {
  shared_block->unregister_Solvers( true );
  delete shared_block;
  }

 if( failures ) {
  std::cout << failures << " checks failed" << std::endl;
  return( 1 );
  }

 std::cout << "All tests passed!!" << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
