/*--------------------------------------------------------------------------*/
/*--------------------- File BinaryKnapsackSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *abstract* class BinaryKnapsackSolver: the raw
 * mirror of the BinaryKnapsackBlock instance with its incremental update
 * under Modification, the normalized (positive, maximisation) core, and the
 * continuous (Dantzig) relaxation of the latter.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Federica Di Pasquale, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <numeric>

#include "BinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolver::set_Block( Block * block )
{
 if( block == f_Block )       // nothing to do
  return;

 Solver::set_Block( block );  // attach to the new Block

 load();                      // load the raw mirror of the instance
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolver::add_Modification( sp_Mod & mod )
{
 if( f_no_Mod )
  return;

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 // if NBModification, reload the raw mirror and clear the queue
 if( std::dynamic_pointer_cast< NBModification >( mod ) ) {
  load();
  v_mod.clear();
  }
 else
  v_mod.push_back( mod );

 f_mod_lock.clear( std::memory_order_release );  // release lock

 }  // end( BinaryKnapsackSolver::add_Modification )

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

bool BinaryKnapsackSolver::update_instance( void )
{
 // drain the queue into a temporary list - - - - - - - - - - - - - - - - - -

 Lst_sp_Mod v_mod_tmp;

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 v_mod_tmp.splice( v_mod_tmp.end() , v_mod );

 bool changed = f_changed || ( ! v_mod_tmp.empty() );
 f_changed = false;

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // apply the changes to the raw mirror - - - - - - - - - - - - - - - - - - -
 // the mirror is raw (original signs and sense) and each change re-reads the
 // CURRENT Block value, so the scan is single-pass and order-insensitive

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( const auto & mod : v_mod_tmp ) {

  // BinaryKnapsackBlockMod - - - - - - - - - - - - - - - - - - - - - - - - -
  const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod.get() );
  if( ! tmod )
   continue;                  // not a physical Modification

  switch( tmod->type() ) {

   case( BinaryKnapsackBlockMod::eChgCapacity ):
    f_Cap = BKB->get_Capacity();
    break;

   case( BinaryKnapsackBlockMod::eChgSense ):
    f_sense = ( BKB->get_objective_sense() == Objective::eMax );
    f_ord_valid = false;
    break;

   default:

    // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - -
    if( const auto rmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >(
                                                                 mod.get() ) )
     for( Index i = rmod->rng().first ; i < rmod->rng().second ; ++i )
      switch( tmod->type() ) {
       case( BinaryKnapsackBlockMod::eChgProfit ):
        v_P[ i ] = BKB->get_Profit( i ); f_ord_valid = false; break;
       case( BinaryKnapsackBlockMod::eChgWeight ):
        v_W[ i ] = BKB->get_Weight( i ); f_ord_valid = false; break;
       case( BinaryKnapsackBlockMod::eChgIntegrality ):
        v_I[ i ] = BKB->get_Integrality( i ); break;
       case( BinaryKnapsackBlockMod::eFixX ):
       case( BinaryKnapsackBlockMod::eUnfixX ):
        v_fxd[ i ] = BKB->get_fxd()[ i ]; break;
       default: break;
       }

    // BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - -
    if( const auto smod = dynamic_cast< BinaryKnapsackBlockSbstMod * >(
                                                                 mod.get() ) )
     for( auto i : smod->nms() )
      switch( tmod->type() ) {
       case( BinaryKnapsackBlockMod::eChgProfit ):
        v_P[ i ] = BKB->get_Profit( i ); f_ord_valid = false; break;
       case( BinaryKnapsackBlockMod::eChgWeight ):
        v_W[ i ] = BKB->get_Weight( i ); f_ord_valid = false; break;
       case( BinaryKnapsackBlockMod::eChgIntegrality ):
        v_I[ i ] = BKB->get_Integrality( i ); break;
       case( BinaryKnapsackBlockMod::eFixX ):
       case( BinaryKnapsackBlockMod::eUnfixX ):
        v_fxd[ i ] = BKB->get_fxd()[ i ]; break;
       default: break;
       }
   }
  }  // end( for( mod ) )

 return( changed );

 }  // end( BinaryKnapsackSolver::update_instance )

/*--------------------------------------------------------------------------*/

bool BinaryKnapsackSolver::pending_changes( void )
{
 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 const bool pending = f_changed || ( ! v_mod.empty() );

 f_mod_lock.clear( std::memory_order_release );  // release lock

 return( pending );
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolver::normalize_instance( void )
{
 c_w.clear();  c_p.clear();  c_orig.clear();  c_comp.clear();
 cc_w.clear(); cc_p.clear(); cc_orig.clear(); cc_comp.clear();
 c_w.reserve( f_N );  c_p.reserve( f_N );
 c_orig.reserve( f_N );  c_comp.reserve( f_N );

 f_x.assign( f_N , 0 );
 f_base = 0;
 f_Cd = f_Cap;

 n_w.resize( f_N );
 n_p.resize( f_N );
 n_comp.resize( f_N );
 n_in.assign( f_N , 0 );

 for( Index i = 0 ; i < f_N ; ++i ) {

  // normalise the profit to a maximisation problem
  const double p = f_sense ? v_P[ i ] : - v_P[ i ];
  const double w = v_W[ i ];

  // per-item normalized data for the cached efficiency order; for items
  // not entering the relaxation (n_in == 0) the values are immaterial
  const bool cmp = ( w < 0 ) && ( p < 0 );
  n_w[ i ] = cmp ? -w : w;
  n_p[ i ] = cmp ? -p : p;
  n_comp[ i ] = cmp;

  // fixed variables: honour the fixing, fold into base profit / capacity
  if( v_fxd[ i ] == 1 ) {
   f_x[ i ] = 0;
   continue;
   }
  if( v_fxd[ i ] == 2 ) {
   f_x[ i ] = 1; f_base += p; f_Cd -= w;
   continue;
   }

  // sign-based pre-fixing
  if( ( w >= 0 ) && ( p <= 0 ) )                 // never profitable: drop
   { f_x[ i ] = 0; continue; }
  if( ( w <= 0 ) && ( p >= 0 ) )                 // always taken: enlarges C
   { f_x[ i ] = 1; f_base += p; f_Cd -= w; continue; }

  if( ( w < 0 ) && ( p < 0 ) ) {                 // complement: tentatively
   f_base += p; f_Cd -= w;                       //  take it, the core toggle
   n_in[ i ] = 1;                                //  undoes it ( x = 1 - y )
   if( v_I[ i ] )
    { c_w.push_back( -w ); c_p.push_back( -p );
      c_orig.push_back( i ); c_comp.push_back( 1 ); }
   else
    { cc_w.push_back( -w ); cc_p.push_back( -p );
      cc_orig.push_back( i ); cc_comp.push_back( 1 ); }
   continue;
   }

  // genuine free item with positive weight and profit
  n_in[ i ] = 1;
  if( v_I[ i ] )
   { c_w.push_back( w ); c_p.push_back( p );
     c_orig.push_back( i ); c_comp.push_back( 0 ); }
  else
   { cc_w.push_back( w ); cc_p.push_back( p );
     cc_orig.push_back( i ); cc_comp.push_back( 0 ); }
  }
 }  // end( BinaryKnapsackSolver::normalize_instance )

/*--------------------------------------------------------------------------*/

double BinaryKnapsackSolver::fractional_relaxation( FracInfo & fi )
{
 // the relaxation makes no difference between integer and continuous items:
 // the greedy fill works on the items entering the relaxation (see n_in) in
 // cached efficiency order, which only has to be (re)computed when profits,
 // weights or the objective sense change (see update_instance() / load()),
 // NOT when the fixings do: re-solves under different fixings are O( n )
 if( ! f_ord_valid ) {
  v_ord.resize( f_N );
  std::iota( v_ord.begin() , v_ord.end() , 0 );
  std::sort( v_ord.begin() , v_ord.end() , [ & ]( int a , int b ) {
   return( n_p[ a ] * n_w[ b ] > n_p[ b ] * n_w[ a ] );
   } );
  f_ord_valid = true;
  }

 fi = FracInfo{ -1 , 1 , false , false , 0 , 0 };

 double obj = f_base;
 double rem = f_Cd;
 double yval = 1;              // 1 before the critical item, then 0

 for( int i : v_ord ) {
  if( ! n_in[ i ] )            // not in the relaxation (fixed or trivial)
   continue;

  const double w = n_w[ i ];
  double y = yval;

  if( ( yval == 1 ) && ( w > rem ) ) {           // the critical item
   y = rem / w;
   fi = FracInfo{ i , n_comp[ i ] ? 1 - y : y , ! v_I[ i ] ,
                  ( bool ) n_comp[ i ] , n_p[ i ] , y };
   yval = 0;
   }

  if( y == 1 )
   { rem -= w; obj += n_p[ i ]; }
  else if( y > 0 )
   { rem = 0; obj += y * n_p[ i ]; }

  f_x[ i ] = n_comp[ i ] ? 1 - y : y;
  }

 return( obj );

 }  // end( BinaryKnapsackSolver::fractional_relaxation )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolver::load( void )
{
 f_changed = true;
 f_ord_valid = false;

 if( ! f_Block ) {
  f_N = 0;
  f_Cap = 0;
  v_P.clear();
  v_W.clear();
  v_I.clear();
  v_fxd.clear();
  return;
  }

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( ! BKB )
  throw( std::invalid_argument( "BinaryKnapsackSolver::load: Block must be "
                                "a BinaryKnapsackBlock" ) );

 // (try to) lock the BinaryKnapsackBlock
 bool owned = BKB->is_owned_by( f_id );
 if( ( ! owned ) && ( ! BKB->read_lock() ) )
  throw( std::runtime_error( "BinaryKnapsackSolver::load: unable to lock "
                             "the Block" ) );

 f_sense = ( BKB->get_objective_sense() == Objective::eMax );
 f_N = BKB->get_NItems();
 f_Cap = BKB->get_Capacity();

 const auto & P = BKB->get_Profits();
 const auto & W = BKB->get_Weights();
 const auto & I = BKB->get_Integrality();
 const auto & FXD = BKB->get_fxd();

 v_P.assign( P.begin() , P.end() );
 v_W.assign( W.begin() , W.end() );
 v_I.assign( I.begin() , I.end() );
 v_fxd.assign( FXD.begin() , FXD.end() );

 if( ! owned )
  BKB->read_unlock();

 }  // end( BinaryKnapsackSolver::load )

/*--------------------------------------------------------------------------*/
/*------------------ End File BinaryKnapsackSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
