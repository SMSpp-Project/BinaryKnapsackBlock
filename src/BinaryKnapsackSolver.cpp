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
    break;

   default:

    // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - -
    if( const auto rmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >(
                                                                 mod.get() ) )
     for( Index i = rmod->rng().first ; i < rmod->rng().second ; ++i )
      switch( tmod->type() ) {
       case( BinaryKnapsackBlockMod::eChgProfit ):
        v_P[ i ] = BKB->get_Profit( i ); break;
       case( BinaryKnapsackBlockMod::eChgWeight ):
        v_W[ i ] = BKB->get_Weight( i ); break;
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
        v_P[ i ] = BKB->get_Profit( i ); break;
       case( BinaryKnapsackBlockMod::eChgWeight ):
        v_W[ i ] = BKB->get_Weight( i ); break;
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

 for( Index i = 0 ; i < f_N ; ++i ) {

  // normalise the profit to a maximisation problem
  const double p = f_sense ? v_P[ i ] : - v_P[ i ];
  const double w = v_W[ i ];

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
   if( v_I[ i ] )                                //  undoes it ( x = 1 - y )
    { c_w.push_back( -w ); c_p.push_back( -p );
      c_orig.push_back( i ); c_comp.push_back( 1 ); }
   else
    { cc_w.push_back( -w ); cc_p.push_back( -p );
      cc_orig.push_back( i ); cc_comp.push_back( 1 ); }
   continue;
   }

  // genuine free item with positive weight and profit
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
 // the greedy fill works on the union of the two cores, [ 0 , m ) being the
 // integer items and [ m , m + mc ) the continuous ones
 const int m = int( c_w.size() );
 const int mc = int( cc_w.size() );

 auto wof = [ & ]( int k ) { return k < m ? c_w[ k ] : cc_w[ k - m ]; };
 auto pof = [ & ]( int k ) { return k < m ? c_p[ k ] : cc_p[ k - m ]; };

 std::vector< int > idx( m + mc );
 std::iota( idx.begin() , idx.end() , 0 );
 std::sort( idx.begin() , idx.end() , [ & ]( int a , int b ) {
  return( pof( a ) * wof( b ) > pof( b ) * wof( a ) );
  } );

 fi = FracInfo{ -1 , 1 , false , false , 0 , 0 };

 double obj = f_base;
 double rem = f_Cd;
 double yval = 1;              // 1 before the critical item, then 0

 for( int k : idx ) {
  const double w = wof( k );
  double y = yval;

  if( ( yval == 1 ) && ( w > rem ) ) {           // the critical item
   y = rem / w;
   const bool cont = ( k >= m );
   const Index orig = cont ? cc_orig[ k - m ] : c_orig[ k ];
   const char comp = cont ? cc_comp[ k - m ] : c_comp[ k ];
   fi = FracInfo{ int( orig ) , comp ? 1 - y : y , cont , ( bool ) comp ,
                  pof( k ) , y };
   yval = 0;
   }

  if( y == 1 )
   { rem -= w; obj += pof( k ); }
  else if( y > 0 )
   { rem = 0; obj += y * pof( k ); }

  if( k < m )
   f_x[ c_orig[ k ] ] = c_comp[ k ] ? 1 - y : y;
  else
   f_x[ cc_orig[ k - m ] ] = cc_comp[ k - m ] ? 1 - y : y;
  }

 return( obj );

 }  // end( BinaryKnapsackSolver::fractional_relaxation )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolver::load( void )
{
 f_changed = true;

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
