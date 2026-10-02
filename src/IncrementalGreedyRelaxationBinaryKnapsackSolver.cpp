/*--------------------------------------------------------------------------*/
/*----- File IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp -----------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* classes
 * IncrementalGreedyChangeBinaryKnapsackSolver and
 * IncrementalGreedyRelaxationBinaryKnapsackSolver, which implement the
 * ChangeSolver concept [see ChangeSolver.h] and the RelaxationSolver one
 * [see RelaxationSolver.h] for solving the continuous relaxation of a
 * Knapsack problem as represented by a BinaryKnapsackBlock.
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Filippo Magi, Federica Di Pasquale, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <numeric>

#include "IncrementalGreedyRelaxationBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register IncrementalGreedyChangeBinaryKnapsackSolver and
// IncrementalGreedyRelaxationBinaryKnapsackSolver to the factory

SMSpp_insert_in_factory_cpp_1( IncrementalGreedyChangeBinaryKnapsackSolver );
SMSpp_insert_in_factory_cpp_1(
			     IncrementalGreedyRelaxationBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*------- METHODS OF IncrementalGreedyChangeBinaryKnapsackSolver ----------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::set_Block( Block * block )
{
 if( block == f_Block )  // nothing to do
  return;

 Solver::set_Block( block );  // attach to the new Block

 load();  // load the Binary Knapsack instance

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::set_Block )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int IncrementalGreedyChangeBinaryKnapsackSolver::compute( bool changedvars )
{
 lock();

 process_outstanding_Modification();

 if( changedData ) {
  obj = - Inf< double >();
  initializeVariables();
  C = presolveCapacity;
  }
 changedData = false;

 // C, P and startingSearchIndex are kept up to date by apply()

 // the fixed items alone exceed the capacity: the problem is empty
 if( presolveCapacity < 0 ) {
  unlock();
  return( kInfeasible );
  }

 // solve the continuous knapsack, resuming from startingSearchIndex
 f_ciVal = 1;  // value of the critical item
 reachTheEnd = true;
 for( Index j = startingSearchIndex ; j < f_N ; ++j ) {
  const Index i = sortedVar[ j ];
  if( skip[ i ] )
   continue;
  if( C - v_W[ i ] < 0 ) {
   f_ci = i;
   reachTheEnd = false;
   f_ciVal = C / v_W[ i ];
   startingSearchIndex = j;
   break;
   }
  C -= v_W[ i ];
  P += v_P[ i ];
  }

 if( reachTheEnd )
  f_ci = sortedVar.back();

 obj = P + ( ( f_ciVal == 0 ) || ( f_ciVal == 1 ) ? 0
	                                            : v_P[ f_ci ] * f_ciVal );

 unlock();

 return( kOK );

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::compute )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::add_Modification(
							     sp_Mod & mod )
{
 if( f_no_Mod )
  return;

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;

 // if NBModification, reload BinaryKnapsack instance and clear modifications
 if( std::dynamic_pointer_cast< NBModification >( mod ) ) {
  load();
  v_mod.clear();
  }
 else
  v_mod.push_back( mod );

 f_mod_lock.clear( std::memory_order_release );  // release lock

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::add_Modification )

/*--------------------------------------------------------------------------*/

Change * IncrementalGreedyChangeBinaryKnapsackSolver::apply( Change * chg ,
							     bool doUndo )
{
 auto CHG = dynamic_cast< BinaryKnapsackBlockChange * >( chg );
 if( ! CHG )
  throw( std::invalid_argument( "IncrementalGreedyChangeBinaryKnapsack"
				"Solver::apply: the Change must be a "
				"BinaryKnapsackBlockChange" ) );

 const auto type = CHG->type();
 if( ( type != BinaryKnapsackBlockChange::eFixX ) &&
     ( type != BinaryKnapsackBlockChange::eUnfixX ) )
  return( CHG->apply( f_Block , doUndo ) );  // not an (un)fix: to the Block

 // the items of the Change, and its data (copied, so that the Change can
 // be applied again)
 auto sbst = dynamic_cast< BinaryKnapsackBlockSbstChange * >( CHG );
 auto rngd = dynamic_cast< BinaryKnapsackBlockRngdChange * >( CHG );
 if( ( ! sbst ) && ( ! rngd ) )
  throw( std::invalid_argument( "IncrementalGreedyChangeBinaryKnapsack"
				"Solver::apply: the Change is neither a "
				"Subset nor a Range one" ) );

 Block::Subset items;
 if( sbst )
  items = sbst->nms();
 else
  for( Index i = rngd->rng().first ; i < rngd->rng().second ; ++i )
   items.push_back( i );

 std::vector< double > data = sbst ? sbst->data() : rngd->data();

 // the Change of a branching carries in its last two positions the residual
 // capacity and the profit of the child [see branch()], and its undo the
 // state of the parent in its last three ones: the position from which
 // compute() resumes, the residual capacity and the profit
 const bool branching = ( items.size() + 2 == data.size() );
 const bool restore = ( items.size() + 3 == data.size() );
 std::vector< double > parent = { double( startingSearchIndex ) , C , P };
 std::vector< double > state;
 if( branching || restore )
  state.assign( data.begin() + items.size() , data.end() );
 data.resize( items.size() );

 if( branching ) {
  C = state[ 0 ];
  P = state[ 1 ];
  }

 // the undo Change, which carries the state of the parent if this is a
 // branching one, so that undoing it brings the parent back
 Change * undo = nullptr;
 if( doUndo ) {
  const auto utype = ( type == BinaryKnapsackBlockChange::eFixX )
                     ? BinaryKnapsackBlockChange::eUnfixX
                     : BinaryKnapsackBlockChange::eFixX;
  std::vector< double > udata( data );
  if( branching )
   udata.insert( udata.end() , parent.begin() , parent.end() );
  if( sbst )
   undo = new BinaryKnapsackBlockSbstChange( utype , std::move( udata ) ,
					     Block::Subset( items ) );
  else
   undo = new BinaryKnapsackBlockRngdChange( utype , std::move( udata ) ,
					     Block::Range( rngd->rng() ) );
  }

 for( Index h = 0 ; h < items.size() ; ++h )
  if( type == BinaryKnapsackBlockChange::eUnfixX )
   unfix_item( items[ h ] , data[ h ] );
  else
   fix_item( items[ h ] , data[ h ] );

 if( restore ) {  // the undo of a branching: the parent comes back
  startingSearchIndex = Index( state[ 0 ] );
  C = state[ 1 ];
  P = state[ 2 ];
  }

 return( undo );

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::apply )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::get_var_solution(
						      Configuration * solc )
{
 auto x = rounded_x();
 static_cast< BinaryKnapsackBlock * >( f_Block )->set_x( x.begin() );

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::get_var_solution )

/*--------------------------------------------------------------------------*/

Solution * IncrementalGreedyChangeBinaryKnapsackSolver::get_Solution(
						      Configuration * solc )
{
 return( new BinaryKnapsackSolution( rounded_x() ) );

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::get_Solution )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::unfix_item( Index i ,
							      double value )
{
 skip[ i ] = false;
 presolveCapacity += value * v_W[ i ];
 startingSearchIndex = std::min( startingSearchIndex , varToSorted[ i ] );

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::unfix_item )

/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::fix_item( Index i ,
							    double value )
{
 if( ( value != 0 ) && ( value != 1 ) )
  throw( std::invalid_argument( "IncrementalGreedyChangeBinaryKnapsack"
				"Solver::fix_item: the value of a fixing "
				"must be 0 or 1" ) );

 v_x[ i ] = complemented[ i ] ? 1 - value : value;
 presolveCapacity -= value * v_W[ i ];
 skip[ i ] = true;
 f_ci = i;
 const Index pos = varToSorted[ i ];

 // takes out of the greedy fill the free items from position a to b
 auto take_out = [ & ]( Index a , Index b ) {
  for( Index j = a ; j <= b ; ++j ) {
   const Index z = sortedVar[ j ];
   if( skip[ z ] )
    continue;
   P -= v_P[ z ];
   C += v_W[ z ];
   }
  };

 if( ( value == 1 ) && ( C < 0 ) ) {
  // the fixing to 1 leaves too little capacity: go back to the first item
  // from which the capacity is enough, and resume the search from there
  for( int j = int( pos ) ; j >= 0 ; --j ) {
   const Index z = sortedVar[ j ];
   if( skip[ z ] )
    continue;
   P -= v_P[ z ];
   C += v_W[ z ];
   if( C >= 0 ) {
    startingSearchIndex = j;
    break;
    }
   }
  }
 else
  if( startingSearchIndex > pos )  // the item is before: resume from it
   startingSearchIndex = pos;
  else
   if( ( value == 1 ) || ( startingSearchIndex < pos ) )
    take_out( startingSearchIndex , pos );

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::fix_item )

/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::initializeVariables( void )
{
 // sort the items by nonincreasing profit / weight ratio
 C = 0;
 P = 0;
 startingSearchIndex = 0;
 sortedVar.resize( f_N );
 varToSorted.resize( f_N );
 std::iota( sortedVar.begin() , sortedVar.end() , 0 );
 std::sort( sortedVar.begin() , sortedVar.end() ,
	    [ & ]( Index a , Index b ) {
	     return( v_P[ a ] / v_W[ a ] > v_P[ b ] / v_W[ b ] );
	     } );
 for( Index i = 0 ; i < f_N ; ++i )
  varToSorted[ sortedVar[ i ] ] = i;

 f_ci = f_N ? sortedVar.front() : 0;
 v_x.assign( f_N , 0 );
 presolveCapacity = f_C;
 skip.assign( f_N , false );

 for( Index i = 0 ; i < f_N ; ++i ) {
  if( v_fxd[ i ] == 1 ) {         // fixed to 0
   v_x[ i ] = 0;
   skip[ i ] = true;
   continue;
   }

  if( v_fxd[ i ] == 2 ) {         // fixed to 1
   v_x[ i ] = 1;
   skip[ i ] = true;
   P += v_P[ i ];
   presolveCapacity -= v_W[ i ];
   continue;
   }

  if( ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ) {  // never convenient
   v_x[ i ] = 0;
   skip[ i ] = true;
   }
  else
   if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ) {  // always convenient
    v_x[ i ] = 1;
    skip[ i ] = true;
    presolveCapacity -= v_W[ i ];
    P += v_P[ i ];
    }
   else
    if( ( v_W[ i ] < 0 ) && ( v_P[ i ] < 0 ) ) {   // to be complemented
     presolveCapacity -= v_W[ i ];
     P += v_P[ i ];
     v_W[ i ] = - v_W[ i ];
     v_P[ i ] = - v_P[ i ];
     complemented[ i ] = true;
     }
    else
     if( complemented[ i ] ) {  // already complemented, hence positive
      presolveCapacity += v_W[ i ];
      P -= v_P[ i ];
      }
  }

 // an item heavier than the residual capacity is not dropped: the
 // continuous relaxation may take a fraction of it, as the critical item

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::initializeVariables )

/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::load( void )
{
 if( ! f_Block ) {
  f_N = 0;
  v_P.clear();
  v_W.clear();
  v_I.clear();
  v_x.clear();
  v_fxd.clear();
  sortedVar.clear();
  varToSorted.clear();
  complemented.clear();
  skip.clear();
  presolveCapacity = 0;
  C = 0;
  P = 0;
  startingSearchIndex = 0;
  return;
  }

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( ! BKB )
  throw( std::invalid_argument( "IncrementalGreedyChangeBinaryKnapsack"
				"Solver::load: the Block must be a "
				"BinaryKnapsackBlock" ) );

 // (try to) lock the BinaryKnapsackBlock
 bool owned = BKB->is_owned_by( f_id );
 if( ( ! owned ) && ( ! BKB->read_lock() ) )
  throw( std::runtime_error( "IncrementalGreedyChangeBinaryKnapsackSolver::"
			     "load: unable to lock the Block" ) );

 // the scalar data: sense of the objective, number of items and capacity
 f_sense = ( BKB->get_objective_sense() == Objective::eMax );
 f_N = BKB->get_NItems();
 f_C = BKB->get_Capacity();

 v_P.resize( f_N );
 v_W.resize( f_N );
 v_I.resize( f_N );
 v_fxd.resize( f_N );
 complemented.resize( f_N );

 const auto & BP = BKB->get_Profits();
 const auto & BW = BKB->get_Weights();
 const auto & FXD = BKB->get_fxd();
 const auto & BI = BKB->get_Integrality();

 for( Index i = 0 ; i < f_N ; ++i ) {
  v_P[ i ] = f_sense ? BP[ i ] : - BP[ i ];
  v_W[ i ] = BW[ i ];
  v_fxd[ i ] = FXD[ i ];
  v_I[ i ] = BI[ i ];
  complemented[ i ] = false;
  if( ( v_W[ i ] < 0 ) && ( v_P[ i ] < 0 ) ) {
   v_W[ i ] = - v_W[ i ];
   v_P[ i ] = - v_P[ i ];
   complemented[ i ] = true;
   }
  }
 changedData = true;

 if( ! owned )
  BKB->read_unlock();

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::load )

/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::
                                     process_outstanding_Modification( void )
{
 // copy v_mod in a temporary list of modifications - - - - - - - - - - - - -

 Lst_sp_Mod v_mod_tmp;

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 // any change in the profits must be processed only AFTER checking the
 // changes of the sense of the objective; hence v_mod_tmp is scanned twice,
 // checking the Modification of the sense (and of the capacity) first

 auto mod = v_mod_tmp.begin();

 while( mod != v_mod_tmp.end() ) {
  changedData = true;

  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >(
							     mod->get() ) )
   switch( tmod->type() ) {
    case( BinaryKnapsackBlockMod::eChgCapacity ):
     f_C = BKB->get_Capacity();
     mod = v_mod_tmp.erase( mod );
     break;
    case( BinaryKnapsackBlockMod::eChgSense ):
     f_sense = ( BKB->get_objective_sense() == Objective::eMax );
     for( Index i = 0 ; i < f_N ; ++i )
      if( complemented[ i ] ) {
       complemented[ i ] = false;
       v_W[ i ] = - v_W[ i ];
       }
      else
       v_P[ i ] = - v_P[ i ];
     mod = v_mod_tmp.erase( mod );
     break;
    default:
     ++mod;
    }
  else
   mod = v_mod_tmp.erase( mod );  // it is not a physical Modification
  }

 // the changes of profits and weights, and the (un)fixings, one item at a
 // time; a complemented item whose profit (or weight) does not change sign
 // goes back to its plain form, which initializeVariables() then
 // complements again if needed
 auto chg_profit = [ & ]( Index i ) {
  const double oldP = v_P[ i ];
  v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );
  if( complemented[ i ] && ( oldP * v_P[ i ] >= 0 ) ) {
   complemented[ i ] = false;
   v_W[ i ] = - v_W[ i ];
   }
  else
   if( complemented[ i ] )
    v_P[ i ] = - v_P[ i ];
  };

 auto chg_weight = [ & ]( Index i ) {
  const double oldW = v_W[ i ];
  v_W[ i ] = BKB->get_Weight( i );
  if( complemented[ i ] && ( oldW * v_W[ i ] >= 0 ) ) {
   complemented[ i ] = false;
   v_P[ i ] = - v_P[ i ];
   }
  else
   if( complemented[ i ] )
    v_W[ i ] = - v_W[ i ];
  };

 auto fix = [ & ]( Index i ) {
  if( BKB->is_fixed( i ) && ( std::abs( BKB->get_x( i ) ) < 1e-6 ) )
   v_fxd[ i ] = 1;
  else
   if( BKB->is_fixed( i ) && ( std::abs( BKB->get_x( i ) - 1 ) < 1e-6 ) )
    v_fxd[ i ] = 2;
  if( complemented[ i ] ) {
   v_P[ i ] = - v_P[ i ];
   v_W[ i ] = - v_W[ i ];
   complemented[ i ] = false;
   }
  };

 auto unfix = [ & ]( Index i ) { v_fxd[ i ] = 0; };

 auto chg_integrality = [ & ]( Index i ) {
  v_I[ i ] = BKB->get_Integrality( i );
  };

 // applies the action of the type of the Modification to each of its items
 auto process = [ & ]( int type , auto && each ) {
  switch( type ) {
   case( BinaryKnapsackBlockMod::eChgProfit ): each( chg_profit ); break;
   case( BinaryKnapsackBlockMod::eFixX ):      each( fix );        break;
   case( BinaryKnapsackBlockMod::eUnfixX ):    each( unfix );      break;
   case( BinaryKnapsackBlockMod::eChgWeight ): each( chg_weight ); break;
   case( BinaryKnapsackBlockMod::eChgIntegrality ):
    each( chg_integrality );
    break;
   }
  };

 for( auto & m : v_mod_tmp ) {
  changedData = true;
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >(
							       m.get() ) )
   process( tmod->type() , [ & ]( auto && f ) {
    for( Index i = tmod->rng().first ; i < tmod->rng().second ; ++i )
     f( i );
    } );
  else
   if( const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >(
							       m.get() ) )
    process( tmod->type() , [ & ]( auto && f ) {
     for( auto i : tmod->nms() )
      f( i );
     } );
  }

 v_mod_tmp.clear();  // clear the temporary list of Modification

 }  // end( IncrementalGreedyChangeBinaryKnapsackSolver::process_outstanding_
    //      Modification )

/*--------------------------------------------------------------------------*/
/*------- METHODS OF IncrementalGreedyRelaxationBinaryKnapsackSolver -------*/
/*--------------------------------------------------------------------------*/

std::vector< Change * >
IncrementalGreedyRelaxationBinaryKnapsackSolver::branch( void )
{
 // the data of each child are the value the critical item is fixed to, the
 // residual capacity after the fixing and the profit of the child, the last
 // two being what apply() takes as the state of the child

 if( presolveCapacity - v_W[ f_ci ] < 0 ) {  // only the child at 0
  std::vector< Change * > branches( 1 );
  branches[ 0 ] = new BinaryKnapsackBlockRngdChange(
			   BinaryKnapsackBlockChange::eFixX ,
			   std::vector< double >{ 0.0 , C , P } ,
			   std::make_pair( f_ci , f_ci + 1 ) );
  return( branches );
  }

 std::vector< Change * > branches( 2 );
 branches[ 0 ] = new BinaryKnapsackBlockRngdChange(
			   BinaryKnapsackBlockChange::eFixX ,
			   std::vector< double >{ 1.0 , C - v_W[ f_ci ] ,
						  P + v_P[ f_ci ] } ,
			   std::make_pair( f_ci , f_ci + 1 ) );
 branches[ 1 ] = new BinaryKnapsackBlockRngdChange(
			   BinaryKnapsackBlockChange::eFixX ,
			   std::vector< double >{ 0.0 , C , P } ,
			   std::make_pair( f_ci , f_ci + 1 ) );
 return( branches );

 }  // end( IncrementalGreedyRelaxationBinaryKnapsackSolver::branch )

/*--------------------------------------------------------------------------*/
/*--- End File IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp ---------*/
/*--------------------------------------------------------------------------*/
