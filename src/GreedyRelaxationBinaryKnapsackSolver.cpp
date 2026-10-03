/*--------------------------------------------------------------------------*/
/*------------- File GreedyRelaxationBinaryKnapsackSolver.cpp --------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class
 * GreedyRelaxationBinaryKnapsackSolver, which implements the Solver concept
 * [see Solver.h] for solving the continuous relaxation of a Knapsack problem
 * as represented by a BinaryKnapsackBlock. The handling of the instance (raw
 * mirror, Modification processing, normalization) and the greedy fill of
 * the solves from scratch are those of the base BinaryKnapsackSolver; the
 * greedy fill kept across the Changes (intIncremental == 1) is implemented
 * here.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Federica Di Pasquale, Filippo Magi,
 *                   Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <numeric>

#include "GreedyRelaxationBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register GreedyRelaxationBinaryKnapsackSolver to the factory

SMSpp_insert_in_factory_cpp_1( GreedyRelaxationBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*------------ METHODS OF GreedyRelaxationBinaryKnapsackSolver -------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int GreedyRelaxationBinaryKnapsackSolver::compute( bool changedvars )
{
 lock();                     // lock the mutex

 if( f_incremental ) {
  const int status = inc_compute();
  unlock();
  return( status );
  }

 update_instance();          // process all the pending modifications

 if( ! f_norm_valid )        // profits / weights / sense changed:
  normalize_instance();      //  raw mirror -> normalized core
 else                        // only the fixings (possibly) changed:
  refresh_fixings();         //  cheap single-pass refresh

 obj = - Inf< double >();

 // the problem is empty iff the residual capacity is negative
 if( f_Cd < 0 ) {
  unlock();
  return( kInfeasible );
  }

 // solve the continuous knapsack
 obj = fractional_relaxation( f_fi );

 // the critical item, for branch(); with no critical item any index will do
 f_ci = f_fi.orig >= 0 ? Index( f_fi.orig ) : ( f_N ? f_N - 1 : 0 );

 unlock();                   // unlock the mutex

 return( kOK );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationBinaryKnapsackSolver::get_var_solution(
                                                      Configuration * solc )
{
 sync_x();
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 BKB->set_x( f_x.begin() );  // write the solution in BinaryKnapsackBlock

 }  // end( GreedyRelaxationBinaryKnapsackSolver::get_var_solution )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

std::vector< Change * > GreedyRelaxationBinaryKnapsackSolver::branch( void )
{
 if( f_incremental )
  return( inc_branch() );

 // reduced-cost (Martello-Toth) pegging, when an incumbent is available: a
 // free integer item whose flip cannot beat it is stuck at its greedy value
 // throughout this subtree [see the branch() doc]; the pegged items are split
 // by the value they are stuck at
 Block::Subset peg0 , peg1;
 // the incumbent is read from the "DoubleProperties" universe of the
 // GlobalInformation, if any [see GlobalInformation.h]
 bool hasIncumbent = false;
 double cutoff = 0;
 if( f_global_information )
  if( auto universe = f_global_information->get_from_Universe< double >(
							"DoubleProperties" ) )
   hasIncumbent = universe->read( "incumbent" , cutoff );
 if( hasIncumbent && ( f_fi.orig >= 0 ) ) {
  // the incumbent in the normalized (maximization) sense; the conservative
  // guard wards off pegging the optimum away by a floating-point whisker
  const double zmax = ( f_sense ? cutoff : - cutoff ) -
                      1e-9 * ( 1 + std::abs( cutoff ) );
  const Index ci = Index( f_fi.orig );
  const double rho = n_p[ ci ] / n_w[ ci ];   // critical efficiency
  for( Index j = 0 ; j < f_N ; ++j ) {
   if( ( ! n_in[ j ] ) || ( j == ci ) || ( ! v_I[ j ] ) )
    continue;                  // not free, critical, or not integer
   const double xj = f_x[ j ];                   // 0 or 1 (original space)
   const char yj = n_comp[ j ] ? char( 1 - xj ) : char( xj );
   const double margin = yj ? n_p[ j ] - rho * n_w[ j ]
                            : rho * n_w[ j ] - n_p[ j ];
   if( obj - margin <= zmax )
    ( xj ? peg1 : peg0 ).push_back( j );
   }
  }

 // the two children fix the critical item to 0 and to 1; each also carries
 // the (subtree-valid) pegged fixings, folded into a GroupChange
 std::vector< Change * > branches( 2 );
 for( int b = 0 ; b < 2 ; ++b ) {
  Change * crit = new BinaryKnapsackBlockRngdChange(
                       BinaryKnapsackBlockChange::eFixX ,
                       std::vector< double >{ double( b ) } ,
                       std::make_pair( f_ci , f_ci + 1 ) );
  if( peg0.empty() && peg1.empty() ) {
   branches[ b ] = crit;
   continue;
   }
  auto grp = new GroupChange();
  grp->add( crit );
  if( ! peg0.empty() )
   grp->add( new BinaryKnapsackBlockSbstChange(
                  BinaryKnapsackBlockChange::eFixX ,
                  std::vector< double >( peg0.size() , 0 ) ,
                  Block::Subset( peg0 ) ) );
  if( ! peg1.empty() )
   grp->add( new BinaryKnapsackBlockSbstChange(
                  BinaryKnapsackBlockChange::eFixX ,
                  std::vector< double >( peg1.size() , 1 ) ,
                  Block::Subset( peg1 ) ) );
  branches[ b ] = grp;
  }

 return( branches );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::branch )

/*--------------------------------------------------------------------------*/

Change * GreedyRelaxationBinaryKnapsackSolver::apply( Change * chg ,
						      bool doUndo )
{
 // a GroupChange (e.g., the composed undo of branching plus separation, see
 // BranchAndXSolver) is applied by decomposing it; the sub-undos are
 // composed back in reverse order
 if( auto grp = dynamic_cast< GroupChange * >( chg ) ) {
  GroupChange * undo = doUndo ? new GroupChange() : nullptr;
  for( auto sub : grp->sub_Changes() )
   if( auto u = apply( sub , doUndo ) ; u && undo )
    undo->add_front( u );
  return( undo );
  }

 auto CHG = dynamic_cast< BinaryKnapsackBlockChange * >( chg );
 if( ! CHG )
  throw( std::invalid_argument( "GreedyRelaxationBinaryKnapsackSolver::"
				"apply: the Change must be a "
				"BinaryKnapsackBlockChange" ) );

 if( ( CHG->type() != BinaryKnapsackBlockChange::eFixX ) &&
     ( CHG->type() != BinaryKnapsackBlockChange::eUnfixX ) )
  return( CHG->apply( f_Block , doUndo ) );  // not an (un)fix: to the Block

 if( f_incremental )
  return( inc_apply( CHG , doUndo ) );

 return( apply_to_mirror( CHG , doUndo ) );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::apply )

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationBinaryKnapsackSolver::sync_x( void )
{
 if( ( ! f_incremental ) || f_x_valid )
  return;

 // the items before the critical one in the efficiency order are taken,
 // those after it are not, the critical one is taken by f_ci_val; the
 // fixed items keep the value v_ix has for them
 double xval = 1;
 for( Index i : v_sorted ) {
  if( v_skip[ i ] ) {
   if( i == f_ci )
    xval = 0;
   continue;
   }
  if( i == f_ci ) {
   v_ix[ i ] = v_icomp[ i ] ? 1 - f_ci_val : f_ci_val;
   xval = 0;
   }
  else
   v_ix[ i ] = v_icomp[ i ] ? 1 - xval : xval;
  }

 f_x = v_ix;
 f_x_valid = true;

 }  // end( GreedyRelaxationBinaryKnapsackSolver::sync_x )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationBinaryKnapsackSolver::inc_initialize( void )
{
 // the data in the maximization sense; a free item with negative weight
 // and profit is complemented, a fixed one is not
 v_ip.resize( f_N );
 v_iw.resize( f_N );
 v_icomp.assign( f_N , false );
 for( Index i = 0 ; i < f_N ; ++i ) {
  v_ip[ i ] = f_sense ? v_P[ i ] : - v_P[ i ];
  v_iw[ i ] = v_W[ i ];
  if( ( ! v_fxd[ i ] ) && ( v_iw[ i ] < 0 ) && ( v_ip[ i ] < 0 ) ) {
   v_ip[ i ] = - v_ip[ i ];
   v_iw[ i ] = - v_iw[ i ];
   v_icomp[ i ] = true;
   }
  }

 // sort the items by nonincreasing profit / weight ratio; an item of weight
 // 0 has ratio + or - infinity (0 for a null profit), so that the
 // comparison is a strict weak ordering
 const auto ratio = [ this ]( Index i ) {
  if( v_iw[ i ] != 0 )
   return( v_ip[ i ] / v_iw[ i ] );
  return( v_ip[ i ] > 0 ? Inf< double >()
	                : ( v_ip[ i ] < 0 ? - Inf< double >() : 0.0 ) );
  };
 v_sorted.resize( f_N );
 v_sorted_pos.resize( f_N );
 std::iota( v_sorted.begin() , v_sorted.end() , 0 );
 std::sort( v_sorted.begin() , v_sorted.end() ,
	    [ & ]( Index a , Index b ) { return( ratio( a ) > ratio( b ) ); } );
 for( Index i = 0 ; i < f_N ; ++i )
  v_sorted_pos[ v_sorted[ i ] ] = i;

 f_ci = f_N ? v_sorted.front() : 0;
 f_ssi = 0;
 f_iP = 0;
 f_free_cap = f_Cap;
 v_ix.assign( f_N , 0 );
 v_skip.assign( f_N , false );

 for( Index i = 0 ; i < f_N ; ++i ) {
  if( v_fxd[ i ] == 1 ) {         // fixed to 0
   v_ix[ i ] = 0;
   v_skip[ i ] = true;
   continue;
   }

  if( v_fxd[ i ] == 2 ) {         // fixed to 1
   v_ix[ i ] = 1;
   v_skip[ i ] = true;
   f_iP += v_ip[ i ];
   f_free_cap -= v_iw[ i ];
   continue;
   }

  if( ( v_iw[ i ] >= 0 ) && ( v_ip[ i ] <= 0 ) ) {  // never convenient
   v_ix[ i ] = 0;
   v_skip[ i ] = true;
   }
  else
   if( ( v_iw[ i ] <= 0 ) && ( v_ip[ i ] >= 0 ) ) {  // always convenient
    v_ix[ i ] = 1;
    v_skip[ i ] = true;
    f_free_cap -= v_iw[ i ];
    f_iP += v_ip[ i ];
    }
   else
    if( v_icomp[ i ] ) {  // complemented: tentatively taken
     f_free_cap += v_iw[ i ];
     f_iP -= v_ip[ i ];
     }
  }

 // an item heavier than the residual capacity is not dropped: the
 // continuous relaxation may take a fraction of it, as the critical item

 f_iC = f_free_cap;
 f_x_valid = false;

 }  // end( GreedyRelaxationBinaryKnapsackSolver::inc_initialize )

/*--------------------------------------------------------------------------*/

int GreedyRelaxationBinaryKnapsackSolver::inc_compute( void )
{
 // any Modification makes the greedy fill start again from the data
 if( update_instance() || ( ! f_inc_valid ) ) {
  obj = - Inf< double >();
  inc_initialize();
  f_inc_valid = true;
  }

 f_x_valid = false;

 // the fixed items alone exceed the capacity: the problem is empty
 if( f_free_cap < 0 )
  return( kInfeasible );

 // solve the continuous knapsack, resuming from f_ssi; f_iC and f_iP are
 // kept up to date by apply()
 f_ci_val = 1;  // value of the critical item
 f_to_end = true;
 for( Index j = f_ssi ; j < f_N ; ++j ) {
  const Index i = v_sorted[ j ];
  if( v_skip[ i ] )
   continue;
  if( f_iC - v_iw[ i ] < 0 ) {
   f_ci = i;
   f_to_end = false;
   f_ci_val = f_iC / v_iw[ i ];
   f_ssi = j;
   break;
   }
  f_iC -= v_iw[ i ];
  f_iP += v_ip[ i ];
  }

 if( f_to_end )
  f_ci = f_N ? v_sorted.back() : 0;

 obj = f_iP + ( ( f_ci_val == 0 ) || ( f_ci_val == 1 ) ? 0
	                                                : v_ip[ f_ci ] * f_ci_val );

 // the critical item, as fractional_relaxation() describes it
 if( f_to_end || ( ! f_N ) )
  f_fi = FracInfo{ -1 , 1 , false , false , 0 , 0 };
 else
  f_fi = FracInfo{ int( f_ci ) ,
		   v_icomp[ f_ci ] ? 1 - f_ci_val : f_ci_val ,
		   ! v_I[ f_ci ] , bool( v_icomp[ f_ci ] ) ,
		   v_ip[ f_ci ] , f_ci_val };

 return( kOK );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::inc_compute )

/*--------------------------------------------------------------------------*/

Change * GreedyRelaxationBinaryKnapsackSolver::inc_apply(
				  BinaryKnapsackBlockChange * chg , bool doUndo )
{
 // the fixings apply to the greedy fill of the current data
 if( update_instance() || ( ! f_inc_valid ) ) {
  inc_initialize();
  f_inc_valid = true;
  }

 f_x_valid = false;

 // the items of the Change, and its data (copied, so that the Change can
 // be applied again)
 const auto type = chg->type();
 auto sbst = dynamic_cast< BinaryKnapsackBlockSbstChange * >( chg );
 auto rngd = dynamic_cast< BinaryKnapsackBlockRngdChange * >( chg );
 if( ( ! sbst ) && ( ! rngd ) )
  throw( std::invalid_argument( "GreedyRelaxationBinaryKnapsackSolver::"
				"inc_apply: the Change is neither a Subset "
				"nor a Range one" ) );

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
 std::vector< double > parent = { double( f_ssi ) , f_iC , f_iP };
 std::vector< double > state;
 if( branching || restore )
  state.assign( data.begin() + items.size() , data.end() );
 data.resize( items.size() );

 if( branching ) {
  f_iC = state[ 0 ];
  f_iP = state[ 1 ];
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
   inc_unfix_item( items[ h ] , data[ h ] );
  else
   inc_fix_item( items[ h ] , data[ h ] );

 if( restore ) {  // the undo of a branching: the parent comes back
  f_ssi = Index( state[ 0 ] );
  f_iC = state[ 1 ];
  f_iP = state[ 2 ];
  }

 return( undo );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::inc_apply )

/*--------------------------------------------------------------------------*/

std::vector< Change * > GreedyRelaxationBinaryKnapsackSolver::inc_branch(
								       void )
{
 // the data of each child are the value the critical item is fixed to, the
 // residual capacity after the fixing and the profit of the child, the last
 // two being what apply() takes as the state of the child

 if( f_free_cap - v_iw[ f_ci ] < 0 ) {  // only the child at 0
  std::vector< Change * > branches( 1 );
  branches[ 0 ] = new BinaryKnapsackBlockRngdChange(
			   BinaryKnapsackBlockChange::eFixX ,
			   std::vector< double >{ 0.0 , f_iC , f_iP } ,
			   std::make_pair( f_ci , f_ci + 1 ) );
  return( branches );
  }

 std::vector< Change * > branches( 2 );
 branches[ 0 ] = new BinaryKnapsackBlockRngdChange(
			   BinaryKnapsackBlockChange::eFixX ,
			   std::vector< double >{ 1.0 , f_iC - v_iw[ f_ci ] ,
						  f_iP + v_ip[ f_ci ] } ,
			   std::make_pair( f_ci , f_ci + 1 ) );
 branches[ 1 ] = new BinaryKnapsackBlockRngdChange(
			   BinaryKnapsackBlockChange::eFixX ,
			   std::vector< double >{ 0.0 , f_iC , f_iP } ,
			   std::make_pair( f_ci , f_ci + 1 ) );
 return( branches );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::inc_branch )

/*--------------------------------------------------------------------------*/

void GreedyRelaxationBinaryKnapsackSolver::inc_unfix_item( Index i ,
							   double value )
{
 v_skip[ i ] = false;
 f_free_cap += value * v_iw[ i ];
 f_ssi = std::min( f_ssi , v_sorted_pos[ i ] );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::inc_unfix_item )

/*--------------------------------------------------------------------------*/

void GreedyRelaxationBinaryKnapsackSolver::inc_fix_item( Index i ,
							 double value )
{
 if( ( value != 0 ) && ( value != 1 ) )
  throw( std::invalid_argument( "GreedyRelaxationBinaryKnapsackSolver::"
				"inc_fix_item: the value of a fixing must "
				"be 0 or 1" ) );

 v_ix[ i ] = v_icomp[ i ] ? 1 - value : value;
 f_free_cap -= value * v_iw[ i ];
 v_skip[ i ] = true;
 f_ci = i;
 const Index pos = v_sorted_pos[ i ];

 // takes out of the greedy fill the free items from position a to b
 auto take_out = [ & ]( Index a , Index b ) {
  for( Index j = a ; j <= b ; ++j ) {
   const Index z = v_sorted[ j ];
   if( v_skip[ z ] )
    continue;
   f_iP -= v_ip[ z ];
   f_iC += v_iw[ z ];
   }
  };

 if( ( value == 1 ) && ( f_iC < 0 ) ) {
  // the fixing to 1 leaves too little capacity: go back to the first item
  // from which the capacity is enough, and resume the search from there
  for( int j = int( pos ) ; j >= 0 ; --j ) {
   const Index z = v_sorted[ j ];
   if( v_skip[ z ] )
    continue;
   f_iP -= v_ip[ z ];
   f_iC += v_iw[ z ];
   if( f_iC >= 0 ) {
    f_ssi = j;
    break;
    }
   }
  }
 else
  if( f_ssi > pos )  // the item is before: resume from it
   f_ssi = pos;
  else
   if( ( value == 1 ) || ( f_ssi < pos ) )
    take_out( f_ssi , pos );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::inc_fix_item )

/*--------------------------------------------------------------------------*/
/*----------- End File GreedyRelaxationBinaryKnapsackSolver.cpp ------------*/
/*--------------------------------------------------------------------------*/
