/*--------------------------------------------------------------------------*/
/*----- File IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp -----------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class
 * IncrementalGreedyRelaxationBinaryKnapsackSolver, the incremental variant of
 * GreedyRelaxationBinaryKnapsackSolver. It reuses, from the base
 * BinaryKnapsackSolver, the raw mirror of the instance and its normalized
 * core, and from GreedyRelaxationBinaryKnapsackSolver the branching, the
 * separation, the bounds and the Solution; it only replaces the way the
 * continuous relaxation is solved and the way an (un)fixing is absorbed, so
 * as to carry the greedy solution along the enumeration tree instead of
 * recomputing it at every node.
 *
 * \author Antonio Frangioni \n
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
 * Copyright &copy by Antonio Frangioni, Filippo Magi, Donato Meoli
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

// register IncrementalGreedyRelaxationBinaryKnapsackSolver to the factory

SMSpp_insert_in_factory_cpp_1(
                       IncrementalGreedyRelaxationBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*--- METHODS OF IncrementalGreedyRelaxationBinaryKnapsackSolver -----------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyRelaxationBinaryKnapsackSolver::rebuild_positions( void )
{
 f_pos.resize( f_N );
 for( Index j = 0 ; j < f_N ; ++j )
  f_pos[ v_ord[ j ] ] = j;
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int IncrementalGreedyRelaxationBinaryKnapsackSolver::compute(
                                                          bool changedvars )
{
 lock();                     // lock the mutex

 // a change in the instance (profits, weights, sense, capacity, fixings via
 // Modification, or a full reload) brings the carried state out of sync with
 // the core, so the relaxation is rebuilt from the top; along the tree, where
 // the (un)fixings come through apply() and keep the state current, the
 // greedy fill is instead resumed from where it stopped
 if( update_instance() || ( ! f_inited ) ) {
  normalize_instance();

  if( ! f_ord_valid )        // (re)build the cached efficiency order
   build_efficiency_order();
  rebuild_positions();

  f_ssi = 0;  f_run_p = 0;  f_run_rem = f_Cd;
  f_inited = true;
  }
 else if( f_ssi == 0 ) {
  // the fill restarts from the top after a move that may have unbalanced the
  // running state (a best-first jump, the composed undo of a separation):
  // the base profit, residual capacity and relaxation membership are taken
  // afresh from the fixings, the only quantities the apply() keeps exact
  refresh_fixings();
  f_run_p = 0;  f_run_rem = f_Cd;
  }

 obj = - Inf< double >();

 // the problem is empty iff the residual capacity is negative
 if( f_Cd < 0 ) { unlock(); return( kInfeasible ); }

 // resume the greedy fractional fill from the saved rank, carrying over the
 // profit and the residual capacity of the items already settled before it
 double rem = f_run_rem;
 double p = f_run_p;
 f_fi = FracInfo{ -1 , 1 , false , false , 0 , 0 };

 for( Index j = f_ssi ; j < f_N ; ++j ) {
  const int i = v_ord[ j ];
  if( ! n_in[ i ] )          // settled (fixed or folded): skip
   continue;

  const double w = n_w[ i ];

  if( w > rem ) {            // the critical item: take it fractionally
   const double y = rem / w;
   f_fi = FracInfo{ i , n_comp[ i ] ? 1 - y : y , ! v_I[ i ] ,
                    ( bool ) n_comp[ i ] , n_p[ i ] , y };
   f_x[ i ] = n_comp[ i ] ? 1 - y : y;

   // the resume point is the critical item, with the running quantities
   // taken just before it, so a child can pick them up from branch()
   f_ssi = j;  f_run_rem = rem;  f_run_p = p;

   obj = f_base + p + y * n_p[ i ];
   f_ci = i;

   // the items past the critical one stay out of the relaxation
   for( Index k = j + 1 ; k < f_N ; ++k ) {
    const int ii = v_ord[ k ];
    if( n_in[ ii ] )
     f_x[ ii ] = n_comp[ ii ] ? 1 : 0;
    }

   f_can_resume = true;
   unlock();
   return( kOK );
   }

  rem -= w;  p += n_p[ i ];  f_x[ i ] = n_comp[ i ] ? 0 : 1;
  }

 // every free item fits: the relaxation is integral
 f_ssi = f_N;  f_run_rem = rem;  f_run_p = p;
 obj = f_base + p;
 f_ci = f_N ? f_N - 1 : 0;

 f_can_resume = true;
 unlock();
 return( kOK );

 }  // end( IncrementalGreedyRelaxationBinaryKnapsackSolver::compute )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

std::vector< Change * >
                  IncrementalGreedyRelaxationBinaryKnapsackSolver::branch()
{
 // branch on the critical item, carrying in each child Change the residual
 // capacity and profit accumulated before it, so that apply() can resume the
 // child relaxation from the parent state instead of from the top
 std::vector< Change * > branches( 2 );

 std::vector< double > zero = { 0 , f_run_rem , f_run_p };
 std::vector< double > one = { 1 , f_run_rem , f_run_p };

 branches[ 0 ] = new BinaryKnapsackBlockRngdChange(
                                   BinaryKnapsackBlockChange::eFixX ,
                                   std::move( zero ) ,
                                   std::make_pair( f_ci , f_ci + 1 ) );

 branches[ 1 ] = new BinaryKnapsackBlockRngdChange(
                                   BinaryKnapsackBlockChange::eFixX ,
                                   std::move( one ) ,
                                   std::make_pair( f_ci , f_ci + 1 ) );

 return( branches );
 }

/*--------------------------------------------------------------------------*/

Change * IncrementalGreedyRelaxationBinaryKnapsackSolver::apply(
                                                Change * chg , bool doUndo )
{
 // a resume is sound only for the first apply() right after a compute(), a
 // single descent step from the node just solved [see f_can_resume]
 const bool resume_ok = f_can_resume;
 f_can_resume = false;

 // a GroupChange (e.g., the composed undo of branching plus separation) is
 // applied by decomposing it; the sub-undos are composed back in reverse
 if( auto grp = dynamic_cast< GroupChange * >( chg ) ) {
  GroupChange * undo = doUndo ? new GroupChange() : nullptr;
  for( auto sub : grp->sub_Changes() )
   if( auto u = apply( sub , doUndo ) ; u && undo )
    undo->add_front( u );
  return( undo );
  }

 auto CHG = dynamic_cast< BinaryKnapsackBlockChange * >( chg );
 if( ! CHG )
  throw( std::invalid_argument(
   "IncrementalGreedyRelaxationBinaryKnapsackSolver::apply: the Change must "
   "be a BinaryKnapsackBlockChange" ) );

 if( ( CHG->type() != BinaryKnapsackBlockChange::eFixX ) &&
     ( CHG->type() != BinaryKnapsackBlockChange::eUnfixX ) )
  return( CHG->apply( f_Block , doUndo ) );  // not an (un)fix: to the Block

 const bool fixing = ( CHG->type() == BinaryKnapsackBlockChange::eFixX );
 const Index n = CHG->num_items();
 const auto & data = CHG->data();

 // branch() Changes carry, past the n fixed values, the residual capacity
 // and profit of the parent at the critical item
 const bool branched = ( data.size() == n + 2 );

 // the undo under the LIFO (climb-the-tree) discipline: (un)fixing the same
 // items, carrying the values needed to revert their fold into the core
 Change * undo = nullptr;
 if( doUndo ) {
  std::vector< double > old_data;
  Block::Subset nms;
  old_data.reserve( n );  nms.reserve( n );
  for( Index k = 0 ; k < n ; ++k ) {
   const Index i = CHG->item( k );
   nms.push_back( i );
   old_data.push_back( fixing ? data[ k ] : ( v_fxd[ i ] == 2 ? 1 : 0 ) );
   }
  undo = new BinaryKnapsackBlockSbstChange(
                     fixing ? BinaryKnapsackBlockChange::eUnfixX
                            : BinaryKnapsackBlockChange::eFixX ,
                     std::move( old_data ) , std::move( nms ) );
  }

 // the normalized core is not in place yet (a parallel-exploration worker
 // replays the path to its subtree before ever solving it): just keep the
 // fixings mirror current, the next compute() builds everything from it
 if( ! f_inited ) {
  for( Index k = 0 ; k < n ; ++k )
   v_fxd[ CHG->item( k ) ] = fixing ? ( data[ k ] > 0.5 ? 2 : 1 ) : 0;
  f_ssi = 0;
  return( undo );
  }

 // fold each item into / out of the core in O(1): a core value of 1 (the
 // item is taken) moves its profit and weight into the base and the residual
 // capacity, a core value of 0 only takes it out of the relaxation
 for( Index k = 0 ; k < n ; ++k ) {
  const Index i = CHG->item( k );
  const double val = data[ k ];
  const bool y = n_comp[ i ] ? ( val < 0.5 ) : ( val > 0.5 );
  if( fixing ) {
   if( y ) { f_base += n_p[ i ];  f_Cd -= n_w[ i ]; }
   n_in[ i ] = 0;  f_x[ i ] = val;  v_fxd[ i ] = ( val > 0.5 ) ? 2 : 1;
   }
  else {
   if( y ) { f_base -= n_p[ i ];  f_Cd += n_w[ i ]; }
   n_in[ i ] = 1;  v_fxd[ i ] = 0;
   }
  }

 // a branch that drops its item from the core (core value 0) leaves the
 // residual capacity untouched, so the greedy fill resumes from the critical
 // item with the carried-over running quantities; any other (un)fixing may
 // move the critical item, so the fill restarts from the top (still cheap:
 // the core was kept up to date above, no re-normalization is needed)
 if( resume_ok && fixing && branched && ( n == 1 ) &&
     ! ( n_comp[ CHG->item( 0 ) ] ? ( data[ 0 ] < 0.5 )
                                  : ( data[ 0 ] > 0.5 ) ) ) {
  f_ssi = f_pos[ CHG->item( 0 ) ];
  f_run_rem = data[ n ];
  f_run_p = data[ n + 1 ];
  }
 else {
  f_ssi = 0;  f_run_p = 0;  f_run_rem = f_Cd;
  }

 return( undo );

 }  // end( IncrementalGreedyRelaxationBinaryKnapsackSolver::apply )

/*--------------------------------------------------------------------------*/
/*--- End File IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp ---------*/
/*--------------------------------------------------------------------------*/
