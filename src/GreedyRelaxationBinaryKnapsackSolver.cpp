/*--------------------------------------------------------------------------*/
/*------------- File GreedyRelaxationBinaryKnapsackSolver.cpp --------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class
 * GreedyRelaxationBinaryKnapsackSolver, which implements the Solver concept
 * [see Solver.h] for solving the continuous relaxation of a Knapsack problem
 * as represented by a BinaryKnapsackBlock.
 * All the instance handling (raw mirror, Modification processing,
 * normalization) and the greedy fractional fill itself live in the base
 * BinaryKnapsackSolver; this class only drives them and exposes the
 * relaxation-specific results (true bounds, branching on the critical item).
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

int GreedyRelaxationBinaryKnapsackSolver::compute( bool changedvars ) {

 lock();                     // lock the mutex

 update_instance();          // process all the pending modifications

 if( ! f_norm_valid )        // profits / weights / sense changed:
  normalize_instance();      //  raw mirror -> normalized core
 else                        // only the fixings (possibly) changed:
  refresh_fixings();         //  cheap single-pass refresh

 obj = - Inf< double >();

 // the problem is empty iff the residual capacity is negative
 if( f_Cd < 0 ) { unlock(); return( kInfeasible ); }

 // solve the continuous knapsack
 obj = fractional_relaxation( f_fi );

 // the critical item, for branch(); with no critical item any index will do
 f_ci = f_fi.orig >= 0 ? Index( f_fi.orig ) : ( f_N ? f_N - 1 : 0 );

 unlock();                   // unlock the mutex

 return( kOK );

 }  // end( GreedyRelaxationBinaryKnapsackSolver::compute() )

/*--------------------------------------------------------------------------*/

std::vector< Change * > GreedyRelaxationBinaryKnapsackSolver::branch() {

 // reduced-cost (Martello-Toth) pegging, when an incumbent is available: a
 // free integer item whose flip cannot beat it is stuck at its greedy value
 // throughout this subtree [see the branch() doc]; the pegged items are split
 // by the value they are stuck at
 Block::Subset peg0 , peg1;
 // the incumbent and the local-fixing permission are read lock-free from
 // the cached atomic cells [see set_global_information()]: a not-finite
 // incumbent means none has been found yet
 const double cutoff = f_incumbent ?
  f_incumbent->load( std::memory_order_relaxed ) :
  std::numeric_limits< double >::quiet_NaN();
 if( std::isfinite( cutoff ) &&
     ( ( ! f_lfa ) || f_lfa->load( std::memory_order_relaxed ) ) &&
     ( f_fi.orig >= 0 ) ) {
  // the incumbent in the normalized (maximisation) sense; the conservative
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
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool GreedyRelaxationBinaryKnapsackSolver::has_true_var_solution( void ) {
 // the rounded greedy solution (see rounded_x()) always exists
 return( true );
 }

/*--------------------------------------------------------------------------*/

void GreedyRelaxationBinaryKnapsackSolver::get_var_solution(
                                             Configuration * solc ) {

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 BKB->set_x( f_x.begin() );  // write the solution in BinaryKnapsackBlock

 }  // end( GreedyRelaxationBinaryKnapsackSolver::get_var_solution )

/*--------------------------------------------------------------------------*/
/*----------- End File GreedyRelaxationBinaryKnapsackSolver.cpp ------------*/
/*--------------------------------------------------------------------------*/
