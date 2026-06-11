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

 // branch on the critical item
 std::vector< Change * > branches( 2 );

 std::vector< double > zero = { 0 };
 std::vector< double > one = { 1 };

 branches[ 0 ] = new BinaryKnapsackBlockRngdChange(
                                   BinaryKnapsackBlockChange::eFixX ,
                                   std::move( zero ) ,
                                   std::make_pair( f_ci , f_ci + 1 ) );

 branches[ 1 ] = new BinaryKnapsackBlockRngdChange(
                                   BinaryKnapsackBlockChange::eFixX ,
                                   std::move( one ) ,
                                   std::make_pair( f_ci , f_ci + 1 ) );

 return branches;
}

/*--------------------------------------------------------------------------*/

Change * GreedyRelaxationBinaryKnapsackSolver::apply( Change * chg ,
                                                      bool doUndo )
{
 return chg->apply( f_Block , doUndo );

 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool GreedyRelaxationBinaryKnapsackSolver::has_true_var_solution( void ) {
 // the rounded greedy solution (see rounded_x()) always exists
 return true;
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
