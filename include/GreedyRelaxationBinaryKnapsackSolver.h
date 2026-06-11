/*--------------------------------------------------------------------------*/
/*-------------- File GreedyRelaxationBinaryKnapsackSolver.h ---------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class
 * GreedyRelaxationBinaryKnapsackSolver, which implements the Solver concept
 * [see Solver.h] for solving the continuous relaxation of a Knapsack problem
 * represented by a BinaryKnapsackBlock.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __GreedyRelaxationBinaryKnapsackSolver
 #define __GreedyRelaxationBinaryKnapsackSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackSolver.h"

#include "ChangeSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------- CLASS GreedyRelaxationBinaryKnapsackSolver ---------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver of the continuous (Dantzig) relaxation of a BinaryKnapsackBlock
/** Thin wrapper around the machinery of BinaryKnapsackSolver: compute()
 * normalizes the instance and solves its continuous relaxation by the greedy
 * fractional fill, which also identifies the *critical* (fractional) item -
 * the natural branching candidate, see branch(). Besides the relaxation
 * value, valid true lower / upper bounds on the original mixed-binary
 * problem are available [see get_true_lb() / get_true_ub()]. */

class GreedyRelaxationBinaryKnapsackSolver : public BinaryKnapsackSolver ,
                                             public RelaxationSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor

 GreedyRelaxationBinaryKnapsackSolver() : BinaryKnapsackSolver() ,
                                          f_ci( 0 ) ,
                                          obj( - Inf< double >() ) {
  f_fi = FracInfo{ -1 , 1 , false , false , 0 , 0 };
  }

/*--------------------------------------------------------------------------*/
 /// destructor

 ~GreedyRelaxationBinaryKnapsackSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving a relaxation of the Binary Knapsack encoded by the current
 * BinaryKnapsackBlock @{ */

 /// solve the continuous relaxation of the BinaryKnapsackBlock

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// return a valid lower bound on the optimal objective function value
 /** For a minimisation problem the relaxation optimum itself is a valid
  * lower bound. For a maximisation problem the bound is the value of a
  * feasible solution of the original problem: the greedy one without the
  * fractional part \f$ \phi \, p_c \f$ of the critical item \f$ c \f$ - or
  * the greedy solution itself when there is no critical item or it is a
  * continuous variable (the solution is then feasible as it is). */

 OFValue get_true_lb( void ) override {
  return( f_sense ? true_feasible_value() : - obj );
  }

/*--------------------------------------------------------------------------*/
 /// return a valid upper bound on the optimal objective function value
 /** Symmetric to get_true_lb(): the relaxation optimum for a maximisation
  * problem, (minus) the value of the feasible greedy solution without the
  * critical fraction for a minimisation one. */

 OFValue get_true_ub( void ) override {
  return( f_sense ? obj : - true_feasible_value() );
  }

/*--------------------------------------------------------------------------*/
 /// tells whether a true solution (a solution of the true original problem
 /// and not of the relaxed one solved by this Solver) is available
 /** Called after compute() this method has to return true if a true solution
  * of the original problem (not the relaxed one solved by this Solver)
  * is available to be read with get_true_var_solution().
  *
  * Once "the first" solution (if ever) has been read, new ones may be
  * produced, if the Solver allows it, by means of new_true_var_solution().*/

 bool has_true_var_solution( void ) override;

/*--------------------------------------------------------------------------*/
 /// write the current true (rounded greedy) solution in the Block

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
  auto sol = rounded_x();
  BKB->set_x( sol.begin() );
  }

/*--------------------------------------------------------------------------*/
 /// valid lower bound on the optimal value of the TRUE problem
 /** For a maximisation problem this is the value of the rounded greedy
  * solution [see rounded_x()], for a minimisation one the value of the
  * relaxation. */

 OFValue get_lb( void ) override { return( get_true_lb() ); }

 /// valid upper bound on the optimal value of the TRUE problem (symmetric)

 OFValue get_ub( void ) override { return( get_true_ub() ); }

/*--------------------------------------------------------------------------*/
 /// write the current true solution in the variables of the Block

 bool new_true_var_solution( void ) { return false; }

/*--------------------------------------------------------------------------*/
 /// write the current solution in the variables of the BinaryKnapsackBlock

 void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// return the value of the (current) solution
 /** Return the value of the current solution. Change the sign according
  * to the sense of the problem (f_sense). */

 OFValue get_var_value() override { return f_sense ? obj : - obj; }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /// branch on the critical item: the two Change fix it to 0 and to 1

 std::vector< Change * > branch() override;

/*--------------------------------------------------------------------------*/
 /// physically construct the true (rounded greedy) Solution

 Solution * get_true_solution( Configuration * solc = nullptr ) override {
  return( new BinaryKnapsackSolution( rounded_x() ) );
  }

/*--------------------------------------------------------------------------*/
 /// apply the Change; (un)fixing Changes are applied INTERNALLY
 /** The (un)fixing Changes produced by branch() - and, in general, any
  * BinaryKnapsackBlockChange of eFixX / eUnfixX type - are applied to the
  * internal mirror of the instance [see
  * BinaryKnapsackSolver::apply_to_mirror()], NOT to the Block: the next
  * compute() then re-solves the relaxation under the new fixings. Any
  * other Change is forwarded to the Block. */

 Change * apply( Change * chg , bool doUndo = false ) override {
  auto CHG = dynamic_cast< BinaryKnapsackBlockChange * >( chg );

  if( ! CHG )
   throw( std::invalid_argument( "GreedyRelaxationBinaryKnapsackSolver::"
          "apply: the Change must be a BinaryKnapsackBlockChange" ) );

  if( ( CHG->type() != BinaryKnapsackBlockChange::eFixX ) &&
      ( CHG->type() != BinaryKnapsackBlockChange::eUnfixX ) )
   return( CHG->apply( f_Block , doUndo ) );  // not an (un)fix: to the Block

  return( apply_to_mirror( CHG , doUndo ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// the value (in the maximisation sense) of the greedy solution made
 /// feasible for the TRUE problem by dropping the critical fraction

 double true_feasible_value( void ) const {
  if( ( f_fi.orig < 0 ) || f_fi.cont )   // feasible as it is
   return( obj );
  return( obj - f_fi.cfrac * f_fi.cp );
  }

/*--------------------------------------------------------------------------*/
 /// the greedy solution rounded to a true-feasible one
 /** Returns the greedy solution with the critical item, if any and not a
  * continuous one, rounded away (to 0, i.e., to 1 in the original space if
  * complemented): the solution whose value true_feasible_value() returns. */

 std::vector< double > rounded_x( void ) const {
  std::vector< double > sol( f_x );
  if( ( f_fi.orig >= 0 ) && ( ! f_fi.cont ) )
   sol[ f_fi.orig ] = f_fi.comp ? 1 : 0;
  return( sol );
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 Index f_ci;        ///< index (in the Block) of the critical item
 FracInfo f_fi;     ///< the critical item of the last compute()
 double obj;        ///< the value of the relaxation (maximisation sense)

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;
 // insert GreedyRelaxationBinaryKnapsackSolver in the factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }; // end( class( GreedyRelaxationBinaryKnapsackSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* GreedyRelaxationBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------ End File GreedyRelaxationBinaryKnapsackSolver.h -------------*/
/*--------------------------------------------------------------------------*/
