/*--------------------------------------------------------------------------*/
/*---------------------- File GreedyRelaxationSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class GreedyRelaxationSolver, which
 * implements the Solver concept [see Solver.h] for solving the continuous
 * relaxation of a Knapsack problem represented by a BinaryKnapsackBlock.
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

#ifndef __GreedyRelaxationSolver
 #define __GreedyRelaxationSolver  
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackBlock.h"

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
/*---------------------- CLASS GreedyRelaxationSolver ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// 
/**   */                       


class GreedyRelaxationSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 @{ */

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor

GreedyRelaxationSolver() : Solver() , f_N( 0 ) , f_C( 0 ) , f_sense( true ) ,
                           f_ci( 0 ) , f_ciVal( 0 ) , obj( -Inf< double >() )  
                           {} 

/*--------------------------------------------------------------------------*/
 /// destructor

~GreedyRelaxationSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// set the (pointer to the) Block that the Solver has to solve

void set_Block( Block * block ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving a relaxation of the Binary Knapsack encoded by the current 
 * BinaryKnapsackBlock @{ */

/// 
/**  */

int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// return a valid lower bound on the optimal objective function value

 OFValue get_true_lb( void ) {
  
  // if it is a minimization problem, the optimal value is a lower bound
  // for the original problem  
  if( ! f_sense )
   return( - obj );
  
  // otherwise it is a maximization problem and the solution without the
  // critical item is a feasible solution and it provides a lower bound for
  // the original problem
  if( f_ciVal == 1 )
   return( obj );

  if( v_W[ f_ci ] < 0 ) {
   return( obj + ( 1 - f_ciVal ) * v_P[ f_ci ] ); 
   }
  
  return( obj - f_ciVal * v_P[ f_ci ] );
  }

/*--------------------------------------------------------------------------*/
 /// return a valid upper bound on the optimal objective function value
 /** */

 OFValue get_true_ub( void ) {
  
  // if it is a maximization problem, the optimal value is un upper bound
  // for the original problem  
  if( f_sense )
   return( obj );
  
  // otherwise it is a minimization problem and the solution without the
  // critical item is a feasible solution and it provides an upper bound for
  // the original problem
  if( f_ciVal == 1 )
   return( - obj );
  
  if( v_W[ f_ci ] < 0 ) {
   return( - obj - ( 1 - f_ciVal ) * v_P[ f_ci ] );
   }

  return( - obj + f_ciVal * v_P[ f_ci ] );
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

bool has_true_var_solution( void );

/*--------------------------------------------------------------------------*/
/// write the current true solution in the variables of the Block

void get_true_var_solution( Configuration * solc = nullptr ) {}

/*--------------------------------------------------------------------------*/
/// write the current true solution in the variables of the Block

bool new_true_var_solution( void ) { return false; }

/*--------------------------------------------------------------------------*/
/// write the current solution in the variables of the BinaryKnapsackBlock

void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// return the value of the (current) solution
/** Return the the value of the current solution. Change the sign according
 * to the sense of the problem (f_sense). */

OFValue get_var_value() override { return f_sense ? obj : - obj; }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

/** GreedyRelaxationSolver::add_Modification() is defined to properly react
 * to NBModification, i.e. the Binary Knapsack instance must be reloaded and
 * the list of modification must be cleared. */
 
void add_Modification( sp_Mod &mod ) override;

/*--------------------------------------------------------------------------*/

std::vector< Change * > branch();

/*--------------------------------------------------------------------------*/

Change * apply( Change * , bool doUndo = false );

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/* data of the Binary Knapsack instance - - - - - - - - - - - - - - - - - - */

 Index f_N;                           ///< the number of Items 
 double f_C;                          ///< the Capacity of the Knapsack
 std::vector< double > v_W;           ///< vector of Weights          
 std::vector< double > v_P;           ///< vector of Profits
 std::vector< unsigned char > v_fxd;  ///< vector saying how the x are fixed
                                      /* < v_fxd[ i ] indicates if x_i is
                                       * fixed, with the following encoding:
                                       * 0 = not fixed , 
                                       * 1 = fixed to 0 , 
                                       * 2 = fixed to 1                     */
 bool f_sense;                        ///< the sense of the objective
 
 Index f_ci;                          ///< Index of the critical item
 double f_ciVal;                      ///< Variable value of the critical  
 double obj;                          ///< the value of the objective
 std::vector< double > v_x;           ///< vector of variables

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/
 
/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// load the Binary Knapsack instance 

 void load();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// process all the pending modifications 

 void process_outstanding_Modification();

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert GreedyRelaxationSolver in the factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }; // end( class( GreedyRelaxationSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* GreedyRelaxationSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File GreedyRelaxationSolver.h -------------------*/
/*--------------------------------------------------------------------------*/

