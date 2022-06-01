/*--------------------------------------------------------------------------*/
/*---------------------- File GreedyRelaxationSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class GreedyRelaxationSolver, which
 * implements the RelaxationSolver concept [see RelaxationSolver.h] and the
 * Solver concept [see Solver.h] for solving the continuous relaxation of a  
 * Knapsack problem represented by a BinaryKnapsackBlock.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Antonio Frangioni
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

#include "RelaxationSolver.h"

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


class GreedyRelaxationSolver : public Solver , public RelaxationSolver {  

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
                           obj( -Inf< double >() ) {} 

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
 
virtual void add_Modification( sp_Mod &mod ) override;

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

