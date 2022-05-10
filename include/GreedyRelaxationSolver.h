/*--------------------------------------------------------------------------*/
/*---------------------- File GreedyRelaxationSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class GreedyRelaxationSolver, which
 * implements the RelaxationSolver concept [see Solver.h] for solving a 
 * relaxation of a Knapsack problem represented by a BinaryKnapsackKBlock.
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


class GreedyRelaxationSolver : public Solver {  // public RelaxationSolver

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
 using c_Index = Block::c_Index;   

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;

 /// tolerance for integrality property of weights- - - - - - - - - - - - - -

 static constexpr double WeightIntegrality = 1e-06;

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor

GreedyRelaxationSolver() : Solver() {} // RelaxationSolver()

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
/** Return the the value of the current solution. */

OFValue get_var_value() override { return f_sense ? obj : - obj; }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the GreedyRelaxationSolver @{ */

/*--------------------------------------------------------------------------*/
/// set the "double" parameters of GreedyRelaxationSolver
/* set the "double" parameters of GreedyRelaxationSolver.
 *  */

void set_par( idx_type par , double value ) override;

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

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

/* data of the Binary Knapsack instance - - - - - - - - - - - - - - - - - - */

 Index f_N;                     ///< the number of Items 
 double f_C;                    ///< the Capacity of the Knapsack
 std::vector< int > v_W;        ///< vector of Weights          
 std::vector< double > v_P;     ///< vector of Profits
 std::vector< bool > v_I;       ///< vector of Integrality (Binary/Continuous)
 bool f_sense;                  ///< the sense of the objective

/* handling of the continuous part- - - - - - - - - - - - - - - - - - - - - */

 Index countCont;        ///< counter for the number of continuous variables
 Subset indexContinuous; ///< indexes of the continuous variables

/* handling of the solution - - - - - - - - - - - - - - - - - - - - - - - - */

 Index besth;            ///< best height of the integer part
 Index lastIndex;        ///< last continuous variables index
 double lastVar;         ///< fraction of solution of lastIndex

 double obj;                            ///< the value of the objective
 std::vector< double > v_x;             ///< vector of variables
 
/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// load the Binary Knapsack instance and perform the preprocessing    

 void load();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// process all the pending modifications 

 void process_outstanding_Modification();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if item i is fixed to 0 because the variable is fixed or the
/// item can be preprocessed (i.e. it has positive weight and negative profit)  

 bool isFixed0( Index i ){
  
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  // if the variable is fixed to 0
  if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) ) < 1e-6 )      
   return( true );                                      

  // if the variable is fixed to 1
  if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) - 1 ) < 1e-6 )      
   return( false );                                     
   
  if( v_W[ i ] >= 0 && v_P[ i ] <= 0 )          // if the item has positive
   return( true );                              // weight and negative profit

  return( false );

 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if item i is fixed to 1 because the variable is fixed or the
/// item can be preprocessed (i.e. it has negative weight and positive profit)  

 bool isFixed1( Index i ){
  
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  // if the variable is fixed to 0
  if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) ) < 1e-6 )      
   return( false );                                     

  // if the variable is fixed to 1
  if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) - 1 ) < 1e-6 )      
   return( true );                                      
   
  if( v_W[ i ] <= 0 && v_P[ i ] >= 0 )          // if the item has negative
   return( true );                              // weight and positive profit

  return( false );
  
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if item i is fixed because the variable is fixed or the item
/// can be preprocessed  

 bool isFixed( Index i ){                       // if the variable is fixed
  return( isFixed0( i ) || isFixed1( i ) );     // to 1 or to 0
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if the item is not fixed and it has negative weight and profit
/// isNeg( i ) is true if the item is not fixed and it has negative weight 
/// and profit [see details in compute()]. 

 bool isNeg( Index i ){
 
  if( isFixed( i ) )                        // if the variable is fixed
   return( false ); 

  if( v_W[ i ] < 0 && v_P[ i ] < 0 )        // if the weight and the profit
   return( true );                          // are negative

  return( false );

 } 
 
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

