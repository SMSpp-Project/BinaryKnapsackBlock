/*--------------------------------------------------------------------------*/
/*---------------------- File DPBinaryKnapsackSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DPBinaryKnapsackSolver
 #define __DPBinaryKnapsackSolver  
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
/*---------------------- CLASS DPBinaryKnapsackSolver ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Dynamic Programming Solver for BinaryKnapsackBlock
/** The DPBinaryKnapsackSolver implements the Solver interface for the Binary
 * Knapsack Problem described by a BinaryKnapsackBlock using the standard 
 * Dynamic Programming approach.
 *
 * The algorithm assumes that weights of the items are integers, otherwise an
 * exception is thrown. Capacity and Profits can be double. */

class DPBinaryKnapsackSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 @{ */

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor

DPBinaryKnapsackSolver() : Solver() , f_NItems( 0 ) , f_C( 0 ) , f_sense( 1 ),
						   obj( -Inf< double >() ) , prp_W( 0 ) , prp_P( 0 ) , 
						   nW( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor

~DPBinaryKnapsackSolver() override = default;


/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// set the (pointer to the) Block that the Solver has to solve

void set_Block( Block * block ) override;

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Binary Knapsack encoded by the current 
 * BinaryKnapsackBlock @{ */

int compute( bool changedvars = true ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

//virtual OFValue get_lb( void ) override;

//virtual OFValue get_ub( void ) override;

/*--------------------------------------------------------------------------*/
/// write the "current" solution in the variables of the BinaryKnapsackBlock

void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// return the value of the (current) solution
/** Return the the value of the current solution. The DP algorithm solves the 
 * maximization problem. If the problem encoded in the BinaryKnapsackBlock is
 * a minimization problem, when the instance is loaded the signs of all 
 * profits are changed and f_sense is set to 0. The sign of obj must change
 * according to f_sense. */

OFValue get_var_value() override { return f_sense ? obj : - obj; }


/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the DPBinaryKnapsackSolver @{ */


/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

/** DPBinaryKnapsackSolver::add_Modification() is defined to properly react
 * to NBModification, i.e. the Binary Knapsack instance must be reloaded and
 * the list of modification must be cleared. */
 
virtual void add_Modification( sp_Mod &mod ) override;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

/* data of the Binary Knapsack instance - - - - - - - - - - - - - - - - - - */

 int f_NItems;                  		///< the number of Items 
 int f_C;                     			///< the Capacity of the Knapsack
 std::vector< int > v_W;        		///< vector of Weights          
 std::vector< double > v_P;     		///< vector of Profits
 bool f_sense;                  		///< the sense of the objective 

 double obj;							///< the value of the objective
 std::vector< bool > v_x;				///< vector of binary variables

/* data of the graph constructed by the DP algorithm- - - - - - - - - - - - */

 std::vector< std::vector< bool > > pred; ///< Matrix of predecessors
 std::vector< double > labels;			///< vector of labels

/* preprocessing data - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int prp_W;
///< total weight of preprocessed items whose variables are set to 1
double prp_P;			
///< total profit of preprocessed items whose variables are set to 1

std::vector< int > pi;
///< indeces of NOT preprocessed items with positive weight
std::vector< int > ni;
///< indeces of NOT preprocessed items with negative weight

int nW; 	///< total weight of NOT preprocessed items with negative weight

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// load the Binary Knapsack instance and perform the preprocessing
 void load();

 /// process all the pending modifications
 void process_outstanding_Modification();

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_h;	// insert DPBinaryKnapsackSolver in the factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
 }; // end( class( DPBinaryKnapsackSolver ) )
}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* DPBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DPBinaryKnapsackSolver.h -------------------*/
/*--------------------------------------------------------------------------*/





