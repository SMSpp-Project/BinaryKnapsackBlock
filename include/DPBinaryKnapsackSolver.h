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
 * exception is thrown. Capacity and Profits can be double. 
 *
 * Even if compute() always solves the maximization problem, the objective
 * sense of the problem encoded in the BinaryKnapsackBlock can be either Min 
 * or Max and all the necessary transformations are automatically handled. 
 *
 * The implemented algorithm also manages the presence of items with negative 
 * weights.                                                                 */                         


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

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;
    
 /// public enum for the algorithmic parameters

 enum dbl_par_type_DPBKSlv{
  dblReopt = dblLastAlgPar,
 };

 /// define struct for data stored in a graph G constructed by the DP algorithm

 typedef struct{
  std::vector< double > lab;
  std::vector< bool > pred;
 } slice;


 typedef struct{
  bool mod;
  Range rng;
  Subset nms;  
 } Modification;

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
                           obj( - Inf< double >() ) , start_item( 0 ) , 
                           reopt( 0 ){
                            
                            G.resize( 1 );          // initialize dummy node 
                            G[ 0 ].lab.resize( 1 ); // in the origin
                            G[ 0 ].lab[ 0 ] = 0;

                            modRng_items.first = + Inf< int >();
                            modRng_items.second = 0;

                            Mod.resize( 6 );

                            for( auto & m : Mod ){
                             m.mod = false; 
                             m.rng.first = + Inf< int >();
                             m.rng.second = 0;
                            }

                            set_par( dblReopt , reopt );
                           
                           }

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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return a valid lower bound on the optimal objective function value

virtual OFValue get_lb( void ) override{ return get_var_value(); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return a valid upper bound on the optimal objective function value

virtual OFValue get_ub( void ) override{ return get_var_value(); }

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

void set_par( idx_type par , double value ) override;

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

 int f_NItems;                          ///< the number of Items 
 int f_C;                               ///< the Capacity of the Knapsack
 std::vector< int > v_W;                ///< vector of Weights          
 std::vector< double > v_P;             ///< vector of Profits
 bool f_sense;                          ///< the sense of the objective 

 double obj;                            ///< the value of the objective
 std::vector< bool > v_x;               ///< vector of binary variables

/* data of the graph constructed by the DP algorithm- - - - - - - - - - - - */

 std::vector< slice > G;                ///< DP Graph 

 int start_item;                        ///< index of the starting item

/* preprocessing and modifications data - - - - - - - - - - - - - - - - - - */

 std::vector< int > items;              ///< preprocessing informations
   
 std::vector< Modification > Mod;       ///< vector of modifications

 Range modRng_items;                    ///< range of modified items
 Subset modSbst_items;                  ///< subset of modified items

/* algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - */

 double reopt;                          ///< reoptimization parameter                      
 
 int step;

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
 /// compute solution

 void compute_var_solution();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// perform preprocessing

 void preprocessing( Range rng );

 void preprocessing( Subset & nms );
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// modifications
 
 void capacity_Modification();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void sense_Modification();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void fixX_Modification( Range rng );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void fixX_Modification( Subset & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void unFixX_Modification( Range rng );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void unFixX_Modification( Subset & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void weight_Modification( Range rng );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void weight_Modification( Subset & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void profit_Modification( Range rng );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void profit_Modification( Subset & nms );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert DPBinaryKnapsackSolver in the factory

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









