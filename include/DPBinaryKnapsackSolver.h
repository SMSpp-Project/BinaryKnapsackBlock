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
 * The implemented algorithm also handles the presence of items with negative 
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

 using Index = Block::Index;
 using c_Index = Block::c_Index;   

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;
    
 /// struct for data stored in a graph G constructed by the DP algorithm

 typedef struct{
  std::vector< double > lab;            ///< vector of labels
  std::vector< bool > pred;             ///< vector of predecessors
 } slice;

 /// public enum for the algorithmic parameters - - - - - - - - - - - - - - - 

 enum dbl_par_type_DPBKSlv{

  dblReopt = dblLastAlgPar,             ///< reoptimization parameter
  
  dblLastDPBKSlvPar

 };



 /// tolerance for integrality property of weights

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

DPBinaryKnapsackSolver() : Solver() , f_NItems( 0 ) , f_C( 0 ) , f_sense( 1 ),
                           obj( - Inf< double >() ) , start_item( 0 ) , 
                           reopt( 0 ){

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
/** Return a valid lower bound on the optimal objective function value.
* get_lb() must be called after compute() has been called. */

OFValue get_lb( void ) override{ return get_var_value(); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return a valid upper bound on the optimal objective function value
/** Return a valid upper bound on the optimal objective function value.
* get_ub() must be called after compute() has been called. */

OFValue get_ub( void ) override{ return get_var_value(); }

/*--------------------------------------------------------------------------*/
/// write the current solution in the variables of the BinaryKnapsackBlock

void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// return the value of the (current) solution
/** Return the the value of the current solution. 
 * Since the implemented DP algorithm always solves the maximization problem,
 * the sign of the value of the solution must change according to the real 
 * sense (f_sense) of the problem encoded in the BinaryKnapsackBlock. */

OFValue get_var_value() override { return f_sense ? obj : - obj; }

/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the DPBinaryKnapsackSolver @{ */

/*--------------------------------------------------------------------------*/
/// set the "double" parameters of DPBinaryKnapsackSolver
/* set the "double" parameters of DPBinaryKnapsackSolver.
 * 
 * The only parameter currently present is:
 * 
 * - dblReopt [0]: Reoptimization parameter. Accepted values in [ 0 , 1 ] 
 *                 dblReopt roughly defines how often the labels computed at
 *                 each iteration of the DP algorithm should be saved.
 *                 
 *                 - 0 if no labels should be saved  
 *                 - 1 if all labels should be saved 
 *                 Intermediate values define a "step" s.t. each time the 
 *                 index i of an item is a multiple of step , the labels
 *                 computed at the corresponding iteration should be saved 
 *                 (unless i is preprocessed).                              */

void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

idx_type get_num_dbl_par( void ) const override{
 return( idx_type( dblLastDPBKSlvPar ) );
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
double get_dbl_par( const idx_type par ) const override{
  
 if( par == dblReopt )
  return( reopt );

 return( Solver::get_dflt_dbl_par( par ) );
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

idx_type dbl_par_str2idx( const std::string & name ) const override {
 
 const auto it = dbl_pars_map.find( name );
 if( it != dbl_pars_map.end() )
  return( it->second );
 
 return( Solver::dbl_par_str2idx( name ) );
}

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

/* algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - */

 double reopt;                          ///< reoptimization parameter                      
 
 int step;                              ///< reoptimization step 


// static fields - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const static std::map< std::string , idx_type > dbl_pars_map;
 ///< the (static const) map for double parameters names

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

  if( BKB->is_fixed( i ) ){
   
   if( BKB->get_x( i ) == 0 )
    return( true );
   
   return( false );
  }
   
  if( v_W[ i ] >= 0 && v_P[ i ] <= 0 ) 
   return( true );  

  return( false );

 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if item i is fixed to 1 because the variable is fixed or the
/// item can be preprocessed (i.e. it has negative weight and positive profit)  

 bool isFixed1( Index i ){
  
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  if( BKB->is_fixed( i ) ){
   
   if( BKB->get_x( i ) == 1 )
    return( true );
   
   return( false );
  }
   
  if( v_W[ i ] <= 0 && v_P[ i ] >= 0 ) 
   return( true );  

  return( false );
  
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if item i is fixed because the variable is fixed or the item
/// can be preprocessed  

 bool isFixed( Index i ){
  return( isFixed0( i ) || isFixed1( i ) );  
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if the item is not fixed and it has negative weight and profit
/// isNeg( i ) is true for the "negative" items processed by the DP algorithm
/// [see details in compute()]. 

 bool isNeg( Index i ){
 
  if( isFixed( i ) )
   return( false ); 

  if( v_W[ i ] < 0 && v_P[ i ] < 0 )
   return( true );

  return( false );

 }
 
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









