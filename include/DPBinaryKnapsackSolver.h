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
 * Knapsack Problem [see BinaryKnapsackBlock.h] using the standard Dynamic
 * Programming approach.
 *
 * The algorithm assumes that the weights of the items are integers (any non
 * integer weight will lead to exception been thrown). Capacity and Profits 
 * can be double.
 *
 * There are no restrictions on the weights and profits sign (both positive 
 * and negative values are allowed), and the objective sense of the 
 * problem can be either Min or Max.                                        */                       


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
    
 /// data structure for the graph G constructed by the DP algorithm - - - - -
 
 typedef struct{
  std::vector< double > lab;            ///< vector of labels
  std::vector< bool > pred;             ///< vector of predecessors
 } slice;

 /// public enum for the algorithmic parameters - - - - - - - - - - - - - - - 

 enum dbl_par_type_DPBKSlv{

  dblReopt = dblLastAlgPar,             ///< reoptimization parameter
  
  dblLastDPBKSlvPar

 };

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

DPBinaryKnapsackSolver() : Solver() , f_N( 0 ) , f_C( 0 ) , 
                           f_sense( true ), obj( - Inf< double >() ) , 
                           start_item( 0 ) , reopt( 0 ){}

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

/// Solve the Binary Knapsack Problem encoded in the BinaryKnapsackBlock
/** Solve the Binary Knapsack Problem encoded in the BinaryKnapsackBlock,
* using the standard Dynamic Programming approach. 
*
* The implemented algorithm processes one item at a time and iteratevily  
* constructs a graph G with the solutions of each sub-problem.
*
* In particular, the graph G is constructed as follows: 
*
* let SP( i , j ) the sub-problem where the objective is to select a subset 
* of the first i items, that maximizes the total profit and whose total 
* weight is equal to j. For each of these sub-problems, let consider a node 
* u_{ i , j } with two outgoing oriented arcs: 
*
*  - A "horizontal" arc: ( u_{ i , j } , u_{ i + 1 , j } ) with "cost" 0                
*                          
*  - A "diagonal" arc: ( u_{ i , j } , u_{ i + 1 , j + w } ) with "cost" p
*                       
* where w and p are the weight and the profit of the ( i + 1 )-th item. 
* The Binary Knapsack problem is equivalent to finding a path in this graph
* that maximizes the total profit, from a "dummy" node u_{ 0 , 0 } to a node 
* u_{ f_N , j } with j <= Capacity of the Knapsack. In a solution path, 
* selecting a "diagonal" arc means that the ( i + 1 )-th item has been 
* selected, whereas selecting a "horizontal" arc is equivalent to discard
* the item.  
* 
* Therefore, G is implemented as a vector with ( f_N + 1 ) entries, one
* for each item + the "dummy" node. Each entry contains a "slice" that is a  
* data structure with two vectors:  G[ i ].lab (or "labels") and G[ i ].pred 
* (or "predecessors").
*
* Labels (of type double) and predecessors (of type bool) are vectors s.t.
*   
*   G[ i ].lab[ j ]  corresponds to node u_{ i , j } and contains the optimal 
*                    value of SP( i , j )
*
*   G[ i ].pred[ j ] corresponds to the last arc that has been selected in 
*                    order to obtain the value in G[ i ].lab[ j ]. 
*                    It is true if it is a "diagonal" arc, i.e. the 
*                    corresponding solution contains the i-th item; it is 
*                    false otherwise.     
*
* At each iteration i, each entry of G[ i + 1 ].lab is computed starting 
* from G[ i ].lab, by comparing the profits of all the possible path reaching
* the corresponding node (choosing the best one), and G[ i ].pred is updated 
* accordingly.
* 
* Eventually G[ f_N ].lab contains the optimal values of the problems 
* containing all the items. The optimal value of the Binary Knapsack problem 
* is the the best value among those in G[ f_N ].lab[ j ] with j less or 
* equal then the Capacity of the Knapsack, and the optimal solution can be 
* reconstructed from the vectors of predecessors.
*
* Note that:
*
*
* - There is no need to store the vectors of labels, since each of them 
*   requires only the labels of the previous iteration to be computed. 
*   Therefore, only two vectors of labels are used: one for the current 
*   labels (currlab) and one for the next labels (nextlab). 
*   However, for reoptimization purposes, some of them are stored in 
*   G[ i ].lab according to the reopt algorithmic parameter.
*
*
* - Some of the items are pre-processed and, therefore, discarded from the
*   computation. An item i is pre-processed if:
*
*       - the corresponding variable is fixed 
*       - it has positive weight and negative profit -> set v_x[ i ] = 0
*       - it has negative weight and positive profit -> set v_x[ i ] = 1
*       - its weight exceeds the "residual" capacity, i.e. the capacity
*         obtained subtrancting the weight of the items that are fixed to 1           
*
*   These conditions are checked at the beginning of each iteration. The last
*   one is directly checked, whereas the first three are checked by calling
*   the methods is_fixed0() and is_fixed1().
*
*
* - The algorithm, as described above, only deals with item with positive 
*   weights. However, items with negative weight can be treated as follows:
*   
*   - if the weight is negative but the profit is positive, then the item is
*     pre-processed.
*
*   - if both weight and profit are negative, the idea is to "select" the 
*     item, and therefore only update the capacity and the profit of the 
*     solution, but discarding it from the computation. Then consider "another"  
*     item with the same weight and profit but of opposite signs (both 
*     positive). At the end of the algorithm, if the added item has been 
*     selected, it neutralizes the effect of the initial selection, i.e. it is
*     equivalent to not select the original item. Conversly, if the added item
*     has not been selected, it is equivalent to select the original item.
*     The method is_Neg() returns true if this transformation must be done.
*
*
* - The implemented algorithm always solves a maximization problem. However,
*   if the problem encoded in the BinaryKnapsackBlock is a minimization one,
*   the signs of all the profits are changed immediately when the instance is
*   loaded and f_sense is set accordingly. The sign of the objective can then
*   be properly changed when needed. 
*                                                                           */

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
 * - dblReopt [0]: Reoptimization parameter. Accepted values in [ 0 , 1 ].
 *                 DPBinaryKnapsackSolver implements a reoptimization 
 *                 technique based on the idea of storing intermediate labels
 *                 computed during the first execution of the DP algorithm, 
 *                 such that, in the next call of compute(), it is possible 
 *                 to re-start from one of these intermediate points. Storing
 *                 labels can be computationally expensive, hence dblReopt   
 *                 defines how many labels have to be stored in G.
 * 
 *                 The possibilities are:
 *                 ----------------------------------------------------------
 *                 # labels to store | Indexes of stored labels
 *                 ----------------------------------------------------------
 *                   1                [ f_N ]
 *                   2                [ (1/2)f_N , f_N ]
 *                   4                [ (1/4)f_N , (2/4)f_N , (3/4)f_N , f_N ]
 *                   8                [ (1/8)f_N , ... , f_N ]
 *                   ...              ...
 *                   2^k              [ (1/k)f_N , ... , f_N ]                
 *                   ...              ...
 *                   f_N              [ 1 , 2 , ... , f_N ]
 *
 *                 hence, the number of possibilies is m = log2( f_N ).
 *                 The correspondence between dblReopt values and one of these
 *                 possibilities is done by dividing the [ 0 , 1 ] interval 
 *                 into m smaller intervals of equal size, and checking to 
 *                 which of these intervals dblReopt belongs. That is, if 
 *                 dblReopt \in k-th interval, then 2^k labels will be stored
 *                 in G. In particular, it follows that:
 *                 
 *                 dblReopt = 0 -> store only G[ f_N ].lab
 *                 dblReopt = 1 -> store G[ i ].lab for all i
 *
 *                 Intermediate dblReopt values are handled by:
 *                 - retrieving the interval k to which they belong
 *                 - defining a step:
 *                      step = f_N / 2^k
 *                   such that at each multiple i of step the corresponding 
 *                   labels are stored in G[ i ].lab
 *                 The step is computed in the compute_step() method.                
 *                                                                         */

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

 int f_N;                               ///< the number of Items 
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

  if( BKB->is_fixed( i ) && BKB->get_x( i ) == 0 )      // if the variable 
   return( true );                                      // is fixed to 0

  if( BKB->is_fixed( i ) && BKB->get_x( i ) == 1 )      // if the variable 
   return( false );                                     // is fixed to 1
   
  if( v_W[ i ] >= 0 && v_P[ i ] <= 0 )          // if the item has positive
   return( true );                              // weight and negative profit

  return( false );

 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// return true if item i is fixed to 1 because the variable is fixed or the
/// item can be preprocessed (i.e. it has negative weight and positive profit)  

 bool isFixed1( Index i ){
  
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  if( BKB->is_fixed( i ) && BKB->get_x( i ) == 0 )      // if the variable 
   return( false );                                     // is fixed to 0

  if( BKB->is_fixed( i ) && BKB->get_x( i ) == 1 )      // if the variable 
   return( true );                                      // is fixed to 1
   
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// compute the step for reoptimization [see set_par() dblReopt for details ]
 
 int compute_step(){

  // find the interval to which reopt belongs
  int k = std::floor( reopt * std::log2( f_N ) );

  // define the step 
  int step = std::floor( f_N / std::exp2( k ) );
  
  return step;
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









