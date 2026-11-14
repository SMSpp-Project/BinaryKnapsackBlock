/*--------------------------------------------------------------------------*/
/*---------------------- File DPBinaryKnapsackSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DPBinaryKnapsackSolver, which
 * implements the Solver concept [see Solver.h] for solving Knapsack problems
 * with linear objective and both integer and continuous variables, as
 * represented by a BinaryKnapsackBlock, using a standard Dynamic
 * Programming approach to deal with the integer variables combined with the
 * exact greedy algorithm to deal with the continuous ones.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Francesca Demelas \n
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Federica Di Pasquale, Francesca Demelas,
 *                      Antonio Frangioni
 */
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
 *  @{ */

 using Index = Block::Index;
 using Subset = Block::Subset;
    
 /// public enum for the algorithmic parameters - - - - - - - - - - - - - - - 

 enum dbl_par_type_DPBKSlv {
  dblReopt = dblLastAlgPar,             ///< reoptimization parameter
  dblLastDPBKSlvPar
  };

 enum int_par_type_DPBKSlv {
  intBoundPrune = intLastAlgPar ,  ///< pruning of the states by bounds
  intLastDPBKSlvPar
  };

 /// tolerance for integrality property of weights- - - - - - - - - - - - - -

 static constexpr double WeightIntegrality = 1e-06;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor

 DPBinaryKnapsackSolver()
     : Solver() , f_N( 0 ) , f_C( 0 ) , f_sense( true ) , besth( 0 ) ,
       obj( -Inf< double >() ) , start_item( 0 ) , reopt( 0 ) , step( 1 ) ,
       f_bound_prune( 0 ) , f_tainted( false ) , f_expanded( 0 ) ,
       f_dominated( 0 ) , f_bounded( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor

 ~DPBinaryKnapsackSolver() override = default;


/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Binary Knapsack encoded by the current
 *  BinaryKnapsackBlock
 *  @{ */

/// Solve the Binary Knapsack Problem encoded in the BinaryKnapsackBlock
/** Solve the Binary Knapsack Problem encoded in the BinaryKnapsackBlock,
* with both integers and continuous variables.
* 
* The standard dynamic programming approach is used to solve the integer part
* of the problem (the problem restricted only to integer variables). The 
* solution of the continuous part is then computed using the greedy algorithm.
*/ 
/*--------------------------------------------------------------------------*/
/*--------------- METHOD FOR SOLVING THE INTEGER KNAPSACK ------------------*/
/*--------------------------------------------------------------------------*/    
/*
* The implemented algorithm processes one item at a time and iteratevily  
* constructs a graph G with the solutions of each sub-problem.
*
* In particular, the graph G is constructed as follows: 
*
* let SP( i , j ) the sub-problem where the objective is to select a subset 
* of the first i items, that maximizes the total profit and whose total 
* weight is equal to j. For each of these sub-problems, let consider a node 
* u_{i,j} with two outgoing oriented arcs: 
*
*  - A "horizontal" arc: ( u_{i,j} , u_{i+1,j} ) with "cost" 0                
*                          
*  - A "diagonal" arc: ( u_{i,j} , u_{i+1,j+w} ) with "cost" p
*                       
* where w and p are the weight and the profit of the ( i + 1 )-th item.
*                      
*                      X 
*                  X  /
*                 /  / X
*        --------/--/-/----------------------O - -  Capacity          
*               O--O-/-O                     O - -  w0 + w1 
* .            /  / /                        .
*           p1/  / O---O                     . - -  w0 + w2
*            /  O-/O---O      ...            . - -  w1
*           /  / /     O                     .
*          O--/-O--O--/O                     . - -  w0 
*      p0 /p1/    /  /                       .
*        /  /    /  /                        
*       *--O----O--O---O                     O - -  0
*          ^    ^  ^   ^                     ^
* items    0    1  2   3      ...            N 
* 
* Basically, starting from a dummy node *, for each new item two possibilites
* are considered: either the next item is selected (diagonal arc with cost p), 
* or it is discarded (horizontal arc with cost 0). 
*  
* The Binary Knapsack problem is equivalent to finding a path in this graph
* from node * to a node u_{ N , j } with j <= Capacity of the Knapsack,
* that maximizes the total profit. To find a Longest path in G, a label is 
* assigned to each node; a label is the optimal value of the sub-problem
* represented by that node and can be computed using dynamic programming.
* To simplify the explanation, let us call slice[ i ] the set of nodes u_{i,j}
* with fixed i (a "vertical slice" in the graph in figure). 
*  
* - Assume all the labels of slice[ i ] have already been computed
* - Each node in slice[ i + 1 ] can be reached either from an horizontal arc
*   or from a diagonal arc
* - The labels of slice[ i + 1 ] can be computed by comparing the two 
*   possibilities (choosing the one with the best total profit).  
* 
* Therefore, G is implemented with two vectors (of vectors) containing
* labels (lab) and predecessors (pred). Each vector has N + 1 entries, one 
* for each item + the dummy node. 
*
* Labels (of type double) and predecessors (of type bool) are s.t.
*   
*   lab[ i ][ j ]    corresponds to node u_{i,j} and contains its label,
*                    that is the optimal value of SP( i , j ).
*
*   pred[ i ][ j ]   corresponds to the last arc that has been selected in 
*                    order to obtain the value in lab[ i ][ j ]. 
*                    It is true if it is a "diagonal" arc, i.e. the 
*                    corresponding solution contains the i-th item; it is 
*                    false otherwise.     
*
* At each iteration i, each entry of lab[ i + 1 ] is computed starting 
* from lab[ i ], by comparing the profits of the two possible path reaching
* the corresponding node (choosing the best one), and pred[ i ] is updated 
* accordingly.
* 
* Eventually the last set of labels lab[ N ] contains the optimal values of 
* the problems containing all the items. The optimal value of the Binary 
* Knapsack problem is the best value among those in lab[ N ][ j ] with j
* less or equal then the Capacity of the Knapsack, and the optimal solution 
* can be reconstructed from the vectors of predecessors.
*
* Note that:
*
*
* - There is no need to store the vectors of labels, since each of them 
*   requires only the labels of the previous iteration to be computed. 
*   Therefore, only two vectors of labels are used: one for the current 
*   labels (currlab) and one for the next labels (nextlab). 
*   However, for reoptimization purposes, some of them are stored in 
*   lab[ i ] according to the reopt algorithmic parameter.
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
*   These conditions are checked in the preprocessing method()
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
*     solution, but discarding it from the computation. Then consider 
*     "another" item with the same weight and profit but of opposite signs 
*     (both positive). At the end of the algorithm, if the added item has 
*     been selected, it neutralizes the effect of the initial selection, i.e. 
*     it is equivalent to not select the original item. Conversly, if the 
*     added item has not been selected, it is equivalent to select the 
*     original item. 
*
*
* - The implemented algorithm always solves a maximization problem. However,
*   if the problem encoded in the BinaryKnapsackBlock is a minimization one,
*   the signs of all the profits are changed immediately when the instance is
*   loaded and f_sense is set accordingly. The sign of the objective can then
*   be properly changed when needed. 
*                                                                           */
/*--------------------------------------------------------------------------*/
/*------------ METHOD FOR SOLVING THE CONTINUOUS KNAPSACK ------------------*/
/*--------------------------------------------------------------------------*/
/*

1) For each height (restricted capacity) consider the best binary solution
   with the correspondent index bestBinIndex
   
2) For each height H, starting from the bigger one (H=B-1):
   -1- compute the solution of the Continuous Knapsack with capacity B-H,
       considering as input also the continuous solution obtained in the 
       previous step:
      *1* sort the components of the continuopus variables
          considering Profits/Weights
      *2* With the new order, start from the first element and
          assign to the variable the maximum value allowed by 
          the residual capacity
   -3- Update/Memorize the optimal objective value obtained in the 
       node lab[H][bestBinIndex[H]]
   -4- update/memorize the complete solution (x_binary,x_continuous)
   
3) For each height H=0,...,B:
   -1- Compare the labels in lab[H][bestBinIndex[H]] in order to find 
       the optimal solution                                                */ 

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// return a valid lower bound on the optimal objective function value
 /** Return a valid lower bound on the optimal objective function value.
  * get_lb() must be called after compute() has been called. */

 OFValue get_lb( void ) override { return( get_var_value() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// return a valid upper bound on the optimal objective function value
 /** Return a valid upper bound on the optimal objective function value.
  * get_ub() must be called after compute() has been called. */

 OFValue get_ub( void ) override { return( get_var_value() ); }

/*--------------------------------------------------------------------------*/
 /// write the current solution in the variables of the BinaryKnapsackBlock

 void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns the current solution as a BinaryKnapsackSolution
 /** Returns the current solution as a BinaryKnapsackSolution [see
  * BinaryKnapsackBlock.h], reconstructed out of the dynamic programming
  * table exactly as get_var_solution() does, but written straight into the
  * Solution rather than into the Variable of the BinaryKnapsackBlock: no
  * abstract representation is therefore required to exist, and the
  * BinaryKnapsackBlock is not written into at all, hence it is not
  * lock()-ed and any number of Solver attached to it can produce their own
  * Solution at the same time. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
/// return the value of the (current) solution
/** Return the value of the current solution.
 * Since the implemented DP algorithm always solves the maximization problem,
 * the sign of the value of the solution must change according to the real 
 * sense (f_sense) of the problem encoded in the BinaryKnapsackBlock. */

 OFValue get_var_value() override { return( f_sense ? obj : - obj ); }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the DPBinaryKnapsackSolver
 *  @{ */

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
 *                 defines how many labels have to be stored in lab.
 * 
 *                 The possibilities are:
 *                 ----------------------------------------------------------
 *                 # labels to store | Indices of stored labels
 *                 ----------------------------------------------------------
 *                   1                [ N ]
 *                   2                [ (1/2)N , N ]
 *                   4                [ (1/4)N , (2/4)N , (3/4)N , N ]
 *                   8                [ (1/8)N , ... , N ]
 *                   ...              ...
 *                   2^k              [ (1/k)N , ... , N ]                
 *                   ...              ...
 *                   N                [ 1 , 2 , ... , N ]
 *
 *                 hence, the number of possibilies is m = log2( N ).
 *                 The correspondence between dblReopt values and one of these
 *                 possibilities is done by dividing the [ 0 , 1 ] interval 
 *                 into m smaller intervals of equal size, and checking to 
 *                 which of these intervals dblReopt belongs. That is, if 
 *                 dblReopt \in k-th interval, then 2^k labels will be stored
 *                 in lab. In particular, it follows that:
 *                 
 *                 dblReopt = 0 -> store only lab[ N ]
 *                 dblReopt = 1 -> store lab[ i ] for all i
 *
 *                 Intermediate dblReopt values are handled by:
 *                 - retrieving the interval k to which they belong
 *                 - defining a step:
 *                      step = N / 2^k
 *                   such that at each multiple i of step the corresponding 
 *                   labels are stored in lab[ i ]               
 *                                                                         */

 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastDPBKSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 double get_dbl_par( const idx_type par ) const override {
  if( par == dblReopt )
   return( reopt );

  return( Solver::get_dflt_dbl_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  if( name == "dblReopt" )
   return( dblReopt );

  return( Solver::dbl_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the "int" parameters of DPBinaryKnapsackSolver
 /** Set the "int" parameters of DPBinaryKnapsackSolver. The only parameter
 * currently present is:
 *
 * - intBoundPrune [0]: if nonzero, a state \f$ ( i , j ) \f$ of the table,
 *   i.e., the best profit \f$ \lambda_{i,j} \f$ of a subset of the first
 *   \f$ i \f$ (integer) items of weight \f$ j \f$, is not expanded when
 *   \f[
 *     \lambda_{i,j} + U_i( C - j ) < z - \epsilon \, ,
 *   \f]
 *   where \f$ U_i( c ) \f$ is the Dantzig bound of the items that are still
 *   to be decided at \f$ i \f$, i.e., the integer items from \f$ i \f$ on
 *   and all the continuous ones, with capacity \f$ c \f$ (the value of the
 *   continuous relaxation of the knapsack of those items), \f$ z \f$ is the
 *   value of the best feasible solution found so far, and \f$ \epsilon =
 *   10^{-9} ( 1 + | z | + \sum_k | p_k | ) \f$ covers the rounding.
 *   \f$ U_i \f$ is admissible, since any completion of the state is a
 *   feasible solution of that relaxation, and therefore a pruned state is
 *   on no optimal path: the optimal value, and the solution, are those of
 *   intBoundPrune = 0 [see dynamic_programming_bound()]. The first
 *   \f$ z \f$ is the greedy solution of the remaining items from the first
 *   state; then, while the states of a column are scanned for the bound,
 *   \f$ z \f$ is raised to \f$ \lambda_{i,j} \f$ plus the profit of the
 *   longest prefix of the items sorted by efficiency that fits in
 *   \f$ C - j \f$, a feasible completion that the bound computes anyway.
 *   The columns are not stored (a pruned column is not one of the full
 *   table), and therefore a solve after one with intBoundPrune restarts
 *   from the first item whatever dblReopt says. */

 void set_par( idx_type par , int value ) override;

 using Solver::set_par;  // keep the other set_par() overloads visible

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastDPBKSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 int get_int_par( const idx_type par ) const override {
  if( par == intBoundPrune )
   return( f_bound_prune );

  return( Solver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 int get_dflt_int_par( const idx_type par ) const override {
  if( par == intBoundPrune )
   return( 0 );

  return( Solver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type int_par_str2idx( const std::string & name ) const override {
  if( name == "intBoundPrune" )
   return( intBoundPrune );

  return( Solver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & int_par_idx2str( const idx_type idx ) const override {
  static const std::string name = "intBoundPrune";
  if( idx == intBoundPrune )
   return( name );

  return( Solver::int_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

/** DPBinaryKnapsackSolver::add_Modification() is defined to properly react
 * to NBModification, i.e. the Binary Knapsack instance must be reloaded and
 * the list of modification must be cleared. */
 
 virtual void add_Modification( sp_Mod &mod ) override;

/*--------------------------------------------------------------------------*/
 /// the kinds of Modification read: the physical ones
 /** The Solver works on the data of the BinaryKnapsackBlock, which the
  * physical Modification describe, and has no use for the abstract ones
  * [see Solver::concerned_by()]. */

 [[nodiscard]] Modification::ModConcern concerned_by( void ) const override {
  return( Modification::eModPhys );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Dynamic Programming to solve the integer knapsack

 virtual void dynamic_programming( Index C );

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the Dynamic Programming of dynamic_programming() with bound pruning
 /** The same table as dynamic_programming(), column by column and with the
  * same dominance between the states of a column, save that a state that
  * passes the dominance is not expanded if its label plus the Dantzig bound
  * of the remaining items is below the incumbent [see intBoundPrune].
  *
  * The result is that of dynamic_programming(). Let \f$ z^* \f$ be the
  * optimal value; a state on an optimal path has
  * \f$ \lambda_{i,j} + U_i( C - j ) \geq z^* \geq z \f$, hence it is
  * never pruned, and so is every state that reaches a node of an optimal
  * path with the label of that node, as it is on an optimal path itself.
  * Therefore the labels of the nodes of the optimal paths, and the arcs
  * that the ties choose at them, are those of the full table. A state that
  * the dominance removes in the full table and not here (because the state
  * that dominates it has been pruned) has a smaller label and a smaller
  * residual capacity than a pruned one, hence its bound is smaller and it
  * is pruned too, and the states that are expanded are a subset of those of
  * the full table. Each column keeps the sorted list of its states, and it
  * is scanned by increasing \f$ j \f$, i.e., decreasing capacity \f$ c = C
  * - j \f$, over the remaining items sorted by decreasing efficiency, so
  * that the break item of \f$ U_i( c ) \f$ only moves backwards: a column
  * costs \f$ O( n + a_i ) \f$, where \f$ a_i \f$ is its number of states,
  * rather than the \f$ O( C ) \f$ of dynamic_programming(), plus the
  * \f$ O( n \log n ) \f$ sorting and the \f$ O( C ) \f$ allocation done
  * once. \p C is the integer residual
  * capacity, as in dynamic_programming(), and \p Cd the residual capacity
  * itself, which the continuous items can fill. */

 void dynamic_programming_bound( Index C , double Cd );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// reconstructs the optimal solution into v_x out of the DP table
 /** Reconstructs the optimal solution out of the dynamic programming table
  * into v_x, which is what both get_var_solution() and get_Solution() need;
  * \p method is the name of the caller, for the exception thrown if
  * compute() has not been called first. */

 void reconstruct_x( const std::string & method );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

/* data of the Binary Knapsack instance - - - - - - - - - - - - - - - - - - */

 Index f_N;                     ///< the number of items 
 double f_C;                    ///< the Capacity of the Knapsack
 bool f_sense;                  ///< the sense of the objective
 std::vector< int > v_W;        ///< vector of Weights          
 std::vector< double > v_P;     ///< vector of Profits
 std::vector< bool > v_I;       ///< vector of Integrality (Binary/Continuous)
 std::vector< unsigned char > v_fxd;  ///< vector saying how the x are fixed
                                      /* < v_fxd[ i ] indicates if x_i is
                                       * fixed, with the following encoding:
                                       * 0 = not fixed , 
                                       * 1 = fixed to 0 , 
                                       * 2 = fixed to 1                     */

/* handling of the continuous part- - - - - - - - - - - - - - - - - - - - - */

 Subset idxCont;               ///< indices of the continuous variables

/* handling of the solution - - - - - - - - - - - - - - - - - - - - - - - - */

 double obj;                   ///< the value of the objective
 std::vector< double > v_x;    ///< vector of variables

 Index besth;                  ///< best height of the integer part

/* data of the graph constructed by the DP algorithm- - - - - - - - - - - - */

 std::vector< std::vector< double > > lab;  ///< labels
 std::vector< std::vector< char > > pred;   ///< predecessors 
 Index start_item;                          ///< index of the starting item

/* preprocessing- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 std::vector< bool > skip;     ///< preprocessed variables to skip

/* algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - */

 double reopt;       ///< reoptimization parameter
 Index step;         ///< step for reoptimization
 int f_bound_prune;  ///< the intBoundPrune parameter
 bool f_tainted;     ///< the last solve has pruned

 /* statistics of the last solve - - - - - - - - - - - - - - - - - - - - - - */

 unsigned long f_expanded;   ///< states expanded by the last compute()
 unsigned long f_dominated;  ///< states removed by the dominance
 unsigned long f_bounded;    ///< states removed by the bound

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
 /// perform the preprocessing and return residual capacity and profit

 std::tuple< double , double > preprocessing();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Greedy algorithm to solve the continuous knapsack

 void greedy_algorithm( double C );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// process all the pending modifications and compute start_item 

 void process_outstanding_Modification();
 
/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert DPBinaryKnapsackSolver in the factory

/*--------------------------------------------------------------------------*/

 }; // end( class( DPBinaryKnapsackSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* DPBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DPBinaryKnapsackSolver.h -------------------*/
/*--------------------------------------------------------------------------*/

