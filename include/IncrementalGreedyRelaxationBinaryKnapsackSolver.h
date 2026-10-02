/*--------------------------------------------------------------------------*/
/*------ File IncrementalGreedyRelaxationBinaryKnapsackSolver.h ------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* classes
 * IncrementalGreedyChangeBinaryKnapsackSolver and
 * IncrementalGreedyRelaxationBinaryKnapsackSolver, which solve the
 * continuous (Dantzig) relaxation of a Binary Knapsack problem represented
 * by a BinaryKnapsackBlock, re-optimizing it after each Change of the Block
 * (typically, the fixing of a variable) rather than solving it from scratch,
 * which makes them suitable as the bounding solver at the nodes of a
 * Branch-and-Bound algorithm.
 *
 * IncrementalGreedyChangeBinaryKnapsackSolver is the ChangeSolver [see
 * ChangeSolver.h] that keeps the relaxation along the Changes;
 * IncrementalGreedyRelaxationBinaryKnapsackSolver extends it with the
 * RelaxationSolver interface, additionally providing bounds and a solution
 * for the original integer problem (when available), as well as the two
 * sub-problems obtained by branching on the critical item.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Antonio Frangioni, Filippo Magi
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __IncrementalGreedyRelaxationBinaryKnapsackSolver
 #define __IncrementalGreedyRelaxationBinaryKnapsackSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackBlock.h"

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
/*------------ CLASS IncrementalGreedyChangeBinaryKnapsackSolver -----------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// ChangeSolver of the continuous relaxation of a BinaryKnapsackBlock
/** Solves the continuous (Dantzig) relaxation of the BinaryKnapsackBlock by
 * the greedy fill along the items sorted by nonincreasing efficiency, and
 * keeps it across the Changes applied to it [see apply()]: a fixing updates
 * the residual capacity and the profit accumulated so far, and moves the
 * position in the efficiency order from which the next compute() resumes
 * the search of the critical item, instead of restarting it from the first
 * item. Items with negative weight and profit are complemented, and the
 * data are normalized to a maximization [see load()]. */

class IncrementalGreedyChangeBinaryKnapsackSolver
 : public virtual ChangeSolver , public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor

 IncrementalGreedyChangeBinaryKnapsackSolver()
  : ChangeSolver() , f_N( 0 ) , f_C( 0 ) , f_sense( true ) , f_ci( 0 ) ,
    f_ciVal( 0 ) , obj( - Inf< double >() ) , changedData( true ) ,
    startingSearchIndex( 0 ) , presolveCapacity( 0 ) , P( 0 ) , C( 0 ) ,
    reachTheEnd( false ) {}

 /// destructor

 ~IncrementalGreedyChangeBinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the continuous relaxation, resuming from the last Changes

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// return a valid lower bound on the optimal objective function value
 /** For a minimization problem the relaxation optimum is a lower bound; for
  * a maximization one, the value of the greedy solution without the
  * critical item, which is feasible, or the relaxation optimum itself if
  * the critical item is a continuous variable (the greedy solution is then
  * feasible as it is). */

 OFValue get_lb( void ) override {
  if( ! f_sense )
   return( - obj );
  if( ( f_ciVal == 1 ) || ( ( f_ci < f_N ) && ( ! v_I[ f_ci ] ) ) )
   return( obj );
  return( obj - f_ciVal * v_P[ f_ci ] );
  }

/*--------------------------------------------------------------------------*/
 /// return a valid upper bound on the optimal objective function value
 /** Symmetric to get_lb(): the relaxation optimum for a maximization
  * problem, the value of the greedy solution without the critical item for
  * a minimization one. */

 OFValue get_ub( void ) override {
  if( f_sense )
   return( obj );
  if( ( f_ciVal == 1 ) || ( ( f_ci < f_N ) && ( ! v_I[ f_ci ] ) ) )
   return( - obj );
  return( - obj + f_ciVal * v_P[ f_ci ] );
  }

/*--------------------------------------------------------------------------*/
 /// return the value of the (current) solution
 /** Returns the value of the current solution, with the sign of the sense
  * of the problem (f_sense). */

 OFValue get_var_value( void ) override { return( f_sense ? obj : - obj ); }

/*--------------------------------------------------------------------------*/
 /// tells whether a solution of the relaxation is available

 bool has_var_solution( void ) override { return( f_ci <= f_N ); }

/*--------------------------------------------------------------------------*/
 /// tells whether the solution of the relaxation is feasible for it

 bool is_var_feasible( void ) override { return( has_var_solution() ); }

/*--------------------------------------------------------------------------*/
 /// write the current solution, the critical item rounded, in the Block

 void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// physically construct the current Solution, the critical item rounded

 Solution * get_Solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// read the solution currently stored in the BinaryKnapsackBlock

 void set_var_blockSolution( void ) {
  auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
  if( ! BKB )
   throw( std::invalid_argument( "IncrementalGreedyChangeBinaryKnapsack"
				 "Solver::set_var_blockSolution: the Block "
				 "is not a BinaryKnapsackBlock" ) );
  BKB->get_x( v_x.begin() );
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// add a Modification to the list of those to be processed
 /** Reacts to a NBModification by reloading the instance and clearing the
  * list of the Modification; any other one is stored, to be processed by
  * the next compute(). */

 void add_Modification( sp_Mod & mod ) override;

/*--------------------------------------------------------------------------*/
 /// apply the Change; (un)fixing Changes are applied internally
 /** The (un)fixing Changes [see
  * IncrementalGreedyRelaxationBinaryKnapsackSolver::branch()] update the
  * state of the greedy fill (residual capacity, accumulated profit and
  * position in the efficiency order) without reaching the Block, and so
  * does the returned undo Change, which for a branching Change carries the
  * state of the parent in its last three data (the position, the residual
  * capacity and the profit), so that undoing a child brings the parent
  * back; any other Change is forwarded to the Block. */

 Change * apply( Change * chg , bool doUndo = false ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// write in v_x the greedy solution, the critical item at f_ciVal

 void update_v_x( void ) {
  double xval = 1;
  for( Index i : sortedVar ) {
   if( skip[ i ] ) {
    if( i == f_ci )
     xval = 0;
    continue;
    }
   if( i == f_ci ) {
    v_x[ i ] = complemented[ i ] ? 1 - f_ciVal : f_ciVal;
    xval = 0;
    }
   else
    v_x[ i ] = complemented[ i ] ? 1 - xval : xval;
   }
  }

/*--------------------------------------------------------------------------*/
 /// the greedy solution with the critical item rounded away, if any
 /** Returns v_x [see update_v_x()] with the critical item, if the greedy
  * fill has one and it is an integer variable, rounded to 0 in the space of
  * the complemented items, i.e., to 1 if it is complemented: the solution
  * of the original problem. */

 std::vector< double > rounded_x( void ) {
  update_v_x();
  std::vector< double > sol( v_x );
  if( ( ! reachTheEnd ) && ( f_ci < f_N ) && v_I[ f_ci ] )
   sol[ f_ci ] = complemented[ f_ci ] ? 1 : 0;
  return( sol );
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 // the data of the Binary Knapsack instance

 Index f_N;                        ///< the number of items
 double f_C;                       ///< the capacity of the knapsack
 std::vector< double > v_W;        ///< the weights
 std::vector< double > v_P;        ///< the profits
 std::vector< bool > v_I;          ///< whether each item is integer

 std::vector< unsigned char > v_fxd;  ///< how the x are fixed
 /**< v_fxd[ i ] says whether x_i is fixed, with the encoding 0 = not fixed,
  * 1 = fixed to 0, 2 = fixed to 1 */

 bool f_sense;                     ///< the sense of the objective

 Index f_ci;                       ///< index of the critical item
 double f_ciVal;                   ///< value of the critical item
 double obj;                       ///< the value of the objective
 std::vector< double > v_x;        ///< the solution

 /// the items sorted by nonincreasing efficiency (profit / weight)
 std::vector< Index > sortedVar;

 /// the position of each item in sortedVar
 std::vector< Index > varToSorted;

 /// the items that compute() skips (fixed, or settled by their signs)
 std::vector< bool > skip;

 /// whether the data have changed (Modification, or first load)
 bool changedData;

 /// whether the item is complemented, i.e., its profit and weight have
 /// changed sign
 std::vector< bool > complemented;

 /// the position in sortedVar from which compute() resumes the search of
 /// the critical item (after a branching, say)
 Index startingSearchIndex;

 /// the capacity left by the fixed items alone, for a quick feasibility
 /// check
 double presolveCapacity;

 /// the profit of the greedy solution without the critical item
 double P;

 /// the residual capacity of the greedy solution without the critical item
 double C;

 /// whether the greedy fill reached the last item, i.e., there is no
 /// critical item
 bool reachTheEnd;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// load the Binary Knapsack instance

 void load( void );

/*--------------------------------------------------------------------------*/
 /// process all the pending Modification

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
 /// initialize the state of the greedy fill after (re)loading the instance

 void initializeVariables( void );

/*--------------------------------------------------------------------------*/
 /// unfix item i, which was fixed to value
 /** Updates the capacity left by the fixed items and moves the position
  * from which compute() resumes back to the item, if it is before. */

 void unfix_item( Index i , double value );

/*--------------------------------------------------------------------------*/
 /// fix item i to value (0 or 1)
 /** Updates the solution, the capacity left by the fixed items, the profit
  * and residual capacity of the greedy fill and the position from which
  * compute() resumes, as the fixing requires. */

 void fix_item( Index i , double value );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( IncrementalGreedyChangeBinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/
/*---------- CLASS IncrementalGreedyRelaxationBinaryKnapsackSolver ---------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// RelaxationSolver of the continuous relaxation of a BinaryKnapsackBlock
/** Extends IncrementalGreedyChangeBinaryKnapsackSolver with the
 * RelaxationSolver interface: the true bounds and solution of the original
 * integer problem, and the branching on the critical item. */

class IncrementalGreedyRelaxationBinaryKnapsackSolver
 : public IncrementalGreedyChangeBinaryKnapsackSolver ,
   public virtual RelaxationSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor

 IncrementalGreedyRelaxationBinaryKnapsackSolver()
  : IncrementalGreedyChangeBinaryKnapsackSolver() {}

 /// destructor

 ~IncrementalGreedyRelaxationBinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// branch on the critical item
 /** The children fix the critical item to 1 and to 0 (only to 0 if it does
  * not fit the capacity left by the fixed items). Besides the value, the
  * data of each child carry, in its last two positions, the residual
  * capacity and the profit of the greedy solution of the child, so that
  * apply() picks up the state of the parent without recomputing it. */

 std::vector< Change * > branch( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// return a valid lower bound on the optimal value of the true problem
 /** As get_lb(). */

 OFValue get_true_lb( void ) override { return( get_lb() ); }

/*--------------------------------------------------------------------------*/
 /// return a valid upper bound on the optimal value of the true problem
 /** As get_ub(). */

 OFValue get_true_ub( void ) override { return( get_ub() ); }

/*--------------------------------------------------------------------------*/
 /// tells whether a true solution (a solution of the true original problem
 /// and not of the relaxed one solved by RelaxationSolver) is available
 /** Called after compute() this method has to return true if a true
  * solution of the original problem (not the relaxed one solved by
  * RelaxationSolver) is available to be read with get_true_var_solution().
  * The greedy solution with the critical item rounded away always is. */

 bool has_true_var_solution( void ) override { return( f_ci <= f_N ); }

/*--------------------------------------------------------------------------*/
 /// write the current true solution in the variables of the Block
 /** The true solution is the greedy one with the critical item rounded
  * away [see rounded_x()]. */

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  auto x = rounded_x();
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
  BKB->lock( this );
  BKB->set_x( x.begin() );
  BKB->unlock( this );
  }

/*--------------------------------------------------------------------------*/
 /// after has_true_var_solution(), whether another true solution exists

 bool new_true_var_solution( void ) override { return( false ); }

/*--------------------------------------------------------------------------*/
 /// physically construct the true Solution, the critical item rounded

 Solution * get_true_solution( Configuration * solc = nullptr ) override {
  return( new BinaryKnapsackSolution( rounded_x() ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( IncrementalGreedyRelaxationBinaryKnapsackSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* IncrementalGreedyRelaxationBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*---- End File IncrementalGreedyRelaxationBinaryKnapsackSolver.h ----------*/
/*--------------------------------------------------------------------------*/
