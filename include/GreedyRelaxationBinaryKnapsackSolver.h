/*--------------------------------------------------------------------------*/
/*-------------- File GreedyRelaxationBinaryKnapsackSolver.h ---------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class
 * GreedyRelaxationBinaryKnapsackSolver, which implements the Solver concept
 * [see Solver.h] for solving the continuous relaxation of a Knapsack problem
 * represented by a BinaryKnapsackBlock, either from scratch at each solve or
 * re-optimizing it along the Changes of a Branch-and-Bound [see
 * intIncremental].
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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Federica Di Pasquale, Filippo Magi,
 *                   Donato Meoli
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

#include "BranchingRule.h"

#include "ChangeSolver.h"

#include "Change.h"

#include "GlobalInformation.h"

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
/** The continuous relaxation is solved by the greedy fractional fill along
 * the items sorted by nonincreasing efficiency, which also identifies the
 * *critical* (fractional) item, the natural branching candidate [see
 * branch()]. Besides the relaxation value, valid true lower and upper bounds
 * on the original mixed-binary problem are available [see get_true_lb() and
 * get_true_ub()], as well as the greedy solution with the critical item
 * rounded away, which is feasible for it.
 *
 * The (un)fixing Changes, as those produced by branch(), are applied to the
 * Solver and not to the Block [see apply()], so that the Solver can serve as
 * the relaxation of the nodes of a Branch-and-Bound (as in BranchAndXSolver).
 * Two ways of doing so are available, chosen by the parameter
 * intIncremental:
 *
 * - 0 (the default): each compute() normalizes the instance (unless only the
 *   fixings have changed) and fills the whole knapsack again, as described
 *   in BinaryKnapsackSolver; branch() also pegs, by reduced cost, the items
 *   that cannot improve on the incumbent [see branch()];
 *
 * - 1: the greedy fill is kept across the Changes: a fixing updates the
 *   residual capacity and the profit accumulated so far, and moves the
 *   position in the efficiency order from which the next compute() resumes
 *   the search of the critical item, instead of restarting it from the
 *   first item; branch() passes to each child the state of its greedy fill,
 *   and the undo of a branching brings back the state of the parent. Any
 *   Modification of the Block, instead, makes the next compute() rebuild
 *   the greedy fill from the data.
 *
 * As a BranchBounder [see BranchingRule.h], with intIncremental == 1 it may
 * also give a valid bound of each child of a branch() without solving it
 * [see intBranchBound]. */

class GreedyRelaxationBinaryKnapsackSolver : public BinaryKnapsackSolver ,
                                             public RelaxationSolver ,
                                             public BranchBounder {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// public enum for the int algorithmic parameters

 enum int_par_type_GRBKSlv {
  intIncremental = intLastAlgPar ,  ///< re-optimize along the Changes
  intBranchBound ,                  ///< the bounds of the children
  intLastGRBKSlvPar                 ///< first allowed new int par for derived
  };

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 using BinaryKnapsackSolver::set_par;  // keep the other set_par() visible

 /// set the int parameters specific of GreedyRelaxationBinaryKnapsackSolver
 /** Set the int parameters specific of GreedyRelaxationBinaryKnapsackSolver:
  *
  * - intIncremental [0]: 0 = each compute() solves the relaxation from
  *   scratch, 1 = the greedy fill is kept across the (un)fixing Changes
  *   applied to the Solver [see the general notes of the class]. Changing
  *   it makes the next compute() start from scratch.
  *
  * - intBranchBound [0]: with 1 and intIncremental == 1, branch_bounds()
  *   gives the bounds of the children of each branch() of Martello and Toth
  *   (1977), which the greedy fill gives at the price of looking for the
  *   item before and after the critical one c: in the normalized
  *   maximization sense, with \f$ P \f$ and \f$ C \f$ the profit and the
  *   residual capacity of the fill without c, \f$ e^- \f$ the efficiency
  *   of the last item taken before c and \f$ e^+ \f$ that of the first
  *   free one after it, \f$ P + p_c - ( w_c - C ) e^- \f$ for the child
  *   that takes c (empty if no item is taken before c) and \f$ P + C e^+
  *   \f$ for the one that leaves it (\f$ P \f$ if no item comes after c),
  *   both no less than the continuous relaxation of the child; with 0, or
  *   with intIncremental == 0, branch_bounds() gives none. */

 void set_par( idx_type par , int value ) override {
  if( par == intIncremental ) {
   if( value != f_incremental ) {
    f_incremental = value;
    f_inc_valid = false;
    }
   }
  else
   if( par == intBranchBound ) {
    if( ( value < 0 ) || ( value > 1 ) )
     throw( std::invalid_argument(
		    "GreedyRelaxationBinaryKnapsackSolver::set_par: "
		    "intBranchBound is neither 0 nor 1" ) );
    f_branch_bound = value;
    }
   else
    BinaryKnapsackSolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastGRBKSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == intIncremental )
   return( f_incremental );
  if( par == intBranchBound )
   return( f_branch_bound );
  return( BinaryKnapsackSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( ( par == intIncremental ) || ( par == intBranchBound ) )
   return( 0 );
  return( BinaryKnapsackSolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intIncremental" )
   return( intIncremental );
  if( name == "intBranchBound" )
   return( intBranchBound );
  return( BinaryKnapsackSolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string name = "intIncremental";
  static const std::string bound = "intBranchBound";
  if( idx == intIncremental )
   return( name );
  if( idx == intBranchBound )
   return( bound );
  return( BinaryKnapsackSolver::int_par_idx2str( idx ) );
  }

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
 /** For a minimization problem the relaxation optimum itself is a valid
  * lower bound. For a maximization problem the bound is the value of a
  * feasible solution of the original problem: the greedy one without the
  * fractional part \f$ \phi \, p_c \f$ of the critical item \f$ c \f$, or
  * the greedy solution itself when there is no critical item or it is a
  * continuous variable (the solution is then feasible as it is). */

 OFValue get_true_lb( void ) override {
  return( f_sense ? true_feasible_value() : - obj );
  }

/*--------------------------------------------------------------------------*/
 /// return a valid upper bound on the optimal objective function value
 /** Symmetric to get_true_lb(): the relaxation optimum for a maximization
  * problem, (minus) the value of the feasible greedy solution without the
  * critical fraction for a minimization one. */

 OFValue get_true_ub( void ) override {
  return( f_sense ? obj : - true_feasible_value() );
  }

/*--------------------------------------------------------------------------*/
 /// tells whether a true solution (a solution of the true original problem
 /// and not of the relaxed one solved by this Solver) is available
 /** Called after compute() this method has to return true if a true solution
  * of the original problem (not the relaxed one solved by this Solver)
  * is available to be read with get_true_var_solution(). The greedy solution
  * with the critical item rounded away always is. */

 bool has_true_var_solution( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
 /// write the current true (rounded greedy) solution in the Block

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
  auto sol = rounded_x();
  BKB->set_x( sol.begin() );
  }

/*--------------------------------------------------------------------------*/
 /// physically construct the Solution of the relaxation
 /** The Solution holds the (possibly fractional) greedy solution, built out
  * of the data structures of the Solver without going through the Block. */

 Solution * get_Solution( Configuration * solc = nullptr ) override {
  sync_x();
  return( new BinaryKnapsackSolution( std::vector< double >( f_x ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// valid lower bound on the optimal value of the true problem
 /** As get_true_lb(). */

 OFValue get_lb( void ) override { return( get_true_lb() ); }

 /// valid upper bound on the optimal value of the true problem
 /** As get_true_ub(). */

 OFValue get_ub( void ) override { return( get_true_ub() ); }

/*--------------------------------------------------------------------------*/
 /// after has_true_var_solution(), whether another true solution exists

 bool new_true_var_solution( void ) override { return( false ); }

/*--------------------------------------------------------------------------*/
 /// write the current solution in the variables of the BinaryKnapsackBlock

 void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// return the value of the (current) solution
 /** Return the value of the current solution, with the sign of the sense
  * of the problem (f_sense). */

 OFValue get_var_value() override { return( f_sense ? obj : - obj ); }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /// branch on the critical item
 /** The two children fix the critical item to 0 and to 1.
  *
  * With intIncremental == 0, each child also carries, when an incumbent is
  * available and local tightenings are sound [see
  * GlobalInformation::str_Incumbent and
  * GlobalInformation::str_LocalFixingAllowed], the reduced-cost
  * (Martello-Toth) pegging valid in its subtree: with \f$ \rho \f$ the
  * efficiency of the critical item and \f$ U \f$ the relaxation value,
  * flipping a free integer item \f$ j \f$ away from its greedy value
  * \f$ \bar{x}_j \f$ cannot yield more than \f$ U - | p_j - \rho w_j | \f$
  * (in the normalized maximization sense), so whenever this does not beat
  * the incumbent \f$ x_j \f$ is fixed at \f$ \bar{x}_j \f$ throughout the
  * subtree; these fixings travel with the fixing of the critical item in a
  * GroupChange.
  *
  * With intIncremental == 1, the data of each child carry, in their last
  * four positions, the state of the greedy fill of the child (the position
  * from which compute() resumes, the residual capacity and the profit) and
  * the epoch of the data it was computed on, so that apply() sets it
  * without recomputing it, whatever Change came before (say, a parallel
  * exploration re-applying the path to a node), as long as the data have
  * not changed since; the value of the fixing is that of the item in the
  * Block, as for any Change, whether or not the fill complements the item;
  * the child that takes the critical item is not produced if the item does
  * not fit the capacity left by the fixed items, unless local reductions
  * are not allowed [see GlobalInformation::str_LocalFixingAllowed], since
  * the child may be non-empty under other data. */

 std::vector< Change * > branch() override;

/*--------------------------------------------------------------------------*/
 /// the bounds of the children of the last branch() [see intBranchBound],
 /// in the sense of the Block, none if intBranchBound == 0

 std::vector< double > branch_bounds( void ) override { return( f_bounds ); }

/*--------------------------------------------------------------------------*/
 /// classify the effect of a Modification on an enumeration tree
 /** The knapsack-specific classification below is not active: classify()
  * is not redefined, and the default of ChangeSolver [see
  * ChangeSolver::classify()] applies. Once enabled, it would classify as
  * follows [see RelaxationSolver::classify()]: profit changes only touch
  * the objective, weight and capacity changes only touch the feasible
  * region, while (un)fixings, integrality and sense changes reshape the
  * branching space (they invalidate the structure of the enumeration tree,
  * not merely its bounds) and so, together with anything structural
  * (NBModification, GroupModification, that update_instance() would not
  * unpack either), conservatively invalidate everything. Non-physical
  * Modification (the abstract-representation mirror of the physical ones)
  * are classified eModNothing, consistently with update_instance() ignoring
  * them: this Solver only listens to the physical representation. */

 /*
 [[nodiscard]] int classify( const sp_Mod & mod ) override {
  const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod.get() );
  if( ! tmod ) {
   if( std::dynamic_pointer_cast< NBModification >( mod ) ||
       std::dynamic_pointer_cast< GroupModification >( mod ) )
    return( eModEverything );
   return( eModNothing );    // abstract-representation mirror: ignored
   }
  switch( tmod->type() ) {
   case( BinaryKnapsackBlockMod::eChgProfit ):
    return( eModObjective );
   case( BinaryKnapsackBlockMod::eChgWeight ):
   case( BinaryKnapsackBlockMod::eChgCapacity ):
    return( eModFeasibility );
   default:    // incl. (un)fixings, integrality and sense changes
    return( eModEverything );
   }
  }
 */

/*--------------------------------------------------------------------------*/
 /// physically construct the true (rounded greedy) Solution

 Solution * get_true_solution( Configuration * solc = nullptr ) override {
  return( new BinaryKnapsackSolution( rounded_x() ) );
  }

/*--------------------------------------------------------------------------*/
 /// apply the Change; (un)fixing Changes are applied to the Solver
 /** The (un)fixing BinaryKnapsackBlockChange, as those produced by branch(),
  * are applied to the Solver and not to the Block: with intIncremental == 0
  * to the internal mirror of the instance [see
  * BinaryKnapsackSolver::apply_to_mirror()], with intIncremental == 1 to
  * the state of the greedy fill, and the next compute() then solves the
  * relaxation under the new fixings. With intIncremental == 1 the undo of
  * a fixing Change carries the state of the parent in its last four
  * data (the position from which compute() resumes, the residual capacity,
  * the profit and the epoch of the data), so that undoing a child brings
  * the parent back, as the branching Change itself carries that of the
  * child [see branch()]; a state of data that have changed since (say, a
  * Change of a tree retained across a Modification) is not taken, and the
  * fill starts again from the first item, out of the fixed items alone, as
  * after any Change that carries no state. A
  * GroupChange is applied by applying its sub-Change, and its undo is the
  * GroupChange of their undo in reverse order; any other Change is
  * forwarded to the Block. */

 Change * apply( Change * chg , bool doUndo = false ) override;

/*--------------------------------------------------------------------------*/
 /// set the GlobalInformation, forgetting the cells read from the previous
 /** The cells of GlobalInformation::str_LocalFixingAllowed and
  * GlobalInformation::str_Incumbent are looked up at the first branch()
  * that finds them, and read through a reference from then on. */

 void set_global_information( GlobalInformation * gi ) override {
  RelaxationSolver::set_global_information( gi );
  f_lfa_cells.reset();
  f_lfa_cell = nullptr;
  f_inc_cells.reset();
  f_inc_cell = nullptr;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// the value (in the maximization sense) of the greedy solution made
 /// feasible for the true problem by dropping the critical fraction

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

 std::vector< double > rounded_x( void ) {
  sync_x();
  std::vector< double > sol( f_x );
  if( ( f_fi.orig >= 0 ) && ( ! f_fi.cont ) )
   sol[ f_fi.orig ] = f_fi.comp ? 1 : 0;
  return( sol );
  }

/*--------------------------------------------------------------------------*/
 /// write in f_x the greedy solution of the last compute()
 /** With intIncremental == 1 the greedy fill does not keep the solution
  * item by item, and this writes it in f_x, the critical item at its
  * fractional value; with intIncremental == 0 f_x is written by compute(),
  * and this does nothing. */

 void sync_x( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 Index f_ci;        ///< index (in the Block) of the critical item
 FracInfo f_fi;     ///< the critical item of the last compute()
 double obj;        ///< the value of the relaxation (maximization sense)

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/
/*- - - - - - - - - - - - - with intIncremental == 1  - - - - - - - - - - -*/

 /// build the greedy fill from the raw mirror of the instance
 /** Normalizes the profits to the maximization sense, complements the free
  * items with negative weight and profit, sets the items fixed or settled
  * by their signs, sorts the others by nonincreasing efficiency and resets
  * the greedy fill to its first item. */

 void inc_initialize( void );

/*--------------------------------------------------------------------------*/
 /// solve the relaxation, resuming the greedy fill from where it was

 int inc_compute( void );

/*--------------------------------------------------------------------------*/
 /// apply an (un)fixing Change to the greedy fill [see apply()]

 Change * inc_apply( BinaryKnapsackBlockChange * chg , bool doUndo );

/*--------------------------------------------------------------------------*/
 /// branch on the critical item, passing the state of each child

 std::vector< Change * > inc_branch( void );

/*--------------------------------------------------------------------------*/
 /// unfix item i, which was fixed to x (in the space of the Block)
 /** Updates the capacity left by the fixed items and moves the position
  * from which compute() resumes back to the item, if it is before. */

 void inc_unfix_item( Index i , double x );

/*--------------------------------------------------------------------------*/
 /// whether the local reductions are allowed
 /** The value of GlobalInformation::str_LocalFixingAllowed, true if there is
  * no such cell. */

 bool local_fixing_allowed( void );

/*--------------------------------------------------------------------------*/
 /// the incumbent, if there is one
 /** Writes in \p value the value of GlobalInformation::str_Incumbent, and
  * returns whether there is such a cell and it is finite. */

 bool incumbent( double & value );

/*--------------------------------------------------------------------------*/
 /// fix item i to x (0 or 1, in the space of the Block)
 /** Updates the capacity left by the fixed items, the profit and residual
  * capacity of the greedy fill and the position from which compute()
  * resumes, as the fixing requires. */

 void inc_fix_item( Index i , double x );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int f_incremental = 0;  ///< the parameter intIncremental
 int f_branch_bound = 0; ///< the parameter intBranchBound

 /// the bounds of the children of the last branch() [see intBranchBound]
 std::vector< double > f_bounds;

 // the cells of the GlobalInformation [see set_global_information()], with
 // the Collection that holds each, so that the reference stays valid

 /// the Collection of GlobalInformation::str_LocalFixingAllowed
 std::shared_ptr< Collection< std::atomic< bool > > > f_lfa_cells;

 /// the cell of GlobalInformation::str_LocalFixingAllowed, nullptr if not
 /// found yet
 std::atomic< bool > * f_lfa_cell = nullptr;

 /// the Collection of GlobalInformation::str_Incumbent
 std::shared_ptr< Collection< std::atomic< double > > > f_inc_cells;

 /// the cell of GlobalInformation::str_Incumbent, nullptr if not found yet
 std::atomic< double > * f_inc_cell = nullptr;

 // the greedy fill kept across the Changes, with intIncremental == 1

 bool f_inc_valid = false;      ///< whether the greedy fill matches the data
 bool f_x_valid = false;        ///< whether f_x holds the last solution

 std::vector< double > v_ip;    ///< profits, normalized and complemented
 std::vector< double > v_iw;    ///< weights, complemented
 std::vector< bool > v_icomp;   ///< whether each item is complemented
 std::vector< double > v_ix;    ///< the solution, the fixed items included

 /// the items sorted by nonincreasing efficiency (profit / weight)
 std::vector< Index > v_sorted;

 /// the position of each item in v_sorted
 std::vector< Index > v_sorted_pos;

 /// the items that compute() skips (fixed, or settled by their signs)
 std::vector< bool > v_skip;

 /// the items settled by their signs: 1 never taken, 2 always taken, 0
 /// neither
 std::vector< char > v_settled;

 /// the position in v_sorted from which compute() resumes the search of
 /// the critical item (after a branching, say)
 Index f_ssi = 0;

 double f_free_cap = 0;  ///< the capacity left by the fixed items alone
 double f_fix_P = 0;     ///< the profit of the fixed items alone

 /// the epoch of the data of the greedy fill, which inc_initialize()
 /// advances; the state a branching Change carries is only taken if it is
 /// of the current epoch [see apply()]
 std::size_t f_inc_epoch = 0;
 double f_iP = 0;        ///< the profit of the fill without the critical item
 double f_iC = 0;        ///< the residual capacity of the same
 double f_ci_val = 1;    ///< the fraction of the critical item that is taken
 bool f_to_end = false;  ///< whether the fill took every item, i.e., there
                         ///< is no critical item

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
