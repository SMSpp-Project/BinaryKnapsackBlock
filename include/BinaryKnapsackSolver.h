/*--------------------------------------------------------------------------*/
/*---------------------- File BinaryKnapsackSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *abstract* class BinaryKnapsackSolver, the common base
 * of the Solver for the (Mixed) Binary Knapsack problem encoded by a
 * BinaryKnapsackBlock. It factors out the two representations every such
 * Solver needs:
 *
 * - the *raw mirror* of the instance (profits, weights, integrality and
 *   fixing status in the original item order and signs), kept incrementally
 *   up to date by processing the BinaryKnapsackBlockMod queue, so that a
 *   re-solve after a Modification never needs to lock and re-traverse the
 *   Block;
 *
 * - the *normalized core*, i.e., the maximisation-sense instance over the
 *   free items only with positive weights and profits: fixed items and
 *   sign-trivial ones are folded into a base profit and a residual capacity,
 *   and an item with negative weight AND profit is complemented
 *   (\f$ x_i = 1 - y_i \f$ with \f$ ( w_i , p_i ) \to ( -w_i , -p_i ) \f$,
 *   the transform being exact on both profit and capacity), split between
 *   the integer and the continuous items.
 *
 * On the normalized core the base also offers the continuous (Dantzig)
 * relaxation
 * \f[
 *  \max \Big\{ \textstyle\sum_i p_i x_i :
 *              \sum_i w_i x_i \le C \,,\; x \in [ 0 , 1 ]^n \Big\} ,
 * \f]
 * solved by the greedy fractional fill in non-increasing efficiency order
 * \f$ p_1 / w_1 \ge p_2 / w_2 \ge \dots \f$, whose optimum takes the items
 * before the *critical* one entirely and the critical one fractionally:
 * this is at once the whole solve of GreedyRelaxationBinaryKnapsackSolver
 * and a valid upper bound for any exact derived Solver.
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

#ifndef __BinaryKnapsackSolver
 #define __BinaryKnapsackSolver
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
/*----------------------- CLASS BinaryKnapsackSolver -----------------------*/
/*--------------------------------------------------------------------------*/
/// abstract base class of the Solver for a BinaryKnapsackBlock
/** Keeps a raw mirror of the BinaryKnapsackBlock instance incrementally up
 * to date under Modification, derives from it the normalized (positive,
 * maximisation, free-items-only) core, and solves the continuous (Dantzig)
 * relaxation of the latter. Derived classes implement compute() on top of
 * these; see the file-level comment. The inheritance from Solver is virtual
 * so that a derived Solver can also implement other Solver-extending
 * interfaces (e.g., the ChangeSolver / RelaxationSolver concepts of the
 * BranchAndBoundSolver module) without duplicating the Solver sub-object. */

class BinaryKnapsackSolver : public virtual Solver {

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

 BinaryKnapsackSolver() : Solver() , f_sense( true ) , f_N( 0 ) , f_Cap( 0 ) ,
                          f_changed( true ) , f_ord_valid( false ) ,
                          f_norm_valid( false ) , f_Cd( 0 ) ,
                          f_base( 0 ) {}

 ~BinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// set the (pointer to the) Block: loads the raw mirror of the instance

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// react to a Modification
 /** An NBModification (the Block has been entirely re-loaded) triggers an
  * immediate reload of the raw mirror and clears the queue; any other
  * Modification is queued, to be applied to the mirror by update_instance()
  * at the next compute(). */

 void add_Modification( sp_Mod & mod ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED TYPES -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// the critical item of the continuous (Dantzig) relaxation

 struct FracInfo {
  int    orig;     ///< original index of the critical item, -1 if none
  double xval;     ///< its value in the original space
  bool   cont;     ///< whether it is a continuous item
  bool   comp;     ///< whether it is a complemented item
  double cp;       ///< its core profit
  double cfrac;    ///< its core fraction
  };

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// apply the queued Modification to the raw mirror
 /** Drains the Modification queue applying each change to the raw mirror in
  * place (no lock on the Block is ever taken). Returns true if anything in
  * the instance changed since the last call - including a reload due to an
  * NBModification or set_Block() - in which case the caller must rebuild
  * whatever it derives from the mirror (normalize_instance() first). */

 bool update_instance( void );

 /// tells whether changes not yet seen by update_instance() are pending

 bool pending_changes( void );

 /// build the normalized core from the raw mirror
 /** Fills f_base, f_Cd, the integer ( c_* ) and continuous ( cc_* ) cores
  * and the fixed / pre-fixed entries of f_x: profits are normalised to the
  * maximisation sense, fixed items and sign-trivial ones (\f$ w_i \ge 0 ,
  * p_i \le 0 \f$ never profitable, \f$ w_i \le 0 , p_i \ge 0 \f$ always
  * taken) are folded into f_base / f_Cd, and an item with \f$ w_i < 0 ,
  * p_i < 0 \f$ is complemented: tentatively taken (folded) and stored with
  * flipped signs and c_comp = 1, so that \f$ x_i = 1 - y_i \f$ with
  * \f$ y_i \f$ the core variable. All core weights and profits are thus
  * positive doubles (no integrality is required here); a f_Cd < 0 on return
  * means the instance is infeasible. */

 void normalize_instance( void );

/*--------------------------------------------------------------------------*/
 /// refresh the normalization after (un)fixing changes only
 /** Cheaper substitute of normalize_instance() valid while #f_norm_valid is
  * true, i.e., when profits, weights and objective sense are unchanged
  * since the last normalize_instance(): re-computes f_x, f_base, f_Cd and
  * the relaxation membership (n_in) in a single pass, leaving the cores of
  * the DP solvers untouched (stale: any solver using them must call
  * normalize_instance() itself, as they all do). */

 void refresh_fixings( void );

/*--------------------------------------------------------------------------*/
 /// apply an (un)fixing Change to the internal mirror of the instance
 /** Applies an (un)fixing BinaryKnapsackBlockChange (its type() must be
  * eFixX or eUnfixX) to the internal mirror of the instance, NOT to the
  * Block: the next compute() of a derived Solver then works under the new
  * fixings without any access to the Block. This is the workhorse of the
  * ChangeSolver concept of the BranchAndXSolver module: the acted-upon
  * items are read via the uniform polymorphic view of the Change [see
  * BinaryKnapsackBlockChange::num_items() / item()], so the concrete
  * (ranged / subset) type of the Change is never needed.
  *
  * With \p doUndo true the returned Change un-does \p chg under a LIFO
  * (climb-the-enumeration-tree) discipline: the undo of fixing is unfixing
  * the same items, and vice-versa with the fixed values captured at call
  * time. */

 Change * apply_to_mirror( const BinaryKnapsackBlockChange * chg ,
                           bool doUndo = false );

 /// solve the continuous (Dantzig) relaxation of the normalized core
 /** Greedy fractional fill of f_Cd over the union of the integer and the
  * continuous core in non-increasing efficiency order: items before the
  * critical one are taken entirely, the critical one fractionally. Writes
  * the complete solution in f_x (in the original space, through the
  * orig / comp maps) and the critical item in @p fi ; returns the optimal
  * value of the relaxation, f_base included, in the maximisation sense.
  * Must be called after normalize_instance(), with f_Cd >= 0. */

 double fractional_relaxation( FracInfo & fi );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /* raw mirror of the instance (original item order, raw signs and sense) - */

 bool   f_sense;                      ///< true = the objective is Max
 Index  f_N;                          ///< number of items
 double f_Cap;                        ///< the Capacity of the Knapsack
 std::vector< double > v_P;           ///< raw profits
 std::vector< double > v_W;           ///< raw weights
 std::vector< char > v_I;             ///< 1 if the variable is integer
 std::vector< unsigned char > v_fxd;  ///< 0 = free, 1 = fixed to 0,
                                      ///< 2 = fixed to 1

 /* normalized core (free items only, positive data, maximisation sense) - -*/

 double f_Cd;                         ///< residual capacity
 double f_base;                       ///< profit of the folded items

 std::vector< double > c_w;           ///< weights of the integer core
 std::vector< double > c_p;           ///< profits of the integer core
 std::vector< Index >  c_orig;        ///< original index of each item
 std::vector< char >   c_comp;        ///< 1 if the item is complemented

 std::vector< double > cc_w;          ///< weights of the continuous core
 std::vector< double > cc_p;          ///< profits of the continuous core
 std::vector< Index >  cc_orig;       ///< original index of each item
 std::vector< char >   cc_comp;       ///< 1 if the item is complemented

 /* solution - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 // cached efficiency order for fractional_relaxation() - - - - - - - - - -
 // the order of the items by efficiency only depends on profits, weights
 // and objective sense, NOT on the fixings: it is computed lazily and kept
 // across (un)fixing changes, which is what makes re-solving the relaxation
 // along an enumeration tree O( n ) per node rather than O( n log n )

 std::vector< int >    v_ord;     ///< ALL the items by efficiency
 bool f_ord_valid;                ///< if v_ord matches profits / weights
 bool f_norm_valid;               ///< if the normalization matches the
                                  ///< current profits / weights / sense
 std::vector< double > n_w;       ///< normalized (complemented) weights
 std::vector< double > n_p;       ///< normalized (complemented) profits
 std::vector< char >   n_in;      ///< 1 if the item enters the relaxation
 std::vector< char >   n_comp;    ///< 1 if the item is complemented

 std::vector< double > f_x;           ///< full solution over ALL items:
                                      ///< fixed / pre-fixed entries set by
                                      ///< normalize_instance(), core ones by
                                      ///< the derived Solver

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// (re)load the raw mirror from the BinaryKnapsackBlock

 void load( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 bool f_changed;   ///< the mirror changed since the last update_instance()

/*--------------------------------------------------------------------------*/

 };  // end( class( BinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File BinaryKnapsackSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
