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
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
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

#include <algorithm>

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
                          f_fix_steps( 0 ) , f_exact( false ) ,
                          f_norm_valid( false ) ,
                          f_relax_data( true ) ,
                          f_Cd( 0 ) , f_base( 0 ) {}

 ~BinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// set the (pointer to the) Block: loads the raw mirror of the instance

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR READING THE SOLUTION -------------------*/
/*--------------------------------------------------------------------------*/

 /// returns the current solution as a BinaryKnapsackSolution
 /** Returns the current solution as a BinaryKnapsackSolution [see
  * BinaryKnapsackBlock.h], built out of the data structures of the Solver
  * rather than by writing it in the Variable of the BinaryKnapsackBlock and
  * having it read back from there: no abstract representation is therefore
  * required to exist, and the BinaryKnapsackBlock is not written into at
  * all, hence it is not lock()-ed and any number of Solver attached to it
  * can produce their own Solution at the same time.
  *
  * The solution is the one over *all* the items, i.e., the very same one
  * that get_var_solution() would write into the Variable; nullptr is
  * returned if none is available. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override {
  if( f_x.empty() )
   return( nullptr );

  return( new BinaryKnapsackSolution( std::vector< double >( f_x ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// the Solution comes from the data of the Solver, not from the Variable

 [[nodiscard]] bool is_get_Solution_physical( void ) const override {
  return( true );
  }

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
 /// the kinds of Modification read: the physical ones
 /** The Solver works on the data of the BinaryKnapsackBlock, which the
  * physical Modification describe, and has no use for the abstract ones
  * [see Solver::concerned_by()]. */

 [[nodiscard]] Modification::ModConcern concerned_by( void ) const override {
  return( Modification::eModPhys );
  }

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
 /// (un)fix an item of the mirror
 /** Sets #v_fxd[ i ] to v (0 free, 1 fixed to 0, 2 fixed to 1). While
  * #f_norm_valid, f_x[ i ], f_base, f_Cd and n_in[ i ] change with it as
  * refresh_fixings() would change them, so that a re-solve under the new
  * fixings needs no pass over all the items; one is made anyway after more
  * than f_N items have been (un)fixed this way [see #f_fix_steps], which
  * bounds the rounding that the updates of f_base and f_Cd pile up, unless
  * these updates are exact [see #f_exact]. Without the data of the
  * relaxation [see #f_relax_data], the normalization is just marked as
  * stale. */

 void set_fixing( Index i , unsigned char v );

/*--------------------------------------------------------------------------*/
 /// (re)build the cached efficiency order of the items
 /** Fills #v_ord with the items sorted by non-increasing efficiency
  * \f$ p_i / w_i \f$ of the normalized instance [see
  * normalize_instance()], which is the order the greedy fill of the
  * continuous relaxation follows. The order only depends on the profits,
  * the weights and the sense of the Objective, NOT on the fixings, hence it
  * survives a re-solve under different fixings and is only rebuilt when
  * #f_ord_valid says so; so is the envelope that max_margin() reads. */

 void build_efficiency_order( void );

/*--------------------------------------------------------------------------*/
 /// the largest margin of an item at the given efficiency
 /** Returns \f$ \max_j | p_j - \rho w_j | \f$ over the items that can
  * enter the relaxation (those of positive normalized weight [see
  * normalize_instance()]), whether or not they are fixed, 0 if there is
  * none, in \f$ O( \log n ) \f$ time: the maximum is the upper envelope of
  * the \f$ 2n \f$ lines \f$ \pm ( p_j - \rho w_j ) \f$, which is built
  * with the efficiency order [see build_efficiency_order()]. */

 double max_margin( double rho ) const;

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
  * the critical item in @p fi and returns the optimal value of the
  * relaxation, f_base included, in the maximisation sense; the solution
  * of the items in the relaxation is written in f_x by fractional_x().
  * The fill resumes from the first position of the efficiency order whose
  * item has entered or left the relaxation since the last call, the
  * critical item being found by a binary search on the part of the fill
  * before it, so that a node of a tree that differs from the last one by
  * a few fixings costs far less than \f$ O( n ) \f$. Must be called
  * after normalize_instance(), with f_Cd >= 0. */

 double fractional_relaxation( FracInfo & fi );

/*--------------------------------------------------------------------------*/
 /// write in f_x the items of the last fractional_relaxation()
 /** fractional_relaxation() leaves the entries of f_x of the items in the
  * relaxation as they were, which this writes, in \f$ O( n ) \f$. */

 void fractional_x( void );

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
 std::vector< double > v_env_m;   ///< slopes, intercepts and breakpoints
 std::vector< double > v_env_a;   ///< of the upper envelope that
 std::vector< double > v_env_x;   ///< max_margin() reads
 std::vector< Index > v_ord_pos;  ///< the position of each item in v_ord

 // the greedy fill of the last fractional_relaxation(), from which the next
 // one resumes: v_fr_w[ k ] and v_fr_p[ k ] are the weight and the profit
 // of the items of the relaxation in the first k positions of v_ord, right
 // for k <= f_fr_valid, which set_fixing() lowers to the position of an
 // item entering or leaving the relaxation

 std::vector< double > v_fr_w;    ///< prefix sums of the weights
 std::vector< double > v_fr_p;    ///< prefix sums of the profits
 Index f_fr_valid = 0;            ///< how many of them are right
 Index f_fr_crit = 0;             ///< the position of the critical item of
                                  ///< the last fill, f_N if none
 double f_fr_y = 1;               ///< and its fraction
 Index f_fix_steps;               ///< items (un)fixed by set_fixing() since
                                  ///< the last pass over all of them
 bool f_exact;                    ///< if profits, weights and capacity are
                                  ///< integer and so small that f_base and
                                  ///< f_Cd are kept with no rounding
 bool f_norm_valid;               ///< if the normalization matches the
                                  ///< current profits / weights / sense
 bool f_relax_data;               ///< if normalize_instance() fills n_w,
                                  ///< n_p, n_comp and n_in (the relaxation
                                  ///< data), false in the Solvers not using
                                  ///< fractional_relaxation()
 std::vector< double > n_w;       ///< normalized (complemented) weights
 std::vector< double > n_p;       ///< normalized (complemented) profits
 std::vector< char >   n_in;      ///< 1 if the item enters the relaxation
 std::vector< char >   n_comp;    ///< 1 if the item is complemented
 std::vector< unsigned char > n_cls;  ///< what a free item is: 0 dropped,
                                      ///< 1 always taken, 2 complemented,
                                      ///< 3 a genuine one

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
 /// the part of update_instance() for a nonempty queue
 /** Called with f_mod_lock held, which it releases. */

 bool drain_Modification( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 bool f_changed;   ///< the mirror changed since the last update_instance()

/*--------------------------------------------------------------------------*/

 };  // end( class( BinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/
/*----------------- INLINE METHODS OF BinaryKnapsackSolver -----------------*/
/*--------------------------------------------------------------------------*/

inline bool BinaryKnapsackSolver::update_instance( void )
{
 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 // nothing queued and nothing reloaded, which is the case of nearly every
 // call [see GreedyRelaxationBinaryKnapsackSolver::apply()]
 if( ( ! f_changed ) && v_mod.empty() ) {
  f_mod_lock.clear( std::memory_order_release );  // release lock
  return( false );
  }

 return( drain_Modification() );

 }  // end( BinaryKnapsackSolver::update_instance )

/*--------------------------------------------------------------------------*/

inline void BinaryKnapsackSolver::set_fixing( Index i , unsigned char v )
{
 const unsigned char o = v_fxd[ i ];
 v_fxd[ i ] = v;
 if( ( o == v ) || ( ! f_norm_valid ) )
  return;

 if( ! f_relax_data ) {
  f_norm_valid = false;
  return;
  }

 if( ! f_exact )
  ++f_fix_steps;

 // under the fixing s, whether the item is (tentatively) taken and whether
 // it is in the relaxation, as normalize_instance() and refresh_fixings()
 // have it: a free item is taken if always taken or complemented
 const unsigned char cls = n_cls[ i ];
 const auto taken = [ cls ]( unsigned char s ) {
  return( s ? ( s == 2 ) : ( ( cls == 1 ) || ( cls == 2 ) ) );
  };

 if( const int d = int( taken( v ) ) - int( taken( o ) ) ) {
  const double p = f_sense ? v_P[ i ] : - v_P[ i ];
  if( d > 0 ) { f_base += p; f_Cd -= v_W[ i ]; }
  else        { f_base -= p; f_Cd += v_W[ i ]; }
  }

 const char in = ( ! v ) && ( cls >= 2 );
 if( ( in != n_in[ i ] ) && f_ord_valid )   // the fill changes from it on
  f_fr_valid = std::min( f_fr_valid , v_ord_pos[ i ] );
 n_in[ i ] = in;
 if( ! in )                 // in the relaxation, f_x[ i ] set by the scan
  f_x[ i ] = taken( v );

 }  // end( BinaryKnapsackSolver::set_fixing )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File BinaryKnapsackSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
