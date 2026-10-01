/*--------------------------------------------------------------------------*/
/*-------------------- File CoreDPBinaryKnapsackSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class CoreDPBinaryKnapsackSolver, an exact
 * Solver for the (pure 0-1) Binary Knapsack encoded by a BinaryKnapsackBlock
 * based on a *core / non-dominated-state* dynamic program rather than the
 * textbook full-table O(n*C) one of DPBinaryKnapsackSolver.
 *
 * RATIONALE. The full-table DP scans the whole n x C table; on hard
 * (correlated, spanner) instances and on large-capacity ones (Pisinger
 * largecoeff R=1e7, Jooken W=1e10) this is hopeless in both time and memory.
 * The state of the art (Pisinger MINKNAP 1997, Martello-Pisinger-Toth COMBO
 * 1999, and the recent RECORD 2026) instead enumerates only the *non-dominated*
 * (profit, weight) states around an efficiency "break item", pruning with LP /
 * dominance / surrogate bounds, and is orders of magnitude faster while keeping
 * a small state set in memory.
 *
 * Why a subclass of BinaryKnapsackSolver and NOT of DPBinaryKnapsackSolver:
 * the latter's dynamic_programming() must fill the full final-slice label
 * profile lab[N][.] and the per-item predecessor graph pred[.][.], on which
 * its continuous-greedy step and solution reconstruction depend.
 * Materialising those is exactly what a core DP avoids; overriding only
 * dynamic_programming() would force rebuilding pred for every item and throw
 * away the memory win. The abstract base instead provides exactly the shared
 * non-algorithmic machinery - the raw mirror of the instance with its
 * incremental update under Modification and the normalized (positive,
 * maximisation) core - on which the core enumeration is built, with ALL the
 * state in the object (no file-scope globals), which also makes the Solver
 * reentrant -- a hard requirement for running many independent knapsack
 * solves concurrently (the coarse-grained parallel avenue via
 * ParallelBundleSolver, the only one that actually pays off).
 *
 * STATUS: implements the break-item core enumeration with dominance (Pareto
 * frontier) and the Dantzig LP bound (per-state pruning + global early-exit
 * ceiling), i.e. the MINKNAP scheme, plus the COMBO / RECORD machinery on top:
 * the Martello-Toth U2 ceiling, the surrogate relaxation with cardinality
 * constraints (fractional bound AND exact recursive solve, which both
 * certifies and raises the incumbent), weight- and profit-divisibility
 * reductions, Dembo-Hammer per-item fixing, fixing-by-dominance, the guarded
 * (read-only-checked) DP extension, item aggregation with multiplicity
 * reduction and binary-decomposed copy expansion, and the initial / pairing
 * heuristic pool (H1-H3, PH, TPH, SSPH, SPH, GCH). Handles the full
 * BinaryKnapsackBlock feature set: both objective senses, fixed variables,
 * negative weights / profits (via complementation), and continuous variables
 * (the integer Pareto frontier is combined with the fractional fill of the
 * continuous part).
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CoreDPBinaryKnapsackSolver
 #define __CoreDPBinaryKnapsackSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------- CLASS CoreDPBinaryKnapsackSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/// core / non-dominated-state Dynamic Programming Solver for the 0-1 Knapsack
/** Exact Solver for the BinaryKnapsackBlock that enumerates only the
 * non-dominated (profit, weight) states of the 0-1 knapsack instead of the
 * full n x C label table. See the file-level comment for the rationale and the
 * relation to DPBinaryKnapsackSolver. */

class CoreDPBinaryKnapsackSolver : public BinaryKnapsackSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// tolerance for the integrality of the weights (must be integers)

 static constexpr double WeightIntegrality = 1e-06;

 /// public enum for the int algorithmic parameters

 enum int_par_type_CDPBKSlv {
  intReopt = intLastAlgPar ,     ///< reoptimization level
  intSurrogate ,                 ///< surrogate relaxation: off / bound / solve
  intSurrTrigger ,               ///< frontier size that triggers it
  intSurrAdapt ,                 ///< adaptive rules of the surrogate
  intDPExtension ,               ///< guarded DP extension on / off
  intDominanceFix ,              ///< fixing by dominance on / off
  intReductionFix ,              ///< Dembo-Hammer reduction fixing on / off
  intPrimalHeur ,                ///< primal heuristics during the enumeration
  intLazyCore ,                  ///< fixed core around the break item
  intLastCDPBKSlvPar             ///< first allowed new int par for derived
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 CoreDPBinaryKnapsackSolver() : BinaryKnapsackSolver() , f_C( 0 ) ,
                                f_obj( - Inf< double >() ) , f_reopt( 0 ) ,
                                f_surrogate( 2 ) , f_surr_trigger( 2000 ) ,
                                f_surr_adapt( 1 ) , f_rec_left( 0 ) ,
                                f_state_cap( -1 ) , f_dp_ext( 1 ) ,
                                f_dom_fix( 2 ) , f_red_fix( 1 ) , f_heur( 1 ) ,
                                f_lazy_core( 1 ) , f_prev_valid( false )
  { f_relax_data = false; }

 ~CoreDPBinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 using BinaryKnapsackSolver::set_par;  // keep the other set_par() visible

 /// set the int parameters of CoreDPBinaryKnapsackSolver
 /** Set the int parameters specific of CoreDPBinaryKnapsackSolver:
  *
  * - intReopt [0]: how much of the previous solve is reused by the next one
  *   after a change of the data. 0 = none, every solve starts from scratch;
  *   1 = the previous optimal solution \f$ \bar{x} \f$, repaired to the new
  *   data (fixings enforced, the items whose sign makes them never or always
  *   profitable set accordingly, and the least efficient taken items dropped
  *   until \f$ \sum_i w_i \bar{x}_i \leq C \f$), seeds the incumbent of the
  *   core enumeration: every state whose bound cannot beat
  *   \f$ \sum_i p_i \bar{x}_i \f$ is pruned from the start, and
  *   \f$ \bar{x} \f$ itself is returned when nothing beats it; 2 = as 1,
  *   and moreover no solve at all when the changes since the last solve
  *   cannot have made any other solution better than \f$ \bar{x} \f$: same
  *   capacity and same core items, of which the taken ones kept their weight
  *   and did not lose profit, and the others did not lose weight nor gain
  *   profit; 3 = as 2, and moreover no solve also when some taken item lost
  *   profit or some other one gained it, if for each such item \f$ k \f$
  *   the Lagrangian bound of the new data (capacity relaxed, multiplier
  *   the efficiency of the break item) with \f$ x_k = 1 - \bar{x}_k \f$
  *   cannot beat \f$ \bar{x} \f$ (any better solution must flip one of
  *   them); after a few failures in a row this test is skipped for a
  *   number of solves that doubles with each further failure. Only the
  *   pure 0-1 case (no continuous variable) is warm started.
  *
  * The following ones switch the components of the core enumeration, all
  * exact (they change the running time, never the optimum), so that the
  * choice can be made per instance:
  *
  * - intSurrogate [2]: the surrogate relaxation with a cardinality
  *   constraint; 0 = never, 1 = as a fractional bound, 2 = solved exactly
  *   by recursion, so that its optimum certifies and possibly raises the
  *   incumbent;
  *
  * - intSurrTrigger [2000]: the size of the Pareto frontier past which the
  *   surrogate relaxation is computed (below it the instance is easy and the
  *   relaxation would only add overhead);
  *
  * - intSurrAdapt [1]: the adaptive rules of the surrogate relaxation, a bit
  *   mask; bit 0 (+1) = the multiplier search starts from the end of its
  *   range where the surrogate reduces to the continuous relaxation, and
  *   stops there when the gradient says that the minimiser is that end (so
  *   that a surrogate adding nothing to the bound costs one or two
  *   evaluations instead of a whole binary search); bit 1 (+2) = the exact
  *   solve of the surrogate (intSurrogate = 2) runs on a budget of twice the
  *   states generated by the enumeration up to that point, falling back to
  *   the fractional surrogate bound when it runs out and retrying once the
  *   enumeration has doubled its states (off by default: a solve that would
  *   certify can need more than the budget, and the retries then cost more
  *   than they save);
  *
  * - intDPExtension [1]: the guarded DP extension, i.e., the read-only check
  *   that skips an item whose expansion could not produce any state;
  *
  * - intDominanceFix [2]: the fixing by dominance of the items that a copy
  *   of an already useless item dominates; 0 = never, 1 = each item is
  *   checked against the last useless item of its side, 2 = every useless
  *   item marks at once all the items still to be processed on its side
  *   that it dominates;
  *
  * - intReductionFix [1]: the Dembo-Hammer fixing of the items whose
  *   reduction bound cannot beat the incumbent;
  *
  * - intPrimalHeur [1]: the primal heuristics run during the enumeration
  *   (pairing, subset sampling and greedy completion);
  *
  * - intLazyCore [1]: the fixed core, on (1) or off (0). When on, the
  *   break item is found by quickselect instead of sorting the instance, and
  *   every item whose Dembo-Hammer bound cannot beat the break solution (or
  *   the given lower bound) is fixed at its break value; the core
  *   enumeration then runs on the remaining items only, or on the whole
  *   instance if they are more than half of it. This spares the sorting of
  *   the whole instance, which dominates the solve of the easy ones. */

 void set_par( idx_type par , int value ) override {
  if( ( par >= intReopt ) && ( par < intLastCDPBKSlvPar ) )
   *int_par_ptr( par ) = value;
  else
   BinaryKnapsackSolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastCDPBKSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( ( par >= intReopt ) && ( par < intLastCDPBKSlvPar ) )
   return( *const_cast< CoreDPBinaryKnapsackSolver * >( this )->int_par_ptr(
                                                                     par ) );
  return( BinaryKnapsackSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  static constexpr int dflt[] = { 0 , 2 , 2000 , 1 , 1 , 2 , 1 , 1 , 1 };
  if( ( par >= intReopt ) && ( par < intLastCDPBKSlvPar ) )
   return( dflt[ par - intReopt ] );
  return( BinaryKnapsackSolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  for( idx_type i = intReopt ; i < intLastCDPBKSlvPar ; ++i )
   if( name == par_names()[ i - intReopt ] )
    return( i );
  return( BinaryKnapsackSolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( ( idx >= intReopt ) && ( idx < intLastCDPBKSlvPar ) )
   return( par_names()[ idx - intReopt ] );
  return( BinaryKnapsackSolver::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the 0-1 knapsack encoded in the BinaryKnapsackBlock

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override { return( get_var_value() ); }

 OFValue get_ub( void ) override { return( get_var_value() ); }

 /// return the value of the (current) solution, in the Block's own sense

 OFValue get_var_value( void ) override { return( f_sense ? f_obj : - f_obj ); }

 /// write the current solution into the variables of the BinaryKnapsackBlock

 void get_var_solution( Configuration * solc = nullptr ) override;

 /// how much of the previous solve the last one reused [see intReopt]
 /** Returns 0 if the last compute() solved the instance from scratch, 1 if
  * it started from the previous solution as the incumbent, 2 if it
  * returned the previous solution with no solve because the changes could
  * not have made it suboptimal, and 3 if it did so because a Lagrangian
  * bound showed it. */

 [[nodiscard]] int get_reopt_outcome( void ) const { return( f_outcome ); }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED TYPES -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// a non-dominated DP state: profit/weight sum plus backtrack info

 struct State {
  long   wsum;       ///< weight sum of the state ( <= residual capacity )
  double psum;       ///< profit sum of the state
  int    parent;     ///< index of the originating state in the previous step
  bool   took;       ///< whether the current item was added to reach this state
  };

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// build the integer-weight cores from the base's normalized instance:
 /// checks the weight integrality, floors the residual capacity into f_C and
 /// pre-fixes the integer items that cannot fit it

 void extract_instance( void );

 /// solve the core and set f_obj / f_x, dispatching to the fast break-item
 /// enumeration (no continuous variables) or to the Pareto-frontier +
 /// fractional-fill scheme (with continuous variables)

 void enumerate_states( void );

 /// break-item core enumeration of the pure 0-1 integer knapsack.
 /** Solves \f$ \max \sum_i p_i x_i \f$ s.t. \f$ \sum_i w_i x_i \le C \f$,
  * \f$ x_i\in\{0,1\} \f$. With the items in non-increasing efficiency order
  * \f$ \tfrac{p_1}{w_1}\ge\tfrac{p_2}{w_2}\ge\dots \f$, the break item
  * \f$ b \f$ is the first with
  * \f$ \sum_{i<b} w_i \le C < \sum_{i\le b} w_i \f$; the break
  * solution takes \f$ x_i=1,\,i<b \f$ (profit \f$ \hat p=\sum_{i<b}p_i \f$,
  * weight \f$ \hat w=\sum_{i<b}w_i \f$, residual \f$ r=C-\hat w \f$). Only the
  * deviations from it are enumerated; a state \f$ s=(s_p,s_w) \f$ is pruned
  * when its LP bound cannot beat the incumbent \f$ z \f$:
  * \f[ B(s)=\begin{cases}
  *   s_p+(C-s_w)\,\dfrac{p_{nr}}{w_{nr}}, & s_w\le C\\[2.2ex]
  *   s_p-(s_w-C)\,\dfrac{p_{nl}}{w_{nl}}, & s_w> C
  *   \end{cases}\qquad\text{prune if } B(s)<z+1, \f]
  * where \f$ nr \f$ / \f$ nl \f$ are the next unprocessed right / left items.
  * Returns the optimal core profit and fills @p in with the chosen items.
  * Virtual so that a derived Solver can bring its own engine for the integer
  * core while keeping all the rest (see RECORDBinaryKnapsackSolver). */

 virtual double solve_integer_core( std::vector< char > & in );

 /// trivial divisibility bound: effective (reduced) capacity.
 /** With \f$ g=\gcd_i w_i \f$ every reachable weight is a multiple of
  * \f$ g \f$, so the usable capacity is
  * \f[ \bar C=\Big\lfloor\tfrac{C}{g}\Big\rfloor\,g\le C. \f]
  * Returned in place of @p C: it is exact for feasibility (all state
  * weights are multiples of \f$ g \f$) and only tightens the bounds; it applies
  * to the pure-integer core, the continuous fill keeps the real \f$ C \f$. */

 long divisibility_capacity( const std::vector< long > & w , int m ,
                             long C ) const;

 /// strong initial heuristic: a good feasible incumbent to seed the pruning.
 /** Returns the best of three O(n) feasible extensions of the break solution
  * \f$ \hat x \f$ (so always \f$ \ge\hat p \f$):
  * \f[ z_{H1}=\hat p+\max\Big\{\textstyle\sum_{i\in F}p_i :
  *     F\subseteq\{b,\dots\}\text{ greedily filled, }
  *     \hat w+\textstyle\sum_{i\in F}w_i\le C\Big\}, \f]
  * \f[ z_{H2}=\hat p-p_{b-1}+\textstyle\sum_{i\in F'}p_i
  *     \quad\text{(drop the least efficient taken item, refill)} , \f]
  * \f[ z_{H3}=\hat p-\textstyle\sum_{i\in D}p_i+p_b
  *     +\textstyle\sum_{i\in F''}p_i
  *     \quad\text{(drop left items until } b \text{ fits, refill)} , \f]
  * with \f$ F',F'' \f$ the greedy refills of the respective residuals. Inputs
  * are in efficiency-sorted core order; @p C is the effective capacity, @p in
  * receives the achieving selection. */

 double initial_heuristic( const std::vector< long > & w ,
                           const std::vector< double > & p , int m , int b ,
                           long wsumb , double psumb , long C ,
                           std::vector< char > & in );

 /// Martello-Toth U2 global upper bound (the break item is wholly in or out).
 /** \f[ U^{2}=\max\Big\{\;
  *   \hat p+\big\lfloor r\,\tfrac{p_{b+1}}{w_{b+1}}\big\rfloor ,
  *   \;\; \hat p+p_b-\big\lceil (w_b-r)\,\tfrac{p_{b-1}}{w_{b-1}}\big\rceil
  *   \;\Big\} , \f]
  * with \f$ \hat p \f$, \f$ r=C-\hat w \f$ as above (the floor/ceil are dropped
  * for non-integer profits). Provably \f$ U^{2}\le U^{D}=\hat p+r\,p_b/w_b \f$
  * (the Dantzig bound), hence always valid. @p C is the effective capacity. */

 double upper_bound_u2( const std::vector< long > & w ,
                        const std::vector< double > & p , int m , int b ,
                        long wsumb , double psumb , long C , bool intp ) const;

 /// surrogate relaxation with a cardinality constraint (COMBO scheme).
 /** Let \f$ N_{\max} \f$ be the largest cardinality that fits \f$ C \f$ (the
  * \f$ N_{\max} \f$ smallest weights sum to \f$ \le C \f$) and \f$ N_{\min} \f$
  * the smallest cardinality whose largest profits exceed \f$ z \f$. When one of
  * them is within 1 of the break cardinality \f$ b \f$ (the cardinality is then
  * essentially forced) the capacity and cardinality constraints are aggregated
  * with a multiplier \f$ \mu \f$ into the single surrogate constraint
  * \f[ \sum_i (w_i+\mu)\,x_i\ \le\ C+\mu\,k, \f]
  * and the surrogate upper bound is \f$ \min_{\mu} U^{SF}(\mu) \f$ with
  * \f$ U^{SF}(\mu) \f$ the fractional-knapsack value of the items
  * \f$ (p_i,\,w_i+\mu) \f$ under capacity \f$ C+\mu k \f$. A dichotomy at
  * \f$ b \f$ (\f$ \sum_i x_i\le b \f$ vs \f$ \ge b+1 \f$) is a complete, always
  * valid case split. Returns \f$ +\infty \f$ (no tightening) for a wide
  * window. */

 double surrogate_bound( const std::vector< long > & w ,
                         const std::vector< double > & p , int m , int b ,
                         long C , double z ) const;

 /// tightest surrogate fractional bound for a fixed cardinality @p card.
 /** \f[ \min_{\,\mathrm{slo}\le\mu\le\mathrm{shi}}\; U^{SF}(\mu),\qquad
  *   U^{SF}(\mu)=\text{fractional knapsack of }(p_i,w_i+\mu)
  *   \text{ under }C+\mu\,\mathrm{card}. \f]
  * Any \f$ \mu \f$ gives a valid bound, so the minimum does too; the minimiser
  * is located by a gradient binary search over the integer multiplier; the
  * minimising \f$ \mu \f$ is returned through @p out_sur when not null. */

 double surrogate_card_bound( const std::vector< long > & w ,
                              const std::vector< double > & p , int m ,
                              int card , long C , long slo , long shi ,
                              long * out_sur = nullptr ) const;

 /// surrogate relaxation with cardinality constraints, solved exactly.
 /** Like surrogate_bound() splits on the forced cardinalities, but each case
  * \f$ k \f$ (with its minimising multiplier \f$ \mu \f$) is improved from the
  * fractional bound to the EXACT optimum of the surrogate subproblem
  * \f[ z^{S}=\max\Big\{\sum_i p_i x_i :\ \sum_i (w_i+\mu)\,x_i \le C+\mu k ,
  *     \ x_i\in\{0,1\}\Big\} \ \ge\ \max\{\,p(x): x\text{ feasible, } |x|
  *     \text{ in the case}\,\} , \f]
  * which being a relaxation is itself a valid - usually certifying - upper
  * bound (the fractional bound is kept when the exact solve is truncated).
  * Moreover, if the optimum \f$ x^\* \f$ has \f$ \sum_i x_i^\*=k \f$ then
  * \f$ \sum_i w_i x_i^\* = \sum_i (w_i+\mu)x_i^\* - \mu k \le C \f$: it is
  * feasible for the original problem and, when \f$ \sum_i p_i x_i^\* > z \f$,
  * raises the incumbent (@p z and @p best_in , sorted-core order). Returns
  * the combined ceiling over the case split, \f$ +\infty \f$ when no
  * cardinality is forced. */

 double surrogate_solve( const std::vector< long > & w ,
                         const std::vector< double > & p , int m , int b ,
                         long C , double & z , std::vector< char > & best_in );

 /// the break-item core enumeration over an explicit instance.
 /** The engine behind solve_integer_core(), reusable on any positive-weight
  * positive-profit 0-1 instance \f$ ( w , p , C ) \f$ — in particular on the
  * surrogate subproblems, which surrogate_solve() solves by recursing here
  * with @p relx set (no surrogate relaxation inside a surrogate solve). The
  * caller's incumbent @p lb floors the bound pruning, making the search only
  * ever work towards a strict improvement: the return is the exact optimum
  * (and @p in the achieving selection, in the items' input order) when it is
  * \f$ >lb \f$, and \f$ -\infty \f$ — a proof that no solution beats @p lb —
  * otherwise. Pass \f$ -\infty \f$ to always obtain the optimum.
  *
  * Identical items are aggregated into one item with a multiplicity
  * \f$ d_i \f$, and an item that exactly doubles another is folded into it
  * (multiplicity reduction):
  * \f[ ( p_{i'} , w_{i'} ) = ( 2 p_i , 2 w_i ) \ \Longrightarrow\
  *     d_i \leftarrow d_i + 2\,d_{i'} , \f]
  * applied transitively in decreasing weight order (chains
  * \f$ 8w \to 4w \to 2w \f$ unwind onto the base item), every solution of the
  * reduced instance mapping back to one of equal value. The \f$ d_i \f$
  * copies are then expanded in binary-decomposed batches
  * \f$ 1 , 2 , 4 , \dots \f$ — so any count in \f$ [ 0 , d_i ] \f$ is
  * reachable in \f$ O( \log d_i ) \f$ merges — and as soon as one batch
  * yields no new non-dominated state the item is closed, since every later
  * translate of the frontier is then a translate of an already-merged one. */

 double core_enumerate( const std::vector< long > & iw ,
                        const std::vector< double > & ip , long C ,
                        std::vector< char > & in , double lb , bool relx );

 /// the core enumeration on the whole of the given instance
 /** The engine behind core_enumerate() (same arguments and result), which
  * calls it either on the whole instance or, with intLazyCore > 0, on a
  * fixed core around the break item (see set_par()). */

 double core_enumerate_full( const std::vector< long > & iw ,
                             const std::vector< double > & ip , long C ,
                             std::vector< char > & in , double lb ,
                             bool relx );

 /// the engine of core_enumerate_full(), with state profits of type PT
 /** core_enumerate_full() calls it with PT = long when all the profits are
  * integers (and their sum is exactly representable as a double), so that
  * the per-state arithmetic of the enumeration is integer, with PT = double
  * otherwise. */

 template< class PT >
 double core_enumerate_engine( const std::vector< long > & iw ,
                               const std::vector< double > & ip , long C ,
                               std::vector< char > & in , double lb ,
                               bool relx );

 /// combine the integer Pareto frontier with the fractional fill of the
 /// continuous variables; returns the optimal core profit, fills @p in (over
 /// the integer core items) and @p cx (over the continuous core items)

 double solve_with_continuous( std::vector< char > & in ,
                               std::vector< double > & cx );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /* the integer-weight view of the base's normalized core - - - - - - - - - */

 long   f_C;                      ///< integer residual capacity = floor(f_Cd)

 /* integer core: free 0-1 items, weights/profits normalised to be positive;
  * "complement" items (originally negative weight and profit) are stored with
  * flipped sign and v_comp = 1, so that x_orig = 1 - ( chosen in the core ) */

 std::vector< long >   v_w;       ///< weights of the integer core items ( > 0 )
 std::vector< double > v_p;       ///< profits of the integer core items ( > 0 )
 std::vector< Index >  v_orig;    ///< original item index of each core item
 std::vector< char >   v_comp;    ///< 1 if the item is complemented

 /* continuous core: free continuous items, same positive normalisation and
  * complement convention as the integer core */

 std::vector< long >   v_cw;      ///< weights of the continuous core items
 std::vector< double > v_cp;      ///< profits of the continuous core items
 std::vector< Index >  v_corig;   ///< original item index of each cont. item
 std::vector< char >   v_ccomp;   ///< 1 if the continuous item is complemented

 /* solution - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 double f_obj;                    ///< optimal value (maximisation sense)

 /* reoptimization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 int f_reopt;                     ///< the intReopt parameter
 int f_surrogate;                 ///< the intSurrogate parameter
 int f_surr_trigger;              ///< the intSurrTrigger parameter
 int f_surr_adapt;                ///< the intSurrAdapt parameter
 long f_rec_left;                 ///< states left to the surrogate recursion
                                  ///< (set before each surrogate solve)
 long f_state_cap;                ///< if >= 0, core_enumerate_full() stops
                                  ///< past this many states and returns the
                                  ///< best solution found (a heuristic)
 int f_dp_ext;                    ///< the intDPExtension parameter
 int f_dom_fix;                   ///< the intDominanceFix parameter
 int f_red_fix;                   ///< the intReductionFix parameter
 int f_heur;                      ///< the intPrimalHeur parameter
 int f_lazy_core;                 ///< the intLazyCore parameter

 /// the field holding the int parameter @p par of this class

 int * int_par_ptr( idx_type par ) {
  int * const f[] = { & f_reopt , & f_surrogate , & f_surr_trigger ,
                      & f_surr_adapt , & f_dp_ext , & f_dom_fix ,
                      & f_red_fix , & f_heur , & f_lazy_core };
  return( f[ par - intReopt ] );
  }

 /// the names of the int parameters of this class, in enum order

 static const std::vector< std::string > & par_names( void ) {
  static const std::vector< std::string > n = { "intReopt" , "intSurrogate" ,
   "intSurrTrigger" , "intSurrAdapt" , "intDPExtension" ,
   "intDominanceFix" , "intReductionFix" , "intPrimalHeur" , "intLazyCore" };
  return( n );
  }
 bool f_prev_valid;               ///< v_prev_x holds the last optimal solution
 std::vector< double > v_prev_x;  ///< last optimal solution (original items)

 /// the previous solution repaired to the current data, as a core incumbent
 /** Repairs v_prev_x to the current mirror (see intReopt) and, if the result
  * is a 0-1 feasible solution, writes in @p in its core image and returns its
  * value in the core (i.e., net of f_base); returns -Inf otherwise. */

 double warm_incumbent( std::vector< char > & in ) const;

 /// the core of the last solve and its optimal solution (intReopt 2)

 long f_last_C = -1;
 std::vector< long > v_last_w;
 std::vector< double > v_last_p;
 std::vector< Index > v_last_orig;
 std::vector< char > v_last_comp;
 std::vector< char > v_last_in;

 /// if the last core solution is still optimal for the current core (see
 /// intReopt 2 and 3), write it in @p in and return its value, else return
 /// -Inf

 double last_still_optimal( std::vector< char > & in ) const;

 /// the multiplier of the Lagrangian bound of last_still_optimal(), < 0 if
 /// it is to be recomputed

 mutable double f_lambda = -1;

 /// the consecutive failures of the bound test of last_still_optimal(), and
 /// the calls still to skip before trying it again

 mutable int f_cert_fail = 0;
 mutable int f_cert_wait = 0;

 /// what get_reopt_outcome() returns

 mutable int f_outcome = 0;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;  // insert CoreDPBinaryKnapsackSolver in the factory

 };  // end( class( CoreDPBinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* CoreDPBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File CoreDPBinaryKnapsackSolver.h -----------------*/
/*--------------------------------------------------------------------------*/
