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

// #include "RelaxationSolver.h"
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
    /** Solves the continuous (Dantzig) relaxation of the BinaryKnapsackBlock
     * by the greedy fill along the items sorted by nonincreasing efficiency,
     * and keeps it across the Changes applied to it [see apply()]: a fixing
     * updates the residual capacity and the profit accumulated so far, and
     * moves the position in the efficiency order from which the next
     * compute() resumes the search of the critical item, instead of
     * restarting it from the first item. Items with negative weight and
     * profit are complemented, and the data are normalized to a maximization
     * [see load()]. */

    class IncrementalGreedyChangeBinaryKnapsackSolver : public virtual ChangeSolver, public Solver
    {

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

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Constructor and Destructor
         *  @{ */

/*--------------------------------------------------------------------------*/
        /// constructor

        IncrementalGreedyChangeBinaryKnapsackSolver() : ChangeSolver(),
                                                        sortedVar(),
                                                        skip(),
                                                        varToSorted(),
                                                        v_x(),
                                                        changedData(true),
                                                        complemented(),
                                                        startingSearchIndex(0),
                                                        presolveCapacity(0),
                                                        C(0),
                                                        P(0),
                                                        reachTheEnd(false),
                                                        f_N(0),
                                                        f_C(0),
                                                        f_sense(true),
                                                        f_ci(0),
                                                        f_ciVal(0),
                                                        obj(-Inf<double>()) {
                                                        };

/*--------------------------------------------------------------------------*/
        /// destructor

        ~IncrementalGreedyChangeBinaryKnapsackSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Other initializations @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        /// set the (pointer to the) Block that the Solver has to solve

        void set_Block(Block *block) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Solving a relaxation of the Binary Knapsack encoded by the
         * current BinaryKnapsackBlock @{ */

        /// solve the continuous relaxation, resuming from the last Changes

        int compute(bool changedvars = true) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Accessing the found solutions (if any)
         *  @{ */

/*--------------------------------------------------------------------------*/
        /// return a valid lower bound on the optimal objective function value
        /** For a minimization problem the relaxation optimum is a lower bound;
         * for a maximization one, the value of the greedy solution without
         * the critical item, which is feasible. */

        OFValue get_lb(void) override
        {

            // if it is a minimization problem, the optimal value is a lower
            // bound for the original problem
            if (!f_sense)
                return (-obj);

            // otherwise it is a maximization problem and the solution without
            // the critical item is a feasible solution and it provides a lower
            // bound for the original problem
            if (f_ciVal == 1)
                return (obj);

            /* 			if (complemented[f_ci] < 0)
                        {
                            return (obj + (1 - f_ciVal) * v_P[f_ci]);
                        }
             */
            return (obj - f_ciVal * v_P[f_ci]);
        }

/*--------------------------------------------------------------------------*/
        /// return a valid upper bound on the optimal objective function value
        /** Symmetric to get_lb(): the relaxation optimum for a maximization
         * problem, the value of the greedy solution without the critical item
         * for a minimization one. */

        OFValue get_ub(void) override
        {

            // if it is a maximization problem, the optimal value is an upper
            // bound for the original problem
            if (f_sense)
                return (obj);

            // otherwise it is a minimization problem and the solution without
            // the critical item is a feasible solution and it provides an
            // upper bound for the original problem
            if (f_ciVal == 1)
                return (-obj);

            /* 			if (complemented[f_ci])
                        {
                            return (-obj - (1 - f_ciVal) * v_P[f_ci]);
                        } */

            return (-obj + f_ciVal * v_P[f_ci]);
        }

/*--------------------------------------------------------------------------*/
        /// return the value of the (current) solution
        /** Returns the value of the current solution, with the sign of the
         * sense of the problem (f_sense). */

        OFValue get_var_value() override { return f_sense ? obj : -obj; }

/*--------------------------------------------------------------------------*/
        /// tells whether a solution of the relaxation is available

        bool has_var_solution(void) override { return (/*f_ci >= 0 && */ f_ci <= f_N); }

        /// tells whether the solution of the relaxation is feasible for it

        bool is_var_feasible(void) override { return has_var_solution(); }

/*--------------------------------------------------------------------------*/
        /// write the current solution in the variables of the Block

        void get_var_solution(Configuration *solc = nullptr) override;

        /// physically construct the Solution, the critical item rounded

        Solution *get_Solution(Configuration *solc) override;

/*--------------------------------------------------------------------------*/
        /// read the solution currently stored in the BinaryKnapsackBlock
        void set_var_blockSolution(void)
        {
            auto BKB = dynamic_cast<BinaryKnapsackBlock *>(this->f_Block);
            if (BKB == nullptr)
                throw(std::invalid_argument("IncrementalGreedyChangeBinaryKnapsackSolver::set_var_blockSolution: the Block is not a BinaryKnapsackBlock"));
            BKB->get_x(v_x.begin());
        }

/*--------------------------------------------------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
        /** @name Changing the data of the model
         *  @{ */

        /// add a Modification to the list of those to be processed
        /** Reacts to a NBModification by reloading the instance and clearing
         * the list of the Modification; any other one is stored, to be
         * processed by the next compute(). */

        void add_Modification(sp_Mod &mod) override;

/*--------------------------------------------------------------------------*/

        /// apply the Change; (un)fixing Changes are applied internally
        /** The (un)fixing Changes [see branch()] update the state of the
         * greedy fill (residual capacity, accumulated profit and position in
         * the efficiency order) without reaching the Block, and so does the
         * returned undo Change; any other Change is forwarded to the Block. */
        Change *apply(Change *, bool doUndo = false) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

    protected:
/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

        // the data of the Binary Knapsack instance

        Index f_N;                        ///< the number of Items
        double f_C;                       ///< the Capacity of the Knapsack
        std::vector<double> v_W;          ///< vector of Weights
        std::vector<double> v_P;          ///< vector of Profits
        std::vector<unsigned char> v_fxd; ///< how the x are fixed
        /**< v_fxd[ i ] says whether x_i is fixed, with the encoding
         * 0 = not fixed, 1 = fixed to 0, 2 = fixed to 1 */
        bool f_sense; ///< the sense of the objective

        Index f_ci;              ///< index of the critical item
        double f_ciVal;          ///< value of the critical item
        double obj;              ///< the value of the objective
        std::vector<double> v_x; ///< vector of variables

        /// the items sorted by nonincreasing efficiency (profit / weight)
        std::vector<Index> sortedVar;
        /// the position of each item in sortedVar
        std::vector<Index> varToSorted;
        /// the items that compute() skips (fixed, or settled by their signs)
        std::vector<bool> skip;
        /// whether the data have changed (Modification, or first load)
        bool changedData = true;
        /// whether the item is complemented, i.e., its profit and weight
        /// have changed sign
        std::vector<bool> complemented;
        /// the position in sortedVar from which compute() resumes the search
        /// of the critical item (after a branching, say)
        Index startingSearchIndex;
        /// the capacity left by the fixed items alone, for a quick
        /// feasibility check
        double presolveCapacity;
        /// the profit of the greedy solution without the critical item
        double P;
        /// the residual capacity of the greedy solution without the critical
        /// item
        double C;
        /// whether the greedy fill reached the last item, i.e., there is no
        /// critical item
        bool reachTheEnd;

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

        /// write in v_x the greedy solution, the critical item at f_ciVal

        void update_v_x()
        {
            double xval = 1;
            for (Index i : sortedVar)
            {
                if (skip[i])
                {
                    if (i == f_ci)
                        xval = 0;
                    continue;
                }

                if (i == f_ci)
                {
                    v_x[i] = complemented[i] ? 1 - f_ciVal : f_ciVal;
                    xval = 0;
                }
                else
                {
                    v_x[i] = complemented[i] ? 1 - xval : xval;
                }
            }
        }

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
        /// process all the pending modifications

        void process_outstanding_Modification();
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        /// initialize the state of the greedy fill after (re)loading the
        /// instance

        void initializeVariables();

        SMSpp_insert_in_factory_h; // insert it in the factory

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

    }; // end( class( IncrementalGreedyChangeBinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/
/*---------- CLASS IncrementalGreedyRelaxationBinaryKnapsackSolver ---------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
    /// RelaxationSolver of the continuous relaxation of a BinaryKnapsackBlock
    /** Extends IncrementalGreedyChangeBinaryKnapsackSolver with the
     * RelaxationSolver interface: the true bounds and solution of the
     * original integer problem, and the branching on the critical item. */

    class IncrementalGreedyRelaxationBinaryKnapsackSolver : public IncrementalGreedyChangeBinaryKnapsackSolver,
                                                            public virtual RelaxationSolver
    {

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

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Constructor and Destructor
         *  @{ */

/*--------------------------------------------------------------------------*/
        /// constructor

        IncrementalGreedyRelaxationBinaryKnapsackSolver() : IncrementalGreedyChangeBinaryKnapsackSolver() {};

/*--------------------------------------------------------------------------*/
        /// destructor

        ~IncrementalGreedyRelaxationBinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
        /// branch on the critical item
        /** The children fix the critical item to 1 and to 0 (only to 0 if it
         * does not fit the capacity left by the fixed items). Besides the
         * value, the data of each child carry, in its last two positions, the
         * residual capacity and the profit of the greedy solution of the
         * child, so that apply() picks up the state of the parent without
         * recomputing it. */

        std::vector<Change *> branch() override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Accessing the found solutions (if any)
         *  @{ */

/*--------------------------------------------------------------------------*/
        /// return a valid lower bound on the optimal value of the true problem
        /** As get_lb(). TODO: the correctness of the true bounds is to be
         * checked. */

        OFValue get_true_lb(void) override
        {

            // if it is a minimization problem, the optimal value is a lower
            // bound for the original problem
            if (!f_sense)
                return (-obj);

            // otherwise it is a maximization problem and the solution without
            // the critical item is a feasible solution and it provides a lower
            // bound for the original problem
            if (f_ciVal == 1)
                return (obj);

            /* 			if (complemented[f_ci])
                        {
                            return (obj + (1 - f_ciVal) * v_P[f_ci]);
                        } */

            return (obj - f_ciVal * v_P[f_ci]);
        }

        /* TODO: whether get_ub() and get_lb(), already defined in the base
         * IncrementalGreedyChangeBinaryKnapsackSolver, should rather be the
         * following ones, the relaxed solution being infeasible:
                 OFValue get_ub(void) override
                {
                    return f_sense ? obj : -obj;
                }

                OFValue get_lb(void) override
                {
                    return !f_sense ? -obj : obj;
                }
        */

/*--------------------------------------------------------------------------*/
        /// return a valid upper bound on the optimal value of the true problem
        /** As get_ub(). TODO: the correctness of the true bounds is to be
         * checked. */

        OFValue get_true_ub(void) override
        {

            // if it is a maximization problem, the optimal value is an upper
            // bound for the original problem
            if (f_sense)
                return (obj);

            // otherwise it is a minimization problem and the solution without
            // the critical item is a feasible solution and it provides an
            // upper bound for the original problem
            if (f_ciVal == 1)
                return (-obj);

            /* 			if (complemented[f_ci])
                        {
                            return (-obj - (1 - f_ciVal) * v_P[f_ci]);
                        } */

            return (-obj + f_ciVal * v_P[f_ci]);
        }

/*--------------------------------------------------------------------------*/
        /// tells whether a true solution (a solution of the true original
        /// problem and not of the relaxed one solved by RelaxationSolver) is
        /// available
        /** Called after compute() this method has to return true if a true
         * solution of the original problem (not the relaxed one solved by
         * RelaxationSolver) is available to be read with
         * get_true_var_solution().
         *
         * Once "the first" solution (if ever) has been read, new ones may be
         * produced, if the Solver allows it, by means of
         * new_true_var_solution().*/

        bool has_true_var_solution(void) override { return (/* f_ci >= 0 &&  */ f_ci <= f_N); }

/*--------------------------------------------------------------------------*/
        /// write the current true solution in the variables of the Block
        /** The true solution is the greedy one with the critical item rounded.
         * TODO: its correctness is to be checked. */

        void get_true_var_solution(Configuration *solc = nullptr) override
        {
            update_v_x();
            BinaryKnapsackBlock *BKB = static_cast<BinaryKnapsackBlock *>(this->f_Block);
            BKB->lock(this);
            for (Index i = 0; i < f_N; i++)
                if (i == f_ci)
                    BKB->set_x(i, complemented[f_ci] ? (v_x[f_ci] == 1 ? 0 : 1) : (v_x[f_ci] == 1 ? 1 : 0));
                else
                    BKB->set_x(i, v_x[i]);
            BKB->unlock(this);
        }

/*--------------------------------------------------------------------------*/
        /// after has_true_var_solution(), whether another true solution exists

        bool new_true_var_solution(void) override
        {
            return false;
        }

/*--------------------------------------------------------------------------*/
        /// physically construct the true Solution, the critical item rounded

        Solution *get_true_solution(Configuration *solc) override
        {
            update_v_x();
            std::vector<double> sol_x(v_x);
            // TODO: whether complemented is needed here, since update_v_x()
            // has already used it
            if (!reachTheEnd)
                sol_x[f_ci] = (complemented[f_ci] && !skip[f_ci]) ? (sol_x[f_ci] == 1 ? 0 : 1) : (sol_x[f_ci] == 1 ? 1 : 0);
            //  debug use only
            //  double value = 0;
            //  for (Index i = 0; i < f_N; ++i)
            //	value += sol_x[i] * v_P[i];
            return new BinaryKnapsackSolution(std::move(sol_x));
        }

/*--------------------------------------------------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

        SMSpp_insert_in_factory_h; // insert it in the factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

    }; // end( class( IncrementalGreedyRelaxationBinaryKnapsackSolver ) )

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* IncrementalGreedyRelaxationBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*---- End File IncrementalGreedyRelaxationBinaryKnapsackSolver.h ----------*/
/*--------------------------------------------------------------------------*/
