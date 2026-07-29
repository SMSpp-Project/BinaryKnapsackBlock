/*--------------------------------------------------------------------------*/
/*---------------------- File IncrementalGreedyRelaxationBinaryKnapsackSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @file IncrementalGreedyRelaxationBinaryKnapsackSolver.h
 * @brief Concrete implementations of incremental greedy-based solvers for the Binary Knapsack Problem.
 *
 * This file defines a specialized hierarchy to solve Knapsack problems using greedy
 * strategies. It resolves the diamond inheritance between heuristic/relaxation
 * logic and greedy-specific behaviors.
 *
 * IncrementalGreedyRelaxationBinaryKnapsackSolver provides the exact solution
 * of the continuous relaxation * of a Binary Knapsack problem, efficiently
 * re-optimized after each Change  * to the Block (typically variable fixings),
 *  making it suitable as the * bounding solver at the nodes of a Branch&Bound algorithm.
 *
 * GreedyRelaxationSolver extends it with the RelaxationSolver interface,
 * additionally providing bounds and a solution for the original integer
 * problem (when available), as well as the generation of the two
 * sub-problems obtained by branching on the critical item.
 *
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
 * 	   	   Dipartimento di Informatica \n
 * 	   	   Universita' di Pisa \n
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
    /*---------------------- CLASS GreedyRelaxationSolver ----------------------*/
    /*--------------------------------------------------------------------------*/
    /*--------------------------- GENERAL NOTES --------------------------------*/
    /*--------------------------------------------------------------------------*/
    ///
    /**   */

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
        /** @name Solving a relaxation of the Binary Knapsack encoded by the current
         * BinaryKnapsackBlock @{ */

        ///
        /**  */

        int compute(bool changedvars = true) override;

        /** @} ---------------------------------------------------------------------*/
        /*---------------------- METHODS FOR READING RESULTS -----------------------*/
        /*--------------------------------------------------------------------------*/
        /** @name Accessing the found solutions (if any)
         *  @{ */

        /*--------------------------------------------------------------------------*/
        /// return a valid lower bound on the optimal objective function value
        /** TODO */

        OFValue get_lb(void) override
        {

            // if it is a minimization problem, the optimal value is a lower bound
            // for the original problem
            if (!f_sense)
                return (-obj);

            // otherwise it is a maximization problem and the solution without the
            // critical item is a feasible solution and it provides a lower bound for
            // the original problem
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
        /** */

        OFValue get_ub(void) override
        {

            // if it is a maximization problem, the optimal value is un upper bound
            // for the original problem
            if (f_sense)
                return (obj);

            // otherwise it is a minimization problem and the solution without the
            // critical item is a feasible solution and it provides an upper bound for
            // the original problem
            if (f_ciVal == 1)
                return (-obj);

            /* 			if (complemented[f_ci])
                        {
                            return (-obj - (1 - f_ciVal) * v_P[f_ci]);
                        } */

            return (-obj + f_ciVal * v_P[f_ci]);
        }
        /// return the value of the (current) solution
        /** Return the the value of the current solution. Change the sign according
         * to the sense of the problem (f_sense). */

        OFValue get_var_value() override { return f_sense ? obj : -obj; }

        /*--------------------------------------------------------------------------*/
        /// tells whether a true solution (a solution of the true original problem
        /// and not of the relaxed one solved by RelaxationSolver) is available
        /** Called after compute() this method has to return true if a true solution
         * of the original problem (not the relaxed one solved by RelaxationSolver)
         * is available to be read with get_true_var_solution().
         *
         * Once "the first" solution (if ever) has been read, new ones may be
         * produced, if the Solver allows it, by means of new_true_var_solution().*/

        bool has_var_solution(void) override { return (/*f_ci >= 0 && */ f_ci <= f_N); }

        bool is_var_feasible(void) override { return has_var_solution(); }

        /*--------------------------------------------------------------------------*/
        /// write the current true solution in the variables of the Block

        void get_var_solution(Configuration *solc = nullptr) override;

        Solution *get_Solution(Configuration *solc) override;

        /*--------------------------------------------------------------------------*/
        /// read the solution currently stored in the BinaryKnapsackBlock
        void set_var_blockSolution(void)
        {
            auto BKB = dynamic_cast<BinaryKnapsackBlock *>(this->f_Block);
            if (BKB == nullptr)
                throw(std::invalid_argument("GreedyRelaxationSolver::set_var_blockSolution(): the Block associated to the Solver is not a BinaryKnapsackBlock!"));
            BKB->get_x(v_x.begin());
        }

        /*--------------------------------------------------------------------------*/

        /** @} ---------------------------------------------------------------------*/
        /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
        /*--------------------------------------------------------------------------*/
        /** @name Changing the data of the model
         *  @{ */

        /** GreedyRelaxationSolver::add_Modification() is defined to properly react
         * to NBModification, i.e. the Binary Knapsack instance must be reloaded and
         * the list of modification must be cleared. */

        void add_Modification(sp_Mod &mod) override;

        /*--------------------------------------------------------------------------*/

        /// Apply the Change, in case of fixing variable it's done internally to the solver
        Change *apply(Change *, bool doUndo = false) override;

        /** @} ---------------------------------------------------------------------*/
        /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
        /*--------------------------------------------------------------------------*/

    protected:
        /* data of the Binary Knapsack instance - - - - - - - - - - - - - - - - - - */

        Index f_N;                        ///< the number of Items
        double f_C;                       ///< the Capacity of the Knapsack
        std::vector<double> v_W;          ///< vector of Weights
        std::vector<double> v_P;          ///< vector of Profits
        std::vector<unsigned char> v_fxd; ///< vector saying how the x are fixed
        /* < v_fxd[ i ] indicates if x_i is
         * fixed, with the following encoding:
         * 0 = not fixed ,
         * 1 = fixed to 0 ,f
         * 2 = fixed to 1                     */
        bool f_sense; ///< the sense of the objective

        Index f_ci;              ///< Index of the critical item
        double f_ciVal;          ///< Variable value of the critical
        double obj;              ///< the value of the objective
        std::vector<double> v_x; ///< vector of variables

        std::vector<Index> sortedVar;   /// indices of variables sorted by profit/weight ratio,
        std::vector<Index> varToSorted; /// < mapping from variable index to position in sortedVar
        std::vector<bool> skip;         ///< vector to mark items to skip during compute
        bool changedData = true;        ///< flag to signal if Data are	 changed, for modification or first load
        std::vector<bool> complemented; /// vector saying if the variable is complemented, then we have swapped the sign for his profit and weight
        Index startingSearchIndex;      ///< index in sortedVar from which to start the search for the critical item, used to speed up the search after branching
        double presolveCapacity;        /// remaining maximal Capacity, consider only fixed variables, used to speed up feasibility check
        double P;                       ///< evaluated profit of the solution without the critical item
        double C;                       ///< evaluated weight of the solution without the critical item
        bool reachTheEnd;               ///< evaluate if in update_v_x the last element is the critical item or not

        /*--------------------------------------------------------------------------*/
        /*---------------------------- PROTECTED FIELDS ----------------------------*/
        /*--------------------------------------------------------------------------*/

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
        /// initiliaze variables for the first time after loading or reloading the instance

        void initializeVariables();

        SMSpp_insert_in_factory_h; // insert GreedyHeuristicSolver in the factory

        /*--------------------------------------------------------------------------*/
        /*--------------------------- PRIVATE FIELDS -------------------------------*/
        /*--------------------------------------------------------------------------*/

        /*--------------------------------------------------------------------------*/
        /*--------------------------------------------------------------------------*/

    }; // end( class( IncrementalGreedyChangeBinaryKnapsackSolver ) )

    /*--------------------------------------------------------------------------*/
    /*----------IncrementalGreedyRelaxationBinaryKnapsackSolver-----------------*/
    /*--------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------*/
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
    /*--------------------------------------------------------------------------*/

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

        std::vector<Change *> branch() override;

        /** @} ---------------------------------------------------------------------*/
        /*---------------------- METHODS FOR READING RESULTS -----------------------*/
        /*--------------------------------------------------------------------------*/
        /** @name Accessing the found solutions (if any)
         *  @{ */

        /*--------------------------------------------------------------------------*/
        /// return a valid lower bound on the optimal objective function value
        /** TODO */
        // TODO check the correctness of true bounds
        OFValue get_true_lb(void) override
        {

            // if it is a minimization problem, the optimal value is a lower bound
            // for the original problem
            if (!f_sense)
                return (-obj);

            // otherwise it is a maximization problem and the solution without the
            // critical item is a feasible solution and it provides a lower bound for
            // the original problem
            if (f_ciVal == 1)
                return (obj);

            /* 			if (complemented[f_ci])
                        {
                            return (obj + (1 - f_ciVal) * v_P[f_ci]);
                        } */

            return (obj - f_ciVal * v_P[f_ci]);
        }

        /*
        Understand if it's needed to update correctly yhr get_ub and get_lb methods,
        because they are already defined in the base class IncrementalGreedyChangeBinaryKnapsackSolver,
        even if the correct version is this one since they are unfeasible
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
        /// return a valid upper bound on the optimal objective function value
        /** */

        // TODO check the correctness of true bounds
        OFValue get_true_ub(void) override
        {

            // if it is a maximization problem, the optimal value is un upper bound
            // for the original problem
            if (f_sense)
                return (obj);

            // otherwise it is a minimization problem and the solution without the
            // critical item is a feasible solution and it provides an upper bound for
            // the original problem
            if (f_ciVal == 1)
                return (-obj);

            /* 			if (complemented[f_ci])
                        {
                            return (-obj - (1 - f_ciVal) * v_P[f_ci]);
                        } */

            return (-obj + f_ciVal * v_P[f_ci]);
        }

        /*--------------------------------------------------------------------------*/
        /// tells whether a true solution (a solution of the true original problem
        /// and not of the relaxed one solved by RelaxationSolver) is available
        /** Called after compute() this method has to return true if a true solution
         * of the original problem (not the relaxed one solved by RelaxationSolver)
         * is available to be read with get_true_var_solution().
         *
         * Once "the first" solution (if ever) has been read, new ones may be
         * produced, if the Solver allows it, by means of new_true_var_solution().*/

        bool has_true_var_solution(void) override { return (/* f_ci >= 0 &&  */ f_ci <= f_N); }

        /*--------------------------------------------------------------------------*/
        /// write the current true solution in the variables of the Block

        void get_true_var_solution(Configuration *solc = nullptr) override
        {
            // TODO check if it's correct
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
        /// after has_true_var_solution(), says whether a there is another true solution available

        bool new_true_var_solution(void) override
        {
            return false;
        }
        Solution *get_true_solution(Configuration *solc) override
        {
            update_v_x();
            std::vector<double> sol_x(v_x);
            // TODO capire se serve complemented dato che l'ho usato anche prima
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

        SMSpp_insert_in_factory_h; // insert GreedyRelaxationSolver in the factory

        /*--------------------------------------------------------------------------*/
        /*--------------------------------------------------------------------------*/

    }; // end( class( GreedyRelaxationSolver ) )

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* GreedyRelaxationSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File GreedyRelaxationSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
