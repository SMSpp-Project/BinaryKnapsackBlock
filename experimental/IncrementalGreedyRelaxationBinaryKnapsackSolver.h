/*--------------------------------------------------------------------------*/
/*------ File IncrementalGreedyRelaxationBinaryKnapsackSolver.h ------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class
 * IncrementalGreedyRelaxationBinaryKnapsackSolver, which solves the very same
 * continuous (Dantzig) relaxation as GreedyRelaxationBinaryKnapsackSolver,
 * and shares its branching, separation, bounds and Solution, but maintains
 * the relaxation *incrementally* along the enumeration tree.
 *
 * Where the base greedy re-scans the whole efficiency order at each node
 * (O(n) per node), this variant remembers where the critical item sat and
 * keeps the running profit and residual capacity up to that point: an
 * (un)fixing then adjusts those running quantities and rewinds the search
 * index just enough, so that the next compute() resumes the greedy fill from
 * there instead of from the top. The per-node work becomes sub-linear (the
 * root solve stays O(n log n) for the sort), which trades raw speed for the
 * finest possible node granularity: it is the configuration that stresses
 * the tree-level parallelism the most, hence the interesting end of the
 * performance range for scalability studies. It is selected, in place of the
 * base greedy, through the inner BlockSolverConfig of the BranchAndXSolver.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Donato Meoli
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

#include "GreedyRelaxationBinaryKnapsackSolver.h"

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
    /*------- CLASS IncrementalGreedyRelaxationBinaryKnapsackSolver -------------*/
    /*--------------------------------------------------------------------------*/
    /*--------------------------- GENERAL NOTES --------------------------------*/
    /*--------------------------------------------------------------------------*/
    /// incremental variant of GreedyRelaxationBinaryKnapsackSolver
    /** Solves the continuous (Dantzig) relaxation of a BinaryKnapsackBlock
     * exactly like GreedyRelaxationBinaryKnapsackSolver, from which it inherits
     * the branching on the critical item, the reduced-cost separation, the true
     * bounds and the Solution. It only changes *how* the relaxation is kept along
     * the enumeration tree: the greedy fill is not redone from scratch at each
     * node but resumed from a saved position in the efficiency order, with the
     * profit and residual capacity accumulated up to there carried over and
     * patched by apply(). The branch Changes carry, alongside the fixed value,
     * the residual capacity and profit of the parent solution, so that a child
     * can pick up the parent state without recomputing it. */

    class IncrementalGreedyRelaxationBinaryKnapsackSolver
        : public GreedyRelaxationBinaryKnapsackSolver
    {

        /*--------------------------------------------------------------------------*/
        /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
        /*--------------------------------------------------------------------------*/

    public:
        /*--------------------------------------------------------------------------*/
        /*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
        /*--------------------------------------------------------------------------*/
        /** @name Constructor and Destructor
         *  @{ */

        /// constructor

        IncrementalGreedyRelaxationBinaryKnapsackSolver()
            : GreedyRelaxationBinaryKnapsackSolver(), f_inited(false),
              f_can_resume(false), f_ssi(0), f_run_p(0), f_run_rem(0) {}

        /// destructor

        ~IncrementalGreedyRelaxationBinaryKnapsackSolver() override = default;

        /** @} ---------------------------------------------------------------------*/
        /*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
        /*--------------------------------------------------------------------------*/
        /** @name Solving a relaxation of the Binary Knapsack
         *  @{ */

        /// solve the relaxation by resuming the greedy fill from the saved state

        int compute(bool changedvars = true) override;

        /** @} ---------------------------------------------------------------------*/
        /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
        /*--------------------------------------------------------------------------*/
        /** @name Changing the data of the model
         *  @{ */

        /// branch on the critical item, carrying the running state to the children
        /** Same two children as the base branch() (the critical item fixed to 0 and
         * to 1), with the residual capacity and profit accumulated up to the
         * critical item packed into each Change, so that apply() can resume the
         * child relaxation from the parent state. */

        std::vector<Change *> branch() override;

        /*--------------------------------------------------------------------------*/
        /// apply the Change, patching the running state for an (un)fixing one
        /** Same contract as the base apply(): (un)fixing Changes act on the
         * internal mirror and any other Change goes to the Block. On top of the
         * mirror, an (un)fixing here also patches the running profit / residual
         * capacity and the resume index, so that the next compute() needs not
         * rescan the items already settled in the parent. */

        Change *apply(Change *chg, bool doUndo = false) override;

        /** @} ---------------------------------------------------------------------*/
        /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
        /*--------------------------------------------------------------------------*/

    protected:
        /*--------------------------------------------------------------------------*/
        /*--------------------------- PROTECTED METHODS ----------------------------*/
        /*--------------------------------------------------------------------------*/

        /// position of every item in the cached efficiency order v_ord
        /** Inverse of v_ord (the order in which the greedy fill visits the items):
         * f_pos[ i ] is the rank of item i, needed to move the resume index to the
         * rank of an (un)fixed item. Rebuilt together with v_ord. */

        void rebuild_positions(void);

        /*--------------------------------------------------------------------------*/
        /*---------------------------- PROTECTED FIELDS ----------------------------*/
        /*--------------------------------------------------------------------------*/

        bool f_inited; ///< whether the running state matches the current core

        /// whether the carried state still describes the items before the resume
        /** The greedy fill can be resumed only when the solution of the items
         * before the resume rank is still the one the running state was built on,
         * which holds for a single descent step right after a compute(). A best-
         * first or breadth-first move walks several nodes between two compute()s,
         * so only the first apply() after a compute() may resume; the others
         * restart the fill from the top, which is correct on any navigation. */

        bool f_can_resume;

        Index f_ssi; ///< resume rank in v_ord: the greedy fill restarts here

        double f_run_p; ///< core profit accumulated by the items before f_ssi

        double f_run_rem; ///< residual capacity left after the items before f_ssi

        std::vector<Index> f_pos; ///< rank of each item in v_ord [see f_pos]

        /*--------------------------------------------------------------------------*/
        /*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
        /*--------------------------------------------------------------------------*/

    private:
        /*--------------------------------------------------------------------------*/
        /*--------------------------- PRIVATE FIELDS -------------------------------*/
        /*--------------------------------------------------------------------------*/

        SMSpp_insert_in_factory_h;
        // insert IncrementalGreedyRelaxationBinaryKnapsackSolver in the factory

        /*--------------------------------------------------------------------------*/
        /*--------------------------------------------------------------------------*/

    }; // end( class( IncrementalGreedyRelaxationBinaryKnapsackSolver ) )

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* IncrementalGreedyRelaxationBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*--- End File IncrementalGreedyRelaxationBinaryKnapsackSolver.h -----------*/
/*--------------------------------------------------------------------------*/
