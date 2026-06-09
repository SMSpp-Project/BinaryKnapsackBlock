/*--------------------------------------------------------------------------*/
/*------------------ File ParallelDPBinaryKnapsackSolver.h -----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class ParallelDPBinaryKnapsackSolver.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ParallelDPBinaryKnapsackSolver
 #define __ParallelDPBinaryKnapsackSolver  
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <memory>

#include "DPBinaryKnapsackSolver.h"

#define PDPBKS_PARALLEL 1
/* If PDPBKS_PARALLEL > 0, the labels of each dynamic programming slice can be
 * computed in parallel; the actual number of workers is controlled at run time
 * by the intMaxThread parameter (0 = use all available cores, 1 = serial,
 * k = k workers). Set to 0 in build setups where FastFlow is not available
 * (e.g. the plain makefiles): the solver then degrades to the serial base. */

#if PDPBKS_PARALLEL
namespace ff { class ParallelFor; }  // forward declaration (FastFlow)
#endif

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
/*------------------ CLASS ParallelDPBinaryKnapsackSolver ------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Parallel Dynamic Programming Solver for BinaryKnapsackBlock
/** Derived from DPBinaryKnapsackSolver, this solver computes the labels of
 * each dynamic programming slice (the inner loop over the residual capacity)
 * in parallel, with a choice of engines (OpenMP, FastFlow, std::thread) and a
 * number of workers controlled by the inherited intMaxThread parameter.
 *
 * Note, however, that this fine-grained inner-loop parallelization does NOT
 * pay off in practice: the inner loop is memory-bandwidth bound (a couple of
 * loads and one store per element, no real arithmetic), and each of the N
 * items needs a fork-join barrier, so adding workers only adds synchronization
 * overhead (see the scaling benchmark in tests/BinaryKnapsackBlock/batch-par).
 * For this reason the default (intWhichParallel < 0) runs the serial base
 * algorithm; the parallel engines are kept opt-in, mostly for benchmarking.
 * Worthwhile parallelism for the binary knapsack lives at a coarser level
 * (many independent knapsack solves), not inside a single dynamic program. */


class ParallelDPBinaryKnapsackSolver : public DPBinaryKnapsackSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 @{ */

 /// public enum for the int algorithmic parameters - - - - - - - - - - - - -

 enum int_par_type_PDPBKSlv {
  intWhichParallel = intLastAlgPar ,  ///< parallel engine selector
  intLastPDPBKSlvPar                  ///< first allowed new int par for derived
  };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor

ParallelDPBinaryKnapsackSolver( Index maxthread = 0 ,
                                int whichparallel = -1 ) :
                                DPBinaryKnapsackSolver() ,
                                f_max_thread( maxthread ) ,
                                f_which_parallel( whichparallel ) {}

/*--------------------------------------------------------------------------*/
 /// destructor

 // defined in the .cpp (not defaulted here) so that the destructor of the
 // pimpl'd ff::ParallelFor is instantiated where the type is complete

~ParallelDPBinaryKnapsackSolver() override;

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 using DPBinaryKnapsackSolver::set_par;  // keep the other set_par() visible

 /// set the int parameters of ParallelDPBinaryKnapsackSolver
 /** Set the int parameters specific of ParallelDPBinaryKnapsackSolver:
  *
  * - intMaxThread [0]: number of workers the parallel dynamic programming
  *   can spawn; 0 means "all the available cores", 1 means serial.
  *
  * - intWhichParallel [-1]: which parallel engine to use; -1 = auto (serial
  *   pruning base for a single worker or a small problem, FastFlow otherwise),
  *   0 = OpenMP, 1 = FastFlow, 2 = std::thread, 3 = OpenMP (push). Mostly
  *   meant to benchmark the different strategies. */

 void set_par( idx_type par , int value ) override {
  if( par == intMaxThread )
   f_max_thread = value;
  else if( par == intWhichParallel )
   f_which_parallel = value;
  else
   DPBinaryKnapsackSolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastPDPBKSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the int parameters of ParallelDPBinaryKnapsackSolver

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == intMaxThread )
   return( f_max_thread );
  if( par == intWhichParallel )
   return( f_which_parallel );
  return( DPBinaryKnapsackSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par == intWhichParallel )
   return( -1 );
  return( DPBinaryKnapsackSolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intWhichParallel" )
   return( intWhichParallel );
  return( DPBinaryKnapsackSolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string par_str = "intWhichParallel";
  if( idx == intWhichParallel )
   return( par_str );
  return( DPBinaryKnapsackSolver::int_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Binary Knapsack encoded by the current 
 * BinaryKnapsackBlock @{ */

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Dynamic Programming to solve the integer knapsack

 void dynamic_programming( Index C ) override;

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 int f_max_thread;      ///< intMaxThread: 0 = all cores, 1 = serial, k workers

 int f_which_parallel;  ///< parallel engine: < 0 = auto (FastFlow), 0 = OpenMP,
                        ///< 1 = FastFlow, 2 = std::thread, 3 = OpenMP (push)

#if PDPBKS_PARALLEL
 /// FastFlow parallel-for engine, created on first parallel use and reused
 /// across re-solves (re-creating it per call would spawn/join threads every
 /// time and tank performance)
 std::unique_ptr< ff::ParallelFor > f_pf;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  
 // insert ParallelDPBinaryKnapsackSolver in the factory

/*--------------------------------------------------------------------------*/

 }; // end( class( ParallelDPBinaryKnapsackSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ParallelDPBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File ParallelDPBinaryKnapsackSolver.h ---------------*/
/*--------------------------------------------------------------------------*/