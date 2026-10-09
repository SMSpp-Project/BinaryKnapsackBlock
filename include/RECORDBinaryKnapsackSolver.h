/*--------------------------------------------------------------------------*/
/*------------------- File RECORDBinaryKnapsackSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class RECORDBinaryKnapsackSolver, an exact
 * Solver for the Binary Knapsack encoded by a BinaryKnapsackBlock that hands
 * the integer core of the instance to RECORD
 * (gitlab.com/renanfernandofranco/record), an external dynamic programming
 * solver for the 0-1 and bounded knapsack problems.
 *
 * Everything but the solve of the integer core is that of
 * CoreDPBinaryKnapsackSolver, from which the class derives: the mirror of
 * the instance and its update under the Modification, the objective sense,
 * the fixed variables, the complementation of the items with negative weight
 * and profit, and the case with continuous variables (which RECORD does not
 * cover, and which is left to the core enumeration). RECORD needs integer
 * profits: an instance whose core has a fractional one is solved by the core
 * enumeration too. With intReopt, the previous solution repaired to the
 * current data is the incumbent RECORD has to beat, and is returned if it
 * finds nothing better.
 *
 * The class is only built if a checkout of RECORD is given at configure time
 * (the RECORD_ROOT CMake variable).
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __RECORDBinaryKnapsackSolver
 #define __RECORDBinaryKnapsackSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CoreDPBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------- CLASS RECORDBinaryKnapsackSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for the BinaryKnapsackBlock based on the external RECORD
/** Exact Solver for the BinaryKnapsackBlock that solves the integer core of
 * the instance with RECORD; see the file-level comment. */

class RECORDBinaryKnapsackSolver : public CoreDPBinaryKnapsackSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 RECORDBinaryKnapsackSolver() : CoreDPBinaryKnapsackSolver() {}

 ~RECORDBinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the integer core with RECORD
 /** Solves the extracted integer core (v_w, v_p, f_C) with RECORD, falling
  * back to the core enumeration of CoreDPBinaryKnapsackSolver if a profit is
  * not integer; returns the optimal core profit and fills @p in. */

 double solve_integer_core( std::vector< char > & in ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;  // insert it in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( RECORDBinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* RECORDBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File RECORDBinaryKnapsackSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
