/*--------------------------------------------------------------------------*/
/*-------------------- File COMBOBinaryKnapsackSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class COMBOBinaryKnapsackSolver, an exact
 * Solver for the Binary Knapsack encoded by a BinaryKnapsackBlock that hands
 * the integer core of the instance to COMBO (S. Martello, D. Pisinger,
 * P. Toth, "Dynamic programming and strong bounds for the 0-1 knapsack
 * problem", Management Science 45, 1999), an external dynamic programming
 * solver for the 0-1 knapsack problem.
 *
 * Everything but the solve of the integer core is that of
 * CoreDPBinaryKnapsackSolver, from which the class derives: the mirror of
 * the instance and its update under the Modification, the objective sense,
 * the fixed variables, the complementation of the items with negative weight
 * and profit, and the case with continuous variables (which COMBO does not
 * cover, and which is left to the core enumeration). COMBO needs integer
 * profits: an instance whose core has a fractional one is solved by the core
 * enumeration too. With intReopt, the previous solution repaired to the
 * current data is the incumbent COMBO has to beat, and is returned if it
 * finds nothing better. COMBO stops on an internal error (e.g., the state
 * space it allocates up front is exhausted), in which case compute()
 * throws.
 *
 * The class is only built if a directory holding combo.c and combo.h, as
 * distributed by their authors, is given at configure time (the COMBO_ROOT
 * CMake variable); the code is for academic or non-commercial use only, and
 * is not distributed with SMS++.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __COMBOBinaryKnapsackSolver
 #define __COMBOBinaryKnapsackSolver
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
/*------------------- CLASS COMBOBinaryKnapsackSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for the BinaryKnapsackBlock based on the external COMBO
/** Exact Solver for the BinaryKnapsackBlock that solves the integer core of
 * the instance with COMBO; see the file-level comment. */

class COMBOBinaryKnapsackSolver : public CoreDPBinaryKnapsackSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 COMBOBinaryKnapsackSolver() : CoreDPBinaryKnapsackSolver() {}

 ~COMBOBinaryKnapsackSolver() override = default;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the integer core with COMBO
 /** Solves the extracted integer core (v_w, v_p, f_C) with COMBO, falling
  * back to the core enumeration of CoreDPBinaryKnapsackSolver if a profit is
  * not integer; returns the optimal core profit and fills @p in. */

 double solve_integer_core( std::vector< char > & in ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;  // insert it in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( COMBOBinaryKnapsackSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* COMBOBinaryKnapsackSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File COMBOBinaryKnapsackSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
