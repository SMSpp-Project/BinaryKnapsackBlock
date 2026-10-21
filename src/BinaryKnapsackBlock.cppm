/*--------------------------------------------------------------------------*/
/*--------------------- File BinaryKnapsackBlock.cppm ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.binaryknapsackblock, which exports
 * every name that the headers of BinaryKnapsackBlock declare at namespace
 * scope, together with the modules it depends on; a translation unit can then
 * write
 *
 *     import smspp.binaryknapsackblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Francesca Demelas \n
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Federica Di Pasquale, Antonio Frangioni,
 *                      Francesca Demelas
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "BinaryKnapsackBlock.h"
#include "BinaryKnapsackSolver.h"
#include "COMBOBinaryKnapsackSolver.h"
#include "COMBOBridge.h"
#include "CoreDPBinaryKnapsackSolver.h"
#include "DPBinaryKnapsackSolver.h"
#include "GreedyRelaxationBinaryKnapsackSolver.h"
#include "RECORDBinaryKnapsackSolver.h"
#include "RECORDBridge.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.binaryknapsackblock;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // BinaryKnapsackBlock.h
 using SMSpp_di_unipi_it::BinaryKnapsackBlock;
 using SMSpp_di_unipi_it::BinaryKnapsackBlockChange;
 using SMSpp_di_unipi_it::BinaryKnapsackBlockMod;
 using SMSpp_di_unipi_it::BinaryKnapsackBlockRngdChange;
 using SMSpp_di_unipi_it::BinaryKnapsackBlockRngdMod;
 using SMSpp_di_unipi_it::BinaryKnapsackBlockSbstChange;
 using SMSpp_di_unipi_it::BinaryKnapsackBlockSbstMod;
 using SMSpp_di_unipi_it::BinaryKnapsackSolution;

 // BinaryKnapsackSolver.h
 using SMSpp_di_unipi_it::BinaryKnapsackSolver;

 // COMBOBinaryKnapsackSolver.h
 using SMSpp_di_unipi_it::COMBOBinaryKnapsackSolver;

 // CoreDPBinaryKnapsackSolver.h
 using SMSpp_di_unipi_it::CoreDPBinaryKnapsackSolver;

 // DPBinaryKnapsackSolver.h
 using SMSpp_di_unipi_it::DPBinaryKnapsackSolver;

 // GreedyRelaxationBinaryKnapsackSolver.h
 using SMSpp_di_unipi_it::GreedyRelaxationBinaryKnapsackSolver;

 // RECORDBinaryKnapsackSolver.h
 using SMSpp_di_unipi_it::RECORDBinaryKnapsackSolver;
}

export namespace SMSpp_di_unipi_it {
 namespace combo_bridge {

  // COMBOBridge.h
  using SMSpp_di_unipi_it::combo_bridge::solve_01;
 }
}

export namespace SMSpp_di_unipi_it {
 namespace record_bridge {

  // RECORDBridge.h
  using SMSpp_di_unipi_it::record_bridge::solve_01;
 }
}

/*--------------------------------------------------------------------------*/
/*------------------- End File BinaryKnapsackBlock.cppm --------------------*/
/*--------------------------------------------------------------------------*/
