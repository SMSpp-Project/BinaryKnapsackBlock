/*--------------------------------------------------------------------------*/
/*--------------------- File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class GreedyRelaxationSolver, which
 * implements the RelaxationSolver concept [see Solver.h] for solving 
 * a relaxation of a  Knapsack problems as represented by a 
 * BinaryKnapsackKBlock.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Antonio Frangioni
 */ 
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "GreedyRelaxationSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using idx_type = ThinComputeInterface::idx_type;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset; 

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register GreedyRelaxationkSolver to the factory

SMSpp_insert_in_factory_cpp_1( GreedyRelaxationSolver );

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF GreedyRelaxationSolver ---------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int GreedyRelaxationSolver::compute( bool changedvars ) {

 for( Index i = 0 ; i < v_I.size() ; ++i )
  v_I[ i ] = false;

 return DPBinaryKnapsackSolver::compute();

 }  // end( GreedyRelaxationSolver::compute() )
 
/*--------------------------------------------------------------------------*/
/*----------------- End File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
