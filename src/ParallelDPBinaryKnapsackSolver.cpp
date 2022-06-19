/*--------------------------------------------------------------------------*/
/*----------------- File ParallelDPBinaryKnapsackSolver.cpp ----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class ParallelDPBinaryKnapsackSolver.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale
 */ 
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ParallelDPBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE AND USING --------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DPBinaryKnapsackSolver to the factory

SMSpp_insert_in_factory_cpp_1( ParallelDPBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS OF ParallelDPBinaryKnapsackSolver -----------------*/
/*--------------------------------------------------------------------------*/

void ParallelDPBinaryKnapsackSolver::dynamic_programming( Index C ) {
 OPENMP_dynamic_programming( C );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// OPENMP

void ParallelDPBinaryKnapsackSolver::OPENMP_dynamic_programming( Index C ) {

 std::vector< double > currlab = lab[ start_item ]; // set of current labels 
 std::vector< double > nextlab;                     // set of next labels  

 for( Index i = start_item ; i < f_N ; ++i ) {      // for each item

  if( i % step == 0  )          // reoptimization: save currlab in lab[ i ]
   lab[ i ] = currlab;             

  if( skip[ i ] || ! v_I[ i ] ) // skip preprocessed and continuous variables   
   continue;                               

  double p = v_P[ i ];          // profit of the current item
  int w = v_W[ i ];             // weight of the current item
  
  if( p < 0 && w < 0 ) {        // if both are negative, change the signs
   p = -p;
   w = -w;
   }                    

  // max size of next labels   
  Index maxnextlab = std::min( Index( currlab.size() + w ) , C + 1 );

  // initialize next labels (with -INF) and allocate precedessors
  nextlab.resize( maxnextlab );
  pred[ i + 1 ].resize( maxnextlab ); 

  // compute nextlab
  #pragma omp parallel for
  for( int j = 0 ; j < nextlab.size() ; ++j ) {

   nextlab[ j ] = -Inf< double >();       // initialize with -INF 

   if( ( j < currlab.size() ) && ( j - w >= 0 ) ) {     // both arcs

    // if both node don't exist
    if( ( currlab[ j ] == -Inf< double >() ) && 
        ( currlab[ j - w ] == -Inf< double >() ) )
     continue; 

    if( currlab[ j ] > currlab[ j - w ] + p ) {         // horizontal
     pred[ i + 1 ][ j ] = false;
     nextlab[ j ] = currlab[ j ];
     }
    else {                                              // diagonal
     pred[ i + 1 ][ j ] = true;
     nextlab[ j ] = currlab[ j - w ] + p;
     }
    continue;   
    }
    
   if( ( j - w >= 0 ) && ( currlab[ j - w ] != -Inf< double >() ) ) {
    pred[ i + 1 ][ j ] = true;
    nextlab[ j ] = currlab[ j - w ] + p;
    }

   if( ( j < currlab.size() ) && ( currlab[ j ] != -Inf< double >() ) ) {
    pred[ i + 1 ][ j ] = false;
    nextlab[ j ] = currlab[ j ];
    }

   }

   std::swap( currlab , nextlab );

  }
 
 // always save last vector of labels (instead of copying just swap)
 std::swap( lab.back() , currlab );

 }

/*--------------------------------------------------------------------------*/
/*------------- End File ParallelDPBinaryKnapsackSolver.cpp ----------------*/
/*--------------------------------------------------------------------------*/