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
using namespace ff;

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

  // define the lambda function to be executed in parallel
  auto f = [ & ]( Index j ) {
   
   nextlab[ j ] = -Inf< double >();       // initialize with -INF 
     
   if( ( j < currlab.size() ) && ( j >= w ) ) {     // both arcs

    // if both node don't exist
    if( ( currlab[ j ] == -Inf< double >() ) && 
        ( currlab[ j - w ] == -Inf< double >() ) )
     return;

     if( currlab[ j ] > currlab[ j - w ] + p ) {         // horizontal
      pred[ i + 1 ][ j ] = false;
      nextlab[ j ] = currlab[ j ];
      }
     else {                                              // diagonal
      pred[ i + 1 ][ j ] = true;
      nextlab[ j ] = currlab[ j - w ] + p;
      }
      return;
    }
    
   if( ( j >= w ) && ( currlab[ j - w ] != -Inf< double >() ) ) {
    pred[ i + 1 ][ j ] = true;
    nextlab[ j ] = currlab[ j - w ] + p;
    }

   if( ( j < currlab.size() ) && ( currlab[ j ] != -Inf< double >() ) ) {
    pred[ i + 1 ][ j ] = false;
    nextlab[ j ] = currlab[ j ];
    } 
  
   }; // end lambda function  


  // compute nextlab in parallel
  switch( WhichParallel ) {

   // seq - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   case( -1 ): {
    
    for( Index j = 0 ; j < nextlab.size() ; ++j )
     f( j );

    break;
   } 
   
   case( 0 ): {
   
   // OpenMP- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
    #pragma omp parallel for num_threads( MaxThread )
    for( Index j = 0 ; j < nextlab.size() ; ++j )
     f( j ); 
    
    break;
    }
   
   case( 1 ): {
    
    // FastFlow - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    // compute chunksize
    int chunksize = ( NChunk == -1 ) ? 0 : nextlab.size() / NChunk;
    if( Schedule )
     chunksize = - chunksize; 
    
    // Run the parallel for
    ParallelFor pf( MaxThread );    
    pf.parallel_for( 0 , nextlab.size() , 1 , chunksize , f );
    
    break;
    }

   case( 2 ): {
    
    // Thread - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    // define the lambda function to be executed in parallel
    auto g = [ & ]( Index start , Index stop ) {
     for( Index j = start ; j < stop ; ++j )
      f( j );
     };

    Index n_threads = MaxThread;
    if( nextlab.size() < n_threads )
     n_threads = nextlab.size();

    Index chunksize = nextlab.size() / n_threads;
    Index last_chunksize = chunksize + nextlab.size() % n_threads;

    // define the vector of threads
    std::vector< std::thread * > VecThreads;

    for( Index t = 0 ; t < n_threads - 1 ; ++t ) {
     VecThreads.push_back( new std::thread( g , t * chunksize , 
                                          ( t + 1 ) * chunksize ) );
     }

    // last thread
    VecThreads.push_back( new std::thread( g , 
                        ( n_threads - 1 ) * chunksize ,
                        ( n_threads - 1 ) * chunksize + last_chunksize ) );

    // join all threads 
    for( auto th : VecThreads )
     th->join(); 

    break;
    }
   
   default:
    throw( std::invalid_argument( " " ) );
    break;
   
   } // end( switch )

  std::swap( currlab , nextlab );

  }

 // always save last vector of labels (instead of copying just swap)
 std::swap( lab.back() , currlab );

 }

/*--------------------------------------------------------------------------*/
/*------------- End File ParallelDPBinaryKnapsackSolver.cpp ----------------*/
/*--------------------------------------------------------------------------*/
