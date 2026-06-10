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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Donato Meoli
 */ 
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ParallelDPBinaryKnapsackSolver.h"

#include <thread>

#if PDPBKS_PARALLEL
 #include <ff/parallel_for.hpp>
#endif

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

// register ParallelDPBinaryKnapsackSolver to the factory

SMSpp_insert_in_factory_cpp_1( ParallelDPBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS OF ParallelDPBinaryKnapsackSolver -----------------*/
/*--------------------------------------------------------------------------*/

// defined here (not defaulted in the header) so that the destructor of the
// pimpl'd ff::ParallelFor is instantiated where the type is complete

ParallelDPBinaryKnapsackSolver::~ParallelDPBinaryKnapsackSolver() = default;

/*--------------------------------------------------------------------------*/

void ParallelDPBinaryKnapsackSolver::dynamic_programming( Index C ) {

#if PDPBKS_PARALLEL
 // an explicit f_which_parallel ( >= 0 ) selects a parallel engine; the auto
 // mode ( < 0, the default ) runs the serial base, because parallelising the
 // inner (capacity) loop does not pay off here: the loop is memory-bandwidth
 // bound and each of the N items incurs a fork-join barrier, so adding workers
 // only adds overhead (see the scaling benchmark tests/.../batch-par). The
 // parallel engines below are kept opt-in for benchmarking and future use.
 const int which = f_which_parallel;
 if( which < 0 ) {
  DPBinaryKnapsackSolver::dynamic_programming( C );
  return;
  }

 // resolve the number of workers: 0 = all the available cores, 1 = serial
 Index nthreads = f_max_thread ? Index( f_max_thread )
                               : Index( std::thread::hardware_concurrency() );
 if( nthreads == 0 )            // hardware_concurrency() may return 0
  nthreads = 1;

 // the ParallelFor (hence its worker threads) is created once per solver
 // instance and reused across re-solves; blocking workers (default) sleep idle
 if( ( which == 1 ) && ( ! f_pf ) )
  f_pf = std::make_unique< ff::ParallelFor >(
          std::max< Index >( 1 , std::thread::hardware_concurrency() ) );

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

   // horizontal arc
   if( j < currlab.size() && currlab[ j ] > nextlab[ j ] ) {
    pred[ i + 1 ][ j ] = false;
    nextlab[ j ] = currlab[ j ];
    }

   // diagonal arc
   if( j >= w && 
      currlab[ j - w ] != - Inf< double >() && 
      currlab[ j - w ] + p > nextlab[ j ] ) {
    pred[ i + 1 ][ j ] = true;
    nextlab[ j ] = currlab[ j - w ] + p;
    }  
   }; // end lambda function  

  // compute nextlab in parallel
  switch( which ) {

   case( 0 ): {

    // OpenMP - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    #pragma omp parallel for num_threads( nthreads )
    for( Index j = 0 ; j < nextlab.size() ; ++j )
     f( j );

    break;
    }
   
   case( 1 ): {
    
    // FastFlow - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // run the parallel for on the persistent worker pool
    f_pf->parallel_for( 0 , nextlab.size() , 1 , 0 , f , nthreads );
    
    break;
    }

   case( 2 ): {

    // std::thread - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // run f() over a contiguous chunk of the slice
    auto g = [ & ]( Index start , Index stop ) {
     for( Index j = start ; j < stop ; ++j )
      f( j );
     };

    Index n_threads = std::min< Index >( nthreads , nextlab.size() );
    if( n_threads < 1 )
     n_threads = 1;

    const Index chunksize = nextlab.size() / n_threads;

    // spawn n_threads - 1 helpers (threads owned by value, joined and
    // destroyed below: no manual new/delete) and run the last, possibly
    // larger, chunk on the current thread
    std::vector< std::thread > threads;
    threads.reserve( n_threads - 1 );

    for( Index t = 0 ; t < n_threads - 1 ; ++t )
     threads.emplace_back( g , t * chunksize , ( t + 1 ) * chunksize );

    g( ( n_threads - 1 ) * chunksize , nextlab.size() );

    for( auto & th : threads )
     th.join();

    break;
    }
   
   case( 3 ): {
    
    // OpenMP in avanti - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // initialize nextlab
    nextlab.assign( maxnextlab , - Inf< double >() );
    
    #pragma omp parallel for num_threads( nthreads )
    for( Index j = 0 ; j < currlab.size() ; ++j ) {
   
     if( currlab[ j ] == -Inf< double >() )   // skip node with label = -inf 
      continue;                          
                          
     if( currlab[ j ] > nextlab[ j ] ) {            // horizontal arc 
      pred[ i + 1 ][ j ] = false;                         
      nextlab[ j ] = currlab[ j ];
      }
     
     }

    Index h = w > maxnextlab ? 0 : maxnextlab - w;
        
    #pragma omp parallel for num_threads( nthreads )
    for( Index j = 0 ; j < h ; ++j ) {

     if( currlab[ j ] == -Inf< double >() )
      continue;  
     
     if( currlab[ j ] + p > nextlab[ j + w ] ) {    // diagonal arc
      pred[ i + 1 ][ j + w ] = true;             
      nextlab[ j + w ] = currlab[ j ] + p;
      }     

     }   

    break;
    }

   case( 4 ): {
    
    // Doppio n avanti- - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // initialize nextlab
    nextlab.assign( maxnextlab , - Inf< double >() );
    
    for( Index j = 0 ; j < currlab.size() ; ++j ) {
   
     if( currlab[ j ] == -Inf< double >() )   // skip node with label = -inf 
      continue;                          
                          
     if( currlab[ j ] > nextlab[ j ] ) {            // horizontal arc 
      pred[ i + 1 ][ j ] = false;                         
      nextlab[ j ] = currlab[ j ];
      }
     
     }

    Index h = w > maxnextlab ? 0 : maxnextlab - w;
        
    for( Index j = 0 ; j < h ; ++j ) {

     if( currlab[ j ] == -Inf< double >() )
      continue;  
     
     if( currlab[ j ] + p > nextlab[ j + w ] ) {    // diagonal arc
      pred[ i + 1 ][ j + w ] = true;             
      nextlab[ j + w ] = currlab[ j ] + p;
      }     

     }   

    break;
    }
   
   default:
    throw( std::invalid_argument( "ParallelDPBinaryKnapsackSolver::"
                       "dynamic_programming: invalid parallel engine" ) );

   } // end( switch )

  std::swap( currlab , nextlab );

  }

 // always save last vector of labels (instead of copying just swap)
 std::swap( lab.back() , currlab );

#else
 // FastFlow not available at build time: fall back to the serial base
 DPBinaryKnapsackSolver::dynamic_programming( C );
#endif

 }

/*--------------------------------------------------------------------------*/
/*------------- End File ParallelDPBinaryKnapsackSolver.cpp ----------------*/
/*--------------------------------------------------------------------------*/
