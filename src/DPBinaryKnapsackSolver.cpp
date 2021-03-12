/*--------------------------------------------------------------------------*/
/*--------------------- File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/

 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DPBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DPBinaryKnapsackSolver to the factory

SMSpp_insert_in_factory_cpp_1( DPBinaryKnapsackSolver );

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF DPBinaryKnapsackSolver ---------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::set_Block( Block * block ){
 
 if( block == f_Block )       // nothing to do  
  return;

 Solver::set_Block( block );  // attach to the new Block

 load();                      // load Binary Knapsack instance

}

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int DPBinaryKnapsackSolver::compute( bool changedvars ){

/* The variable i is the first item processed by the DP algorithm. If compute
 is called for the first time then i = 0 and the algorithm is entirely 
 executed. Otherwise process_outstanding_Modification(), after processing all 
 the modifications, returns the index of the item from which to restart the 
 algorithm (if it is possible according to the reopt parameter).            */
 
 int i = process_outstanding_Modification();

// check if the problem is empty - - - - - - - - - - - - - - - - - - - - - - -

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( BKB->is_empty() ){
  obj = - Inf< double >();
  return( kInfeasible );
 }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

int N = pi.size() + ni.size();      // number of items (not preprocessed)

if( N == 0 ){                       // if N == 0 all the items have been
 obj = prp_P;                       // preprocessed, no need of DP
 return( kOK );                     // return
}

if( i == N )                        // if i == N it is assumed that everything 
 return( kOK );                     // has already been done in 
                                    // process_outstanding_Modification()

int C = f_C - prp_W;                // new capacity after the preprocessing

G.resize( N + 1 );                  // resize the DP graph

int maxcurrlab;                     // max current height of the graph
int maxnextlab;                     // max next height of the graph    
int size;                           // size of the current set of labels

std::vector< double > currlab;      // set of current labels
std::vector< double > nextlab;      // set of next labels

// Initialization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


std::vector< int >::iterator it;    // iterator over ni[] and pi[]

if( i == 0 ){                       // if i == 0 the DP algorithm is entirely
 G[ 0 ].minLab = 0;                 // executed, initialize dummy node in the 
 maxcurrlab = 0;                    // origin
 size = 1;                         
 currlab.resize( size );
 currlab[ 0 ] = 0;

 it = ni.begin();                   // and start from the beginning   
}
else{                               // otherwise initialize with data stored 
                                    // in G[ i ]
 size = G[ i ].lab.size();         
 maxcurrlab = size + G[ i ].minLab - 1;
 currlab.resize( size );
 std::copy( G[ i ].lab.begin() , G[ i ].lab.end() , currlab.begin() );

 // and start from the correct item
 it = i >= ni.size() ? pi.begin() + i - ni.size() : ni.begin() + i;
}

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
/*
The algorithm processes first the items with negative weight and second the 
items with positive weights.

At each iteration it computes the maximum and the minimum "height" of the 
constructed graph in order to allocate only the necessary memory to store the 
current labels.

According to the reopt parameter, labels can be stored in G[ i ].lab so that
they will be used for reoptimization.
                                                                            */

// DP Algorithm- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for( ; i < N ; i++ , it++ ){

 // if all the items with negative weights have been processed 
 // start process the items with positive weights
 if( it == ni.end() )
  it = pi.begin();

  double p = v_P[ *it ];        // profit of the current item
  int w = v_W[ *it ];           // weight of the current item  

  // compute max and min "height" of the graph and resize accordingly
  // and if w > 0 check if the maximum capacity is reached
  // Also store minLab in G[ i ].minLab
  maxnextlab = w > 0 ? std::min( maxcurrlab + w , C + nW ) : maxcurrlab;
  G[ i + 1 ].minLab = std::min( G[ i ].minLab , G[ i ].minLab + w );

  int nextsize = maxnextlab - G[ i + 1 ].minLab + 1;

  // Initialize labels 
  nextlab.resize( nextsize );
  std::fill( nextlab.begin() , nextlab.end() , -Inf< double >() );
  
  // Allocate predecessors
  G[ i + 1 ].pred.resize( nextsize );
  
  // Initialize the best label among those of this slice 
  double bestlab = -Inf< double >();

  for( int j = 0 ; j < size ; j++ ){
   
   // if the node does not exists (currlab = -Inf) 
   // or if it does not have a better label continue
   if( currlab[ j ] <= bestlab )
    continue;

   // update bestlab
   bestlab = currlab[ j ];

   // horizontal and diagonal indeces (are shifted by minLab)
   int hi = j + G[ i ].minLab - G[ i + 1 ].minLab;
   int di = hi + w;

   // horizontal arc
   if( currlab[ j ] > nextlab[ hi ] ){
    G[ i + 1 ].pred[ hi ] = false;
    nextlab[ hi ] = currlab[ j ];
   }

   // check if the maximum capacity is reached
   if( ( w > 0 ) && ( di > C + nW ) )
    continue;

   // diagonal arc 
   if( currlab[ j ] + p > nextlab[ di ] ){
    G[ i + 1 ].pred[ di ] = true;
    nextlab[ di ] = currlab[ j ] + p;
   }

  }

 // if reopt save currlab
 if( reopt ){
  G[ i ].lab.resize( size );
  std::copy( currlab.begin() , currlab.end() , G[ i ].lab.begin() );
 }

  
 size = nextsize;
 maxcurrlab = maxnextlab;
 std::swap( currlab , nextlab );

}

// find optimal value- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// save last labels
G[ N ].lab.resize( size );
std::copy( currlab.begin() , currlab.end() , G[ N ].lab.begin() );

// find optimal value
int maxh = std::min( size , C + nW + 1 );
int besth;

obj = -Inf< double >();      // Initialize objective value
for( int i = 0 ; i < maxh ; i++ ){
 if( G[ N ].lab[ i ] > obj ){
  obj = G[ N ].lab[ i ];
  besth = i;
 }
}

obj += prp_P;               // add the profit coming from the preprocessing

// reconstruct the optimal solution - - - - - - - - - - - - - - - - - - - - - 

// if there are no items with positive weights start from ni.end() - 1
it = pi.size() ? pi.end() - 1 : ni.end() - 1;

for( int i = N ; i > 0 ; i-- , it-- ){

 if( G[ i ].pred[ besth ] ){
  v_x[ *it ] = 1;
  besth += G[ i ].minLab - G[ i - 1 ].minLab - v_W[ *it ];
 }
 else
  besth += G[ i ].minLab - G[ i - 1 ].minLab;

 if( it == pi.begin() )         // if all the items with positive weights have
  it = ni.end();                // been processed, start process the items 
                                // with negative weights
} 

return( kOK );

}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::get_var_solution( Configuration * solc ){
 
 // check if there is an available solution
 if( ! v_x.size() )
  throw( std::invalid_argument( "No available solution" ) );

 // check if the BinaryKnapsackBlock has the correct number of variables, 
 // that is it uses the correct abstract representation
 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( v_x.size() != BKB->get_VarSize() )
  throw( std::invalid_argument( "Incompatible variables size" ) );

 // write the current solution in the BinaryKnapsackBlock
 BKB->set_x( v_x );

}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::add_Modification( sp_Mod &mod ){
  
  if( f_no_Mod )
   return;

  // try to acquire lock, spin on failure
  while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

  const auto tmod = std::dynamic_pointer_cast< NBModification >( mod );
  
  if( tmod ){
   load();        // re-load the Binary Knapsack instance
   v_mod.clear(); 
  } 
  else 
   v_mod.push_back( mod );

  f_mod_lock.clear( std::memory_order_release );  // release lock

}

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 void DPBinaryKnapsackSolver::load(){

  if( f_Block ){

   auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
   if( ! BKB )
    throw( std::invalid_argument( "Block must be a BinaryKnapsackBlock" ) );

   // try to read
   if( ! BKB->read_lock() )
    throw(std::logic_error("cannot acquire read_lock on BinariKnapsackBlock"));

   // load Binary Knapsack instance - - - - - - - - - - - - - - - - - - - - - 
   f_sense = BKB->get_objective_sense();    // get the sense of the objective
   f_NItems = BKB->get_NItems();            // get the number of items
   f_C = std::floor( BKB->get_Capacity() ); // get the Capacity

   // get profits
   v_P.resize( f_NItems );
   const std::vector< double > & P = BKB->get_Profits();

   // if the sense is minimization change the sign of the profits      
   for( int i = 0 ; i < P.size() ; i++ )
    v_P[ i ] = f_sense ? P[ i ] : - P[ i ];

   // prepare vector of weights
   v_W.resize( f_NItems );
   const std::vector< double > & W = BKB->get_Weights();

   // load weights and compute nW, prp_P and prp_W
   nW = 0;
   prp_P = 0;
   prp_W = 0;
   
   for( int i = 0 ; i < W.size() ; i++ ){
   
    // check that weights are integers
    v_W[ i ] = std::round( W[ i ] );

    if( std::abs( v_W[ i ] - W[ i ] ) > 1e-06 )
     throw( std::invalid_argument( "Weights must be integers" ) );

    // if the variables is fixed to 1 update prp_P and prp_W
    if( BKB->is_fixed( i ) ){
     if( BKB->get_x( i ) ){
      prp_P += v_P[ i ];
      prp_W += v_W[ i ];
     }
     continue; 
    }

    // items with negative weights involved in the DP algorithm
    if( ( v_W[ i ] < 0 ) && ( v_P[ i ] < 0 ) )
     nW -= v_W[ i ];
    
    
    // items with negative weights NOT involved in the DP algorithm
    if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){
     prp_P += v_P[ i ];
     prp_W += v_W[ i ]; 
    }
   }

   BKB->read_unlock();

 // PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 /*
 The following items are preprocessed:
  
  - Items i with negative weight and positive profit => set x[ i ] = 1
 
  - Items i with positive weight and negative profit => set x[ i ] = 0
 
  - Items with a weight greater than C + nW, where C = f_C - prp_W 

 Actually the variables corresponding to all the remaining items are also
 initialized to 0, but the second set of preprocessed items are discarded
 from the computation                                                      */
 
 int C = f_C - prp_W;
 v_x.resize( f_NItems );

 for( int i = 0 ; i < f_NItems ; i++ ){

  // if the variable is fixed then update v_x and continue (discard it)
  if( BKB->is_fixed( i ) ){
   v_x[ i ] = BKB->get_x( i );
   continue;
  }

  // items to be picked
  if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){
   v_x[ i ] = 1;
   continue;
  }

  // items not to be picked
  if( ( v_W[ i ] > C + nW ) || ( ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ) ){
   v_x[ i ] = 0;
   continue; 
  }

  // store the indeces of the items with positive weights that are involved 
  // in the DP algorithm in pi[] and the indeces of the items with negative 
  // weights in ni[]

  if( v_W[ i ] > 0 )    // items with positive weights
   pi.push_back( i ); 
  else                  // items with negative weights
   ni.push_back( i ); 
  
  v_x[ i ] = 0;         // initialize all variables to 0

 }

 } // end if( f_Block )

} // end( DPBinaryKnapsackSolver::load() )

/*--------------------------------------------------------------------------*/

int DPBinaryKnapsackSolver::process_outstanding_Modification(){

 while( ! v_mod.empty() ){
  
  auto mod = v_mod.front();

  // BinaryKnapsackBlockMod - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod.get() );

   if( tmod ){

    switch( tmod->type() ){
     case( BinaryKnapsackBlockMod::eChgCapacity ): break;
     case( BinaryKnapsackBlockMod::eChgSense ): break;  
    }

   }
  }   

  // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >( mod.get() );

   if( tmod ){

    switch( tmod->type() ){
     case( BinaryKnapsackBlockMod::eChgWeight ): break;
     case( BinaryKnapsackBlockMod::eChgProfit ): break;
     case( BinaryKnapsackBlockMod::eFixX ): break;
     case( BinaryKnapsackBlockMod::eUnfixX ): break;         
    }

   }
  } 

  // BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >( mod.get() ); 

   if( tmod ){

    switch( tmod->type() ){
     case( BinaryKnapsackBlockMod::eChgWeight ): break;
     case( BinaryKnapsackBlockMod::eChgProfit ): break;
     case( BinaryKnapsackBlockMod::eFixX ): break;
     case( BinaryKnapsackBlockMod::eUnfixX ): break;         
    }

   }
  }

 } // end( while( ! v_mod.empty() ) )

return 0;

}// end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/


