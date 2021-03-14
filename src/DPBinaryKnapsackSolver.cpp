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

// process all the pending modifications and compute the first item from 
// which to restart the DP algorithm (start_item)

process_outstanding_Modification();

// check if the problem is empty - - - - - - - - - - - - - - - - - - - - - - -

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( BKB->is_empty() ){
  obj = - Inf< double >();
  return( kInfeasible );
 }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

/* The variable i is the first item processed by the DP algorithm. If compute
 is called for the first time then i = 0 and the algorithm is entirely 
 executed. Otherwise process_outstanding_Modification(), after processing all 
 the modifications, computes the index of the item from which to restart the 
 algorithm (if it is possible according to the reopt parameter).            */

int i = start_item;

int N = items.size();               // number of items (not preprocessed)

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

std::vector< double > currlab;      // set of current labels
std::vector< double > nextlab;      // set of next labels

std::vector< int >::iterator it;    // iterator over items[]

// Initialization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if( i == 0 ){                       // if i == 0 the DP algorithm is entirely
 maxcurrlab = 0;                    // executed, initialize dummy node in the                                        
 currlab.resize( maxcurrlab + 1 );  // origin
 currlab[ 0 ] = 0;

 it = items.begin();                // and start from the beginning   
}
else{                         // otherwise initialize with data stored in G[ i ]
 maxcurrlab = G[ i ].lab.size() - 1;
 currlab.resize( maxcurrlab + 1 );
 std::copy( G[ i ].lab.begin() , G[ i ].lab.end() , currlab.begin() );

 // and start from the correct item
 it = items.begin() + i;
}

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
/*
At each iteration the algorithm computes the maximum height of the constructed 
graph in order to allocate only the necessary memory to store the current 
labels.

According to the reopt parameter, labels can be stored in G[ i ].lab so that
they will be used for reoptimization.
                                                                            */

// DP Algorithm- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for( ; i < N ; i++ , it++ ){

  double p = v_P[ *it ];        // profit of the current item
  int w = v_W[ *it ];           // weight of the current item  

  // max "height" of the graph checking if the maximum capacity is reached
  maxnextlab = std::min( maxcurrlab + w , C );

  // Initialize labels 
  nextlab.resize( maxnextlab + 1 );
  std::fill( nextlab.begin() , nextlab.end() , - Inf< double >() );
  
  // Allocate predecessors
  G[ i + 1 ].pred.resize( maxnextlab + 1 );
  
  // Initialize the best label among those of this slice 
  double bestlab = -Inf< double >();

  for( int j = 0 ; j <= maxcurrlab ; j++ ){
   
   // if the node does not exists (currlab = -Inf) 
   // or if it does not have a better label continue
   if( currlab[ j ] <= bestlab )
    continue;

   // update bestlab
   bestlab = currlab[ j ];

   // horizontal arc
   if( currlab[ j ] > nextlab[ j ] ){
    G[ i + 1 ].pred[ j ] = false;
    nextlab[ j ] = currlab[ j ];
   }

   // check if the maximum capacity is reached
   if( j + w > C )
    continue;

   // diagonal arc 
   if( currlab[ j ] + p > nextlab[ j + w ] ){
    G[ i + 1 ].pred[ j + w ] = true;
    nextlab[ j + w ] = currlab[ j ] + p;
   }

  }

 // if reopt save currlab
 if( reopt ){
  G[ i ].lab.resize( maxcurrlab + 1 );
  std::copy( currlab.begin() , currlab.end() , G[ i ].lab.begin() );
 }

 maxcurrlab = maxnextlab;
 std::swap( currlab , nextlab );

}

// always save last labels
G[ N ].lab.resize( maxcurrlab + 1 );
std::copy( currlab.begin() , currlab.end() , G[ N ].lab.begin() );

// find optimal value- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

int maxh = std::min( maxcurrlab + 1 , C + 1 );
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

it = items.end() - 1;

for( int i = N ; i > 0 ; i-- , it-- ){

 if( G[ i ].pred[ besth ] ){
  v_x[ *it ] = 1;
  besth -= v_W[ *it ];
 }

}

// x = 1 - x for the items in ni
for( auto i : ni )
 v_x[ i ] = 1 - v_x[ i ];    

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

   // load weights and update prp_W (total negative weight) and prp_P
   prp_W = 0;
   prp_P = 0;

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

    // items with negative weights 
    if( v_W[ i ] <= 0 )
     prp_W += v_W[ i ]; 
    
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
   prp_P += v_P[ i ];
   continue;
  }

  // items not to be picked
  if( ( v_W[ i ] > C ) || ( ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ) ){
   v_x[ i ] = 0;
   continue; 
  }

  // if the item has negative weight (and negative profit) update prp_P and
  // change the sign of v_W[ i ] and v_P[ i ]
  if( v_W[ i ] <= 0 ){
   prp_P += v_P[ i ];

   v_W[ i ] = - v_W[ i ];
   v_P[ i ] = - v_P[ i ];

   ni.push_back( i ); 
  }

  // store the indeces of the items that are involved in the DP algorithm

  items.push_back( i );
  
  v_x[ i ] = 0;                         // initialize all variables to 0

 }

 start_item = 0;                        // starting item 

 obj = -Inf< double >();                // objective value

 } // end if( f_Block )

} // end( DPBinaryKnapsackSolver::load() )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::process_outstanding_Modification(){

 if( v_mod.empty() )
  return;

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );

 // item is the index from which to restart the DP algorithm
 int item = Inf< int >();

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

start_item = 0;

}// end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/


