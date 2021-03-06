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
 
 //process_outstanding_Modification();

// check if the problem is empty - - - - - - - - - - - - - - - - - - - - - - -

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( BKB->is_empty() ){
  obj = - Inf< double >();
  return( kInfeasible );
 }

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
/*
The algorithm processes first the items with negative weight and second the 
items with positive weights.
Three phases are performed:

 - First phase: the items to be processed are the ones with negative weight  

 - Second phase: the constructed graph has not yet reached the maximum 
                  capacity. Avoid to check if the maximum capacity is reached
                  and avoid to allocate unnecessary memory for the matrix of 
                  predecessors  

 - Third phase: the graph has already reached the maximum capacity and
                the items to be processed are the ones with positive weight 
                                                                            */
                          
// initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int N = pi.size() + ni.size();  // number of items (not preprocessed)

if( N == 0 ){                   // if N == 0 all the items have been
 obj = prp_P;                   // preprocessed, no need of DP
 return( kOK );                 // return
}

int C = f_C - prp_W;            // new capacity after the preprocessing
pred.resize( N );               // Matrix of predecessors

// the only two necessary sets of labels 
std::vector< double > currlab( nW + 1 , -Inf< double >() );
std::vector< double > nextlab( nW + 1 );

currlab[ nW ] = 0;              // start DP from nW

// index over all items 
int i = 0;

auto nii = ni.begin();      // iterator over the items with negative weight
int nSize = ni.size();      // number of items with negative weight


// first phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

int mincurrlab = nW;        // minimum label

for( ; i < nSize ; i++ , nii++ ){

 pred[ i ].resize( nW + 1 );
  
 // Initialize nextlab
 std::fill( nextlab.begin() , nextlab.end() , -Inf< double >() );

 double p = v_P[ *nii ];        // profit of the current item
 int w = v_W[ *nii ];           // weight of the current item

 for( int j = mincurrlab ; j <= nW ; j++ ){

  // check if the node exists 
  if( std::isinf( currlab[ j ] ) ) 
   continue;

  // horizontal arc
  if( currlab[ j ] > nextlab[ j ] ){
   pred[ i ][ j ] = 0;
   nextlab[ j ] = currlab[ j ]; 
  }

  // diagonal arc 
  if( currlab[ j ] + p > nextlab[ j + w ] ){
   pred[ i ][ j + w ] = 1;
   nextlab[ j + w ] = currlab[ j ] + p;
  }

 }
  
 std::swap( currlab , nextlab );
 mincurrlab += w;
}
// end first phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

auto pii = pi.begin();      // iterator over the items with positive weight

int maxcurrlab = nW;        // maximum label

currlab.resize( C + nW + 1 );
nextlab.resize( C + nW + 1 );

// second phase- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for( ; i < N ; i++ , pii++ ){
 
 double p = v_P[ *pii ];      // profit of the current item
 int w = v_W[ *pii ];         // weight of the current item

 int maxnextlab = maxcurrlab + w;
 if( maxnextlab > C + nW )
  break;

 pred[ i ].resize( maxnextlab + 1 );

 // Initialize nextlab
 std::fill( nextlab.begin() , nextlab.end() , -Inf< double >() ); 

 double bestlab = -Inf< double >();

 for( int j = 0 ; j <= maxcurrlab ; j++ ){

  // check if the node has a better label and if it exists 
  if( currlab[ j ] <= bestlab ) 
   continue;

  bestlab = currlab[ j ];

  // horizontal arc
  if( currlab[ j ] > nextlab[ j ] ){
   pred[ i ][ j ] = 0;
   nextlab[ j ] = currlab[ j ]; 
  }

  // diagonal arc 
  if( currlab[ j ] + p > nextlab[ j + w ] ){
   pred[ i ][ j + w ] = 1;
   nextlab[ j + w ] = currlab[ j ] + p;
  }

 }
  
 std::swap( currlab , nextlab );
 maxcurrlab = maxnextlab;
}
// end second phase- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// this part has not been initialized before
std::fill( currlab.begin() + maxcurrlab + 1 , currlab.end() , -Inf< double >() );

// third phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for( ; i < N ; i++ , pii++ ){

 pred[ i ].resize( C + nW + 1 );
 // Initialize nextlab
 std::fill( nextlab.begin() , nextlab.end() , -Inf< double >() );

 double p = v_P[ *pii ];      // profit of the current item
 int w = v_W[ *pii ];         // weight of the current item

 double bestlab = -Inf< double >();

 for( int j = 0 ; j <= C + nW ; j++ ){

  // check if the node has a better label and if it exists 
  if( currlab[ j ] <= bestlab ) 
   continue;

  bestlab = currlab[ j ];

  // horizontal arc
  if( currlab[ j ] > nextlab[ j ] ){
   pred[ i ][ j ] = 0;
   nextlab[ j ] = currlab[ j ]; 
  }

  // check if the maximum capacity is reached
  if( j + w > C + nW )
   continue;

  // diagonal arc 
  if( currlab[ j ] + p > nextlab[ j + w ] ){
   pred[ i ][ j + w ] = 1;
   nextlab[ j + w ] = currlab[ j ] + p;
  }

 }
  
 std::swap( currlab , nextlab );
}
// end third phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// find optimal value- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

labels.resize( currlab.size() );  // save currlab in labels 

int bestj;

// initialize objective
obj = -Inf< double >();

for( int i = 0 ; i <= C + nW ; i++ ){
 if( currlab[ i ] > obj ){
  obj = currlab[ i ];
  bestj = i;
 }
 labels[ i ] = currlab[ i ];
}

obj +=  prp_P;        // add the profit coming from the preprocessing 

// reconstruct the optimal solution - - - - - - - - - - - - - - - - - - - - - 

pii = pi.end() - 1;

for( i = N - 1 ; i >= nSize ; i-- , pii-- ){  // positive part
 if( pred[ i ][ bestj ] ){
  v_x[ *pii ] = 1;
  bestj -= v_W[ *pii ];
 }
}

nii = ni.end() - 1;

for( ; i >= 0 ; i-- , nii-- ){                // negative part
 if( pred[ i ][ bestj ] ){
  v_x[ *nii ] = 1;
  bestj -= v_W[ *nii ];
 }
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

    // if the variables is fixed update prp_P and prp_W
    if( BKB->is_fixed( i ) ){
     if( BKB->get_x( i ) ){
      prp_P += v_P[ i ];
      prp_W += v_W[ i ];
     }
     continue; 
    }

    // these are the items with negative weights involved in the DP algorithm
    if( ( v_W[ i ] < 0 ) && ( v_P[ i ] < 0 ) ){
     nW -= v_W[ i ];
    }
    
    // these are the items with positive weights involved in the DP algorithm
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
 from the computation.                                                     */
 
 int C = f_C - prp_W;
 int count_pos = 0;
 int count_neg = 0;

 v_x.resize( f_NItems );
 pi.resize( f_NItems );
 ni.resize( f_NItems );

 auto pii = pi.begin();
 auto nii = ni.begin();

 for( int i = 0 ; i < f_NItems ; i++ ){

  // if the variable is fixed then update v_x and discard it
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

  if( v_W[ i ] > 0 ){    // items with positive weights
   *pii++ = i;
   count_pos++; 
  }
  else{                  // items with negative weights
   *nii++ = i;
   count_neg++; 
  }

  v_x[ i ] = 0;          // initialize all variables to 0

 }

 pi.resize( count_pos );
 ni.resize( count_neg );

 } // end if( f_Block )

} // end load()

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::process_outstanding_Modification(){

} 

/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
