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

using Index = Block::Index;

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
 
 int status = 0;

 process_outstanding_Modification();

 // PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 //
 // The following items are preprocessed:
 // 
 //  - Items i with negative weight and positive profit => set x[ i ] = 1
 //
 //  - Items i with positive weight and negative profit => set x[ i ] = 0
 //
 //  - Items with a weight greater than C - \sum w_ , where w_ are the 
 //    negative weights  
 //
 // Actually the variables corresponding to all the remaining items are also
 // initialized to 0, but the second set of preprocessed items are discarded
 // from the computation.   
 //
 // Since the third condition requires to know in advance the sum of all the
 // negative weights, two for loops are performed. 
 // 
 // For the items with negative weights that are picked and so discarded from
 // the computation, the Capacity is updated accordingly.
 // The variable nW instead will contain the sum of all negative weights (in 
 // absolute value) of the items that will be involved in the computation.


 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_x.resize( f_NItems );  // vector of solution
 int C = f_C;             // Capacity
 int N = f_NItems;        // Number of items involved in the computation
 double P = 0;            // Profit
 int count_neg = 0;       // Number of items with negative weight
 int nW = 0;              // Sum of negative weights

 // First for loop: preprocess only items to be picked
 for( Index i = 0 ; i < f_NItems ; i++ ){
  
  if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){ 
   v_x[ i ] = 1;   // pick the item 
   C -= v_W[ i ];  // update Capacity 
   P += v_P[ i ];  // update Profit
   N--;            // update number of Items
   continue;
  
  }

  // otherwise count items with negative weights and update nW  
  if( v_W[ i ] < 0 ){
   count_neg++;
   nW -= v_W[ i ];
  }

  v_x[ i ] = 0;        // and initialize the variable to 0 
}

 // items vector stores the indeces of the items not discarded
 std::vector< Index > items( N );

 // the vector is filled in reverse: in the last part are stored the items
 // with negative weight; since they are exacty count_neg, the positive items
 // are stored starting from .end() - 1 - count_neg
 auto pii = items.end() - 1 - count_neg;
 auto nii = items.end() - 1;

 // Second for loop: preprocess items not to be picked. Althouth the variable
 // is already 0, it is necessary to discard them from the computation
 for( Index i = 0 ; i < f_NItems ; i++ ){

  // already preprocessed items
  if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) )
    continue; 

  // items still to be preprocessed 
  if( ( v_W[ i ] > C + nW ) || ( ( v_W[ i ] >= 0 ) && ( v_P[ i ] < 0 ) ) ){
   N--;         // variable already 0 but discard 
   continue;    // the item from the computation
  }
  
  // store indeces of items that are not discarded
  if( v_W[ i ] < 0 )
   *nii-- = i;
  else
   *pii-- = i;
 }

// resize items vector because some of the items with positive weights
// could have been discarded in the second for loop
 items.erase( items.begin() , items.begin() + items.size() - N );

// END PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//
// The algorithm processes first the items with positive weight and second the 
// items with negative weights (as they are positioned in the items vector) 
// Three phases are performed:
//
//  - First phase: the constructed graph has not yet reached the maximum 
//                 capacity. Avoid to check if the maximum capacity is reached
//                 and avoid to allocate unnecessary memory for the matrix of 
//                 predecessors  
//
//  - Second phase: the graph has already reached the maximum capacity but 
//                  the items to be processed are still the ones with potive 
//                  weight
//
//  - Third phase: the items to be processed are the ones with negative weight


// Matrix of predecessors
std::vector< std::vector< bool > > pred( N );  // the last is for obj

// the only two necessary sets of labels 
std::vector< double > currlab( C + nW + 1 );
std::vector< double > nextlab( C + nW + 1 );

currlab[ 0 ] = 0;
int maxcurrlab = 0;   // maximum meaningful label

Index i = 0;

// First phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for( ; i < N - count_neg ; i++ ){
 int w = v_W[ items[ i ] ];    // weight of the current item
 int maxnextlab = maxcurrlab + w;
 if( maxnextlab > C + nW )
  break;

 pred[ i ].resize( maxnextlab + 1 );
 
 // Initialize nextlab
 std::fill( nextlab.begin() , nextlab.begin() + maxnextlab + 1 , 
            -Inf< double >() );
 
 double p = v_P[ items[ i ] ]; // profit of the current item

 for( Index j = 0 ; j <= maxcurrlab ; j++ ){

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
 maxcurrlab = maxnextlab;
} 
// end first phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// this part of the label vector (if any) may have never been
// initialized before, but it is checked in the next part
std::fill( currlab.begin() + maxcurrlab + 1  , 
           currlab.end() , -Inf< double >() );

// Second phase- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for( ; i < N - count_neg ; i++ ){
 
 pred[ i ].resize( C + nW  + 1 );

 // Initialize nextlab
 std::fill( nextlab.begin() , nextlab.end() , -Inf< double >() );
 
 int w = v_W[ items[ i ] ];     // weight of the current item
 double p = v_P[ items[ i ] ];  // profit of the current item

 for( int j = 0 ; j <= C + nW ; j++ ){

  // check if the node exists
  if( std::isinf( currlab[ j ] ) ) 
   continue; 

  // horizontal arc
  if( currlab[ j ] > nextlab[ j ] ){
   pred[ i ][ j ] = 0;
   nextlab[ j ] = currlab[ j ];  
  }

  // check if the maximum capacity is reached
  if( j + w > C + nW )
   continue;

  // diagonal arc
  if( currlab[ j ] + p  > nextlab[ j + w ] ){
   pred[ i ][ j + w ] = 1;
   nextlab[ j + w ] = currlab[ j ]  + p; 
  }
 }

 std::swap( currlab , nextlab );
}
// end second phase  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// Third phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
for( ; i < N ; i++ ){
 
 pred[ i ].resize( C + nW  + 1 );

 // Initialize nextlab
 std::fill( nextlab.begin() , nextlab.end() , -Inf< double >() );
 
 int w = v_W[ items[ i ] ];     // weight of the current item
 double p = v_P[ items[ i ] ];  // profit of the current item

 for( int j = 0 ; j <= C + nW ; j++ ){

  // check if the node exists
  if( std::isinf( currlab[ j ] ) ) 
   continue; 

  // horizontal arc
  if( currlab[ j ] > nextlab[ j ] ){
   pred[ i ][ j ] = 0;
   nextlab[ j ] = currlab[ j ];  
  }

  // check if a value < 0 is reached
  if( j + w < 0 )
   continue;

  // diagonal arc
  if( currlab[ j ] + p  > nextlab[ j + w ] ){
   pred[ i ][ j + w ] = 1;
   nextlab[ j + w ] = currlab[ j ]  + p; 
  }
 }

 std::swap( currlab , nextlab );
}
// end third phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// find optimal value- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int bestj;

for( int i = 0 ; i <= C ; i++ ){
 if( currlab[ i ] > obj ){
  obj = currlab[ i ];
  bestj = i;
  }
 }

// add the profit coming from the preprocessing 
obj +=  P;

// reconstruct the optimal solution- - - - - - - - - - - - - - - - - - - - - -
for( int i = N - 1 ; i >= 0 ; i-- ){
 if( pred[ i ][ bestj ] ){
  v_x[ items[ i ] ] = 1;
  bestj -= v_W[ items[ i ] ];
 }
}

 return( status );
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
   f_NItems = BKB->get_NItems();
   f_C = static_cast< int >( BKB->get_Capacity() );

   // get weights and cast to int
   v_W.resize( f_NItems );
   const std::vector< double > & W = BKB->get_Weights();

   std::transform( W.begin() , W.end() , v_W.begin() , 
        []( double w ) -> int { return static_cast< int >(w); } );

   // get profits
   v_P.resize( f_NItems );
   const std::vector< double > & P = BKB->get_Profits();
   std::copy( P.begin() , P.end() , v_P.begin() );

   BKB->read_unlock();

  }
 
 }

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::process_outstanding_Modification(){

} 

/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
