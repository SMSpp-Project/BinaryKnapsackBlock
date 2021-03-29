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

using idx_type = ThinComputeInterface::idx_type;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset; 

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

void mergeRange( Range & rng1 , c_Range rng2 ){

 // Merge two ranges rng1 and rng2 and put the result in rng1 

 if( rng2.first < rng1.first )       
   rng1.first = rng2.first;       
                      
 if( rng2.second > rng1.second )
  rng1.second = rng2.second;  

}

/*--------------------------------------------------------------------------*/

void mergeSubset( Subset & nms1 , c_Subset & nms2 ){

 // Merge two Subsets nms1 and nms2 and put the result in nms1

 Subset dest;               

 std::set_union( nms1.begin(), nms1.end(),
                 nms2.begin(), nms2.end(),                 
                 std::back_inserter(dest) );

 std::swap( dest , nms1 );
}

/*--------------------------------------------------------------------------*/

void SbstRngdiff( Subset & nms , c_Range rng ){

 // given a Subset nms and a Range rng, remove from nms the elements in rng 

 Subset result;

 for( auto i : nms ){
  if( i < rng.first || i >= rng.second )
   result.push_back( i ); 
 }

 std::swap( result , nms );

}


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

auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 
if( BKB->is_empty() ){
 obj = - Inf< double >();
 return( kInfeasible );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

/* start_item is the first item processed by the DP algorithm. If compute is
 called for the first time then the algorithm is entirely executed.
 Otherwise process_outstanding_Modification(), after processing all the
 modifications, computes the index of the item from which to restart the 
 algorithm.                                                                 */

if( start_item == + Inf< int >() )  // if start_item == Inf it is assumed that 
 return( kOK );                     // everything has already been done in 
                                    // process_outstanding_Modification()

int i = start_item;                 // otherwise start from start_item

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int C = f_C;                // compute new capacity after the preprocessing
double prp_P = 0;           // compute the profit coming from the preprocessing
                
for( int i = 0 ; i < f_NItems ; i++ ){

 if( items[ i ] == 1 || items[ i ] == - 1 ){    // if i has been selected
  C -= v_W[ i ];                                // or if it has negative
  prp_P += v_P[ i ];                            // weight and profit
 }

}

int maxcurrlab;                     // max current height of the graph
int maxnextlab;                     // max next height of the graph    

std::vector< double > currlab;      // set of current labels
std::vector< double > nextlab;      // set of next labels

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Initialize maxcurrlab and currlab with data stored in G[ i ].lab

maxcurrlab = G[ i ].lab.size() - 1;         
maxcurrlab = std::min( maxcurrlab , C );    // Capacity may have been changed

currlab.resize( maxcurrlab + 1 ); 
std::copy( G[ i ].lab.begin() , G[ i ].lab.begin() + maxcurrlab + 1 , 
           currlab.begin() );

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
/*
At each iteration the algorithm computes the maximum height of the constructed 
graph in order to allocate only the necessary memory to store the current 
labels.

According to the reopt parameter, labels can be stored in G[ i ].lab so that
they will be used for reoptimization.
                                                                            */

// DP Algorithm- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for( ; i < f_NItems ; i++ ){

  if( items[ i ] == 0 || items[ i ] == 1 )      // if i has been preprocessed
    continue;                                   // continue

  double p = items[ i ] != -1 ? v_P[ i ] : - v_P[ i ];  // profit of i   
  int w = items[ i ] != -1 ? v_W[ i ] : - v_W[ i ];     // weight of i                     

  // max height of the graph checking if the maximum capacity is reached
  maxnextlab = std::min( maxcurrlab + w , C );

  // Initialize labels 
  nextlab.resize( maxnextlab + 1 );
  std::fill( nextlab.begin() , nextlab.end() , - Inf< double >() );

  // Allocate predecessors
  G[ i + 1 ].pred.resize( maxnextlab + 1 );

  // Initialize the best label among those of this slice 
  double bestlab = -Inf< double >();

  for( int j = 0 ; j <= maxcurrlab ; j++ ){
                                                
   if( currlab[ j ] <= bestlab )    // if the node does not exists 
    continue;                       // (currlab = -Inf) or if it does  
                                    // not have a better label continue                 
 
   bestlab = currlab[ j ];          // update bestlab

   if( currlab[ j ] > nextlab[ j ] ){           // horizontal arc
    G[ i + 1 ].pred[ j ] = false;               
    nextlab[ j ] = currlab[ j ];
   }

   if( j + w > C )                              // check if the maximum 
    continue;                                   // capacity is reached

   if( currlab[ j ] + p > nextlab[ j + w ] ){   // diagonal arc 
    G[ i + 1 ].pred[ j + w ] = true;            
    nextlab[ j + w ] = currlab[ j ] + p;
   }

  }

 // reoptimization- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( i % step == 0 ) && ( i != start_item ) ){ 
  G[ i ].lab.resize( maxcurrlab + 1 ); 
  std::copy( currlab.begin() , currlab.end() , G[ i ].lab.begin() );
 }

 // swaps
 maxcurrlab = maxnextlab;
 std::swap( currlab , nextlab );

} // end( DP algorithm )- - - - - - - - - - - - - - - - - - - - - - - - - - - 

// always save last labels in G[ f_NItems ].lab

G[ f_NItems ].lab.resize( maxcurrlab + 1 );
std::copy( currlab.begin() , currlab.end() , G[ f_NItems ].lab.begin() );

// find optimal value - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

int maxh = std::min( maxcurrlab , C );
int besth;

obj = - Inf< double >();                    // Initialize objective value

for( int i = 0 ; i <= maxh ; i++ ){ 
 if( G[ f_NItems ].lab[ i ] > obj ){
  obj = G[ f_NItems ].lab[ i ]; 
  besth = i;
 }
}

G[ f_NItems ].lab.resize( besth + 1 ); 

obj += prp_P;             // add the profit coming from the preprocessing

start_item = + Inf< int >();

v_x.clear();              // clear previous solution (if any)

return( kOK );

}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::get_var_solution( Configuration * solc ){

 // check if the BinaryKnapsackBlock uses the correct abstract representation
 
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 
 if( BKB->get_NItems() != BKB->get_VarSize() )
  throw( std::invalid_argument( "Incompatible variables size" ) );

 if( start_item != + Inf< int >() )
  throw( std::invalid_argument( "compute() must be called first" ) );

 // reconstruct the optimal solution - - - - - - - - - - - - - - - - - - - -

 compute_var_solution();

 // write the current solution in the BinaryKnapsackBlock

 BKB->set_x( v_x );

}
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::set_par( idx_type par , double value ){
 
 switch( par ){
 
  case( dblReopt ):{

   if( value < 0 || value > 1 )
    throw(std::invalid_argument("dblReopt parameter must be in [ 0 , 1 ]"));

   /* The variable step defines how often the labels should be saved. 
    Each time i % step == 0 currlab is saved in G[ i ].lab. */
   reopt = value;
   
   if( reopt < 1e-06 )
    step = + Inf< int >();  // never save currlab
   else
    step = reopt > 0.5 ? 1 : std::floor( 1 / reopt );
   
   start_item = 0;          // restart from 0 in the next call of compute()
   break;
  }  
 
 } // end( switch( par ) )

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
  load();                   // the Binary Knapsack instance must be re-loaded
  v_mod.clear();
 }
 else
  v_mod.push_back( mod );

 f_mod_lock.clear( std::memory_order_release );  // release lock
 
} // end( add_Modification() )

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

   v_P.resize( f_NItems );                  
   const auto & P = BKB->get_Profits();     // get profits and if the sense
   for( int i = 0 ; i < P.size() ; i++ )    // is minimization change the sign 
    v_P[ i ] = f_sense ? P[ i ] : - P[ i ]; // of the profits 

   v_W.resize( f_NItems );                  // prepare vector of weights
   const auto & W = BKB->get_Weights();

   for( int i = 0 ; i < W.size() ; i++ ){   // load weights and check 
    v_W[ i ] = std::round( W[ i ] );        // that they are integers

    if( std::abs( v_W[ i ] - W[ i ] ) > 1e-06 )
     throw( std::invalid_argument( "Weights must be integers" ) );

   }

   BKB->read_unlock();

  // end load Binary Knapsack instance- - - - - - - - - - - - - - - - - - - -
 
  items.resize( f_NItems );           // resize vector of items
 
  G.resize( f_NItems + 1 );           // resize Graph
  
  preprocessing( std::make_pair( 0 , f_NItems ) );    // preprocessing

  start_item = 0;

  obj = -Inf< double >();

  v_x.clear();                        // clear previous solution (if any)

 } // end if( f_Block )

} // end( DPBinaryKnapsackSolver::load() )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::compute_var_solution(){

if( ! v_x.empty() )       // if the solution has already been computed
 return;                  // return

// reconstruct the optimal solution- - - - - - - - - - - - - - - - - - - - - -

 int besth = G[ f_NItems ].lab.size() - 1;

 v_x.resize( f_NItems );

 for( int i = f_NItems - 1 ; i >= 0 ; i-- ){

  if( items[ i ] == 0 || items[ i ] == 1 ){             // preprocessed items
   v_x[ i ] = items[ i ];
   continue;
  }

  if( G[ i + 1 ].pred[ besth ] ){
   v_x[ i ] = 1;
   besth -= items[ i ] != -1 ? v_W[ i ] : - v_W[ i ];  
  }
  else
   v_x[ i ] = 0;

  if( items[ i ] == -1 )            // items with negative weight and profit
   v_x[ i ] = ! v_x[ i ];

 }
  
}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::preprocessing( Range rng ){

 // PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 /*
 The following items are preprocessed:
  
  - Items i with negative weight and positive profit => set x[ i ] = 1
  - Items i with positive weight and negative profit => set x[ i ] = 0
  - Fixed items are preprocessed according to their current value

  The vector items is filled as follows:

   - items[ i ] = -1    if i has negative weight and negative profit
   - items[ i ] = 0     if i is preprocessed and x[ i ] = 0
   - items[ i ] = 1     if i is preprocessed and x[ i ] = 1
   - items[ i ] = 2     otherwise

 The DP algorithm discards the preprocessed items and changes the sign of the
 weights and the profits of the items for which items[ i ] = -1
                                                                           */
  
 if( rng.second <= rng.first )    // nothing to do
  return;   

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( int i = rng.first ; i < rng.second ; i++ ){
   
  if( BKB->is_fixed( i ) ){                         // fixed items
   items[ i ] = BKB->get_x( i );
   continue;
  }

  if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){     // items to select
   items[ i ] = 1;
   continue;
  }

  if(  ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ){    // items not to select
   items[ i ] = 0;
   continue; 
  }

  if( v_W[ i ] <= 0 ){              // if the item has negative weight 
   items[ i ] = -1;                 // and negative profit set items to -1
   continue; 
  }

  items[ i ] = 2;                   // otherwise
   
 }

} // end( preprocessing( Range ) )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::preprocessing( Subset & nms ){

 // PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 /*
 The following items are preprocessed:
  
  - Items i with negative weight and positive profit => set x[ i ] = 1
  - Items i with positive weight and negative profit => set x[ i ] = 0
  - Fixed items are preprocessed according to their current value

  The vector items is filled as follows:

   - items[ i ] = -1    if i has negative weight and negative profit
   - items[ i ] = 0     if i is preprocessed and x[ i ] = 0
   - items[ i ] = 1     if i is preprocessed and x[ i ] = 1
   - items[ i ] = 2     otherwise

 The DP algorithm discards the preprocessed items and changes the sign of the
 weights and the profits of the items for which items[ i ] = -1
                                                                           */

 if( nms.empty() )    // nothing to do
  return;

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  for( auto i : nms ){
 
   if( BKB->is_fixed( i ) ){                        // fixed items
    items[ i ] = BKB->get_x( i );
    continue;
   }

   if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){    // items to select
    items[ i ] = 1;
    continue;
   }

   if(  ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ){   // items not to select
    items[ i ] = 0;
    continue; 
   }

   if( v_W[ i ] <= 0 ){             // if the item has negative weight 
    items[ i ] = -1;                // and negative profit set items to -1
    continue; 
   }

   items[ i ] = 2;                  // otherwise
   
  }

 } // end( preprocessing( Subset ) )


/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::process_outstanding_Modification(){

/* start_item != + Inf means that: 
 
 - compute is called for the first time
 - or a Nuclear modification has been issued 
 - or the reopt parameter has been modified

 In all cases it is not possible to reoptimize
                                                                            */

if( start_item != + Inf< int >() ){    
 if( ! v_mod.empty() ){             // if there are any modifications 
  load();                           // re-load the instance
  v_mod.clear();                    // and clear v_mod  
 }                                   
 return;                            // in any case return                              
}                                    

// otherwise copy v_mod in a temporary list of modifications - - - - - - - - -

 Lst_sp_Mod v_mod_tmp;              // temporary list of modifications

 // try to acquire lock, spin on failure

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );       // copy v_mod in v_mod_tmp

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 /* Since different type of Modifications require to redo the same preprocessing
  phase, all the modified items are stored in modRng_items and modSbst_items
  and the preprocessing is done only once at the end.                        */


 for( auto mod : v_mod_tmp ){

  // BinaryKnapsackBlockMod - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod.get() );
   if( tmod ){ 
    switch( tmod->type() ){

     case( BinaryKnapsackBlockMod::eChgCapacity ):  // Change Capacity
      capacity_Modification(); break;                      

     case( BinaryKnapsackBlockMod::eChgSense ):     // Change the sense of the     
      sense_Modification(); break;                  // objective                     
                                                    
    }                                               
   }                                                

  // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >( mod.get() );

   if( tmod ){ 
    switch( tmod->type() ){
    
     case( BinaryKnapsackBlockMod::eChgProfit ):    // change Profits       
      profit_Modification( tmod->rng() ); break; 

     case( BinaryKnapsackBlockMod::eFixX ):         // fix x        
      fixX_Modification( tmod->rng() ); break; 

     case( BinaryKnapsackBlockMod::eUnfixX ):       // unfix x    
      unFixX_Modification( tmod->rng() ); break; 

     case( BinaryKnapsackBlockMod::eChgWeight ):    // change Weights
      weight_Modification( tmod->rng() ); break;  
    
    }
   }
  } 

  // BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >( mod.get() ); 

   if( tmod ){ 
    switch( tmod->type() ){
     
     case( BinaryKnapsackBlockMod::eChgProfit ):    // change Profits        
      profit_Modification( tmod->nms() ); break; 

     case( BinaryKnapsackBlockMod::eFixX ):         // fix x              
      fixX_Modification( tmod->nms() ); break; 

     case( BinaryKnapsackBlockMod::eUnfixX ):       // unfix x       
      unFixX_Modification( tmod->nms() ); break; 

     case( BinaryKnapsackBlockMod::eChgWeight ):    // change Weights
      weight_Modification( tmod->nms() ); break; 
    
    }
   }
  }

 } // end( for(  ) )

v_mod_tmp.clear();              // clear the temporary list of modification

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// compute start_item- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

/* Each modification updates start_item with the first modified item. At this
 point it is necessary to retrieve the first item smaller than start_item whose 
 labels have been previously stored in G[ i ].lab.

 The labels corresponding to i have been stored in G[ i ].lab if i % step == 0 
 and i wasn't preprocessed. Therefore this computation must be done before 
 ro-doing the preprocessing.                                                 */

 if( start_item != + Inf< int >() ){
  
  int i = ( start_item / step ) * step;

  while( i > 0 && ( i % step != 0 || items[ i ] == 0 || items[ i ] == 1 ) ) 
   i--;
  
  start_item = i;
 }

// re-do preprocessing - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

preprocessing( modRng_items );

SbstRngdiff( modSbst_items , modRng_items ); 

preprocessing( modSbst_items );

// re-initialize - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

modRng_items.first = + Inf< int >();
modRng_items.second = 0;

modSbst_items.clear();

}// end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::capacity_Modification(){
 
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 int nC = std::floor( BKB->get_Capacity() );     // get new Capacity

 /* If the new capacity is smaller than the previous one, restart the DP
  algorithm from f_NItems, i.e. only recompute the optimal value. Otherwise
  restart from the beginning.                                               */

 start_item = nC > f_C ? 0 : std::min( f_NItems , start_item );

 f_C = nC;                    // update the Capacity

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::sense_Modification(){ 

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 
 f_sense = BKB->get_objective_sense();  // update f_sense

 for( int i = 0 ; i < f_NItems ; i++ )  // change the sign of all profits
  v_P[ i ] = - v_P[ i ];

 modRng_items.first = 0;                // preprocessing must be re-done
 modRng_items.second = f_NItems;

 start_item = 0;                        // restart from the beginning 

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::fixX_Modification( Range rng ){

 mergeRange( modRng_items , rng ); 

 start_item = std::min( start_item , int( rng.first ) ); 

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::fixX_Modification( c_Subset & nms ){

 mergeSubset( modSbst_items , nms );

 start_item = std::min( start_item , int( nms[ 0 ] ) );  

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::unFixX_Modification( Range rng ){
 
 mergeRange( modRng_items , rng );

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( int i = rng.first ; i < rng.second ; i++ ){

  if( v_W[ i ] > 0 && items[ i ] == 1 )
   start_item = 0;

  if( v_W[ i ] < 0 && items[ i ] == 0 )
   start_item = 0;      
    
 }

 start_item = std::min( start_item , int( rng.first ) );

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::unFixX_Modification( c_Subset & nms ){

 mergeSubset( modSbst_items , nms );

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( auto i : nms ){

  if( v_W[ i ] > 0 && items[ i ] == 1 )
   start_item = 0;

  if( v_W[ i ] < 0 && items[ i ] == 0 )
   start_item = 0;      
    
 }


 start_item = std::min( start_item , int( nms[ 0 ] ) );     

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::weight_Modification( Range rng ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( int i = rng.first ; i < rng.second ; i++ ){
  
  int nw = std::round( BKB->get_Weight( i ) ); 
  
  if( std::abs( nw - BKB->get_Weight( i ) ) > 1e06 )
   throw( std::invalid_argument( "Weights must be integers!" ) );
   
  if( nw < 0 && nw < v_W[ i ] )         // it is not possible to re-optimize
   start_item = 0;  

  v_W[ i ] = nw;          // update weight

 }

 mergeRange( modRng_items , rng );

 start_item = std::min( start_item , int( rng.first ) );

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::weight_Modification( c_Subset & nms ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( auto i : nms ){
  
  int nw = std::round( BKB->get_Weight( i ) );
  
  if( std::abs( nw - BKB->get_Weight( i ) ) > 1e06 )
   throw( std::invalid_argument( "Weights must be integers" ) );
   
  if( nw < 0 && nw < v_W[ i ] )         // it is not possible to re-optimize
   start_item = 0;    

  v_W[ i ] = nw;

 }

 mergeSubset( modSbst_items , nms ) ;

 start_item = std::min( start_item , int( nms[ 0 ] ) );   

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::profit_Modification( Range rng ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( int i = rng.first ; i < rng.second ; i++ )
  v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );

 mergeRange( modRng_items , rng );

 start_item = std::min( start_item , int( rng.first ) );

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::profit_Modification( c_Subset & nms ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( auto i : nms )
  v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );
 
 mergeSubset( modSbst_items , nms );

 start_item = std::min( start_item , int( nms[ 0 ] ) );  

}


/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/

















