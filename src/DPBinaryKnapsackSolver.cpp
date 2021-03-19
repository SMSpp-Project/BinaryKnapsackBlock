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
 called for the first time then i = 0 and the algorithm is entirely executed.
 Otherwise process_outstanding_Modification(), after processing all the
 modifications, computes the index of the item from which to restart the 
 algorithm (if it is possible according to the reopt parameter).            */

if( start_item == Inf< int >() )    // if start_item == Inf it is assumed that 
 return( kOK );                     // everything has already been done in 
                                    // process_outstanding_Modification()

int i = start_item;                 // otherwise start from start_item

// reoptimization- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/* The variable step defines how often the labels should be saved. 
 Each time i % step == 0 currlab is saved in G[ i ].lab. The frequency depends 
 on the reopt parameter.                                                     */ 

int step = f_NItems;

if( reopt != 0 )
 step = reopt > 0.5 ? 1 : std::floor( 1 / reopt );

while( i % step != 0 )              // start from the first item whose
 i--;                               // corresponding label is stored in G[ i ]
                                    // or from 0

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int C = f_C;                // compute new capacity after the preprocessing
double prp_P = 0;           // compute the profit coming from the preprocessing
                
for( int i = 0 ; i < f_NItems ; i++ ){

 if( items[ i ] == 1 ){             // if i has been selected
  C -= v_W[ i ];
  prp_P += v_P[ i ];
 }

 if( items[ i ] == -1 ){            // if i has negative weight and profit
  C -= v_W[ i ];
  prp_P += v_P[ i ];

  v_W[ i ] = - v_W[ i ];            // change the sign of the items with
  v_P[ i ] = - v_P[ i ];            // negative weight and negative profit
 }

} 

int maxcurrlab;                     // max current height of the graph
int maxnextlab;                     // max next height of the graph    

std::vector< double > currlab;      // set of current labels
std::vector< double > nextlab;      // set of next labels

// Initialization- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

maxcurrlab = G[ i ].lab.size() - 1;    
currlab.resize( maxcurrlab + 1 );
std::copy( G[ i ].lab.begin() , G[ i ].lab.end() , currlab.begin() );


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

  double p = v_P[ i ];        // profit of the current item
  int w = v_W[ i ];           // weight of the current item

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

 // reoptimization- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( i % step == 0 ) && ( i != start_item ) ){
  G[ i ].lab.resize( maxcurrlab + 1 );
  std::copy( currlab.begin() , currlab.end() , G[ i ].lab.begin() );
 }

 // swap
 maxcurrlab = maxnextlab;
 std::swap( currlab , nextlab );

} // end( DP algorithm )- - - - - - - - - - - - - - - - - - - - - - - - - - - 

// always save last labels in G[ f_NItems ].lab

G[ f_NItems ].lab.resize( maxcurrlab + 1 );
std::copy( currlab.begin() , currlab.end() , G[ f_NItems ].lab.begin() );

// find optimal value - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

int maxh = std::min( maxcurrlab , C );

obj = - Inf< double >();                    // Initialize objective value

for( int i = 0 ; i <= maxh ; i++ ){ 
 if( G[ f_NItems ].lab[ i ] > obj ){
  obj = G[ f_NItems ].lab[ i ]; 
  besth = i;
 }
}

obj += prp_P;               // add the profit coming from the preprocessing

start_item = + Inf< int >();

for( int i = 0 ; i < f_NItems ; i++ ){      // change again the signs of the 
 if( items[ i ] == -1 ){                    // items with negative weight
  v_W[ i ] = - v_W[ i ];                    // and negative profit
  v_P[ i ] = - v_P[ i ];   
 }  
}     

HasSol = false;

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
  throw( std::invalid_argument( "No available solution" ) );

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
   
   reopt = value;

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
  reload = true;             // the Binary Knapsack instance must be re-loaded
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

   // get profits
   v_P.resize( f_NItems );
   const std::vector< double > & P = BKB->get_Profits();

   // if the sense is minimization change the sign of the profits      
   for( int i = 0 ; i < P.size() ; i++ )
    v_P[ i ] = f_sense ? P[ i ] : - P[ i ];

   // prepare vector of weights
   v_W.resize( f_NItems );
   const std::vector< double > & W = BKB->get_Weights();

   // load weights 
   for( int i = 0 ; i < W.size() ; i++ ){
   
    // check that weights are integers
    v_W[ i ] = std::round( W[ i ] );

    if( std::abs( v_W[ i ] - W[ i ] ) > 1e-06 )
     throw( std::invalid_argument( "Weights must be integers" ) );

   }

   BKB->read_unlock();

 v_x.resize( f_NItems );            // resize vector of variables
 
 items.resize( f_NItems );          // resize vector of items
 
 G.resize( f_NItems + 1 );          // resize Graph

 prp_nms.resize( f_NItems );
 
 preprocessing( std::make_pair( 0 , f_NItems ) );  // perform the preprocessing

 reload = false;

 } // end if( f_Block )

} // end( DPBinaryKnapsackSolver::load() )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::compute_var_solution(){

// reconstruct the optimal solution- - - - - - - - - - - - - - - - - - - - - -

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

 HasSol = true;    

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::preprocessing( Range rng ){

 // PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 /*
 The following items are preprocessed:
  
  - Items i with negative weight and positive profit => set x[ i ] = 1
 
  - Items i with positive weight and negative profit => set x[ i ] = 0

                                                                           */

  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  for( int i = rng.first ; i < rng.second ; i++ ){
   
   if( BKB->is_fixed( i ) ){
    items[ i ] = BKB->get_x( i );
    continue;
   }

   // items to be picked
   if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){
    items[ i ] = 1;
    continue;
   }

   // items not to be picked
   if(  ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ){
    items[ i ] = 0;
    continue; 
   }

   if( v_W[ i ] <= 0 ){             // if the item has negative weight 
    items[ i ] = -1;                // and negative profit set items to -1
    continue; 
   }

   items[ i ] = -Inf< int >();
   
  }
   
  start_item = rng.first;

 } // end( preprocessing( Range ) )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::preprocessing( Subset & nms ){

 // PREPROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 /*
 The following items are preprocessed:
  
  - Items i with negative weight and positive profit => set x[ i ] = 1
 
  - Items i with positive weight and negative profit => set x[ i ] = 0

                                                                           */

  auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

  for( auto i : nms ){
   
   if( BKB->is_fixed( i ) ){
    items[ i ] = BKB->get_x( i );
    continue;
   }

   // items to be picked
   if( ( v_W[ i ] <= 0 ) && ( v_P[ i ] >= 0 ) ){
    items[ i ] = 1;
    continue;
   }

   // items not to be picked
   if(  ( v_W[ i ] >= 0 ) && ( v_P[ i ] <= 0 ) ){
    items[ i ] = 0;
    continue; 
   }

   if( v_W[ i ] <= 0 ){             // if the item has negative weight 
    items[ i ] = -1;                // and negative profit set items to -1
    continue; 
   }

   items[ i ] = -Inf< int >();
   
  }
   
  start_item = nms[ 0 ];

 } // end( preprocessing( Subset ) )


/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::process_outstanding_Modification(){

 if( reload ){              // if a nuclear modification has 
  load();                   // been issued previously
  v_mod.clear();            // reload the instance
 }

 if( start_item != + Inf< int >() ){    // if compute is called for the first 
                                        // time, clear all the modifications 
  if( ! v_mod.empty() ){                // (if any) and reload the instance
   load();
   v_mod.clear();   
  }

  return;
 }

 // otherwise copy v_mod in a temporary list of modifications - - - - - - - -

 Lst_sp_Mod v_mod_tmp;      // temporary list of modifications

 // try to acquire lock, spin on failure

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );       // copy v_mod in v_mod_tmp

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 prp_rng.first = f_NItems;
 prp_rng.second = 0;

 std::fill( prp_nms.begin() , prp_nms.end() , 0 );
 
 for( auto mod : v_mod_tmp ){
  
  // BinaryKnapsackBlockMod - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod.get() );

   if( tmod ){ 

    switch( tmod->type() ){

     case( BinaryKnapsackBlockMod::eChgCapacity ):      // change Capacity
      capacity_Modification(); break;                   
     
     case( BinaryKnapsackBlockMod::eChgSense ):         // change sense 
      sense_Modification(); break;                      // of the objective
     
    } 

   } 

  } // end( for )


 for( auto mod : v_mod_tmp ){

  // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >( mod.get() );

   if( tmod ){ 

    switch( tmod->type() ){

     case( BinaryKnapsackBlockMod::eChgWeight ):        // change weights 
      weight_Modification( tmod->rng() ); 
     break;

     case( BinaryKnapsackBlockMod::eChgProfit ):        // change profits
      profit_Modification( tmod->rng() ); 
     break;

     case( BinaryKnapsackBlockMod::eFixX ):             // fix x
      fixX_Modification( tmod->rng() ); 
     break;
     
     case( BinaryKnapsackBlockMod::eUnfixX ):           // unfix x
      unFixX_Modification( tmod->rng() ); 
     break;          
    
    }

   }
  } 

  // BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >( mod.get() ); 

   if( tmod ){

    switch( tmod->type() ){

     case( BinaryKnapsackBlockMod::eChgWeight ):        // change weights 
      weight_Modification( std::move( tmod->nms() ) ); 
     break; 

     case( BinaryKnapsackBlockMod::eChgProfit ):        // change profits
      profit_Modification( std::move( tmod->nms() ) ); 
     break;

     case( BinaryKnapsackBlockMod::eFixX ):             // fix x
      fixX_Modification( std::move( tmod->nms() ) ); 
     break;
     
     case( BinaryKnapsackBlockMod::eUnfixX ):           // unfix x
      unFixX_Modification( std::move( tmod->nms() ) ); 
     break;          
    
    }

   }
  }

 } // end( for(  ) )

v_mod_tmp.clear();

if( prp_rng.first < prp_rng.second )
 preprocessing( prp_rng );

Subset nms;

for( int i = 0 ; i < prp_rng.first ; i++ ){
 if( prp_nms[ i ] )
  nms.push_back( i );   
}

for( int i = prp_rng.second ; i < f_NItems ; i++ ){
 if( prp_nms[ i ] )
  nms.push_back( i );   
} 

if( nms.size() )
 preprocessing( nms );


}// end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::capacity_Modification(){
 
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 int nC = std::floor( BKB->get_Capacity() );     // get new Capacity

 if( nC > f_C )             // if the new Capacity is > f_C
  start_item = 0;           // resolve from the beginning
 else{                             
                                
  int C = nC;               // otherwise only recompute the 
  int prp_P = 0;            // optimal value

  for( int i = 0 ; i < f_NItems ; i++ ){
   if( items[ i ] == 1 || items[ i ] == -1 ){
    C -= v_W[ i ];
    prp_P += v_P[ i ];
   }
  }
  
  int maxh = std::min( besth , C );
   
  obj = - Inf< int >();

  for( int i = 0 ; i <= maxh ; i++ ){
   if( G[ f_NItems ].lab[ i ] > obj ){
    obj = G[ f_NItems ].lab[ i ];
    besth = i;    
   }
  }

  obj += prp_P; 

  HasSol = false;           // the solution must be recomputed (if necessary)   
 }

 f_C = nC;                  // update the Capacity

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::sense_Modification(){ 

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 
 f_sense = BKB->get_objective_sense();

 for( int i = 0 ; i < f_NItems ; i++ )
  v_P[ i ] = - v_P[ i ]; 
 

 prp_rng.first = 0;
 prp_rng.second = f_NItems;

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::fixX_Modification( Range rng ){

 prp_rng.first = std::min( prp_rng.first , rng.first );
 prp_rng.second = std::max( prp_rng.second , rng.second );  

}


/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::fixX_Modification( c_Subset && nms ){
 
 for( auto i : nms )
  prp_nms[ i ] = 1;  

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::unFixX_Modification( Range rng ){
 
 prp_rng.first = std::min( prp_rng.first , rng.first );
 prp_rng.second = std::max( prp_rng.second , rng.second );  

}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::unFixX_Modification( c_Subset && nms ){
 
 for( auto i : nms )
  prp_nms[ i ] = 1;  
    
}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::weight_Modification( Range rng ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( int i = rng.first ; i < rng.second ; i++ ){
  
  int nw = std::round( BKB->get_Weight( i ) );
  
  if( std::abs( nw - BKB->get_Weight( i ) ) > 1e06 )
   throw( std::invalid_argument( "Weights must be integers" ) );
   
  v_W[ i ] = nw;

 }

 prp_rng.first = std::min( prp_rng.first , rng.first );
 prp_rng.second = std::max( prp_rng.second , rng.second );  
 
}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::weight_Modification( c_Subset && nms ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( auto i : nms ){
  
  int nw = std::round( BKB->get_Weight( i ) );
  
  if( std::abs( nw - BKB->get_Weight( i ) ) > 1e06 )
   throw( std::invalid_argument( "Weights must be integers" ) );
   
  v_W[ i ] = nw;

  prp_nms[ i ] = 1;

 }
 
}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::profit_Modification( Range rng ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( int i = rng.first ; i < rng.second ; i++ )
  v_P[ i ] = BKB->get_Profit( i );

 prp_rng.first = std::min( prp_rng.first , rng.first );
 prp_rng.second = std::max( prp_rng.second , rng.second );     
    
}

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::profit_Modification( c_Subset && nms ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 for( auto i : nms ){
  v_P[ i ] = BKB->get_Profit( i );
  prp_nms[ i ] = 1;
 }

}


/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/


