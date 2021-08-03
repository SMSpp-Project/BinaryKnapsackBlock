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

// define and initialize here the map for double parameters names
const std::map< std::string , DPBinaryKnapsackSolver::idx_type >
 DPBinaryKnapsackSolver::dbl_pars_map = {
 
 { "dblReopt" , DPBinaryKnapsackSolver::dblReopt }
 
 };

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

// process all the pending modifications and compute start_item 

process_outstanding_Modification(); 

// check if the problem is empty - - - - - - - - - - - - - - - - - - - - - - -

auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );
 
if( BKB->is_empty() ){ 
 obj = - Inf< double >();
 return( kInfeasible );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// start_item is the first item processed by the DP algorithm. If compute() is
// called for the first time, start_item = 0 and the algorithm is entirely 
// executed. Otherwise it is computed in process_outstanding_Modification()                                                               

if( start_item == + Inf< int >() )  // if start_item == Inf it is assumed that 
 return( kOK );                     // everything has already been done in 
                                    // process_outstanding_Modification()

int i = start_item;                 // otherwise start from start_item

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// As some of the items may be pre-processed, it is necessary to update
// the capacity by subtracting the weight of the items that are fixed to 1, 
// also saving their profits which need to be added to the objective

int C = f_C;                            // "residual" capacity
double prp_P = 0;                       // "residual" profit
             

for( int i = 0 ; i < f_N ; i++ ){       // also the items with both negative
 if( isFixed1( i ) || isNeg( i ) ){     // weight and profit are treated as if
  C -= v_W[ i ];                        // they were fixed to 1         
  prp_P += v_P[ i ];                                  
 }                          
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int maxcurrlab;                         // max current height of the graph
int maxnextlab;                         // max next height of the graph    

std::vector< double > currlab;          // set of current labels
std::vector< double > nextlab;          // set of next labels
                     
G[ 0 ].lab.resize( 1 );                 // initialize dummy node              
G[ 0 ].lab[ 0 ] = 0;                    // in the origin

// each time i % step == 0, currlab is saved in G[ i ].lab
int step;
 
if( reopt == 1 )
 step = 1;
else
 step = f_N * std::exp2( - std::log2( f_N ) * reopt ); 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Initialize maxcurrlab and currlab with data stored in G[ i ].lab

currlab = G[ i ].lab;

if( C + 1 < currlab.size() )            // Capacity may have been changed
 currlab.resize( C + 1 );

maxcurrlab = currlab.size() - 1;

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// At each iteration the algorithm computes the maximum height of the 
// constructed graph in order to allocate only the necessary memory to store  
// the next labels.
//
// According to the reopt parameter, labels can be saved in G[ i ].lab so that
// they will be used for reoptimization.
                                                                            

// DP Algorithm- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for( ; i < f_N ; i++ ){

 // reoptimization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 if( i % step == 0  ) 
  G[ i ].lab = currlab;             // save labels in G[ i ].lab
 else 
  G[ i ].lab.clear();               // else clear previous labels (if any) 
 
 // pre-processing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( isFixed( i ) )                 // if i is fixed (or can be pre-processed) 
  continue;                         // continue       

 double p = v_P[ i ];               // profit of the current item   
 int w = v_W[ i ] ;                 // weight of the current item

 if( isNeg( i ) ){                  // if the item has negative 
  p = -p;                           // profit and weight
  w = -w;                           // change both signs
 }

 if( w > C )                        // if the weight exceeds the capacity
  continue;                         // continue

                    
 // max height of the graph
 maxnextlab = std::min( maxcurrlab + w , C );

 // Initialize labels 
 nextlab.resize( maxnextlab + 1 );
 std::fill( nextlab.begin() , nextlab.end() , - Inf< double >() );

 // Allocate predecessors
 G[ i + 1 ].pred.resize( maxnextlab + 1 );

 // Initialize the best label among those of this slice 
 double bestlab = -Inf< double >();

 for( int j = 0 ; j <= maxcurrlab ; j++ ){
                                                
  if( currlab[ j ] <= bestlab )     // if the node does not exists 
   continue;                        // (currlab = -Inf) or if it does  
                                    // not have a better label continue                 
 
  bestlab = currlab[ j ];                       // update bestlab

  if( currlab[ j ] > nextlab[ j ] ){            // horizontal arc
   G[ i + 1 ].pred[ j ] = false;               
   nextlab[ j ] = currlab[ j ];
  }

  if( j + w > C )                               // check if the maximum 
   continue;                                    // capacity is reached

  if( currlab[ j ] + p > nextlab[ j + w ] ){    // diagonal arc 
   G[ i + 1 ].pred[ j + w ] = true;            
   nextlab[ j + w ] = currlab[ j ] + p;
  }

 }

 // swaps
 maxcurrlab = maxnextlab;
 std::swap( currlab , nextlab );

} // end( DP algorithm )- - - - - - - - - - - - - - - - - - - - - - - - - - - 

// always save last labels in G[ f_N ].lab

std::swap( G[ f_N ].lab , currlab );

// find optimal value - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

int besth;

obj = - Inf< double >();                    // Initialize objective value

for( int i = 0 ; i <= maxcurrlab ; i++ ){ 
 if( G[ f_N ].lab[ i ] > obj ){
  obj = G[ f_N ].lab[ i ]; 
  besth = i;
 }
}

G[ f_N ].lab.resize( besth + 1 ); 

obj += prp_P;             // add the profit coming from the preprocessing

start_item = + Inf< int >();

v_x.clear();              // clear previous solution (if any)

return( kOK );

}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::get_var_solution( Configuration * solc ){

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 // check if the BinaryKnapsackBlock uses the correct abstract representation
 if( BKB->get_NItems() != BKB->get_VarSize() )
  throw( std::invalid_argument( "Incompatible variables size" ) );

 // check if compute has been called before
 if( start_item != + Inf< int >() )
  throw( std::invalid_argument( "compute() must be called first" ) );

 // check if the solution has already been computed
 if( ! v_x.empty() ){       
  BKB->set_x( v_x );        // write the solution in the BinaryKnapsackBlock 
  return;                     
 }         

 // reconstruct the optimal solution - - - - - - - - - - - - - - - - - - - - -
 // for each item check if it is fixed (because the corresponding variable is
 // fixed or the item has been preprocessed) or if its weight exceeds the 
 // total capacity and updates the variables accordingly. Reconstruct the 
 // rest of the solution from G[ i ].pred. For the "negative items" (with 
 // negative weight and profit) change x with 1 - x

 int besth = G[ f_N ].lab.size() - 1;

 v_x.resize( f_N );

 for( int i = f_N - 1 ; i >= 0 ; i-- ){

  if( isFixed0( i ) ){                          // items fixed to 0
   v_x[ i ] = 0;
   continue;    
  }

  if( isFixed1( i ) ){                          // items fixed to 1
   v_x[ i ] = 1;
   continue;    
  }  

  int w = isNeg( i ) ? - v_W[ i ] : v_W[ i ];   // weight of the current item
  
  if( w > besth ){                              // if w > besth either the 
   v_x[ i ] = isNeg( i ) ? 1 : 0;               // weight exceeds the capacity 
   continue;                                    // or it has surely not been
  }                                             // selected

  if( G[ i + 1 ].pred[ besth ] ){
   v_x[ i ] = 1; 
   besth -= w;  
  }
  else
   v_x[ i ] = 0;

  if( isNeg( i ) )                  // items with negative weight and profit
   v_x[ i ] = ! v_x[ i ];

 }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
 BKB->set_x( v_x );         // write the solution in the BinaryKnapsackBlock

}
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::set_par( idx_type par , double value ){
 
 switch( par ){
 
  // Reoptimization parameter- - - - - - - - - - - - - - - - - - - - - - - - - 

  case( dblReopt ):{

   if( value < 0 || value > 1 )
    throw(std::invalid_argument("dblReopt parameter must be in [ 0 , 1 ]"));

   reopt = value;
   
   start_item = 0;          // restart from 0 in the next call of compute()

   break;
  }  
 
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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

   // read the block
   bool owned = BKB->is_owned_by( f_id );
      
   if( !owned && !BKB->read_lock() )
    throw std::runtime_error( "Unable to read the Block" );
   
   // load Binary Knapsack instance - - - - - - - - - - - - - - - - - - - - - 

   f_sense = BKB->get_objective_sense();    // get the sense of the objective

   f_N = BKB->get_NItems();                 // get the number of items

   f_C = std::floor( BKB->get_Capacity() ); // get the Capacity

   v_P.resize( f_N );                  

   const auto & P = BKB->get_Profits();     // get profits  
                                             
   for( int i = 0 ; i < P.size() ; i++ )    // if the sense is minimization 
    v_P[ i ] = f_sense ? P[ i ] : - P[ i ]; // change the sign of the profits 

   v_W.resize( f_N );                       // prepare vector of weights

   const auto & W = BKB->get_Weights();

   for( int i = 0 ; i < W.size() ; i++ ){   // load weights and check 
    v_W[ i ] = std::round( W[ i ] );        // that they are integers

    if( std::abs( v_W[ i ] - W[ i ] ) > WeightIntegrality )
     throw( std::invalid_argument( "Weights must be integers" ) );

   }

  if( !owned )
   BKB->read_unlock();

  // end load Binary Knapsack instance- - - - - - - - - - - - - - - - - - - -
 
  G.clear();                            // clear previous graph (if any)  

  G.resize( f_N + 1 );                  // resize Graph
  
  start_item = 0;                       // (re-)start from the beginning

  v_x.clear();                          // clear previous solution (if any)

 } // end if( f_Block )

} // end( DPBinaryKnapsackSolver::load() )

/*--------------------------------------------------------------------------*/


void DPBinaryKnapsackSolver::process_outstanding_Modification(){
                                 
// copy v_mod in a temporary list of modifications - - - - - - - - - - - - - - 

 Lst_sp_Mod v_mod_tmp;              // temporary list of modifications

 // try to acquire lock, spin on failure

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );       // copy v_mod in v_mod_tmp

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

// Any changes in Profits must be processed only AFTER checking the changes 
// on the sense of the Objective. Hence v_mod_tmp is scanned twice, checking
// modifications on Capacity and Objective sense first 

 auto mod = v_mod_tmp.begin(); 

 while( mod != v_mod_tmp.end() ){

  // BinaryKnapsackBlockMod- - - - - - - - - - - - - - - - - - - - - - - - - -
  
  const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod->get() );

   if( tmod ){ 
    switch( tmod->type() ){ 

    // Change Capacity - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // If the new capacity is smaller than the previous one, restart the DP
    // algorithm from f_N, i.e. only recompute the optimal value
    // Otherwise restart from the first item   
    
     case( BinaryKnapsackBlockMod::eChgCapacity ): {
     
      int nC = std::floor( BKB->get_Capacity() );   // get new Capacity

      start_item = nC > f_C ? 0 : std::min( f_N , start_item );

      f_C = nC;                                     // update the Capacity
      
      mod = v_mod_tmp.erase( mod );
      break;
     }                    

    // Change Objective Sense- - - - - - - - - - - - - - - - - - - - - - - - - 
    // Change the sign of all profits and restart from the first item

     case( BinaryKnapsackBlockMod::eChgSense ):{        
      
      f_sense = BKB->get_objective_sense();  // update f_sense

      for( int i = 0 ; i < f_N ; i++ )       // change the sign of all profits
       v_P[ i ] = - v_P[ i ];

      start_item = 0;                        // restart from the beginning 
      
      mod = v_mod_tmp.erase( mod );
      break;
     }                                      
     
     default: mod++;
    
    }                                              
   }
   else
    mod = v_mod_tmp.erase( mod );   // it is not a physical modification
 }

for( auto mod : v_mod_tmp ){                                                

  // BinaryKnapsackBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >( mod.get() );

   if( tmod ){ 
    switch( tmod->type() ){
    
     // change Profits - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
     // update modified profits according to f_sense and update start_item

     case( BinaryKnapsackBlockMod::eChgProfit ):{          
     
      for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ )
       v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );

      start_item = std::min( start_item , int( tmod->rng().first ) );

      break;
     } 
    
    // fix x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    // update start item with the first fixed item

     case( BinaryKnapsackBlockMod::eFixX ):               
      { start_item = std::min( start_item , int( tmod->rng().first ) ); break; }
    
    // unfix x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // unfix a variable x could increase the residual capacity making
    // reoptimization not exploitable 

     case( BinaryKnapsackBlockMod::eUnfixX ): start_item = 0; break;         

    // change Weights- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // update modified weights checking the integrality property. 
    // check also if a new weight increases the residual capacity making
    // reoptimization not exploitable 

     case( BinaryKnapsackBlockMod::eChgWeight ):{    
       
      for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ){
  
       int nw = std::round( BKB->get_Weight( i ) );     // new weight
  
       if( std::abs( nw - BKB->get_Weight( i ) ) > WeightIntegrality )
        throw( std::invalid_argument( "Weights must be integers!" ) );
   
       if( nw < 0 && nw < v_W[ i ] )        // it is not possible 
        start_item = 0;                     // to re-optimize

       v_W[ i ] = nw;                       // update weight

      }

      start_item = std::min( start_item , int( tmod->rng().first ) );

      break;

     }  
    
    }
   }
  } 

  // BinaryKnapsackBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - 
  {
   const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >( mod.get() ); 

   if( tmod ){ 
    switch( tmod->type() ){
    
     // change Profits - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
     // update modified profits according to f_sense and update start_item 

     case( BinaryKnapsackBlockMod::eChgProfit ):{          
     
      for( auto i : tmod->nms() )
       v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );

      start_item = std::min( start_item , int( tmod->nms()[ 0 ] ) );

      break;
     } 
    
    // fix x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    // update start item with the first fixed item 

     case( BinaryKnapsackBlockMod::eFixX ):                
      { start_item = std::min( start_item , int( tmod->nms()[ 0 ] ) ); break; }
    
    // unfix x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // unfix a variable x could increase the residual capacity making
    // reoptimization not exploitable 

     case( BinaryKnapsackBlockMod::eUnfixX ): start_item = 0; break;

    // change Weights- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // update modified weights checking the integrality property. 
    // check also if a new weight increases the residual capacity making
    // reoptimization not exploitable 

     case( BinaryKnapsackBlockMod::eChgWeight ):{    
       
      for( auto i : tmod->nms() ){
  
       int nw = std::round( BKB->get_Weight( i ) );     // new weight
  
       if( std::abs( nw - BKB->get_Weight( i ) ) > WeightIntegrality )
        throw( std::invalid_argument( "Weights must be integers!" ) );
   
       if( nw < 0 && nw < v_W[ i ] )        // it is not possible 
        start_item = 0;                     // to re-optimize

       v_W[ i ] = nw;                       // update weight

      }

      start_item = std::min( start_item , int( tmod->nms()[ 0 ] ) );

      break;

     }  
    
    }
   }
  }

 } // end( for(  ) )

v_mod_tmp.clear();              // clear the temporary list of modification

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// compute start_item- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// Each modification updated start_item with the first modified item. At this
// point it is necessary to retrieve the first item smaller than start_item  
// whose corresponding labels have been previously stored in G[ i ].lab. 

 int step;
 
 if( reopt == 1 )
  step = 1;
 else
  step = f_N * std::exp2( - std::log2( f_N ) * reopt ); 

 if( start_item != + Inf< int >() )
  start_item = ( start_item / step ) * step;


}// end( process_outstanding_Modification() )


/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/

















