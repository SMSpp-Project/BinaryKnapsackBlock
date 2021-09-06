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

// if all the variables are continuous skip the resolution of the integer part

if( countCont == f_N )
 start_item = f_N;

indexContinuous.clear();
for( Index i = 0 ; i < f_N ; i++ ){
 if( ! v_I[ i ] )
  indexContinuous.push_back( i ); 
 }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// As some of the items may be pre-processed, it is necessary to update
// the capacity by subtracting the weight of the items that are fixed to 1, 
// also saving their profits which need to be added to the objective

double C = f_C;                         // "residual" capacity
double prp_P = 0;                       // "residual" profit
             

for( int i = 0 ; i < f_N ; i++ ){       // also the items with both negative
 if(( isFixed1( i ) || isNeg( i ) ) ){  // weight and profit are treated as if
  C -= v_W[ i ];                        // they were fixed to 1         
  prp_P += v_P[ i ];                                  
 }                          
}

int iC = std::floor( C );               // integer Capacity

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int maxcurrlab;                         // max current height of the graph
int maxnextlab;                         // max next height of the graph    

std::vector< double > currlab;          // set of current labels
std::vector< double > nextlab;          // set of next labels
                     
G[ 0 ].lab.resize( 1 );                 // initialize dummy node              
G[ 0 ].lab[ 0 ] = 0;                    // in the origin

// each time i % step == 0, currlab is saved in G[ i ].lab
int step = compute_step();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Initialize maxcurrlab and currlab with data stored in G[ start_item ].lab 

currlab = G[ start_item ].lab;

if( iC + 1 < currlab.size() )            // Capacity may have been changed
 currlab.resize( iC + 1 );

maxcurrlab = currlab.size() - 1;

// DYNAMIC PROGRAMMING - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// At each iteration the algorithm computes the maximum height of the 
// constructed graph in order to allocate only the necessary memory to store  
// the next labels.
//
// According to the reopt parameter, labels can be saved in G[ i ].lab so that
// they will be used for reoptimization.
                                                                            

// DP Algorithm- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for( Index i = start_item ; i < f_N ; i++ ){
 
 // reoptimization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 if( i % step == 0  ) 
  G[ i ].lab = currlab;             // save labels in G[ i ].lab
 else 
  G[ i ].lab.clear();               // else clear previous labels (if any) 
 
 // pre-processing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( isFixed( i ) || ! v_I[ i ] )   // if i is fixed or continuous 
  continue;                         // continue       

 double p = v_P[ i ];               // profit of the current item   
 int w = v_W[ i ] ;                 // weight of the current item

 if( isNeg( i ) ){                  // if the item has negative 
  p = -p;                           // profit and weight
  w = -w;                           // change both signs
 }

 if( w > iC )                       // if the weight exceeds the capacity
  continue;                         // continue

 // max height of the graph 
 maxnextlab = std::min( maxcurrlab + w , iC );
 
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
}

// always save last labels in G[ f_N ].lab
std::swap( G[ f_N ].lab , currlab );

//end integer part - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// sort continuous variable using as weight the proportion of profits/weights
// if( !is_sorted ){
   sort( indexContinuous.begin(), indexContinuous.end(), [&]( const int & a, const int & b ){     
            return ( v_P[ a ] / v_W[ a ] > v_P[ b ] / v_W[ b ] );} );
//   is_sorted = true;
//  }


// compute the objective calculating also the contribution of the continuous variables

obj = -Inf< double >();

besth = 0;
lastIndex = 0;               // index of the last element considered in the continuous knapsack
int tempLastIndex = 0;       // the same but this one is used only in the resolution, the previous one to compute the solution  
lastVar = 0;                 // value of the last variable considered in the knapsack 
double boundX = 1;           // current bound for the variable, 0<=boundX<=1  
double contProfit = 0;       // cumulative profits related to the continuous solution 
double contWeight = 0;       // cumulative wieghts related to the continuous solution
double rC = C - maxcurrlab;  // residual capacity

// compute the continuous part - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//if(maxcurrlab==0)
// obj=G[ f_N ].lab[ 0 ];

// we start from the higher heigth, that have lower residual capacity
for( int i = maxcurrlab ; i >= 0 ; i-- ){

//starting from the index analyzed during the previous iteration (height) 
 for( int j = tempLastIndex ; j < indexContinuous.size() ; j++ ){
 
 // if the item is fixed to zero or one we skip the current iteration
   if( isFixed( indexContinuous[ j ] ) ){
   tempLastIndex++;
   continue;
   }
  
 double p = v_P[ indexContinuous[ j ] ];               // profit of the current item   
 int w = v_W[ indexContinuous[ j ] ] ;                 // weight of the current item
 
 if( isNeg(  indexContinuous[ j ]  ) ){                  // if the item has negative 
  p = -p;                                                // profit and weight
  w = -w;                                                // change both signs
 }
 
  // if the item does not fit entirely in the knapsack
  if( w * boundX + contWeight > rC ){

   // take only a fraction and break
   boundX -= ( rC - contWeight ) / ( w );  
   contProfit += ( rC - contWeight ) / ( w) * p ;
   contWeight = rC; // / ( v_W[ indexContinuous[ j ] ] );
   break; 
  }else{
  // otherwise take the whole item and consider the following item
   contProfit += p * boundX;
   contWeight += w * boundX;
   boundX = 1;
   tempLastIndex++;
  }

  }
 // check if it is the best solution and update obj
 if( G[ f_N ].lab[ i ] + contProfit > obj ){
  obj = G[ f_N ].lab[ i ] + contProfit;
  besth = i;
  lastVar = 1 - boundX;  
  lastIndex = tempLastIndex;
 }
 rC += 1;        // residual capacity

}

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

 v_x.resize( f_N );

 // Integer part - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( int i = f_N - 1 ; i >= 0 ; i-- ){

  if( ! v_I[ i ] )
    continue;

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
 
 // Continuous part - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( int i = 0 ; i < indexContinuous.size() ; i++ ){

  if( isFixed0( indexContinuous[ i ] ) ){        // items fixed to 0
   v_x[ indexContinuous[ i ] ] = 0;
   continue;    
  }

  if( isFixed1( indexContinuous[ i ] ) ){        // items fixed to 1
   v_x[ indexContinuous[ i ] ] = 1;
   continue;    
  } 
  

  if( i < lastIndex )                               // if the item preceed the last item considered
    v_x[ indexContinuous[ i ] ] = 1;                // we fix the variable to one
  else if( i == lastIndex )                         // if the item is exactly the last item considered
    v_x[ indexContinuous[ lastIndex ] ] = lastVar;  // the variable has value lastVar \in (0,1]
  else                                              // otherwise we fix the variable to zero
    v_x[ indexContinuous[ i ] ] = 0;
    
  // items with negative weight and profit
  if( isNeg( indexContinuous[ i ] ) )               
   v_x[ indexContinuous[ i ] ] = 1 - v_x[ indexContinuous[ i ] ];

  } // end( continuous part )
   
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

   f_C = BKB->get_Capacity(); // get the Capacity

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
   v_I.resize( f_N );    
   
   const auto & I = BKB->get_Integrality();

   for( int i = 0 ; i < I.size() ; i++ )  // load weights and check 
    v_I[ i ] = ( bool ) I[ i ] ;          // that they are booleans
 
   countCont = BKB->get_N_Cont_Items();

  if( !owned )
   BKB->read_unlock();

  indexContinuous.clear();
  for( Index i = 0 ; i < f_N ; i++ ){
   if( ! v_I[ i ] )
    indexContinuous.push_back( i ); 
  }

  is_sorted = false;

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
     
      double nC = BKB->get_Capacity();   // get new Capacity

      start_item = nC > f_C ? 0 : std::min( f_N , start_item );

      f_C = nC;                          // update the Capacity
      
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
  
       double nw = std::round( BKB->get_Weight( i ) );     // new weight
  
       if( std::abs( nw - BKB->get_Weight( i ) ) > WeightIntegrality )
        throw( std::invalid_argument( "Weights must be integers!" ) );
   
       if( nw < 0 && nw < v_W[ i ] )        // it is not possible 
        start_item = 0;                     // to re-optimize

       v_W[ i ] = nw;                       // update weight

      }

      start_item = std::min( start_item , int( tmod->rng().first ) );

      break;

     }
     case( BinaryKnapsackBlockMod::eChgIntegrality ):{
     
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ){
     
       bool ni = BKB->get_Integrality( i );                  // new integrality
       
       v_I[ i ] = ni;                                        // update 
     
     }
    
     start_item = 0;
     
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
  
       double nw = std::round( BKB->get_Weight( i ) );     // new weight
  
       if( std::abs( nw - BKB->get_Weight( i ) ) > WeightIntegrality )
        throw( std::invalid_argument( "Weights must be integers!" ) );
   
       if( nw < 0 && nw < v_W[ i ] )        // it is not possible 
        start_item = 0;                     // to re-optimize

       v_W[ i ] = nw;                       // update weight

      }

      start_item = std::min( start_item , int( tmod->nms()[ 0 ] ) );

      break;

     }
     
     
     case( BinaryKnapsackBlockMod::eChgIntegrality ):{
     
     for( auto i : tmod->nms() ){
      
       bool ni = BKB->get_Integrality( i );   //new integrality
       
       v_I[ i ] = ni;                         //update integrality
     
     }
     
     start_item = 0;                         //no-reoptimization
     
     break;
     
     }   
    
    }
   }
  }

 } // end( for(  ) )

v_mod_tmp.clear();              // clear the temporary list of modification


// TODO: count while managing the modification
// count again the number of continuous variables
countCont = std::count( v_I.begin() , v_I.end() , false );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// compute start_item- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// Each modification updated start_item with the first modified item. At this
// point it is necessary to retrieve the first item smaller than start_item  
// whose corresponding labels have been previously stored in G[ i ].lab. 

 int step = compute_step();

 if( start_item != + Inf< int >() )
  start_item = ( start_item / step ) * step;


}// end( process_outstanding_Modification() )


/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
