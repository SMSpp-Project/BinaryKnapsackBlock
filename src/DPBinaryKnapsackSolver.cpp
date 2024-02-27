/*--------------------------------------------------------------------------*/
/*--------------------- File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class DPBinaryKnapsackSolver, which
 * implements the Solver concept [see Solver.h] for solving Knapsack problems
 * with linear objective and both integer and continuous variables, as
 * represented by a BinaryKnapsackBlock, using a standard Dynamic Programming
 * approach to deal with the integer variables combined with the exact greedy
 * algorithm to deal with the continuous ones.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Francesco Demelas \n
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Federica Di Pasquale, Francesco Demelas,
 *                      Antonio Frangioni
 */ 
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DPBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE AND USING --------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using idx_type = ThinComputeInterface::idx_type;

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::set_Block( Block * block ) {
 
 if( block == f_Block )       // nothing to do        
  return;

 Solver::set_Block( block );  // attach to the new Block

 load();                      // load Binary Knapsack instance
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/// Solve the Binary Knapsack problem

int DPBinaryKnapsackSolver::compute( bool changedvars ) {
 
 lock();                     // lock the mutex

 // Process modifications and compute the first item to process (start_item)
 process_outstanding_Modification();

 // check start_item
 if( start_item == Inf< int >() ) {   // INF means everything altready done
  unlock();                           // unlock the mutex and
  return( kOK );                      // return
  }
 
 // apply preprocessing and compute residual capacity C and residual profit P
 auto [ C , P ] = preprocessing();           

 // check if the problem is empty (iff residual capacity < 0)
 if( C < 0 ) {
  obj = -Inf< double >();
  unlock();
  return( kInfeasible );
  }

 // compute step for reoptimization
 Index k = std::floor( reopt * std::log2( f_N ) );

 if( k == 0 )
  step = f_N;
 else  
  step = std::floor( f_N / std::exp2( k ) ); 

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 dynamic_programming( std::floor( C ) );    // solve the integer knapsack
 
 greedy_algorithm( C );                     // solve the continuous knapsack

 obj += P;                        // add residual profit (from preprocessing)
 
 // End - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 start_item = Inf< int >(); // to avoid solving again the same instance

 Return_OK:
 
 unlock();                  // unlock the mutex     
 
 return( kOK );

 }  // end( DPBinaryKnapsackSolver::compute )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// perform the preprocessing and return residual Capacity and Profit

std::tuple< double , double > DPBinaryKnapsackSolver::preprocessing() {
  
 // Initialize variables to skip during the solution algorithms
 skip.assign( f_N , false );

 // residual capacity and residual profit
 double C = f_C;
 double P = 0;

 // for each item
 for( Index i = 0 ; i < f_N ; ++i ) {

  // if the variable is fixed to 0 by the user
  if( v_fxd[ i ] == 1 ) {
   v_x[ i ] = 0;
   skip[ i ] = true;
   continue;
   }

  // if the variable is fixed to 1 by the user
  if( v_fxd[ i ] == 2 ) {
   v_x[ i ] = 1;
   skip[ i ] = true;
   C -= v_W[ i ];         // update Capacity
   P += v_P[ i ];         // update total Profit
   continue;
   }

  // if the item has positive weight and negative profit
  if( v_W[ i ] >= 0 && v_P[ i ] <= 0 ) {
   v_x[ i ] = 0;
   skip[ i ] = true;
   continue; 
   }

  // if the item has negative weight and positive profit
  if( v_W[ i ] <= 0 && v_P[ i ] >= 0 ) {
   v_x[ i ] = 1;
   skip[ i ] = true;
   C -= v_W[ i ];         // update Capacity
   P += v_P[ i ];         // update total Profit
   continue;
   }
   
  // if the item has both negative weight and profit
  if( v_W[ i ] < 0 && v_P[ i ] < 0 ) {
   C -= v_W[ i ];         // update Capacity
   P += v_P[ i ];         // update total Profit
   }

  }

 // check if some item has a weight > C (check only for integer items)
 for( Index i = 0 ; i < f_N ; ++i ) {
  
  if( skip[ i ] || ( ! v_I[ i ] ) )
   continue;

  if( v_W[ i ] > C ) {
   skip[ i ] = true;
   v_x[ i ] = 0;
   }
 
  }

 // if the problem is not empty
 if( C >= 0 ) {

  // new integer residual capacity
  Index iC = std::floor( C );

  // resize lab[ start_item ] according to the new iC
  if( iC + 1 < lab[ start_item ].size() )
   lab[ start_item ].resize( iC + 1 );
  }

 return { C , P };

 } // end( DPBinaryKnapsackSolver::preprocessing )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Dynamic Programming

void DPBinaryKnapsackSolver::dynamic_programming( Index C ) {

 std::vector< double > currlab = lab[ start_item ]; // set of current labels 
 std::vector< double > nextlab;                     // set of next labels  

 for( Index i = start_item ; i < f_N ; ++i ) {      // for each item

  if( i % step == 0  )          // reoptimization: save currlab in lab[ i ]
   lab[ i ] = currlab;             

  if( skip[ i ] || ( ! v_I[ i ] ) ) // skip preprocessed and continuous variables
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
  nextlab.assign( maxnextlab , -Inf< double >() );
  pred[ i + 1 ].resize( maxnextlab ); 

  // initialize the best label among ( i , j ) nodes with fixed i
  double bestlab = -Inf< double >();
  
  // compute nextlab
  for( Index j = 0 ; j < currlab.size() ; ++j ) {

   if( currlab[ j ] <= bestlab )        // skip node with label = -inf or
    continue;                           // with a worse label than bestlab
   
   bestlab = currlab[ j ];                        // update bestlab
                    
   if( currlab[ j ] > nextlab[ j ] ) {            // horizontal arc
    pred[ i + 1 ][ j ] = false;                         
    nextlab[ j ] = currlab[ j ];
    }

   if( j + w > C )                                // check capacity limit                               
    continue;                                    

   if( currlab[ j ] + p > nextlab[ j + w ] ) {    // diagonal arc
    pred[ i + 1 ][ j + w ] = true;             
    nextlab[ j + w ] = currlab[ j ] + p;
    }
  
   }

  std::swap( currlab , nextlab );

  }

 // always save last vector of labels (instead of copying just swap)
 std::swap( lab.back() , currlab );

 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// Greedy algorithm

void DPBinaryKnapsackSolver::greedy_algorithm( double C ) {

 // sort continuous variable in order of profits/weights
 sort( idxCont.begin() , idxCont.end() , [ & ]( auto a , auto b ) { 
       return( v_P[ a ] / v_W[ a ] > v_P[ b ] / v_W[ b ] );
      } );
 
 // last vector of labels contains (partial) solutions of the integer part
 const auto & lastlab = lab.back();
 
 // compute continuous contribution for each label of lastlab and eventually
 // select the best one. Each label lastlab[ h ] corresponds to a (partial)
 // solution with total weight = h, hence the continuous part is computed by
 // considering only the residual capacity rC - h, starting from the maximum
 // h and increasing rC by 1 each time the next (lower) label is considered.  
 
 Index maxh = lastlab.size() - 1;   // maximum h
 C -= maxh;                         // residual capacity
 
 // Initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 besth = 0;               // lastlab[ besth ] will be optimal 
 Index lastIndex = 0;     // index of the critical variable           
 Index tempLastIndex = 0; // temporary lastIndex
 double lastVar = 0;      // value of the critical variable
 double boundX = 1;       // current bound for the critical variable,
                          // 0 <= boundX <=1

 double contProfit = 0;   // cumulative profits of the continuous solution 
 double contWeight = 0;   // cumulative weights of the continuous solution
 
 obj = -Inf< double >(); // initialize objective with -INF

 // start from the highest height, that has lower residual capacity
 for( int i = maxh ; i >= 0 ; --i ) {
  
  // starting from the index analyzed during the previous iteration
  for( Index j = tempLastIndex ; j < idxCont.size() ; ++j ) {
   
   if( skip[ idxCont[ j ] ] ) {       // skip preprocessed variables
    tempLastIndex++;
    continue;
    }

   double p = v_P[ idxCont[ j ] ];    // profit of the current item    
   int w = v_W[ idxCont[ j ] ];       // weight of the current item    
   if( p < 0 && w < 0 ) {             // if both are negative  
    p = -p;                           // change the signs                         
    w = -w;                         
    }
 
   // if the item does not fit entirely in the knapsack
   if( w * boundX + contWeight > C ) {             
    // take only a fraction, update profit and weight and break
    boundX -= ( C - contWeight ) / ( w );
    contProfit += ( C - contWeight ) / ( w ) * p ;
    contWeight = C;
    break;
    }
   else {  
    // otherwise take the whole item and consider the following item
    contProfit += p * boundX;
    contWeight += w * boundX;
    boundX = 1;
    tempLastIndex++;
    }
   }

  // check if it is the best solution and update obj 
  if( lastlab[ i ] + contProfit > obj ) {
   obj = lastlab[ i ] + contProfit;
   besth = i;
   lastVar = 1 - boundX;  
   lastIndex = tempLastIndex;
   }

  C += 1;        // update residual capacity
  }

 // update the solution
 for( Index i = 0 ; i < idxCont.size() ; ++i ) {
  
  Index item = idxCont[ i ];

  if( skip[ item ] )
   continue;

  if( i < lastIndex )          // it is before the critical variable
    v_x[ item ] = 1;           // we fix the variable to one
  else
   if( i == lastIndex )        // it is exactly the critical variable
    v_x[ item ] = lastVar;     // has value lastVar
   else                        // it is after the critical variable
    v_x[ item ] = 0;           // we fix the variable to zero

  // items with negative weight and profit
  if( v_W[ item ] < 0 && v_P[ item ] < 0 )               
   v_x[ item ] = 1 - v_x[ item ];

  }

 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::get_var_solution( Configuration * solc ) {
 
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 // check if compute has been called before
 if( start_item != Inf< int >() )
  throw( std::invalid_argument( "compute() must be called first" ) );

 // reconstruct the optimal solution - - - - - - - - - - - - - - - - - - - - -
 // for each item check if it is fixed (because the corresponding variable is
 // fixed or the item has been preprocessed) or if its weight exceeds the 
 // total capacity and updates the variables accordingly. Reconstruct the 
 // rest of the solution from pred[ i ]. For the "negative items" (with 
 // negative weight and profit) change x with 1 - x

 // Integer part - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 int h = besth;
 
 for( int i = f_N - 1 ; i >= 0 ; i-- ) {
  
  if( skip[ i ] || ( ! v_I[ i ] ) )  // skip preprocessed and continuous variables
   continue;

  int w = v_W[ i ];                    // weight of the current item
  if( v_W[ i ] < 0 && v_P[ i ] < 0 )
    w = -w;

  if( pred[ i + 1 ][ h ] ) {
   v_x[ i ] = 1;
   h -= w;
   }
  else
   v_x[ i ] = 0;

  if( v_W[ i ] < 0 && v_P[ i ] < 0 )               
   v_x[ i ] = 1 - v_x[ i ];
  }  
 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 BKB->set_x( v_x.begin() );  // write the solution in BinaryKnapsackBlock

 }  // end( DPBinaryKnapsackSolver::get_var_solution )
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::set_par( idx_type par , double value ) {
 
 // Reoptimization parameter - - - - - - - - - - - - - - - - - - - - - - - -  
 if( par == dblReopt ) {
  
  if( ( value < 0 ) || ( value > 1 ) )
   throw( std::invalid_argument("dblReopt parameter must be in [ 0 , 1 ]") );
  
  // update reopt
  reopt = value;
    
  // restart from 0 in the next call of compute()
  start_item = 0;  
  
  return;
  }

 Solver::set_par( par , value );
 }


/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::add_Modification( sp_Mod & mod ) {

 if( f_no_Mod )
  return;

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 // if NBModification, reload BinaryKnapsack instance and clear modifications
 if( const auto tmod = std::dynamic_pointer_cast< NBModification >( mod ) ) {
  load();                   
  v_mod.clear();
  }
 else
  v_mod.push_back( mod );

 // release lock
 f_mod_lock.clear( std::memory_order_release );  

 }  // end( DPBinaryKnapsackSolver::add_Modification )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::load( void ) {
 
 if( ! f_Block ) {  // detaching the DPBinaryKnapsackSolver
  f_N = 0;
  v_P.clear();
  v_W.clear();
  v_I.clear();
  idxCont.clear();
  lab.clear();
  pred.clear();
  v_x.clear();
  v_fxd.clear();
  skip.clear();
  return;
  }

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( f_Block );
 if( ! BKB )
  throw( std::invalid_argument( "Block must be a BinaryKnapsackBlock" ) );

 // (try to) lock the BinaryKnapsackBlock
 bool owned = BKB->is_owned_by( f_id );
 if( ( ! owned ) && ( ! BKB->read_lock() ) )
  throw( std::runtime_error( "Unable to lock the Block" ) );
   
 // load Binary Knapsack instance - - - - - - - - - - - - - - - - - - - - - 
 
 // get scalar data: sense of the objective, number of items and capacity
 f_sense = ( BKB->get_objective_sense() == Objective::eMax );
 f_N = BKB->get_NItems();
 f_C = BKB->get_Capacity();

 // data for handling the continuous part
 idxCont.clear();  

 // resize all 
 v_P.resize( f_N );
 v_W.resize( f_N );
 v_I.resize( f_N );
 v_fxd.resize( f_N );

 // get references to profits, weights and integrality values
 const auto & P = BKB->get_Profits();
 const auto & W = BKB->get_Weights();
 const auto & I = BKB->get_Integrality();

 // read and store the data of each item
 for( Index i = 0 ; i < f_N ; ++i ) {

  // Profits: if the sense is minimization change the sign of each profit
  v_P[ i ] = f_sense ? P[ i ] : - P[ i ];

  // Weights: check if they are integers
  v_W[ i ] = std::round( W[ i ] );
  if( std::abs( v_W[ i ] - W[ i ] ) > WeightIntegrality )
   throw( std::invalid_argument( "Weights must be integers" ) );

  // Integrality: store indices of continuous variables
  v_I[ i ] = ( bool ) I[ i ];
  if( ! v_I[ i ] )
   idxCont.push_back( i );

  // check if the variables are fixed
  if( BKB->is_fixed( i ) )
   v_fxd[ i ] = std::abs( BKB->get_x( i ) ) < 1e-6 ? 1 : 2;
  else
   v_fxd[ i ] = 0;

  }

 if( ! owned )
  BKB->read_unlock();

 // end load Binary Knapsack instance - - - - - - - - - - - - - - - - - - - - 
 
 // clear previous data and initialize labels with a dummy node in the origin
 lab.clear();
 lab.resize( f_N  + 1 );
 lab[ 0 ] = { 0 };

 pred.clear();
 pred.resize( f_N  + 1 );

 // (re-)start from the first item
 start_item = 0;                       
 
 // compute step for reoptimization
 Index k = std::floor( reopt * std::log2( f_N ) );
 step = std::floor( f_N / std::exp2( k ) );

 // Initialize solution 
 v_x.resize( f_N );                          

 } // end( DPBinaryKnapsackSolver::load )

/*--------------------------------------------------------------------------*/

void DPBinaryKnapsackSolver::process_outstanding_Modification( void ) {

 // copy v_mod in a temporary list of modifications - - - - - - - - - - - - - 

 Lst_sp_Mod v_mod_tmp;              // temporary list of modifications

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );       // copy v_mod in v_mod_tmp

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 // Any changes in Profits must be processed only AFTER checking the changes 
 // on the sense of the Objective. Hence v_mod_tmp is scanned twice, checking
 // modifications on Objective sense (and also Capacity) first 

 auto mod = v_mod_tmp.begin(); 

 while( mod != v_mod_tmp.end() ) {
  
  // BinaryKnapsackBlockMod - - - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod->get() ) ) {
   switch( tmod->type() ) {
    // Change Capacity- - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

   // Change Objective Sense - - - - - - - - - - - - - - - - - - - - - - - -
   // Change the sign of all profits and restart from the first item
   case( BinaryKnapsackBlockMod::eChgSense ):
    f_sense = ( BKB->get_objective_sense() == Objective::eMax );
    // update f_sense

    for( Index i = 0 ; i < f_N ; ++i )     // change the sign of all profits
     v_P[ i ] = - v_P[ i ];

    start_item = 0;                        // restart from the beginning 
      
    mod = v_mod_tmp.erase( mod );
    break;

    default: mod++;
    }
   }
  else
   mod = v_mod_tmp.erase( mod );   // it is not a physical modification
 
  } // end( while )

 for( auto mod : v_mod_tmp ) {
  
  // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >(
							      mod.get() ) ) { 
   switch( tmod->type() ) {
    // change Profits - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // update modified profits according to f_sense and update start_item

    case( BinaryKnapsackBlockMod::eChgProfit ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ )
      v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );

     start_item = std::min( start_item , tmod->rng().first );
     break;

    // fix x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // update start item with the first fixed item
    case( BinaryKnapsackBlockMod::eFixX ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ) {
      // check if the variable is fixed to 0 or to 1
      if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) ) < 1e-6 )
       v_fxd[ i ] = 1;
      else if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) - 1 ) < 1e-6 )
       v_fxd[ i ] = 2; 
      }

     start_item = std::min( start_item , tmod->rng().first );
     break;
    
    // unfix x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // unfix a variable x could increase the residual capacity making
    // reoptimization not exploitable 
    case( BinaryKnapsackBlockMod::eUnfixX ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ )
      v_fxd[ i ] = 0;

     start_item = 0;
     break;

    // change Weights- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // update modified weights checking the integrality property. 
    // check also if a new weight increases the residual capacity making
    // reoptimization not exploitable 
    case( BinaryKnapsackBlockMod::eChgWeight ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ) {  
      double nw = std::round( BKB->get_Weight( i ) );     // new weight
 
      if( std::abs( nw - BKB->get_Weight( i ) ) > WeightIntegrality )
       throw( std::invalid_argument( "Weights must be integers!" ) );

      if( nw < v_W[ i ] )                  // it is not possible 
       start_item = 0;                     // to re-optimize

      v_W[ i ] = nw;                       // update weight
      }

     start_item = std::min( start_item , tmod->rng().first );
     break;
     
    // change Integrality - - - - - - - - - - - - - - - - - - - - - - - - - -
    // update modified integrality vector
    case( BinaryKnapsackBlockMod::eChgIntegrality ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ) {
      bool ni = BKB->get_Integrality( i );  // new integrality 
      v_I[ i ] = ni;                       // update
      }

     start_item = 0;
     break;
    }
   }

  // BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >(
							      mod.get() ) ) {
   switch( tmod->type() ) {
    // change Profits - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // update modified profits according to f_sense and update start_item 
    case( BinaryKnapsackBlockMod::eChgProfit ):
     for( auto i : tmod->nms() )
      v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );

     start_item = std::min( start_item , tmod->nms()[ 0 ] );
     break;
    
    // fix x- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // update start item with the first fixed item 
    case( BinaryKnapsackBlockMod::eFixX ):
     for( auto i : tmod->nms() ) {
      // check if the variable is fixed to 0 or to 1
      if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) ) < 1e-6 )
       v_fxd[ i ] = 1;
      else if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) - 1 ) < 1e-6 )
       v_fxd[ i ] = 2; 
      }

     start_item = std::min( start_item , tmod->nms()[ 0 ] );
     break;

    // unfix x- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // unfix a variable x could increase the residual capacity making
    // reoptimization not exploitable 
    case( BinaryKnapsackBlockMod::eUnfixX ):
     for( auto i : tmod->nms() )
      v_fxd[ i ] = 0;
     
     start_item = 0;
     break;

    // change Weights- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // update modified weights checking the integrality property. 
    // check also if a new weight increases the residual capacity making
    // reoptimization not exploitable 
    case( BinaryKnapsackBlockMod::eChgWeight ):
     for( auto i : tmod->nms() ) {
      double nw = std::round( BKB->get_Weight( i ) );     // new weight
  
      if( std::abs( nw - BKB->get_Weight( i ) ) > WeightIntegrality )
       throw( std::invalid_argument( "Weights must be integers!" ) );

      if( ( nw < 0 ) && ( nw < v_W[ i ] ) )  // it is not possible 
       start_item = 0;                       // to re-optimize

      v_W[ i ] = nw;                         // update weight
      }

     start_item = std::min( start_item , tmod->nms()[ 0 ] );
     break;
     
    // change Integrality- - - - - - - - - - - - - - - - - - - - - - - - - - 
    // update modified integrality vector
    case( BinaryKnapsackBlockMod::eChgIntegrality ):
     for( auto i : tmod->nms() ) {
      bool ni = BKB->get_Integrality( i );  // new integrality
      v_I[ i ] = ni;                        // update integrality
      }

     start_item = 0;                         // no-reoptimization     
     break;
    }
   }
  }  // end( for(  ) )

 v_mod_tmp.clear();  // clear the temporary list of Modification

 // store again indices of continuous variables
 idxCont.clear();
 for( Index i = 0 ; i < f_N ; ++i ) {
  if( ! v_I[ i ] ) {
   idxCont.push_back( i );
   } 
  }

 // compute start_item- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // Each Modification updated start_item with the first modified item. At this
 // point it is necessary to retrieve the first item smaller than start_item  
 // whose corresponding labels have been previously stored in lab[ i ]. 

 // first, compute step for reoptimization
 int k = std::floor( reopt * std::log2( f_N ) );
 int step = std::floor( f_N / std::exp2( k ) );

 if( start_item != Inf< int >() ) {
  
  start_item = ( start_item / step ) * step;
  
  // if all the variables are continuous
  if( idxCont.size() == f_N ) {
   start_item = f_N;
   lab.back() = { 0 };
   }
  }

 }  // end( DPBinaryKnapsackSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*----------------- End File DPBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
