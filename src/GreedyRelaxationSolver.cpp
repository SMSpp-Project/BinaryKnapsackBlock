/*--------------------------------------------------------------------------*/
/*--------------------- File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class GreedyRelaxationSolver, which
 * implements the Solver concept [see Solver.h] for solving the continuous
 * relaxation of a Knapsack problem as represented by a BinaryKnapsackBlock.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Federica Di Pasquale, Donato Meoli
 */ 
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "GreedyRelaxationSolver.h"

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

// register GreedyRelaxationSolver to the factory

SMSpp_insert_in_factory_cpp_1( GreedyRelaxationSolver );

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF GreedyRelaxationSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::set_Block( Block * block ) {
 
 if( block == f_Block )       // nothing to do        
  return;

 Solver::set_Block( block );  // attach to the new Block

 load();                      // load Binary Knapsack instance
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int GreedyRelaxationSolver::compute( bool changedvars ) {
 
 lock();                     // lock the mutex

 // Process all the pending modifications 
 process_outstanding_Modification();

 // Initialize objective and solution
 obj = - Inf< double >(); 
 v_x.assign( f_N , 0 );

 // Preprocessing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 double C = f_C;    // residual capacity
 double P = 0;      // residual profit

 // item to skip
 std::vector< bool > skip( f_N , false );

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

 // check if the problem is empty (iff residual capacity < 0)
 if( C < 0 ) { unlock(); return( kInfeasible ); } 
                 
 // solve the continuous knapsack - - - - - - - - - - - - - - - - - - - - - -
 
 // indeces of items to sort in order of profit/weight
 Subset idx( f_N );
 std::iota( idx.begin() , idx.end() , 0 );

 // sort indeces in order of profit/weight
 sort( idx.begin() , idx.end() , [ & ]( auto a , auto b ) { 
       return ( v_P[ a ] / v_W[ a ] > v_P[ b ] / v_W[ b ] ); 
      } );

 double tot_W = 0;              // initiliaze total weight
 double tot_P = P;              // initialize total profit

 f_ci = idx.back();             // index of the critical item
 f_ciVal = 1;                   // variable value of the critical item
 
 for( Index i : idx ) {

  if( skip[ i ] )               // skip preprocessed items
   continue;
  
  double w = v_W[ i ];          // weight of the current item
  double p = v_P[ i ];          // profit of the current item
  if( w < 0 && p < 0 ) {
   w = -w;
   p = -p;
   }  
  
  if( tot_W + w > C ) {
   // take only a fraction
   f_ci = i;
   f_ciVal = ( C - tot_W ) / w; 
   tot_W = C;
   tot_P += f_ciVal * p;
   break; 
   }

  tot_W += w;                   // updata total weight
  tot_P += p;                   // update total profit
 }  

 // update solution - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 double xval = 1;
 for( Index i : idx ) {
  
  if( skip[ i ] )
   continue;
  
  if( i == f_ci ) {
   v_x[ i ] = f_ciVal;
   xval = 0;
   }
  else
   v_x[ i ] = xval; 
 
  if( v_W[ i ] < 0 && v_P[ i ] < 0 ) {
   v_x[ i ] = 1 - v_x[ i ];
   if( i == f_ci )
    f_ciVal = 1 - f_ciVal;
  }

 }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 obj = tot_P;               // update objective value

 unlock();                  // unlock the mutex
 
 return( kOK );

 }  // end( GreedyRelaxationSolver::compute() )

/*--------------------------------------------------------------------------*/

std::vector< Change * > GreedyRelaxationSolver::branch() {

 // branch on the critical item
 std::vector< Change * > branches( 2 );
 
 std::vector< double > zero = { 0 };
 std::vector< double > one = { 1 };
 
 branches[ 0 ] = new BinaryKnapsackBlockRngdChange( 
                                   BinaryKnapsackBlockChange::eFixX , 
                                   std::move( zero ) , 
                                   std::make_pair( f_ci , f_ci + 1 ) );

 branches[ 1 ] = new BinaryKnapsackBlockRngdChange( 
                                   BinaryKnapsackBlockChange::eFixX , 
                                   std::move( one ) , 
                                   std::make_pair( f_ci , f_ci + 1 ) );

 return branches;
}

/*--------------------------------------------------------------------------*/

Change * GreedyRelaxationSolver::apply( Change * chg , bool doUndo ) 
{
 return chg->apply( f_Block , doUndo );

 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool GreedyRelaxationSolver::has_true_var_solution( void ) { return false; }

/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::get_var_solution( Configuration * solc ) {
 
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block ); 
  
 BKB->set_x( v_x.begin() );  // write the solution in BinaryKnapsackBlock

 }  // end( GreedyRelaxationSolver::get_var_solution )

 /*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::add_Modification( sp_Mod & mod ) {

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

 }  // end( add_Modification() )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::load() {

 if( ! f_Block ) {
  f_N = 0;
  v_P.clear();
  v_W.clear();
  v_x.clear();
  v_fxd.clear();
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

 // resize all 
 v_P.resize( f_N );
 v_W.resize( f_N );
 v_fxd.resize( f_N );

 // get references to profits, weights and integrality values
 const auto & P = BKB->get_Profits();
 const auto & W = BKB->get_Weights();
 const auto & FXD = BKB->get_fxd();

 for( Index i = 0 ; i < f_N ; ++i ) {
  v_P[ i ] = f_sense ? P[ i ] : - P[ i ];  
  v_W[ i ] = W[ i ];
  v_fxd[ i ] = FXD[ i ];
  } 
 
 if( ! owned )
  BKB->read_unlock();

} // end( GreedyRelaxationSolver::load() )

/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::process_outstanding_Modification() {

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
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod->get()
                                  ) ){ 
   switch( tmod->type() ) {

    case( BinaryKnapsackBlockMod::eChgCapacity ): {
     f_C = BKB->get_Capacity();    
     mod = v_mod_tmp.erase( mod );
     break;
     }

    case( BinaryKnapsackBlockMod::eChgSense ):
     f_sense = ( BKB->get_objective_sense() == Objective::eMax );
     for( Index i = 0 ; i < f_N ; ++i )     
      v_P[ i ] = - v_P[ i ];
     mod = v_mod_tmp.erase( mod );
     break;

     default: mod++;
    }
   }
  else
   mod = v_mod_tmp.erase( mod );   // it is not a physical modification
 
  } // end( while() )

 for( auto mod : v_mod_tmp ) {
  
  // BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockRngdMod * >(
                                  mod.get() ) ) { 
   switch( tmod->type() ) {
    
    case( BinaryKnapsackBlockMod::eChgProfit ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ )
      v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );
     break;

    case( BinaryKnapsackBlockMod::eFixX ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ) {
      // check if the variable is fixed to 0 or to 1
      if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) ) < 1e-6 )
       v_fxd[ i ] = 1;
      else if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) - 1 ) < 1e-6 )
       v_fxd[ i ] = 2; 
      }
     break;

    case( BinaryKnapsackBlockMod::eUnfixX ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ )
      v_fxd[ i ] = 0;
     break;

    case( BinaryKnapsackBlockMod::eChgWeight ):
     for( Index i = tmod->rng().first ; i < tmod->rng().second ; i++ ) {  
      v_W[ i ] =  BKB->get_Weight( i );     
      }
     break;
    }
   }

  // BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = dynamic_cast< BinaryKnapsackBlockSbstMod * >(
                                  mod.get() ) ) {
   switch( tmod->type() ) {

    case( BinaryKnapsackBlockMod::eChgProfit ):
     for( auto i : tmod->nms() )
      v_P[ i ] = f_sense ? BKB->get_Profit( i ) : - BKB->get_Profit( i );
      break;

    case( BinaryKnapsackBlockMod::eFixX ):
     for( auto i : tmod->nms() ) {
      // check if the variable is fixed to 0 or to 1
      if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) ) < 1e-6 )
       v_fxd[ i ] = 1;
      else if( BKB->is_fixed( i ) && std::abs( BKB->get_x( i ) - 1 ) < 1e-6 )
       v_fxd[ i ] = 2; 
      }
     break;

    case( BinaryKnapsackBlockMod::eUnfixX ):
     for( auto i : tmod->nms() )
      v_fxd[ i ] = 0;

     break;

    case( BinaryKnapsackBlockMod::eChgWeight ):
     for( auto i : tmod->nms() ) {
      v_W[ i ] = BKB->get_Weight( i );     // new weight
      }
     break;
    }
   }
  }  // end( for(  ) )

 v_mod_tmp.clear();  // clear the temporary list of Modification

 }  // end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/
/*----------------- End File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
