/*--------------------------------------------------------------------------*/
/*--------------------- File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class GreedyRelaxationSolver, which
 * implements the RelaxationSolver concept [see Solver.h] for solving 
 * a relaxation of a  Knapsack problems as represented by a 
 * BinaryKnapsackKBlock.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Antonio Frangioni
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

// register GreedyRelaxationkSolver to the factory

SMSpp_insert_in_factory_cpp_1( GreedyRelaxationSolver );

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF GreedyRelaxationSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::set_Block( Block * block )
{
 if( block == f_Block )       // nothing to do  
  return;

 Solver::set_Block( block );  // attach to the new Block

 load();                      // load Binary Knapsack instance
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int GreedyRelaxationSolver::compute( bool changedvars )
{
 lock();  // lock the mutex

 // process all the pending modifications and compute start_item 

 process_outstanding_Modification(); 

 Return_OK:
 unlock();  // unlock the mutex

 return( kOK );

 }  // end( GreedyRelaxationSolver::compute() )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::get_var_solution( Configuration * solc )
{
 auto BKB = static_cast< BinaryKnapsackBlock * >( f_Block );

 // check if the solution has already been computed
 if( ! v_x.empty() ) {
  BKB->set_x( v_x.begin() );  // write the solution in BinaryKnapsackBlock 
  return;                     
  }

 // reconstruct the optimal solution - - - - - - - - - - - - - - - - - - - - -

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 BKB->set_x( v_x.begin() );  // write the solution in BinaryKnapsackBlock

 }  // end( GreedyRelaxationSolver::get_var_solution )
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::set_par( idx_type par , double value ) {
 Solver::set_par( par , value );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::add_Modification( sp_Mod & mod )
{
 if( f_no_Mod )
  return;

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 if( const auto tmod = std::dynamic_pointer_cast< NBModification >( mod ) ) {
  load();                   // the Binary Knapsack instance must be re-loaded
  v_mod.clear();
  }
 else
  v_mod.push_back( mod );

 f_mod_lock.clear( std::memory_order_release );  // release lock

 }  // end( add_Modification() )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::load( void ) {
 if( ! f_Block ) {  // detaching the GreedyRelaxationSolver
  f_N = 0;
  v_P.clear();
  v_W.clear();
  v_I.clear();
  indexContinuous.clear();
  v_x.clear();
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
 
 // get the sense of the objective
 f_sense = ( BKB->get_objective_sense() == Objective::eMax );
 
 f_N = BKB->get_NItems();    // get the number of items

 f_C = BKB->get_Capacity();  // get the Capacity

 v_P.resize( f_N );     

 const auto & P = BKB->get_Profits();      // get profits  
                                             
 for( Index i = 0 ; i < P.size() ; ++i )   // if the sense is minimization 
  v_P[ i ] = f_sense ? P[ i ] : - P[ i ];  // change the sign of the profits 

 v_W.resize( f_N );                        // prepare vector of weights

 const auto & W = BKB->get_Weights();

 for( Index i = 0 ; i < W.size() ; ++i ) {  // load weights and check 
  v_W[ i ] = std::round( W[ i ] );          // that they are integers
  if( std::abs( v_W[ i ] - W[ i ] ) > WeightIntegrality )
   throw( std::invalid_argument( "Weights must be integers" ) );
  }

 v_I.resize( f_N );    

 const auto & I = BKB->get_Integrality();
 for( Index i = 0 ; i < I.size() ; ++i )  // load values and check 
  v_I[ i ] = ( bool ) I[ i ] ;            // that they are booleans

 countCont = BKB->get_N_Cont_Items();

 if( ! owned )
  BKB->read_unlock();

 indexContinuous.clear();
 for( Index i = 0 ; i < f_N ; ++i )
  if( ! v_I[ i ] )
   indexContinuous.push_back( i ); 

 // end load Binary Knapsack instance- - - - - - - - - - - - - - - - - - - -

 v_x.clear();                          // clear previous solution (if any)

 } // end( GreedyRelaxationSolver::load() )

/*--------------------------------------------------------------------------*/

void GreedyRelaxationSolver::process_outstanding_Modification( void )
{
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

 v_mod_tmp.clear();  // clear the temporary list of Modification

 }  // end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/
/*----------------- End File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
