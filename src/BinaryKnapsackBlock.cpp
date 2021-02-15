/*--------------------------------------------------------------------------*/
/*---------------------- File BinaryKnapsackBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector<T> & vec , c_Subset & nms ,
		      typename std::vector<T>::const_iterator cmp ,
		      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another

template< typename T >
static void copyidx( std::vector<T> & vec , c_Subset & nms ,
		     typename std::vector<T>::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BinaryKnapsackBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( BinaryKnapsackBlock );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF BinaryKnapsackBlock ----------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( Index n , double Capacity , 
               std::vector<double> & Weights , std::vector<double> & Profits )
{
 
 // sanity checks 

 if( Weights.size() != n )
  throw( std::invalid_argument( "Vector of Weights of the wrong size" ) );

 if( Profits.size() != n )
  throw( std::invalid_argument( "Vector of Profits of the wrong size" ) );

 // erase previous instance, if any

 if( get_NItems() )
  guts_of_destructor();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
 NItems = n;
 C = Capacity;

 W.resize( NItems );  
 P.resize( NItems ); 

 std::copy( Weights.begin() , Weights.end() , W.begin() );
 std::copy( Profits.begin() , Profits.end() , P.begin() );

 generate_abstract_variables();

 // reset conditional bounds
 f_cond_lower = -Inf<double>();
 f_cond_upper = Inf<double>();

 // Modification

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

} // end( BinaryKnapsackBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::deserialize( const netCDF::NcGroup & group ){

 // erase previous instance, if any
 if( get_NItems() )
  guts_of_destructor();

// read problem data 
 netCDF::NcDim ni = group.getDim( "NItems" );
 if( ni.isNull() )
  throw( std::logic_error( "NItems dimension is required" ) );
 NItems = ni.getSize();

 netCDF::NcDim c = group.getDim( "Capacity" );
 if( c.isNull() )
  throw( std::logic_error( "Capacity is required" ) );
 C = c.getSize();

 netCDF::NcVar w = group.getVar( "Weights" );
 if( w.isNull() )
  throw( std::logic_error( "Weights are required" ) );
 
 W.resize( NItems );
 w.getVar( W.data() );

 netCDF::NcVar p = group.getVar( "Profits" );
 if( p.isNull() )
  throw( std::logic_error( "Profits are required" ) );
 
 P.resize( NItems );
 p.getVar( P.data() );

 generate_abstract_variables();

 // reset conditional bounds
 f_cond_lower = -Inf<double>();
 f_cond_upper = Inf<double>();

 // call the method of Block
 // inside this the NBModification, the "nuclear option",  is issued

 Block::deserialize( group );

} // end( BinaryKnapsackBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_abstract_variables( Configuration * stvv )
{
 if( AR & HasVar )  // the variables are there already
  return;           // nothing to do
  
 x.resize( get_NItems() );
   
  for( auto & var : x )
   var.set_type( ColVariable::kBinary , eNoBlck );
   
 add_static_variable( x );

 AR |= HasVar;
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( AR & HasCns )  // the constraint is there already
  return;           // nothing to do

 LinearFunction::v_coeff_pair w( get_NItems() );
 for( Index i = 0 ; i < get_NItems() ; i++ ){
  w[ i ].first = &x[i];
  w[ i ].second = W[i];
 }

 cnst.set_function( new LinearFunction( std::move( w ) , 0 ) , eNoBlck );
 cnst.set_rhs( C );
 cnst.set_lhs( - Inf<double>() );
 
 add_static_constraint( cnst );

 AR |= HasCns;
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_objective( Configuration * objc ){
 
 if( AR & HasObj )  // the objective is there already
  return;           // nothing to do 

 LinearFunction::v_coeff_pair p( get_NItems() );
  for( Index i = 0 ; i < get_NItems() ; i++ ){
   p[ i ].first = &x[i];
   p[ i ].second = P[i];
  }

 LinearFunction * obj = new LinearFunction( std::move( p ) , 0  );
 f.set_function( obj , eNoMod );
 f.set_sense( FRealObjective::eMax , eNoMod );
 set_objective( & f , eNoMod );
 
 AR |= HasObj;
}

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool BinaryKnapsackBlock::is_feasible( bool useabstract , 
                                      Configuration * fsbc ){
 if( useabstract ) {
  // do it using the abstract representation
  if( ! ( AR & HasCns ) )
   throw( std::logic_error( "Constraint required for is_feasible( true , )") );
   
  return( cnst.rel_viol() > 0 ? false : true );
}

// do it using the physical representation
 double tot_weight = 0; 
 for( Index i = 0 ; i < W.size() ; i++ )
  tot_weight += W[ i ] * x[ i ].get_value();

 return( tot_weight <= C ? true : false );    
}

/*--------------------------------------------------------------------------*/

bool BinaryKnapsackBlock::is_empty( bool useabstract , Configuration * optc ){
   
 if( C >= 0 )       // if the Capacity is positive then x = 0 is a feasible 
  return( false );  // solution and the problem is not empty

 double neg_weights = 0;

 for( Index i = 0 ; i < W.size() ; i++ ){
  if( W[ i ] < 0 ){
   neg_weights += W[ i ];
   if( neg_weights <= C )
    return( false );
  }
 }
 return( true );
}

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * BinaryKnapsackBlock::get_Solution( Configuration *solc , 
                                             bool emptys ){
 auto * sol = new BinaryKnapsackSolution(); 

 sol->v_x.resize( get_NItems() );

 if( ! emptys )
  sol->read( this );

 return( sol ); 
}

/*--------------------------------------------------------------------------*/

bool BinaryKnapsackBlock::get_x( c_Index item ){
 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );
 return x[ item ].get_value(); 
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::get_x( std::vector<int> & xSol , Range rng ){
 
 rng.second = std::min( rng.second , get_NItems() );

 auto xSoli = xSol.begin();
 for( Index i = rng.first ; i < rng.second ; i++ )
  ( * xSoli++ ) = x[ i ].get_value();

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::get_x( std::vector<int> & xSol , c_Subset & nms ){
 
 auto xSoli = xSol.begin();
 for( auto i : nms ){
  if( i >= get_NItems() )
   throw( std::invalid_argument( "invalid item" ) );
  ( * xSoli++ ) = x[ i ].get_value();
 }
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( c_Index item , bool value ){
 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );
 x[ item ].set_value( value ); 
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( std::vector<int> & xSol , Range rng ){
 
 rng.second = std::min( rng.second , get_NItems() );

 auto xSoli = xSol.begin();
 for( Index i = rng.first ; i < rng.second ; i++ )
  x[ i ].set_value( * xSoli++ );

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( std::vector<int> & xSol , c_Subset & nms ){
 
 auto xSoli = xSol.begin();
 for( auto i : nms ){
  if( i >= get_NItems() )
   throw( std::invalid_argument( "invalid item" ) );
  x[ i ].set_value( * xSoli++ );
 }

}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::add_Modification( sp_Mod mod , ChnlName chnl ){
 
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*------ METHODS FOR LOADING, PRINTING & SAVING THE BinaryKnapsackBlock ----*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::serialize( netCDF::NcGroup & group ) const {
 
 // call the method of Block
 Block::serialize( group );

// BinaryKnapsackBlock data 
 
 netCDF::NcDim ni = group.addDim( "NItems" , get_NItems() );
 
 group.addDim( "Capacity" , get_Capacity() );
 
 ( group.addVar( "Weights" , netCDF::NcDouble() , ni ) ).putVar( W.data() );
 
 ( group.addVar( "Profits" , netCDF::NcDouble() , ni ) ).putVar( P.data() );

}// end( BinaryKnapsackBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( const bool value, c_Index item , 
                                 c_ModParam issueMod, c_ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid index item" ) );

 if( ( x[ item ].is_fixed() ) && ( x[ item ].get_value() == value ) )// nothing 
  return;                                                            // to do

 // reset conditional bounds
 f_cond_lower = -Inf<double>();
 f_cond_upper = Inf<double>();

if( not_dry_run( issueAMod ) ){
 x[ item ].set_value( value );
 x[ item ].is_fixed( true , issueAMod ); 
}
 
 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
          BinaryKnapsackBlockMod::eFixX , std::make_pair( item , item + 1 ) ) , 
          Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( c_Index item , c_ModParam issueMod, 
                                   c_ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid index item" ) );

 if( !( x[ item ].is_fixed() ) ) // already unfixed
  return;                        // nothing to do

 // reset conditional bounds
 f_cond_lower = -Inf<double>();
 f_cond_upper = Inf<double>();
 
 if( not_dry_run( issueAMod ) )
  x[ item ].is_fixed( false , issueAMod ); 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
          BinaryKnapsackBlockMod::eUnfixX , std::make_pair( item , item + 1 ) ) , 
          Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weight( const double NWeight , c_Index item , 
                          c_ModParam issueMod , c_ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );

 if( W[ item ] == NWeight )  // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  LinearFunction *lf = dynamic_cast<LinearFunction*>( cnst.get_function() );
  lf->modify_coefficient( item , NWeight , issueAMod ); // AR
  W[ item ] = NWeight;                                  // PR
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation
  W[ item ] = NWeight; 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
      BinaryKnapsackBlockMod::eChgWeight , std::make_pair( item , item + 1 ) ), 
      Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::chg_weight )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weights( const std::vector<double>::iterator NWeight,
                                       Range rng , 
                                       c_ModParam issueMod ,
                                       c_ModParam issueAMod ){

 rng.second = std::min( rng.second , get_NItems() );
 if( rng.second <= rng.first )  // nothing to change
  return;   

 if( std::equal( NWeight , NWeight + ( rng.second - rng.first ) ,
		 W.begin() + rng.first ) )
  return;  // nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  
  // physical representation
  std::copy( NWeight , NWeight + ( rng.second - rng.first ) ,
             W.begin() + rng.first );
  
  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>( cnst.get_function() );
  lf->modify_coefficients( std::vector<double>( NWeight , 
                  NWeight + ( rng.second - rng.first ) ) , rng , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  std::copy( NWeight , NWeight + ( rng.second - rng.first ) ,
             W.begin() + rng.first );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                                    BinaryKnapsackBlockMod::eChgWeight , rng ), 
                                    Observer::par2chnl( issueMod ) );

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weights( const std::vector<double>::iterator NWeight,
                                       Subset && nms , const bool ordered ,  
                                       c_ModParam issueMod ,
                                       c_ModParam issueAMod ){

 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( W , nms , NWeight , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  
  // physical representation
  copyidx( W , nms , NWeight );
  
  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>( cnst.get_function() );
  lf->modify_coefficients( std::vector<double>( NWeight , 
                  NWeight + nms.size() ) , 
                  Subset( nms ) , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  copyidx( W , nms , NWeight );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ){ 
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(this,
                      BinaryKnapsackBlockMod::eChgWeight , std::move( nms ) ), 
                      Observer::par2chnl( issueMod ) );
 }

} 

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profit( const double NProfit , c_Index item , 
                          c_ModParam issueMod , c_ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );

 if( P[ item ] == NProfit ) // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  P[ item ] = NProfit; 

  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>( f.get_function() );
  lf->modify_coefficient( item , NProfit , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation 
  P[ item ] = NProfit;
 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
      BinaryKnapsackBlockMod::eChgProfit , std::make_pair( item , item + 1 ) ), 
      Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::chg_profit )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profits( const std::vector<double>::iterator NProfit ,
                                       Range rng , 
                                       c_ModParam issueMod ,
                                       c_ModParam issueAMod ){

 rng.second = std::min( rng.second , get_NItems() );
 if( rng.second <= rng.first )  // nothing to change
  return;   

 if( std::equal( NProfit , NProfit + ( rng.second - rng.first ) ,
		 P.begin() + rng.first ) )
  return;  // nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  std::copy( NProfit , NProfit + ( rng.second - rng.first ) ,
             P.begin() + rng.first );
  
  // abstract representation  
  LinearFunction *lf = dynamic_cast<LinearFunction*>( f.get_function() );
  lf->modify_coefficients( std::vector<double>( NProfit , 
                    NProfit + ( rng.second - rng.first ) ) , rng , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){ // otherwise only physical representation 
  std::copy( NProfit , NProfit + ( rng.second - rng.first ) ,
             P.begin() + rng.first );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                                    BinaryKnapsackBlockMod::eChgProfit , rng ), 
                                    Observer::par2chnl( issueMod ) );

} 

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profits( const std::vector<double>::iterator NProfit,
                                       Subset && nms , const bool ordered ,  
                                       c_ModParam issueMod ,
                                       c_ModParam issueAMod ){

 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( P , nms , NProfit , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  copyidx( P , nms , NProfit );
  
  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>( f.get_function() );
  lf->modify_coefficients( std::vector<double>( NProfit , 
                  NProfit + nms.size() ) , 
                  Subset( nms ) , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  copyidx( P , nms , NProfit );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ){ 
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(this,
                      BinaryKnapsackBlockMod::eChgProfit , std::move( nms ) ), 
                      Observer::par2chnl( issueMod ) );
 }
} 

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_capacity( const double NC , 
                          c_ModParam issueMod , c_ModParam issueAMod ){

 if( C == NC ) // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  C = NC; 
  cnst.set_rhs( NC , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation 
  C = NC;
 

 if( issue_pmod( issueMod ) ) // issue physical Modification 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockMod >( this ,
                           BinaryKnapsackBlockMod::eChgCapacity ) , 
                           Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::chg_capacity )


/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::print( std::ostream & output ) const {
 
 output << "BinaryKnapsackBlock\n";
 output << "Number of items: " << NItems << std::endl; 
 output << "Capacity: " << C << std::endl;
 
 std::cout << "\tWeights\tProfits\n";
 for( int i = 0 ; i < NItems ; i++ )
  std::cout << "Item " << i << "\t" << W[i] << "\t" << P[i] << std::endl;

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( std::istream & input ){ 

// erase previous instance, if any

 if( get_NItems() )
  guts_of_destructor();

// read problem data
if( !( input >> NItems ) )
 throw( std::invalid_argument( "error reading number of items" ) );

if( !( input >> C ) )
 throw( std::invalid_argument( "error reading Capacity" ) );

W.resize( NItems );
P.resize( NItems );

for( Index i = 0 ; i < NItems ; i++ ){
 if( !( input >> W[i] ) )
  throw( std::invalid_argument( "error reading Weights" ) );
}

for( Index i = 0 ; i < NItems ; i++ ){
 if( !( input >> P[i] ) )
  throw( std::invalid_argument( "error reading Profits" ) );
}

 generate_abstract_variables();

 // Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

}

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::guts_of_destructor(){

 // clear the constraint
 cnst.clear();
 // clear the objective function
 f.clear();
 // clear all variables
 x.clear();

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // explicitly reset Constraint and Variables
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 AR = 0;
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::guts_of_add_Modification(p_Mod mod , ChnlName chnl){

// C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
{
 const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod );
 if( tmod ){ 
// the change may concern the objective function or the constraint
// check which one it is
  auto lf = dynamic_cast< LinearFunction * >( tmod->function() );
  if( ! lf )
   throw( std::invalid_argument( "invalid Modification to Linear Function" ) ); 

  int range_size = tmod->range().second - tmod->range().first;

  // if the AR of the objective exists
  if( AR & HasObj ){ 
   // get the objective function
   auto lfo = dynamic_cast< LinearFunction * >( f.get_function() ); 
   // and check if the modification is on the objective  
    if( lf == lfo ){ 
     // vector of new profits 
     std::vector<double> new_profits( range_size );  
     std::vector<double>::iterator npi = new_profits.begin(); 
     // get the values of the new profits from coefficients of lf
     for( Index i = tmod->range().first ; i < tmod->range().second ; i++ )
      ( * npi++ ) = lf->get_coefficient( i );

     chg_profits( new_profits.begin() , tmod->range() ,
                  make_par( eNoBlck , chnl ) , eDryRun );
     return;
    }   
  }

  // if the AR of the constraint exists
  if( AR & HasCns ){ 
   // get the constraint 
   auto lfc = dynamic_cast< LinearFunction * >( cnst.get_function() ); 
   // and check if the modification is on the constraint  
    if( lf == lfc ){ 
     // vector of new weights 
     std::vector<double> new_weights( range_size ); 
     std::vector<double>::iterator nwi = new_weights.begin();
     // get the values of the new weights from coefficients of lf
     for( Index i = tmod->range().first ; i < tmod->range().second ; i++ )
      ( * nwi++ ) = lf->get_coefficient( i );

     chg_weights( new_weights.begin() , tmod->range() ,
                  make_par( eNoBlck , chnl ) , eDryRun );
     return;
    }   
  }

 throw( std::invalid_argument( "illegal Modification to Linear Function" ) );
}
}

// C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
{
 const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod );
 if( tmod ){
// the change may concern the objective function or the constraint
// check which one it is
  auto lf = dynamic_cast< LinearFunction * >( tmod->function() );
  if( ! lf )
   throw( std::invalid_argument( "invalid Modification to Linear Function" ) ); 

  // if the AR of the objective exists
  if( AR & HasObj ){
   // get the objective function
   auto lfo = dynamic_cast< LinearFunction * >( f.get_function() ); 
   // and check if the modification is on the objective  
    if( lf == lfo ){
     // vector of new profits 
     std::vector<double> new_profits( tmod->subset().size() );
     std::vector<double>::iterator npi = new_profits.begin();   
     // get the values of the new profits from coefficients of lf
     for( auto i : tmod->subset() )
      ( * npi++ ) = lf->get_coefficient( i );

     chg_profits( new_profits.begin() , Subset( tmod->subset() ) , true , 
                  make_par( eNoBlck , chnl ) , eDryRun );
     return;
    }   
  }

  // if the AR of the constraint exists
  if( AR & HasCns ){
   // get the constraint 
   auto lfc = dynamic_cast< LinearFunction * >( cnst.get_function() ); 
   // and check if the modification is on the constraint  
    if( lf == lfc ){
     // vector of new weights 
     std::vector<double> new_weights( tmod->subset().size() );  
     std::vector<double>::iterator nwi = new_weights.begin();
     // get the values of the new weights from coefficients of lf
     for( auto i : tmod->subset() )
      ( * nwi++ ) = lf->get_coefficient( i );

     chg_weights( new_weights.begin() , Subset( tmod->subset() ) , true , 
                  make_par( eNoBlck , chnl ) , eDryRun );
     return;
    }   
  }

 throw( std::invalid_argument( "illegal Modification to Linear Function" ) );
}
}

// RowConstraintMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
{
 const auto tmod = dynamic_cast< RowConstraintMod * >( mod );
 if( tmod ){
  if( !( AR & HasCns ) )
   throw( std::invalid_argument("modification to non-constructed constraint"));
  
   auto cp = dynamic_cast< FRowConstraint * >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) ); 

  if( tmod->type() == RowConstraintMod::eChgRHS ){
   chg_capacity( cp->get_rhs(), make_par( eNoBlck , chnl ) , eDryRun );
   return; 
  }
  throw( std::invalid_argument( "illegal Modification to Constraint" ) );
 }
}

// VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
{
 const auto tmod = dynamic_cast< VariableMod * >( mod );
 if( tmod ){ 
  auto xi = dynamic_cast< ColVariable * const >( tmod->variable() );
  if( ! xi )
   throw( std::logic_error( "Modification to wrong type of Variable" ) );
  
  int i = p2i_x( xi );
  if( xi->is_fixed() )
   unfix_x( i , make_par( eNoBlck , chnl ) , eDryRun );
  else 
   fix_x( xi->get_value() , i , make_par( eNoBlck , chnl ) , eDryRun );
  
  return;
}

}

throw(std::invalid_argument( "unsupported Modification to BinaryKnapsackBlock"));
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::compute_conditional_bounds(){

 f_cond_lower = 0;
 f_cond_upper = 0;
 
 for( Index i = 0 ; i < get_NItems() ; i++ ){
  // items contained in the optimal solution 
  if( ( W[ i ] < 0 ) && ( P[ i ] > 0 ) ){ 
   f_cond_lower += P[ i ];                
   f_cond_upper += P[ i ];
   continue;
  }
  // items not contained in the optimal solution
  if( ( W[ i ] > 0 ) && ( P[ i ] < 0 ) )
   continue;
  // remaining items
  if( P[ i ] > 0 )
   f_cond_upper += P[ i ];
  else 
   f_cond_lower += P[ i ];
 }
 
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF BinaryKnapsackSolution --------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::deserialize( const netCDF::NcGroup & group ){
 
 netCDF::NcDim ni = group.getDim( "n" );
 if( ni.isNull() )
  v_x.clear();
 else{
  netCDF::NcVar bx = group.getVar( "x" );
  if( bx.isNull() )
   v_x.clear();
  else{
   v_x.resize( ni.getSize() );
   bx.getVar( v_x.data() );
  }
 }

} // end( BinaryKnapsackSolution::deserialize )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::read( const Block * const block ){

 auto BKB = dynamic_cast< const BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "block is not a BinaryKnapsackBlock" ));

 v_x.resize( BKB->get_NItems() );

 auto vxi = v_x.begin();

 for( auto & xi : BKB->x ) 
  *( vxi++ ) = xi.get_value(); 

} // end( BinaryKnapsackSolution::read )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::write( Block * const block ){
 
 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "block is not a BinaryKnapsackBlock" ));

 // write binary variables - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < BKB->get_NItems() )
   throw( std::invalid_argument( "incompatible variables size" ) );

  auto vxi = v_x.begin();

  for( auto & xi : BKB->x )
   xi.set_value( *( vxi++ ) );
 }
} // end( BinaryKnapsackSolution::write )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::serialize( netCDF::NcGroup & group ){
 if( ! v_x.empty() ) {
  netCDF::NcDim ni = group.addDim( "n" , v_x.size() ); 
( group.addVar( "x" , netCDF::NcInt() , ni ) ).putVar( v_x.data() );
 }
} // end( BinaryKnapsackSolution::serialize )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::print( std::ostream & output ){
  for(Index i = 0 ; i < v_x.size() ; i++ ){
    output << "x" << i << ": " << v_x[i] << " ";
  }
}

/*--------------------------------------------------------------------------*/

BinaryKnapsackSolution * BinaryKnapsackSolution::scale( double factor ) const
{
 auto * sol = BinaryKnapsackSolution::clone( true );

 if( ! v_x.empty() )
  for( Index i = 0 ; i < v_x.size() ; i++ )
   sol->v_x[ i ] = v_x[ i ] * factor;

 return( sol );
} // end( BinaryKnapsackSolution::scale )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::sum( const Solution * solution, 
                                  double multiplier ){
 auto BKB = dynamic_cast< const BinaryKnapsackSolution * >( solution );
 if( ! BKB )
  throw( std::invalid_argument("solution is not a BinaryKnapsackSolution") );
 
 if( ! v_x.empty() ) {
  if( v_x.size() != BKB->v_x.size() )
   throw( std::invalid_argument( "incompatible variables size" ) );
  
  for( Index i = 0 ; i < v_x.size() ; i++ )
   v_x[ i ] = BKB->v_x[ i ] * multiplier;
  }
} // end( BinaryKnapsackSolution::sum )

/*--------------------------------------------------------------------------*/

BinaryKnapsackSolution * BinaryKnapsackSolution::clone( bool empty ) const{
 auto * sol = new BinaryKnapsackSolution();
 
 if( empty ){
  if( ! v_x.empty() )
   sol->v_x.resize( v_x.size() );
 }
 else 
  sol->v_x = v_x;
 
 return( sol );
} // end( BinaryKnapsackSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File BinaryKnapsackBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/