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

// register BinaryKnapsackSolution to the Solution factory

SMSpp_insert_in_factory_cpp_1( BinaryKnapsackSolution );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF BinaryKnapsackBlock ----------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( Index n , double Capacity , 
                                const std::vector<double> & Weights , 
                                const std::vector<double> & Profits,
                                const std::vector<bool> & Integrality )
{
 
 // sanity checks 

 if( Weights.size() != n )
  throw( std::invalid_argument( "Vector of Weights of the wrong size" ) );

 if( Profits.size() != n )
  throw( std::invalid_argument( "Vector of Profits of the wrong size" ) );
 
 if( Integrality.size() != n && Integrality.size()!=0 )
   throw( std::invalid_argument( "Vector of Integrality of the wrong size" ) );
 
 // copy vectors and call load( , , && , && )
 std::vector<double> W( n );
 std::vector<double> P( n );
 std::vector<bool>   I( n );

 std::copy( Weights.begin() , Weights.end() , W.begin() );
 std::copy( Profits.begin() , Profits.end() , P.begin() );
 
 
 if(Integrality.size()==0){
 load( n , Capacity , std::move( W ) , std::move( P ));
 }else{
 std::copy( Integrality.begin() , Integrality.end() , I.begin() );
 load( n , Capacity , std::move( W ) , std::move( P ) ,std::move( I ));
 }
 
} // end( BinaryKnapsackBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( Index n , double Capacity , 
                                std::vector<double> && Weights , 
                                std::vector<double> && Profits ,
                                std::vector<bool> && Integrality  )
{
 
 // sanity checks 

 if( Weights.size() != n )
  throw( std::invalid_argument( "Vector of Weights of the wrong size" ) );

 if( Profits.size() != n )
  throw( std::invalid_argument( "Vector of Profits of the wrong size" ) );

 std::cout<<Integrality.size();
 if( Integrality.size() != n && Integrality.size()!=0 )
   throw( std::invalid_argument( "Vector of Integrality of the wrong size" ) );

 // erase previous instance, if any

 if( get_NItems() )
  guts_of_destructor();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
 f_N = n;
 f_C = Capacity;

 v_W = std::move( Weights );  
 v_P = std::move( Profits ); 
 v_I = std::move( Integrality );  
 
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
 f_N = ni.getSize();

 netCDF::NcDim c = group.getDim( "Capacity" );
 if( c.isNull() )
  throw( std::logic_error( "Capacity is required" ) );
 f_C = c.getSize();

 netCDF::NcVar w = group.getVar( "Weights" );
 if( w.isNull() )
  throw( std::logic_error( "Weights are required" ) );
 
 v_W.resize( f_N );
 w.getVar( v_W.data() );

 netCDF::NcVar p = group.getVar( "Profits" );
 if( p.isNull() )
  throw( std::logic_error( "Profits are required" ) );
 
 v_P.resize( f_N );
 p.getVar( v_P.data() );
 
 netCDF::NcVar i = group.getVar( "Integrality" );
 if( i.isNull() )
  v_I.clear();
 else {
  v_I.resize( f_N );
  std::vector< unsigned char > tmpI( v_I.size() );
 i.getVar( tmpI.data() );
 for( Index i = 0 ; i < v_I.size() ; ++i )
  v_I[ i ] = ( tmpI[ i ] != 0 );
  }
 
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
 
 if( ! get_NItems() )
  return;

 v_x.resize( get_NItems() );
   
  countCont=0;
  if( v_I.empty()){
  for(int i=0; i< get_NItems(); i++ )
     v_x[i].set_type( ColVariable::kBinary , eNoBlck );
  }
  else{
    for(int i=0; i< get_NItems(); i++ ){
     if(v_I[i]==true)
       v_x[i].set_type( ColVariable::kBinary , eNoBlck );
     else{
       v_x[i].set_type( ColVariable::kPosUnitary , eNoBlck );
       countCont++;  
      }
    }
  }
 add_static_variable( v_x );

 AR |= HasVar;
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( AR & HasCns )  // the constraint is there already
  return;           // nothing to do

 if( ! get_NItems() )
  return;

 v_cnst.resize( 1 );

 LinearFunction::v_coeff_pair w( get_NItems() );
 for( Index i = 0 ; i < get_NItems() ; i++ ){
  w[ i ].first = &v_x[ i ];
  w[ i ].second = v_W[ i ];
 }

 v_cnst[ 0 ].set_function( new LinearFunction( std::move( w ) , 0 ) , eNoBlck );
 v_cnst[ 0 ].set_rhs( f_C );
 v_cnst[ 0 ].set_lhs( - Inf<double>() );
 
 add_static_constraint( v_cnst[ 0 ] );

 AR |= HasCns;
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_objective( Configuration * objc ){
 
 if( AR & HasObj )  // the objective is there already
  return;           // nothing to do 

 if( ! get_NItems() )
  return;

 LinearFunction::v_coeff_pair p( get_NItems() );
  for( Index i = 0 ; i < get_NItems() ; i++ ){
   p[ i ].first = &v_x[ i ];
   p[ i ].second = v_P[ i ];
  }

 LinearFunction * obj = new LinearFunction( std::move( p ) , 0  );
 f_obj.set_function( obj , eNoMod );
 f_obj.set_sense( f_sense , eNoMod );
 set_objective( & f_obj , eNoMod );
 
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
   
   for( auto & xi : v_x )
     if( ! xi.is_feasible() )
       return( false );

  return( v_cnst[ 0 ].rel_viol() > 0 ? false : true );
}

// do it using the physical representation
 double tot_weight = 0; 
 for( Index i = 0 ; i < v_W.size() ; i++ )
  tot_weight += v_W[ i ] * get_x( i ); 

 return( tot_weight <= f_C ? true : false );    
}

/*--------------------------------------------------------------------------*/

bool BinaryKnapsackBlock::is_empty( bool useabstract , Configuration * optc ){
 
 // check if there are fixed variables and compute the residual capacity
 
 double C = f_C;

 for( Index i = 0 ; i < f_N ; i++ ){
  if( v_x[ i ].is_fixed()  && v_x[ i ].get_value() )
   C -= v_W[ i ]; 
 }

 if( C >= 0 )                
  return( false );    
 

 double neg_weights = 0;

 for( Index i = 0 ; i < f_N ; i++ ){
  
  if( is_fixed( i ) )
   continue; 

  if( v_W[ i ] < 0 ){
   neg_weights += v_W[ i ];
   
   if( neg_weights <= C )
    return( false );
   
  }
 }

 return( true );
}

/* -------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/

 Block * BinaryKnapsackBlock::get_R3_Block( Configuration *r3bc , 
                                            Block * base , Block * father ){
 
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 BinaryKnapsackBlock * BKB;
 if( base ){
  BKB = dynamic_cast< BinaryKnapsackBlock * >( base );
  if( ! BKB )
   throw( std::invalid_argument( "base is not a BinaryKnapsackBlock" ) );
 } else
    BKB = new BinaryKnapsackBlock( father );

  BKB->load( f_N , f_C , v_W , v_P );
  BKB->set_objective_sense( f_sense );

  return( BKB );

 }// end( BinaryKnapsackBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::map_back_solution( Block *R3B , Configuration *r3bc, 
                        Configuration *solc ){

BinaryKnapsackBlock * BKB = dynamic_cast< BinaryKnapsackBlock * >( R3B );
if( ! BKB )
 throw( std::invalid_argument( "R3B is not a BinaryKnapsackBlock" ) );

// check if the size of the variables are equal. Better not to use NItems
// to allow different formulations.
if( BKB->get_VarSize() != get_VarSize() )
 throw( std::invalid_argument( "incompatible variables size" ) );

// copy solution
std::vector< double > xSol( get_VarSize() );
BKB->get_x( xSol );
set_x( xSol );
 
}// end( BinaryKnapsackBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::map_forward_solution( Block * R3B , 
                              Configuration *r3bc , Configuration *solc ){

BinaryKnapsackBlock * BKB = dynamic_cast< BinaryKnapsackBlock * >( R3B );
if( ! BKB )
 throw( std::invalid_argument( "R3B is not a BinaryKnapsackBlock" ) );

// check if the size of the variables are equal. Better not to use NItems
// to allow different formulations.
if( BKB->get_VarSize() != get_VarSize() )
 throw( std::invalid_argument( "incompatible variables size" ) );

// copy solution
std::vector< double > xSol( get_VarSize() );
get_x( xSol );
BKB->set_x( xSol );
 
}// end( BinaryKnapsackBlock::map_forward_solution )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * BinaryKnapsackBlock::get_Solution( Configuration *solc , 
                                              bool emptys ){
 auto * sol = new BinaryKnapsackSolution(); 

 sol->v_x.resize( get_VarSize() );

 if( ! emptys )
  sol->read( this );

 return( sol ); 
}

/*--------------------------------------------------------------------------*/

double BinaryKnapsackBlock::get_x( Index i ){
 if( i >= get_VarSize() )
  throw( std::invalid_argument( "invalid item" ) );
 return ( v_x[ i ].get_value());// < BinaryTol ? false : true ); 
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::get_x( doubleVec & xSol , Range rng ){
 
 rng.second = std::min( rng.second , get_VarSize() );

 auto xSoli = xSol.begin();
 for( Index i = rng.first ; i < rng.second ; i++ )
  ( * xSoli++ ) = v_x[ i ].get_value();// < BinaryTol ? false : true;

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::get_x( doubleVec & xSol , c_Subset & nms ){
 
 auto xSoli = xSol.begin();
 for( auto i : nms ){
  if( i >= get_VarSize() )
   throw( std::invalid_argument( "invalid item" ) );
  ( * xSoli++ ) = v_x[ i ].get_value();// < BinaryTol ? false : true;
 }
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( Index i , double value ){
 if( i >= get_VarSize() )
  throw( std::invalid_argument( "invalid item" ) );
 v_x[ i ].set_value( value ); 
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( c_doubleVec & xSol , Range rng ){
 
 rng.second = std::min( rng.second , get_NItems() );

 auto xSoli = xSol.begin();
 for( Index i = rng.first ; i < rng.second ; i++ )
  v_x[ i ].set_value( * xSoli++ );

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( c_doubleVec & xSol , c_Subset & nms ){
 
 auto xSoli = xSol.begin();
 for( auto i : nms ){
  if( i >= get_VarSize() )
   throw( std::invalid_argument( "invalid item" ) );
  v_x[ i ].set_value( * xSoli++ );
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
 
 ( group.addVar( "Weights" , netCDF::NcDouble() , ni ) ).putVar( v_W.data() );
 
 ( group.addVar( "Profits" , netCDF::NcDouble() , ni ) ).putVar( v_P.data() );
 
 if( ! v_I.empty() ){
   std::vector< int > tempI;
   for(int i=0; i<v_I.size(); i++){
     tempI[i]=(int) v_I[i]; //.data();
   }
   ( group.addVar( "Integrality", netCDF::NcInt(), ni ) ).putVar( tempI.data() );
 }
}// end( BinaryKnapsackBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( bool value, Index i , 
                                 ModParam issueMod, ModParam issueAMod ){

 if( i >= get_VarSize() )
  throw( std::invalid_argument( "invalid index item" ) );

 // if i is already fixed return
 if( v_x[ i ].is_fixed() )
  return;      

 // reset conditional bounds
 f_cond_lower = -Inf<double>();
 f_cond_upper = Inf<double>();


if( not_dry_run( issueAMod ) ){
 v_x[ i ].set_value( value ); 
 v_x[ i ].is_fixed( true , issueAMod ); 
}
 
 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                BinaryKnapsackBlockMod::eFixX , std::make_pair( i , i + 1 ) ) , 
                Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( c_boolVec & value , Range rng , 
                                 ModParam issueMod, ModParam issueAMod ){

 rng.second = std::min( rng.second , get_VarSize() );
 if( rng.second <= rng.first )  // nothing to do
  return;

 auto vi = value.begin();

 for( Index i = rng.first ; i < rng.second ; i++ , vi++ )
  fix_x( * vi , i , eDryRun , issueAMod );

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                           BinaryKnapsackBlockMod::eFixX , rng ) , 
                           Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( c_boolVec & value , Subset && nms , 
                                 ModParam issueMod, ModParam issueAMod ){

 if( nms.empty() )
  return;

 auto vi = value.begin();

 for( auto i : nms )
  fix_x( * vi++ , i , eDryRun , issueAMod );

 std::sort( nms.begin() , nms.end() );

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(this,
                           BinaryKnapsackBlockMod::eFixX , std::move( nms ) ) , 
                           Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( Index i , ModParam issueMod, 
                                   ModParam issueAMod ){

 if( i >= get_VarSize() )
  throw( std::invalid_argument( "invalid index item" ) );

 if( !( v_x[ i ].is_fixed() ) ) // already unfixed
  return;                        // nothing to do

 // reset conditional bounds
 f_cond_lower = -Inf<double>();
 f_cond_upper = Inf<double>();
 
 if( not_dry_run( issueAMod ) )
  v_x[ i ].is_fixed( false , issueAMod ); 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
            BinaryKnapsackBlockMod::eUnfixX , std::make_pair( i , i + 1 ) ) , 
            Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( Range rng , ModParam issueMod, 
                   ModParam issueAMod ){


 rng.second = std::min( rng.second , get_VarSize() );
 if( rng.second <= rng.first )  // nothing to do
  return;

 for( Index i = rng.first ; i < rng.second ; i++ )
  unfix_x( i , eDryRun , issueAMod );

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                           BinaryKnapsackBlockMod::eUnfixX , rng ) , 
                           Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( Subset && nms , ModParam issueMod, 
                   ModParam issueAMod ){

 if( nms.empty() )
  return;

 for( auto i : nms )
  unfix_x( i , eDryRun , issueAMod );

 std::sort( nms.begin() , nms.end() );

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(this,
                        BinaryKnapsackBlockMod::eUnfixX , std::move( nms ) ) , 
                        Observer::par2chnl( issueMod ) );

} // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weight( double NWeight , Index item , 
                          ModParam issueMod , ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );

 if( v_W[ item ] == NWeight )  // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  LinearFunction *lf = dynamic_cast<LinearFunction*>(v_cnst[0].get_function());
  lf->modify_coefficient( item , NWeight , issueAMod ); // AR
  v_W[ item ] = NWeight;                                // PR
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation
  v_W[ item ] = NWeight; 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
      BinaryKnapsackBlockMod::eChgWeight , std::make_pair( item , item + 1 ) ), 
      Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::chg_weight )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weights( const dblVec_it NWeight,
                                       Range rng , 
                                       ModParam issueMod ,
                                       ModParam issueAMod ){

 rng.second = std::min( rng.second , get_NItems() );
 if( rng.second <= rng.first )  // nothing to change
  return;   

 if( std::equal( NWeight , NWeight + ( rng.second - rng.first ) ,
   v_W.begin() + rng.first ) )
  return;  // nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  
  // physical representation
  std::copy( NWeight , NWeight + ( rng.second - rng.first ) ,
             v_W.begin() + rng.first );
  
  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>(v_cnst[0].get_function());
  lf->modify_coefficients( std::vector<double>( NWeight , 
                  NWeight + ( rng.second - rng.first ) ) , rng , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  std::copy( NWeight , NWeight + ( rng.second - rng.first ) ,
             v_W.begin() + rng.first );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                                    BinaryKnapsackBlockMod::eChgWeight , rng ), 
                                    Observer::par2chnl( issueMod ) );

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weights( const dblVec_it NWeight,
                                       Subset && nms , bool ordered ,  
                                       ModParam issueMod ,
                                       ModParam issueAMod ){

 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( v_W , nms , NWeight , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  
  // physical representation
  copyidx( v_W , nms , NWeight );
  
  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>(v_cnst[0].get_function());
  lf->modify_coefficients( std::vector<double>( NWeight , 
                  NWeight + nms.size() ) , 
                  Subset( nms ) , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  copyidx( v_W , nms , NWeight );
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

void BinaryKnapsackBlock::chg_profit( double NProfit , Index item , 
                          ModParam issueMod , ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );

 if( v_P[ item ] == NProfit ) // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  v_P[ item ] = NProfit; 

  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>( f_obj.get_function() );
  lf->modify_coefficient( item , NProfit , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation 
  v_P[ item ] = NProfit;
 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
      BinaryKnapsackBlockMod::eChgProfit , std::make_pair( item , item + 1 ) ), 
      Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::chg_profit )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profits( const dblVec_it NProfit ,
                                       Range rng , 
                                       ModParam issueMod ,
                                       ModParam issueAMod ){

 rng.second = std::min( rng.second , get_NItems() );
 if( rng.second <= rng.first )  // nothing to change
  return;   

 if( std::equal( NProfit , NProfit + ( rng.second - rng.first ) ,
             v_P.begin() + rng.first ) )
  return;  // nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  std::copy( NProfit , NProfit + ( rng.second - rng.first ) ,
             v_P.begin() + rng.first );
  
  // abstract representation  
  LinearFunction *lf = dynamic_cast<LinearFunction*>( f_obj.get_function() );
  lf->modify_coefficients( std::vector<double>( NProfit , 
                    NProfit + ( rng.second - rng.first ) ) , rng , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){ // otherwise only physical representation 
  std::copy( NProfit , NProfit + ( rng.second - rng.first ) ,
             v_P.begin() + rng.first );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                                    BinaryKnapsackBlockMod::eChgProfit , rng ), 
                                    Observer::par2chnl( issueMod ) );

} 

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profits( const dblVec_it NProfit,
                                       Subset && nms , bool ordered ,  
                                       ModParam issueMod ,
                                       ModParam issueAMod ){

 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( v_P , nms , NProfit , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  copyidx( v_P , nms , NProfit );
  
  // abstract representation
  LinearFunction *lf = dynamic_cast<LinearFunction*>( f_obj.get_function() );
  lf->modify_coefficients( std::vector<double>( NProfit , 
                           NProfit + nms.size() ) , 
                           Subset( nms ) , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  copyidx( v_P , nms , NProfit );
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

void BinaryKnapsackBlock::chg_integrality( bool NIntegrality , Index item , 
                          ModParam issueMod , ModParam issueAMod ){

 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );
  
 if(NIntegrality == false && v_I.empty()){
   v_I.resize( get_NItems() , true );
   countCont=0;
 }
 
 if( v_I[ item ] == NIntegrality ) // nothing to do
  return;


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  v_I[ item ] = NIntegrality; 

  // abstract representation
  if(NIntegrality==true){
    if(v_x[item].get_type()==ColVariable::kPosUnitary){
      countCont--;
      v_x[item].set_type( ColVariable::kBinary , eNoBlck );
    }
  }else{
    if(v_x[item].get_type()==ColVariable::kBinary){
      countCont++;
      v_x[item].set_type( ColVariable::kPosUnitary , eNoBlck );
      }
    }
  } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation 
  v_I[ item ] = NIntegrality;
 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
      BinaryKnapsackBlockMod::eChgIntegrality , std::make_pair( item , item + 1 ) ), 
      Observer::par2chnl( issueMod ) );

   if( ! countCont )
         v_I.clear();

 } // end( BinaryKnapsackBlock::chg_integrality )

/*--------------------------------------------------------------------------*/


void BinaryKnapsackBlock::chg_integrality( const boolVec_it NIntegrality ,
                                       Range rng , 
                                       ModParam issueMod ,
                                       ModParam issueAMod ){

 rng.second = std::min( rng.second , get_NItems() );
 
 if(v_I.empty()){
   v_I.resize( get_NItems() , true );
   countCont=0;
 }
 if( rng.second <= rng.first )  // nothing to change
  return;   

 if( std::equal(  NIntegrality , NIntegrality + ( rng.second - rng.first ) ,
             v_I.begin() + rng.first ) )
  return;  // nothing changes, avoid issuing the Modification

 
 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  std::copy( NIntegrality , NIntegrality + ( rng.second - rng.first ) ,
             v_I.begin() + rng.first );
  
  // abstract representation
  //ColVaruable *vx = dynamic_cast<ColVariable*>( v_x );
     for(int i=0;i<=rng.second-rng.first;i++){
        if(NIntegrality[i]==true){
          if(v_x[rng.first+i].get_type()==ColVariable::kPosUnitary){
            countCont--;
            v_x[rng.first+i].set_type( ColVariable::kBinary , eNoBlck );
          }
        }else{
          if(v_x[rng.first+i].get_type()==ColVariable::kBinary){
            countCont++;
            v_x[rng.first+i].set_type( ColVariable::kPosUnitary , eNoBlck );
          }
       }
     }
 } 
 else if( not_dry_run( issueMod ) ){ // otherwise only physical representation 
  std::copy( NIntegrality , NIntegrality + ( rng.second - rng.first ) ,
             v_I.begin() + rng.first );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(this,
                                    BinaryKnapsackBlockMod::eChgIntegrality , rng ), 
                                    Observer::par2chnl( issueMod ) );
                                    
                                     if( ! countCont )
         v_I.clear();

} 

/*--------------------------------------------------------------------------*/


void BinaryKnapsackBlock::chg_integrality( const boolVec_it NIntegrality,
                                       Subset && nms , bool ordered , 
                                       ModParam issueMod ,
                                       ModParam issueAMod ){

 if( nms.empty() )  // nothing to change
  return;            
 
 if( v_I.empty()){
   v_I.resize( get_NItems() , true );
   countCont=0;
 }
 
 if( is_equal( v_I , nms , NIntegrality , get_NItems() ) )
   return;  // actually nothing changes, avoid issuing the Modification


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  
  // physical representation
  copyidx( v_I , nms , NIntegrality );
  
 // abstract representation
  for(int i=0;i<=nms.size();i++){
    if(NIntegrality[i]==true){
      if(v_x[nms[i]].get_type()==ColVariable::kPosUnitary){
      countCont--;
      v_x[nms[i]].set_type( ColVariable::kBinary , eNoBlck );
      }
    }else{
    if(v_x[nms[i]].get_type()==ColVariable::kBinary){
      countCont++;
      v_x[nms[i]].set_type( ColVariable::kPosUnitary , eNoBlck );
      }
    } 
  }
 } 
 else if( not_dry_run( issueMod ) ){
  // otherwise change only physical representation 
  copyidx( v_I , nms , NIntegrality );
 }

 // issue physical Modification 
 if( issue_pmod( issueMod ) ){ 
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );


  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(this,
                      BinaryKnapsackBlockMod::eChgIntegrality , std::move( nms ) ), 
                      Observer::par2chnl( issueMod ) );
 }
 
  if( ! countCont )
         v_I.clear();
} 
/*-----------------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_capacity( double NC , 
                          ModParam issueMod , ModParam issueAMod ){

 if( f_C == NC )    // nothing to do
  return;


 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ){
  f_C = NC; 
  v_cnst[0].set_rhs( NC , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation 
  f_C = NC;
 

 if( issue_pmod( issueMod ) ) // issue physical Modification
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockMod >( this ,
                           BinaryKnapsackBlockMod::eChgCapacity ) , 
                           Observer::par2chnl( issueMod ) );


 } // end( BinaryKnapsackBlock::chg_capacity )

 /*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_objective_sense( bool sense , ModParam issueMod , 
                                               ModParam issueAMod ){

 if( f_sense == sense ) // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();


 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ){
  f_sense = sense; 
  f_obj.set_sense( sense , issueAMod );
 } 
 else if( not_dry_run( issueMod ) ) // otherwise only physical representation 
  f_sense = sense;
 

 if( issue_pmod( issueMod ) ) // issue physical Modification 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockMod >( this ,
                           BinaryKnapsackBlockMod::eChgSense ) , 
                           Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::set_objective_sense )


/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::print( std::ostream & output ) const {
 
 output << "BinaryKnapsackBlock\n";
 output << "Number of items: " << f_N << std::endl; 
 output << "Capacity: " << f_C << std::endl;
 
 output << "\tWeights\tProfits\n";
 for( Index i = 0 ; i < f_N ; i++ )
  output << "Item " << i << "\t" << v_W[i] << "\t" << v_P[i] << std::endl;

}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( std::istream & input ){ 

// erase previous instance, if any

 if( get_NItems() )
  guts_of_destructor();

// read problem data
if( !( input >> f_N ) )
 throw( std::invalid_argument( "error reading number of items" ) );

if( !( input >> f_C ) )
 throw( std::invalid_argument( "error reading Capacity" ) );

v_W.resize( f_N );
v_P.resize( f_N );
v_I.resize( f_N );

for( Index i = 0 ; i < f_N ; i++ ){
 if( !( input >> v_W[i] ) )
  throw( std::invalid_argument( "error reading Weights" ) );
}

for( Index i = 0 ; i < f_N ; i++ ){
 if( !( input >> v_P[i] ) )
  throw( std::invalid_argument( "error reading Profits" ) );
}

for( Index i = 0 ; i < f_N ; i++ ){
 if( !( (bool) input >> v_I[i] ) )
  throw( std::invalid_argument( "error reading Integrality Constraints" ) );
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
 for( auto & c : v_cnst )
  c.clear();
 v_cnst.clear();
 // clear the objective function
 f_obj.clear();
 // clear all variables
 v_x.clear();

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

//change only this for modification?? add Integrality modification
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
   auto lfo = dynamic_cast< LinearFunction * >( f_obj.get_function() ); 
   // and check if the modification is on the objective  
    if( lf == lfo ){ 
     // vector of new profits 
     std::vector<double> new_profits( range_size );  
     dblVec_it npi = new_profits.begin(); 
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
   auto lfc = dynamic_cast< LinearFunction * >( v_cnst[0].get_function() ); 
   // and check if the modification is on the constraint  
    if( lf == lfc ){ 
     // vector of new weights 
     std::vector<double> new_weights( range_size ); 
     dblVec_it nwi = new_weights.begin();
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
   auto lfo = dynamic_cast< LinearFunction * >( f_obj.get_function() ); 
   // and check if the modification is on the objective  
    if( lf == lfo ){
     // vector of new profits 
     std::vector<double> new_profits( tmod->subset().size() );
     dblVec_it npi = new_profits.begin();   
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
   auto lfc = dynamic_cast< LinearFunction * >( v_cnst[0].get_function() ); 
   // and check if the modification is on the constraint  
    if( lf == lfc ){
     // vector of new weights 
     std::vector<double> new_weights( tmod->subset().size() );  
     dblVec_it nwi = new_weights.begin();
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

// ObjectiveMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
{
 const auto tmod = dynamic_cast< ObjectiveMod * >( mod );
 if( tmod ){ 
  if( !( AR & HasObj ) ) // check if the objective exists
   throw(std::invalid_argument("Modification to non-constructed objective"));

   auto obj = dynamic_cast< Objective * >( & f_obj );
   if( tmod->of() == obj ){ 
    bool sense = tmod->type() == Objective::eMax ? 1 : 0;
    set_objective_sense( sense , make_par( eNoBlck , chnl ) , eDryRun );
    return;
   }
   throw( std::invalid_argument("Modification to the wrong objective") );
  }
}

throw(std::invalid_argument( "unsupported Modification to BinaryKnapsackBlock"));
}

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::compute_conditional_bounds(){

f_cond_lower = 0;
f_cond_upper = 0;

for( Index i = 0 ; i < f_N ; i++ ){
 
 double w = v_W[ i ];                           // weight of the current item
 double p = f_sense ? v_P[ i ] : -v_P[ i ];     // profit of the current item
 

 // items contained in the optimal solution
 if( w <= 0 && p >= 0 ){
  f_cond_lower += p;                
  f_cond_upper += p;  
  continue;
 }

 // items not contained in the optimal solution
 if( w >= 0 && p <= 0 )
  continue;

 // remaining items
 if( p > 0 )
  f_cond_upper += p;
 else 
  f_cond_lower += p;

}

// if it is a minimization problem
if( ! f_sense ){
 f_cond_lower = -f_cond_lower;
 f_cond_upper = -f_cond_upper;
 std::swap( f_cond_lower , f_cond_upper );  
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

 for( auto & xi : BKB->v_x ) 
  *( vxi++ ) = static_cast< double >( xi.get_value() ); 

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

  for( auto & xi : BKB->v_x )
   xi.set_value( static_cast< bool >( *( vxi++ ) ) );
 }
} // end( BinaryKnapsackSolution::write )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::serialize( netCDF::NcGroup & group ) const {
 if( ! v_x.empty() ) {
  netCDF::NcDim ni = group.addDim( "n" , v_x.size() ); 
( group.addVar( "x" , netCDF::NcDouble() , ni ) ).putVar( v_x.data() );
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
