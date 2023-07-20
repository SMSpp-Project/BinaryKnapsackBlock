/*--------------------------------------------------------------------------*/
/*---------------------- File BinaryKnapsackBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class BinaryKnapsackKBlock, which
 * implements the Block concept for a Knapsack problem with linear objective
 * and both integer and continuous variables.
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
 * \copyright &copy; by Federica Di Pasquale, Antonio Frangioni,
 *                      Francesco Demelas
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
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
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector< T > & vec , c_Subset & nms ,
                      typename std::vector< T >::const_iterator cmp ,
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
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
                     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/

static LinearFunction * LF( Function * f )
{
 return( static_cast< LinearFunction * >( f ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BinaryKnapsackBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( BinaryKnapsackBlock );

// register BinaryKnapsackSolution to the Solution factory
SMSpp_insert_in_factory_cpp_1( BinaryKnapsackSolution );

// register BinaryKnapsackChange to the Change factory
SMSpp_insert_in_factory_cpp_1( BinaryKnapsackBlockChange );

// register BinaryKnapsackRngdChange to the Change factory
SMSpp_insert_in_factory_cpp_1( BinaryKnapsackBlockRngdChange );

// register BinaryKnapsackChange to the Change factory
SMSpp_insert_in_factory_cpp_1( BinaryKnapsackBlockSbstChange );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF BinaryKnapsackBlock ----------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( Index n , double Capacity , 
                                const std::vector< double > & Weights , 
                                const std::vector< double > & Profits,
                                const std::vector< bool > & Integrality )
{
 // sanity checks 

 if( Weights.size() != n )
  throw( std::invalid_argument( "Vector of Weights of the wrong size" ) );

 if( Profits.size() != n )
  throw( std::invalid_argument( "Vector of Profits of the wrong size" ) );
 
 if( ( ! Integrality.empty() ) && ( Integrality.size() != n ) )
  throw( std::invalid_argument( "Vector of Integrality of the wrong size" ) );
 
 // call load( , , && , && , && ) on newly constructed copies
 load( n , Capacity , std::vector< double >( Weights ) ,
       std::vector< double >( Profits ) ,
       std::vector< bool >( Integrality ) );

 } // end( BinaryKnapsackBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( Index n , double Capacity , 
                                std::vector< double > && Weights , 
                                std::vector< double > && Profits ,
                                std::vector< bool > && Integrality )
{
 // sanity checks 

 if( Weights.size() != n )
  throw( std::invalid_argument( "Vector of Weights of the wrong size" ) );

 if( Profits.size() != n )
  throw( std::invalid_argument( "Vector of Profits of the wrong size" ) );

 if( ( ! Integrality.empty() ) && ( Integrality.size() != n ) )
   throw( std::invalid_argument( "Vector of Integrality of the wrong size" ) );

 // erase previous instance, if any
 if( get_NItems() )
  guts_of_destructor();

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 f_C = Capacity;
 v_W = std::move( Weights );  
 v_P = std::move( Profits );

 if( Integrality.empty() ) 
  v_I.assign( n , true );
 else 
  v_I = std::move( Integrality );

 countCont = std::count( v_I.begin() , v_I.end() , false );

 v_fxd.assign( n , 0 ); // all the variables are not fixed
 
 generate_abstract_variables();

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 // Modification
 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( BinaryKnapsackBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any
 if( get_NItems() )
  guts_of_destructor();

 // read problem data
 Index n;
 if( ! ( input >> eatcomments >> n ) )
  throw( std::invalid_argument( "error reading number of items" ) );

 if( ! ( input >> eatcomments >> f_C ) )
  throw( std::invalid_argument( "error reading Capacity" ) );

 v_W.resize( n );
 v_P.resize( n );
 v_I.assign( n , true );
 v_fxd.assign( n , 0 ); // all the variables are not fixed         

 for( Index i = 0 ; i < get_NItems() ; ++i )
  if( ! ( input >> eatcomments >> v_W[ i ] ) )
   throw( std::invalid_argument( "error reading Weights" ) );

 for( Index i = 0 ; i < get_NItems() ; ++i )
  if( ! ( input >> eatcomments >> v_P[ i ] ) )
   throw( std::invalid_argument( "error reading Profits" ) );

 input >> eatcomments;
 if( ! input.eof() )
  for( Index i = 0 ; i < get_NItems() ; ++i ) {
   input >> eatcomments;
   if( ! ( ( bool ) input >> v_I[ i ] ) )
    throw( std::invalid_argument( "error reading Integrality Constraints" ) );
   }

 generate_abstract_variables();

 // Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( BinaryKnapsackBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::deserialize( const netCDF::NcGroup & group )
{
 // erase previous instance, if any
 if( get_NItems() )
  guts_of_destructor();

 // read problem data 
 netCDF::NcDim ni = group.getDim( "NItems" );
 if( ni.isNull() )
  throw( std::logic_error( "NItems dimension is required" ) );
 Index n = ni.getSize();

 netCDF::NcDim c = group.getDim( "Capacity" );
 if( c.isNull() )
  throw( std::logic_error( "Capacity is required" ) );
 f_C = c.getSize();

 netCDF::NcVar w = group.getVar( "Weights" );
 if( w.isNull() )
  throw( std::logic_error( "Weights are required" ) );
 
 v_W.resize( n );
 w.getVar( v_W.data() );

 netCDF::NcVar p = group.getVar( "Profits" );
 if( p.isNull() )
  throw( std::logic_error( "Profits are required" ) );
 
 v_P.resize( n );
 p.getVar( v_P.data() );
 
 netCDF::NcVar i = group.getVar( "Integrality" );

 v_I.resize( n );
 if( i.isNull() )
  std::fill( v_I.begin() , v_I.end() , true );
 else {
  std::vector< unsigned char > tmpI( v_I.size() );
  i.getVar( tmpI.data() );
  for( Index i = 0 ; i < v_I.size() ; ++i )
   v_I[ i ] = ( tmpI[ i ] != 0 );
  }
 
 countCont = std::count( v_I.begin() , v_I.end() , false );

 v_fxd.assign( n , 0 ); // all the variables are not fixed  
 
 generate_abstract_variables();

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 // call the method of Block
 // inside this the NBModification, the "nuclear option",  is issued
 Block::deserialize( group );

 }  // end( BinaryKnapsackBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_abstract_variables( Configuration * stvv )
{
 if( AR & HasVar )  // the variables are there already
  return;           // nothing to do
 
 if( ! get_NItems() )
  return;

 v_x.resize( get_NItems() );
 
 for( Index i = 0 ; i < get_NItems() ; ++i ) {
  if( v_I[ i ] )
   v_x[ i ].set_type( ColVariable::kBinary , eNoBlck );
  else
   v_x[ i ].set_type( ColVariable::kPosUnitary , eNoBlck ); 
  }

 add_static_variable( v_x );

 AR |= HasVar;

 } // end( BinaryKnapsackBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_abstract_constraints(Configuration * stcc)
{
 if( AR & HasCns )  // the constraint is there already
  return;           // nothing to do

 if( ! get_NItems() )
  return;

 LinearFunction::v_coeff_pair w( get_NItems() );
 for( Index i = 0 ; i < get_NItems() ; i++ ) {
  w[ i ].first = &v_x[ i ];
  w[ i ].second = v_W[ i ];
  }

 f_cnst.set_function( new LinearFunction( std::move( w ) , 0 ) , eNoBlck );
 f_cnst.set_rhs( f_C );
 f_cnst.set_lhs( -Inf< double >() );
 
 add_static_constraint( f_cnst );

  AR |= HasCns;

 } // end( BinaryKnapsackBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::generate_objective( Configuration * objc )
{
 if( AR & HasObj )  // the objective is there already
  return;           // nothing to do 

 if( ! get_NItems() )
  return;

 LinearFunction::v_coeff_pair p( get_NItems() );
 for( Index i = 0 ; i < get_NItems() ; i++ ) {
  p[ i ].first = &v_x[ i ];
  p[ i ].second = v_P[ i ];
  }

 LinearFunction * obj = new LinearFunction( std::move( p ) , 0  );
 f_obj.set_function( obj , eNoMod );
 f_obj.set_sense( f_sense ? Objective::eMax : Objective::eMin , eNoMod );
 set_objective( & f_obj , eNoMod );
 
 AR |= HasObj;

 } // end( BinaryKnapsackBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool BinaryKnapsackBlock::is_feasible( bool useabstract , 
                                       Configuration * fsbc )
{
 // first check bound constraints, for which only one way is available
 for( auto & xi : v_x )
  if( ! xi.is_feasible() )
   return( false );

 // now check the knapsack constraint(s)

 if( useabstract && ( AR & HasCns ) ) {
  // do it using the abstract representation, if there is any

  // ensure the knapsack constraint(s) is compute()-d, as this is
  // necessary for abs_viol() to work
  if( auto ret = f_cnst.compute() ;
      ( ret <= FRowConstraint::kUnEval ) || ( ret > FRowConstraint::kOK ) )
   return( false );  // failure to compute() means no proof it is feasible

  // TODO: some epsilon might be required
  return( f_cnst.abs_viol() == 0 );
  }

 // do it using the physical representation
 double tot_weight = 0; 
 for( Index i = 0 ; i < v_W.size() ; ++i )    // check fixed variables?
  tot_weight += v_W[ i ] * get_x( i ); 

 return( tot_weight <= f_C );

 } // end( BinaryKnapsackBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool BinaryKnapsackBlock::is_empty( bool useabstract , Configuration * optc )
{
 // check if there are fixed variables and compute the residual capacity
 
 double C = f_C;
 for( Index i = 0 ; i < get_NItems() ; ++i )
  if( v_x[ i ].is_fixed()  && v_x[ i ].get_value() )
   C -= v_W[ i ]; 

 if( C >= 0 )
  return( false );    
 
 double neg_weights = 0;

 for( Index i = 0 ; i < get_NItems() ; ++i ) {
  if( is_fixed( i ) )
   continue; 

  if( v_W[ i ] < 0 ) {
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

Block * BinaryKnapsackBlock::get_R3_Block( Configuration * r3bc ,
             Block * base , Block * father )
{
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 BinaryKnapsackBlock * BKB;
 if( base ) {
  BKB = dynamic_cast< BinaryKnapsackBlock * >( base );
  if( ! BKB )
   throw( std::invalid_argument( "base is not a BinaryKnapsackBlock" ) );
  }
 else
  BKB = new BinaryKnapsackBlock( father );

 BKB->load( get_NItems() , f_C , v_W , v_P );
 BKB->set_objective_sense( f_sense );

 return( BKB );

 }  // end( BinaryKnapsackBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::map_back_solution( Block * R3B , 
                                             Configuration * r3bc , 
                                             Configuration * solc )
{
 BinaryKnapsackBlock * BKB = dynamic_cast< BinaryKnapsackBlock * >( R3B );
 if( ! BKB )
  throw( std::invalid_argument( "R3B is not a BinaryKnapsackBlock" ) );

 // check if the size of the variables are equal. Better not to use NItems
 // to allow different formulations.
 if( BKB->v_x.size() != v_x.size() )
  throw( std::invalid_argument( "incompatible variables size" ) );

 // copy solution
 doubleVec xSol( v_x.size() );
 BKB->get_x( xSol.begin() );
 set_x( xSol.begin() );
 
 }  // end( BinaryKnapsackBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::map_forward_solution( Block * R3B , 
                                                Configuration * r3bc , 
                                                Configuration * solc )
{
 // fantastically dirty trick
 BinaryKnapsackBlock * BKB = dynamic_cast< BinaryKnapsackBlock * >( R3B );
 if( ! BKB )
  throw( std::invalid_argument( "R3B is not a BinaryKnapsackBlock" ) );

 BKB->map_back_solution( this , r3bc , solc );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * BinaryKnapsackBlock::get_Solution( Configuration * solc , 
                                              bool emptys )
{
 auto * sol = new BinaryKnapsackSolution(); 

 sol->v_x.resize( v_x.size() );

 if( ! emptys )
  sol->read( this );

 return( sol ); 
 }

/*--------------------------------------------------------------------------*/

double BinaryKnapsackBlock::get_x( Index i ) const
{
 if( i >= v_x.size() )
  throw( std::invalid_argument( "BinaryKnapsackBlock::get_x: invalid item"
        ) );
 return( v_x[ i ].get_value());
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::get_x( dblVec_it xSol , Range rng ) const
{
 rng.second = std::min( rng.second , Index( v_x.size() ) );

 for( Index i = rng.first ; i < rng.second ; ++i )
  *(xSol++) = v_x[ i ].get_value();
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::get_x( dblVec_it xSol , c_Subset & nms ) const
{ 
 for( auto i : nms ) {
  if( i >= v_x.size() )
   throw( std::invalid_argument( "BinaryKnapsackBlock::get_x: invalid item"
         ) );
  *(xSol++) = v_x[ i ].get_value();
  }
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( Index i , double value )
{
 if( i >= v_x.size() )
  throw( std::invalid_argument( "BinaryKnapsackBlock::set_x: invalid item"
        ) );
 v_x[ i ].set_value( value ); 
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( c_dblVec_it xSol , Range rng )
{ 
 rng.second = std::min( rng.second , get_NItems() );

 for( Index i = rng.first ; i < rng.second ; i++ )
  v_x[ i ].set_value( *(xSol++) );
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_x( c_dblVec_it xSol , c_Subset & nms )
{ 
 for( auto i : nms ) {
  if( i >= v_x.size() )
   throw( std::invalid_argument( "BinaryKnapsackBlock::set_x: invalid item"
         ) );
  v_x[ i ].set_value( *(xSol++) );
  }
 }

/*--------------------------------------------------------------------------*/

 bool BinaryKnapsackBlock::is_fixed( Index i ) const { 
  
  if( i >= v_fxd.size() )
    throw( std::invalid_argument( "Invalid index of the item" ) );
  
  if( v_fxd[ i ] == 0 )
   return( false );

  return( true );
 }


/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
 }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*------ METHODS FOR LOADING, PRINTING & SAVING THE BinaryKnapsackBlock ----*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::print( std::ostream & output , char vlvl ) const
{ 
 output << "BinaryKnapsackBlock" << std::endl;
 output << "Number of items: " << get_NItems() << std::endl; 
 output << "Capacity: " << f_C << std::endl;
 output << "Sense of the objective: ";
 if( f_sense )
  output << "Max" << std::endl;
 else
  output << "Min" << std::endl;

 if( vlvl != 'C' )
  return;

 output << "\tWeights\tProfits\tBinary\tFixed\tFixed to\n";
 for( Index i = 0 ; i < get_NItems() ; ++i ) {
  output << "Item " << i << "\t" << v_W[ i ] << "\t" << v_P[ i ] << "\t";
  output << v_I[ i ] << "\t";
  if( v_fxd[ i ] == 0 )
   output << "0" << "\t" << "-" << std::endl;
  else {
   if( v_fxd[ i ] == 1 )
    output << "1" << "\t" << "0" << std::endl;
   else
    output << "1" << "\t" << "1" << std::endl;  
  } 
 }

 } // end( BinaryKnapsackBlock::print )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::serialize( netCDF::NcGroup & group ) const
{ 
 // call the method of Block
 Block::serialize( group );

 // BinaryKnapsackBlock data 
 netCDF::NcDim ni = group.addDim( "NItems" , get_NItems() );
 
 group.addDim( "Capacity" , get_Capacity() );
 
 ( group.addVar( "Weights" , netCDF::NcDouble() , ni ) ).putVar( v_W.data() );
 
 ( group.addVar( "Profits" , netCDF::NcDouble() , ni ) ).putVar( v_P.data() );
 
 if( countCont ) {
   std::vector< int > tempI;
   for( Index i = 0 ; i < v_I.size() ; ++i )
    tempI[ i ] = ( int ) v_I[ i ]; 

   ( group.addVar( "Integrality" ,
       netCDF::NcInt(), ni ) ).putVar( tempI.data() );
  }
 }  // end( BinaryKnapsackBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( bool value , Index i , 
                                 ModParam issueMod , ModParam issueAMod )
{
 if( i >= v_x.size() )
  throw( std::invalid_argument( "BinaryKnapsackBlock::fix_x: invalid item"
        ) );

 // check if it is already fixed
 if( ( v_x[ i ].is_fixed() ) && ( v_fxd[ i ] != 0 ) )  
  return;                                               

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 if( not_dry_run( issueAMod ) ) {
  v_x[ i ].set_value( value ); 
  v_x[ i ].is_fixed( true , un_ModBlock( issueAMod ) );
  v_fxd[ i ] = value ? 2 : 1; 
  }
 else
  if( not_dry_run( issueMod ) )
   v_fxd[ i ] = value ? 2 : 1;
 
 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
  this , BinaryKnapsackBlockMod::eFixX , std::make_pair( i , i + 1 ) ) , 
         Observer::par2chnl( issueMod ) );

 } // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( c_boolVec_it value , Range rng , 
                                 ModParam issueMod, ModParam issueAMod )
{
 rng.second = std::min( rng.second , Index( v_x.size() ) );
 if( rng.second <= rng.first )  // nothing to do
  return;

 Index i = rng.first;
 for(  ; i < rng.second ; ++i )
  if( ( ! v_x[ i ].is_fixed() ) && ( v_fxd[ i ] == 0 ) )
   break;      

 if( i == rng.second )  // all fixed already
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 // TODO: use a GroupModification
 for( i = rng.first ; i < rng.second ; i++ ) {
  double val = *(value++); // new value

  if( ( ! v_x[ i ].is_fixed() ) && ( v_fxd[ i ] == 0 ) ) {
    if( not_dry_run( issueAMod ) ) {
     v_x[ i ].set_value( val );
     v_x[ i ].is_fixed( true , un_ModBlock( issueAMod ) );
     v_fxd[ i ] = val ? 2 : 1;
    }
    else
     if( not_dry_run( issueMod ) )
      v_fxd[ i ] = val ? 2 : 1; 
  }
 }

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
                           this , BinaryKnapsackBlockMod::eFixX , rng ) , 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::fix_x( c_boolVec_it value ,Subset && nms , 
                                 bool ordered , ModParam issueMod , 
                                 ModParam issueAMod )
{
 if( nms.empty() )
  return;

 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= v_x.size() )
  throw( std::invalid_argument( "BinaryKnapsackBlock::fix_x: invalid item"
        ) );
 bool done = true;
 for( auto i : nms )
  if( ( ! v_x[ i ].is_fixed() ) && ( v_fxd[ i ] == 0 ) ) {
   done = false;
   break;      
   }

 if( done )  // all fixed already
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 // TODO: use a GroupModification
 for( auto i : nms ) {
  double val = *(value++); // new value

  if( ( ! v_x[ i ].is_fixed() ) && ( v_fxd[ i ] == 0 ) ) {
    if( not_dry_run( issueAMod ) ) {
     v_x[ i ].set_value( val );
     v_x[ i ].is_fixed( true , un_ModBlock( issueAMod ) );
     v_fxd[ i ] = val ? 2 : 1;
    }
    else
     if( not_dry_run( issueMod ) )
      v_fxd[ i ] = val ? 2 : 1; 
   }
  }

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(
                           this , BinaryKnapsackBlockMod::eFixX , 
                           std::move( nms ) ) , 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::fix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( Index i , ModParam issueMod ,
                                   ModParam issueAMod )
{
 if( i >= v_x.size() )
  throw( std::invalid_argument( "BinaryKnapsackBlock::unfix_x: invalid item"
        ) );

 if( ( ! v_x[ i ].is_fixed() ) && ( v_fxd[ i ] == 0 ) )   // already unfixed
  return;                                                 // nothing to do

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();
 
 if( not_dry_run( issueAMod ) ) {
  v_x[ i ].is_fixed( false , un_ModBlock( issueAMod ) );
  v_fxd[ i ] = 0;
  }
 else
  if( not_dry_run( issueMod ) )
   v_fxd[ i ] = 0; 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
     this , BinaryKnapsackBlockMod::eUnfixX , std::make_pair( i , i + 1 ) ) , 
     Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( Range rng , ModParam issueMod ,
                                   ModParam issueAMod )
{
 rng.second = std::min( rng.second , Index( v_x.size() ) );
 if( rng.second <= rng.first )  // nothing to do
  return;

 Index i = rng.first;
 for(  ; i < rng.second ; ++i )
  if( v_x[ i ].is_fixed() && ( v_fxd[ i ] != 0 ) )
   break;      

 if( i == rng.second )  // all unfixed already
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 // TODO: use a GroupModification
 if( not_dry_run( issueAMod ) )
  for( i = rng.first ; i < rng.second ; ++i ) {
   v_x[ i ].is_fixed( false , un_ModBlock( issueAMod ) );
   v_fxd[ i ] = 0; 
  }
 else
  if( not_dry_run( issueMod ) ) {
   for( i = rng.first ; i < rng.second ; ++i )
    v_fxd[ i ] = 0;
   }

 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
                           this , BinaryKnapsackBlockMod::eUnfixX , rng ) , 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::unfix_x( Subset && nms , bool ordered ,
                                   ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )
  return;

 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= v_x.size() )
  throw( std::invalid_argument( "BinaryKnapsackBlock::unfix_x: invalid item"
        ) );

 bool done = true;
 for( auto i : nms )
  if( v_x[ i ].is_fixed() && ( v_fxd[ i ] != 0 ) ) {
   done = false;
   break;      
   }

 if( done )  // all unfixed already
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = Inf< double >();

 // TODO: use a GroupModification
 if( not_dry_run( issueAMod ) ) {
  for( auto i : nms ) {
   v_x[ i ].is_fixed( false , un_ModBlock( issueAMod ) );
   v_fxd[ i ] = 0;
   }
  }
 else
  if( not_dry_run( issueMod ) ) {
   for( auto i : nms )
    v_fxd[ i ] = 0;
   } 
 
 // issue physical Modification
 if( issue_pmod( issueMod ) )  
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >( 
                           this , BinaryKnapsackBlockMod::eUnfixX , 
                           std::move( nms ) ) , 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::unfix_x )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weight( double NWeight , Index item , 
                                      ModParam issueMod , ModParam issueAMod )
{
 if( item >= get_NItems() )
  throw( std::invalid_argument( "invalid item" ) );

 if( v_W[ item ] == NWeight )  // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ) {
  LF( f_cnst.get_function() )->modify_coefficient( item , NWeight ,
               un_ModBlock( issueAMod ) );
  v_W[ item ] = NWeight;
  }
 else
  if( not_dry_run( issueMod ) ) // otherwise only physical representation
   v_W[ item ] = NWeight; 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
                           this , BinaryKnapsackBlockMod::eChgWeight ,
                           std::make_pair( item , item + 1 ) ), 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_weight )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weights( c_dblVec_it NWeight , 
                                       Range rng , 
                                       ModParam issueMod , 
                                       ModParam issueAMod )
{
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
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ) {
  // physical representation
  std::copy( NWeight , NWeight + ( rng.second - rng.first ) ,
             v_W.begin() + rng.first );

  // abstract representation
  LF( f_cnst.get_function() )->modify_coefficients(
    doubleVec( NWeight , NWeight + ( rng.second - rng.first ) ) ,
    rng , un_ModBlock( issueAMod ) );
  }
 else
  if( not_dry_run( issueMod ) )
   // otherwise change only physical representation 
   std::copy( NWeight , NWeight + ( rng.second - rng.first ) ,
        v_W.begin() + rng.first );

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
        this , BinaryKnapsackBlockMod::eChgWeight , rng ) ,
         Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_weights( Range ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_weights( c_dblVec_it NWeight,
                                       Subset && nms , bool ordered ,  
                                       ModParam issueMod ,
                                       ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( v_W , nms , NWeight , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ) {
  // physical representation
  copyidx( v_W , nms , NWeight );
  
  // abstract representation
  LF( f_cnst.get_function() )->modify_coefficients(
             doubleVec( NWeight , NWeight + nms.size() ) , 
             Subset( nms ) , ordered ,
             un_ModBlock( issueAMod ) );
  } 
 else
  if( not_dry_run( issueMod ) )
   // otherwise change only physical representation 
   copyidx( v_W , nms , NWeight );

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) {
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(
       this , BinaryKnapsackBlockMod::eChgWeight , std::move( nms ) ) ,
         Observer::par2chnl( issueMod ) );
  }
 }  // end( BinaryKnapsackBlock::chg_weights( Subset ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profit( double NProfit , Index item , 
                                      ModParam issueMod , ModParam issueAMod )
{
 if( item >= get_NItems() )
  throw( std::invalid_argument(
        "BinaryKnapsackBlock::chg_profit: invalid item" ) );

 if( v_P[ item ] == NProfit )  // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // physical representation
  v_P[ item ] = NProfit; 

  // abstract representation
  LF( f_obj.get_function() )->modify_coefficient( item , NProfit ,
              un_ModBlock( issueAMod ) );
  }
 else
  if( not_dry_run( issueMod ) ) // otherwise only physical representation 
   v_P[ item ] = NProfit;

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
                           this , BinaryKnapsackBlockMod::eChgProfit ,
                           std::make_pair( item , item + 1 ) ), 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_profit )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profits( c_dblVec_it NProfit , Range rng , 
                                       ModParam issueMod ,
                                       ModParam issueAMod )
{
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
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  
  // physical representation
  std::copy( NProfit , NProfit + ( rng.second - rng.first ) ,
             v_P.begin() + rng.first );
  
  // abstract representation  
  LF( f_obj.get_function() )->modify_coefficients(
    doubleVec( NProfit , NProfit + ( rng.second - rng.first ) ) ,
    rng , un_ModBlock( issueAMod ) );
  }
 else
  if( not_dry_run( issueMod ) )  // otherwise only physical representation 
   std::copy( NProfit , NProfit + ( rng.second - rng.first ) ,
        v_P.begin() + rng.first );

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
          this , BinaryKnapsackBlockMod::eChgProfit , rng ) ,
         Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_profits( Range ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_profits( c_dblVec_it NProfit ,
                                       Subset && nms , bool ordered ,  
                                       ModParam issueMod ,
                                       ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( v_P , nms , NProfit , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // physical representation
  copyidx( v_P , nms , NProfit );
  
  // abstract representation
  LF( f_obj.get_function() )->modify_coefficients(
       doubleVec( NProfit , NProfit + nms.size() ) , 
       Subset( nms ) , ordered , un_ModBlock( issueAMod ) );
  } 
 else
  if( not_dry_run( issueMod ) )
   // otherwise change only physical representation 
   copyidx( v_P , nms , NProfit );

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) {
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(
        this , BinaryKnapsackBlockMod::eChgProfit , std::move( nms ) ) ,
         Observer::par2chnl( issueMod ) );
  }
 }  // end( BinaryKnapsackBlock::chg_profits( Subset ) )

/*--------------------------------------------------------------------------*/ 

void BinaryKnapsackBlock::chg_integrality( bool NIntegrality , Index item , 
                                           ModParam issueMod ,
                                           ModParam issueAMod )
{
 if( item >= get_NItems() )
  throw( std::invalid_argument(
       "BinaryKnapsackBlock::chg_integrality: invalid item" ) );
 
 if( v_I[ item ] == NIntegrality )  // nothing to do
  return;

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // physical representation
  v_I[ item ] = NIntegrality; 

  // abstract representation
  if( NIntegrality ) {
    if( v_x[ item ].get_type() == ColVariable::kPosUnitary ) {
     v_x[ item ].set_type( ColVariable::kBinary ,
         un_ModBlock( issueAMod ) );
     countCont--;
     }
   }
  else
   if( v_x[ item ].get_type() == ColVariable::kBinary ) {
    v_x[ item ].set_type( ColVariable::kPosUnitary ,
        un_ModBlock( issueAMod ) );
    countCont++;
    }
  } 
 else
  if( not_dry_run( issueMod ) )  // otherwise only physical representation 
  v_I[ item ] = NIntegrality;
 

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
                           this , BinaryKnapsackBlockMod::eChgIntegrality ,
                           std::make_pair( item , item + 1 ) ), 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_integrality )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_integrality( c_boolVec_it NIntegrality ,
                                           Range rng , ModParam issueMod ,
                                                       ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NItems() );

 if( rng.second <= rng.first )  // nothing to change
  return;   

 if( std::equal(  NIntegrality , NIntegrality + ( rng.second - rng.first ) ,
                  v_I.begin() + rng.first ) )
  return;  // nothing changes, avoid issuing the Modification

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // physical representation
  std::copy( NIntegrality , NIntegrality + ( rng.second - rng.first ) ,
             v_I.begin() + rng.first );
  
  // abstract representation
  for( Index i = 0 ; i < rng.second - rng.first ; ++i ) {
   if( NIntegrality[ i ] ) {
    if( v_x[ rng.first + i ].get_type() == ColVariable::kPosUnitary ) {
     v_x[ rng.first + i ].set_type( ColVariable::kBinary ,
            un_ModBlock( issueAMod ) );
     countCont--;
     }
    }
   else
    if( v_x[ rng.first + i ].get_type() == ColVariable::kBinary ) {
     v_x[ rng.first + i ].set_type( ColVariable::kPosUnitary ,
            un_ModBlock( issueAMod ) );
     countCont++;
     }
   }
  }
 else
  if( not_dry_run( issueMod ) )  // otherwise only physical representation 
  std::copy( NIntegrality , NIntegrality + ( rng.second - rng.first ) ,
             v_I.begin() + rng.first );

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockRngdMod >(
                           this , BinaryKnapsackBlockMod::eChgIntegrality , 
                           rng ) , Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_integrality( Range ) )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_integrality( c_boolVec_it NIntegrality,
                                           Subset && nms , bool ordered , 
                                           ModParam issueMod ,
                                           ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;            

 if( is_equal( v_I , nms , NIntegrality , get_NItems() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // physical representation
  copyidx( v_I , nms , NIntegrality );
  
 // abstract representation
  for(Index i = 0 ; i < nms.size() ; i++ ) {
   if( NIntegrality[ i ] ) {
    if( v_x[ nms[ i ] ].get_type() == ColVariable::kPosUnitary ) {
     v_x[ nms[ i ] ].set_type( ColVariable::kBinary ,
             un_ModBlock( issueAMod ) );
     countCont--;
     }
    }
   else
    if( v_x[ nms[ i ] ].get_type() == ColVariable::kBinary) {
     v_x[ nms[ i ] ].set_type( ColVariable::kPosUnitary ,
             un_ModBlock( issueAMod ) );
     countCont++;
     }
   }
  } 
 else
  if( not_dry_run( issueMod ) )
   // otherwise change only physical representation 
   copyidx( v_I , nms , NIntegrality );

 // issue physical Modification 
 if( issue_pmod( issueMod ) ) {
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< BinaryKnapsackBlockSbstMod >(
    this, BinaryKnapsackBlockMod::eChgIntegrality , std::move( nms ) ), 
         Observer::par2chnl( issueMod ) );
  }
 }  // end( BinaryKnapsackBlock::chg_integrality( Subset ) )

/*---------------------------------------------------------------------------*/

void BinaryKnapsackBlock::chg_capacity( double NC , ModParam issueMod ,
                      ModParam issueAMod )
{
 if( f_C == NC )  // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasCns ) ) {
  f_C = NC; 
  f_cnst.set_rhs( NC , un_ModBlock( issueAMod ) );
  } 
 else
  if( not_dry_run( issueMod ) )  // otherwise only physical representation 
   f_C = NC;

 if( issue_pmod( issueMod ) )  // issue physical Modification
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockMod >( this ,
                           BinaryKnapsackBlockMod::eChgCapacity ) , 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::chg_capacity )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::set_objective_sense( bool sense , ModParam issueMod , 
                                               ModParam issueAMod )
{
 if( f_sense == sense )  // nothing to do
  return;

 // reset conditional bounds
 f_cond_lower = -Inf< double >();
 f_cond_upper = +Inf< double >();

 // change both physical and abstract representation (if it exists)
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  f_sense = sense; 
  f_obj.set_sense( sense ? Objective::eMax : Objective::eMin ,
       un_ModBlock( issueAMod ) );
  } 
 else
  if( not_dry_run( issueMod ) )  // otherwise only physical representation 
   f_sense = sense;
 
 if( issue_pmod( issueMod ) ) // issue physical Modification 
  Block::add_Modification( std::make_shared< BinaryKnapsackBlockMod >( this ,
                           BinaryKnapsackBlockMod::eChgSense ) , 
                           Observer::par2chnl( issueMod ) );

 }  // end( BinaryKnapsackBlock::set_objective_sense )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::guts_of_destructor( void )
{
 f_cnst.clear();  // clear the constraint
 f_obj.clear();   // clear the objective function
 v_x.clear();     // clear all variables

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

void BinaryKnapsackBlock::guts_of_add_Modification( c_p_Mod mod ,
                ChnlName chnl )
{
 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinRngd * >( mod ) ) { 
  // the change may concern the objective function or the constraint
  // check which one it is
  auto lf = dynamic_cast< LinearFunction * >( tmod->function() );
  if( ! lf )
   throw( std::invalid_argument( "invalid Modification to Linear Function" ) ); 

  int range_size = tmod->range().second - tmod->range().first;

  if( AR & HasObj ) {  // if the AR of the objective exists
   // get the objective function
   auto lfo = dynamic_cast< LinearFunction * >( f_obj.get_function() ); 
   // and check if the modification is on the objective  
   if( lf == lfo ) {
    // vector of new profits 
    doubleVec new_profits( range_size );  
    dblVec_it npi = new_profits.begin(); 
    // get the values of the new profits from coefficients of lf
    for( Index i = tmod->range().first ; i < tmod->range().second ; i++ )
     ( * npi++ ) = lf->get_coefficient( i );

    chg_profits( new_profits.begin() , tmod->range() ,
     make_par( eNoBlck , chnl ) , eDryRun );
    return;
    }   
   }

  if( AR & HasCns ) {  // if the AR of the constraint exists
   // get the constraint 
   auto lfc = dynamic_cast< LinearFunction * >( f_cnst.get_function() ); 
   // and check if the modification is on the constraint  
   if( lf == lfc ) {
    // vector of new weights 
    doubleVec new_weights( range_size ); 
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

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinSbst * >( mod ) ) {
  // the change may concern the objective function or the constraint
  // check which one it is
  auto lf = dynamic_cast< LinearFunction * >( tmod->function() );
  if( ! lf )
   throw( std::invalid_argument( "invalid Modification to Linear Function" ) ); 

  if( AR & HasObj ) {  // if the AR of the objective exists
   // get the objective function
   auto lfo = dynamic_cast< LinearFunction * >( f_obj.get_function() ); 
   // and check if the modification is on the objective  
    if( lf == lfo ) {
     // vector of new profits 
     doubleVec new_profits( tmod->subset().size() );
     dblVec_it npi = new_profits.begin();   
     // get the values of the new profits from coefficients of lf
     for( auto i : tmod->subset() )
      ( * npi++ ) = lf->get_coefficient( i );

     chg_profits( new_profits.begin() , Subset( tmod->subset() ) , true , 
                  make_par( eNoBlck , chnl ) , eDryRun );
     return;
     }
   }

  if( AR & HasCns ) {  // if the AR of the constraint exists
   // get the constraint 
   auto lfc = dynamic_cast< LinearFunction * >( f_cnst.get_function() ); 
   // and check if the modification is on the constraint  
    if( lf == lfc ) {
     // vector of new weights 
     doubleVec new_weights( tmod->subset().size() );  
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

 // RowConstraintMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const RowConstraintMod * >( mod ) ) {
  if( ! ( AR & HasCns ) )
   throw( std::invalid_argument( "Modification to non-constructed constraint "
         ) );
 
   auto cp = dynamic_cast< FRowConstraint * >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) ); 

  if( tmod->type() == RowConstraintMod::eChgRHS ) {
   chg_capacity( cp->get_rhs(), make_par( eNoBlck , chnl ) , eDryRun );
   return; 
   }

  throw( std::invalid_argument( "illegal Modification to Constraint" ) );
  }

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 if( auto tmod = dynamic_cast< const VariableMod * >( mod ) ) {
  auto xi = dynamic_cast< ColVariable * const >( tmod->variable() );
  if( ! xi )
   throw( std::logic_error( "Modification to wrong type of Variable" ) );

  int i = p2i_x( xi );

  auto new_state = tmod->new_state();  // get new state of the variable
  if( ( ! ColVariable::is_unitary( new_state ) ) ||
      ( ! ColVariable::is_positive( new_state ) ) )
   throw( std::invalid_argument( "invalid ColVariable Modification" ) );

  auto old_state = tmod->old_state();  // get old state of the variable

  // if the "state of fixing" has changed- - - - - - - - - - - - - - - - - - -
  if( Variable::is_fixed( new_state ) != Variable::is_fixed( old_state ) ) {
   
   if( Variable::is_fixed( new_state ) ) 
    fix_x( xi->get_value() , i , make_par( eNoBlck , chnl ) , eDryRun );
   else
    unfix_x( i , make_par( eNoBlck , chnl ) , eDryRun );
  
   }
  
  // if the Integrality has changed- - - - - - - - - - - - - - - - - - - - - - 
  if( ColVariable::is_integer( new_state ) !=
      ColVariable::is_integer( old_state ) )
   chg_integrality( ColVariable::is_integer( new_state ) , i ,
        make_par( eNoBlck , chnl ) , eDryRun );

  return;
  }

 // ObjectiveMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const ObjectiveMod * >( mod ) ) {
  if( ! ( AR & HasObj ) )  // check if the objective exists
   throw( std::invalid_argument( "Modification to non-constructed objective" )
    );

  auto obj = dynamic_cast< Objective * >( & f_obj );
  if( tmod->of() == obj ) {
   bool sense = tmod->type() == Objective::eMax ? 1 : 0;
   set_objective_sense( sense , make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  throw( std::invalid_argument( "Modification to the wrong objective" ) );
  }

 throw( std::invalid_argument( 
                        "Unsupported Modification to BinaryKnapsackBlock" ) );

}  // end( BinaryKnapsackBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackBlock::compute_conditional_bounds( void )
{
 f_cond_lower = 0;
 f_cond_upper = 0;

 for( Index i = 0 ; i < get_NItems() ; i++ ) { 
  double w = v_W[ i ];                        // weight of the current item
  double p = f_sense ? v_P[ i ] : -v_P[ i ];  // profit of the current item

  // items contained in the optimal solution
  if( ( w <= 0 ) && ( p >= 0 ) ) {
   f_cond_lower += p;                
   f_cond_upper += p;  
   continue;
   }

  // items not contained in the optimal solution
  if( ( w >= 0 ) && ( p <= 0 ) )
   continue;

  // remaining items
  if( p > 0 )
   f_cond_upper += p;
  else 
   f_cond_lower += p;
  }

 // if it is a minimization problem
 if( ! f_sense ) {
  f_cond_lower = -f_cond_lower;
  f_cond_upper = -f_cond_upper;
  std::swap( f_cond_lower , f_cond_upper );  
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF BinaryKnapsackSolution --------------------*/
/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::deserialize( const netCDF::NcGroup & group )
{ 
 netCDF::NcDim ni = group.getDim( "n" );
 if( ni.isNull() )
  v_x.clear();
 else {
  netCDF::NcVar bx = group.getVar( "x" );
  if( bx.isNull() )
   v_x.clear();
  else {
   v_x.resize( ni.getSize() );
   bx.getVar( v_x.data() );
   }
  }
 }  // end( BinaryKnapsackSolution::deserialize )

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::read( const Block * const block )
{
 auto BKB = dynamic_cast< const BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "block is not a BinaryKnapsackBlock" ));

 v_x.resize( BKB->get_NItems() );

 auto vxi = v_x.begin();

 for( auto & xi : BKB->v_x ) 
  *(vxi++) = static_cast< double >( xi.get_value() );
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::write( Block * const block )
{
 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "block is not a BinaryKnapsackBlock" ) );

 if( ! v_x.empty() ) {
  if( v_x.size() < BKB->get_NItems() )
   throw( std::invalid_argument( "incompatible variables size" ) );

  auto vxi = v_x.begin();
  for( auto & xi : BKB->v_x )
   xi.set_value( *(vxi++) );
  }
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::serialize( netCDF::NcGroup & group ) const
{
 if( ! v_x.empty() ) {

  // always call the method of the base class first
  Solution::serialize( group );
  
  netCDF::NcDim ni = group.addDim( "n" , v_x.size() ); 
  ( group.addVar( "x" , netCDF::NcDouble() , ni ) ).putVar( v_x.data() );
  }
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::print( std::ostream & output )
{
 for( Index i = 0 ; i < v_x.size() ; ++i )
  output << "x" << i << ": " << v_x[ i ] << std::endl;
 }

/*--------------------------------------------------------------------------*/

BinaryKnapsackSolution * BinaryKnapsackSolution::scale( double factor ) const
{
 auto * sol = BinaryKnapsackSolution::clone( true );

 if( ! v_x.empty() )
  for( Index i = 0 ; i < v_x.size() ; i++ )
   sol->v_x[ i ] = v_x[ i ] * factor;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void BinaryKnapsackSolution::sum( const Solution * solution , 
                                  double multiplier )
{
 auto BKB = dynamic_cast< const BinaryKnapsackSolution * >( solution );
 if( ! BKB )
  throw( std::invalid_argument( "solution is not a BinaryKnapsackSolution" )
   );

 if( ! v_x.empty() ) {
  if( v_x.size() != BKB->v_x.size() )
   throw( std::invalid_argument( "incompatible variables size" ) );
  
  for( Index i = 0 ; i < v_x.size() ; i++ )
   v_x[ i ] += BKB->v_x[ i ] * multiplier;
  }
 }

/*--------------------------------------------------------------------------*/

BinaryKnapsackSolution * BinaryKnapsackSolution::clone( bool empty ) const
{
 auto * sol = new BinaryKnapsackSolution();
 
 if( empty ) {
  if( ! v_x.empty() )
   sol->v_x.resize( v_x.size() );
  }
 else 
  sol->v_x = v_x;
 
 return( sol );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF BinaryKnapsackBlockChange-------------------*/
/*--------------------------------------------------------------------------*/

Change * BinaryKnapsackBlockChange::apply( Block * block , 
                                           bool doUndo , 
                                           ModParam issueMod , 
                                           ModParam issueAMod ) {
 
 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "Block must be a BinaryKnapsackBlock" ) );

 switch( f_type ) {
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eEmpty:
   throw( std::invalid_argument( "Empty BinaryKnapsackBlockChange" ) );
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgSense: {

   Change * undoChg = nullptr;

   if( doUndo ) {
    bool old_sense = ( BKB->get_objective_sense() ==  Objective::eMax );
    undoChg = new BinaryKnapsackBlockChange( eChgSense , 
                                            { double( old_sense ) } ); 
    }
   
   // get new objective sense
   bool new_sense = v_data[ 0 ];

   // Apply the Change to BinaryKnapsackBlock
   BKB->set_objective_sense( new_sense , issueMod , issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgCapacity: {

   Change * undoChg = nullptr;

   if( doUndo ) {
    double old_C = BKB->get_Capacity();
    undoChg = new BinaryKnapsackBlockChange( eChgCapacity , { old_C } );
   }

   double new_C = v_data[ 0 ];

   // Apply the Change to BinaryKnapsackBlock
   BKB->chg_capacity( new_C , issueMod , issueAMod );
   
   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  default:
   throw( std::invalid_argument( 
          "BinaryKnapsackBlockChange type not supported" ) );
  }
 }

/*--------------------------------------------------------------------------*/

Change * BinaryKnapsackBlockRngdChange::apply( Block * block ,
                                               bool doUndo ,
                                               ModParam issueMod ,
                                               ModParam issueAMod ) {
 
 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "Block must be a BinaryKnapsackBlock" ) );

 switch( f_type ) {
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eEmpty:
   throw( std::invalid_argument( "Empty BinaryKnapsackBlockChange" ) );
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgWeight: {

   // check data size and range 
   if( v_data.size() != f_rng.second - f_rng.first )
    throw( std::invalid_argument( "v_data.size() != range size" ) );

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< double > old_data( f_rng.second - f_rng.first ); 
    BKB->get_Weights( old_data.begin() , f_rng );
    undoChg = new BinaryKnapsackBlockRngdChange( eChgWeight , 
                                                 std::move( old_data ) , 
                                                 f_rng );
    }

   // Change data
   BKB->chg_weights( v_data.begin() , f_rng , issueMod , issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgProfit: {
  
   // check data size and range 
   if( v_data.size() != f_rng.second - f_rng.first )
    throw( std::invalid_argument( "v_data.size() != range size" ) );

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< double > old_data( f_rng.second - f_rng.first );
    BKB->get_Profits( old_data.begin() , f_rng );
    undoChg = new BinaryKnapsackBlockRngdChange( eChgProfit , 
                                                 std::move( old_data ) , 
                                                 f_rng );
    }
  
   // Change data
   BKB->chg_profits( v_data.begin() , f_rng , issueMod , issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgIntegrality: {
  
   // check data size and range 
   if( v_data.size() != f_rng.second - f_rng.first )
    throw( std::invalid_argument( "v_data.size() != range size" ) );

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< bool > old_int( f_rng.second - f_rng.first );
    BKB->get_Integrality( old_int.begin() , f_rng );
    std::vector< double > old_data( old_int.begin() , old_int.end() );
    undoChg = new BinaryKnapsackBlockRngdChange( eChgIntegrality , 
                                                 std::move( old_data ) , 
                                                 f_rng );
    }
  
   // Change data (cast to bool first)
   std::vector< bool > new_int( v_data.begin() , v_data.end() );
   BKB->chg_integrality( new_int.begin() , f_rng , issueMod , issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eFixX: {

   // check data size and range 
   if( v_data.size() != f_rng.second - f_rng.first )
    throw( std::invalid_argument( "v_data.size() != range size" ) );

   // if some of the variables are already fixed, throw an exception. This is
   // needed to provide an undo change 
   for( Index i = f_rng.first ; i < f_rng.second ; ++i ) {
    if( BKB->is_fixed( i ) )
     throw( std::invalid_argument( "variable " + std::to_string( i ) + 
                                   "already fixed" ) ); 
    }

   Change * undoChg = nullptr;

   if( doUndo )
    undoChg = new BinaryKnapsackBlockRngdChange( eUnfixX , {} , f_rng );
   
   // Change data
   std::vector< bool > new_x( v_data.begin() , v_data.end() ); 
   BKB->fix_x( new_x.begin() , f_rng , issueMod , issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eUnfixX: {

   // if some of the variables are not fixed, throw an exception. This is
   // needed to provide an undo change 
   for( Index i = f_rng.first ; i < f_rng.second ; ++i ) {
    if( ! BKB->is_fixed( i ) )
     throw( std::invalid_argument( "variable " + std::to_string( i ) + 
                                   " is not fixed" ) ); 
    }

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< double > old_data( f_rng.second - f_rng.first );
    BKB->get_x( old_data.begin() , f_rng );
    undoChg = new BinaryKnapsackBlockRngdChange( eFixX , 
                                                 std::move( old_data ) ,
                                                 f_rng );
    }


   // Change data
   BKB->unfix_x( f_rng );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  default:
   throw( std::invalid_argument( 
          "BinaryKnapsackBlockChange type not supported" ) );
  }

 }

/*--------------------------------------------------------------------------*/

Change * BinaryKnapsackBlockSbstChange::apply( Block * block ,
                                               bool doUndo ,
                                               ModParam issueMod ,
                                               ModParam issueAMod ) {

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >( block );
 if( ! BKB )
  throw( std::invalid_argument( "Block must be a BinaryKnapsackBlock" ) );

 switch( f_type ) {
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eEmpty:
   throw( std::invalid_argument( "Empty BinaryKnapsackBlockChange" ) );
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgWeight: {

   // check data size and subset size 
   if( v_data.size() != v_nms.size() )
    throw( std::invalid_argument( "v_data.size() != v_nms.size()" ) );

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< double > old_data( v_nms.size() );
    BKB->get_Weights( old_data.begin() , Subset( v_nms ) );
    undoChg = new BinaryKnapsackBlockSbstChange( eChgWeight , 
                                                 std::move( old_data ) , 
                                                 Subset( v_nms ) );
    }
   
   // Change data
   BKB->chg_weights( v_data.begin() , 
                     Subset( v_nms ) , 
                     issueMod , 
                     issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgProfit: {

   // check data size and subset size 
   if( v_data.size() != v_nms.size() )
    throw( std::invalid_argument( "v_data.size() != v_nms.size()" ) );

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< double > old_data( v_nms.size() );
    BKB->get_Profits( old_data.begin() , Subset( v_nms ) );
    undoChg = new BinaryKnapsackBlockSbstChange( eChgProfit , 
                                                 std::move( old_data ) ,
                                                 Subset( v_nms ) );
    }

   // Change data
   BKB->chg_profits( v_data.begin() , 
                     Subset( v_nms ) , 
                     issueMod , 
                     issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eChgIntegrality: {

   // check data size and subset size 
   if( v_data.size() != v_nms.size() )
    throw( std::invalid_argument( "v_data.size() != v_nms.size()" ) );
   
   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< bool > old_int( v_nms.size() ); 
    BKB->get_Integrality( old_int.begin() , Subset( v_nms ) );
    std::vector< double > old_data( old_int.begin() , old_int.end() );
    undoChg = new BinaryKnapsackBlockSbstChange( eChgIntegrality , 
                                                 std::move( old_data ) ,
                                                 Subset( v_nms ) );
    }

   // Change data (cast to bool first)
   std::vector< bool > new_int( v_data.begin() , v_data.end() );
   BKB->chg_integrality( new_int.begin() , 
                         Subset( v_nms ) , 
                         issueMod , 
                         issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eFixX: {

   // check data size and subset size 
   if( v_data.size() != v_nms.size() )
    throw( std::invalid_argument( "v_data.size() != v_nms.size()" ) );

   // if some of the variables are already fixed, throw an exception. This is
   // needed to provide an undo change 
   for( Index i : v_nms ) {
    if( BKB->is_fixed( i ) )
     throw( std::invalid_argument( "variable " + std::to_string( i ) + 
                                   "already fixed" ) ); 
    }

   Change * undoChg = nullptr;

   if( doUndo )
    undoChg = new BinaryKnapsackBlockSbstChange( eUnfixX , {} , 
                                                 Subset( v_nms ) );
   
   // Change data
   std::vector< bool > new_x( v_data.begin() , v_data.end() ); 
   BKB->fix_x( new_x.begin() , Subset( v_nms ) , issueMod , issueAMod );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case eUnfixX: {

   // if some of the variables are not fixed, throw an exception. This is
   // needed to provide an undo change 
   for( Index i : v_nms ) {
    if( ! BKB->is_fixed( i ) )
     throw( std::invalid_argument( "variable " + std::to_string( i ) + 
                                   " is not fixed" ) ); 
    }

   Change * undoChg = nullptr;

   if( doUndo ) {
    std::vector< double > old_data( v_nms.size() );
    BKB->get_x( old_data.begin() , Subset( v_nms ) );
    undoChg = new BinaryKnapsackBlockSbstChange( eFixX , 
                                                 std::move( old_data ) ,
                                                 Subset( v_nms ) );
    }


   // Change data
   BKB->unfix_x( Subset( v_nms ) );

   return( undoChg );
   }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  default:
   throw( std::invalid_argument( 
          "BinaryKnapsackBlockChange type not supported" ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File BinaryKnapsackBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/

