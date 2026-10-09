/*--------------------------------------------------------------------------*/
/*----------------------------- File bk2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing BinaryKnapsackBlock netCDF files out of
 * textual ones. It creates one BinaryKnapsackBlock, loads it from a text file
 * in either the Pisinger/Jooken benchmark format or the native one, and
 * serializes it to a netCDF file.
 *
 * The sense of the objective can be reversed while doing so: a knapsack is
 * written as a maximization problem, but whoever uses it as a component of a
 * decomposition of a minimization problem needs it as a minimization one, and
 * the two are the same problem provided that the profits change sign.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>

#include <BinaryKnapsackBlock.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// Custom terminate function to print the exception message

static void smspp_terminate( void ) {
 std::cerr << "Uncaught exception in executing SMS++:\n";
 try {
  std::rethrow_exception( std::current_exception() );
 }
 catch( const std::exception & e ) {
  std::cerr << "\tException type: " << typeid( e ).name() << "\n";
  std::cerr << "\tException message: " << e.what() << "\n";
 } catch( ... ) {
  std::cerr << "\tUnknown exception" << std::endl;
 }
 std::abort();
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // override the default terminate handler to print the exception message
 std::set_terminate( smspp_terminate );

 char frmt = 'P';
 bool minimize = false;

 switch( argc ) {
  case( 4 ): minimize = ( std::string( argv[ 3 ] ) == "min" );
  case( 3 ): frmt = argv[ 2 ][ 0 ];
  case( 2 ): break;
  default:
   std::cerr << "Usage: " << argv[ 0 ] << " knapsack_in [ format [ sense ] ]"
	     << std::endl
	     << "        format P = Pisinger/Jooken, N = native [P]"
	     << std::endl
	     << "        sense min reverses it, changing the sign of the"
	     << std::endl
	     << "        profits, so that the problem is the same one [max]"
	     << std::endl
	     << "        if knapsack_in ends in .txt the suffix is removed"
	     << std::endl
	     << "        and knapsack_in[-suffix][-min].nc4 is created"
	     << std::endl;
   return( 1 );
  }

 std::ifstream ProbFile( argv[ 1 ] );
 if( ! ProbFile.is_open() ) {
  std::cerr << "Error: cannot open file " << argv[ 1 ] << std::endl;
  return( 1 );
  }

 auto BKB = dynamic_cast< BinaryKnapsackBlock * >(
                                 Block::new_Block( "BinaryKnapsackBlock" ) );
 if( ! BKB ) {
  std::cerr << "Failed to initialize BinaryKnapsackBlock" << std::endl;
  return( 1 );
  }

 BKB->load( ProbFile , frmt );
 ProbFile.close();

 // the same problem written the other way round: whoever minimizes wants the
 // profits with the other sign, the optimal solution being the same one
 if( minimize ) {
  auto n = BKB->get_NItems();
  std::vector< double > P( n );
  for( decltype( n ) i = 0 ; i < n ; ++i )
   P[ i ] = - BKB->get_Profit( i );

  BKB->chg_profits( P.begin() , Block::Range( 0 , n ) );
  BKB->set_objective_sense( false );
  }

 std::string name( argv[ 1 ] );
 if( name.size() > 4 ) {
  std::string sffx = name.substr( name.size() - 4 , 4 );
  std::string txt( ".txt" );

  if( std::equal( sffx.begin() , sffx.end() , txt.begin() ,
		  []( auto a , auto b ) {
		   return( std::tolower( a ) == std::tolower( b ) ); } ) )
   name.erase( name.size() - 4 , 4 );
  }

 if( minimize )
  name.append( "-min" );
 name.append( ".nc4" );

 netCDF::NcFile f( name , netCDF::NcFile::replace );

 f.putAtt( "SMS++_file_type" , netCDF::NcInt() , eBlockFile );

 BKB->Block::serialize( f , eBlockFile );

 delete BKB;

 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File bk2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
