/*--------------------------------------------------------------------------*/
/*---------------------------- File bk2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing instance files, be them netCDF ones or text
 * ones, of the (Binary) Knapsack problem by reading them into the
 * BinaryKnapsackBlock (in whatever format it supports) and either
 * print()-ing to a text file (in whatever format this is supported) or
 * serialize()-ing into a netCDF file.
 *
 * The Pisinger/Jooken knapsack benchmark format is read with frmt_in == 'P'
 * [see BinaryKnapsackBlock::load( std::istream & , char )].
 *
 * While this main() is written for BinaryKnapsackBlock, in fact it does not
 * even include BinaryKnapsackBlock.h: by having the name used in the factory
 * call (new_Block), say, be provided in the command line it would work with
 * any other kind of :Block.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include "Block.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Custom terminate function to print the exception message

void smspp_terminate( void ) {
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
 std::abort(); // or exit(1)
}

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // override the default terminate handler to print the exception message
 std::set_terminate( smspp_terminate );

 if( argc < 4 ) {
  std::cerr << "Usage: " << argv[ 0 ]
	    << " file_in frmt_in file_out [frmt_out]" << std::endl
	    << "        if file_in ends in .nc4 a netCDF file is assumed"
	    << std::endl
	    << "        else frmt is the parameter of load()"
	    << std::endl
	    << "        ('P' = Pisinger/Jooken knapsack benchmark format)"
	    << std::endl
	    << "        if file_out ends in .nc4 a netCDF file is created"
	    << std::endl
	    << "        else frmt_out (if any) is the parameter of print()"
	    << std::endl;
  return( 1 );
  }

 char frmt_in = argv[ 2 ][ 0 ];
 char frmt_out = argc > 4 ? argv[ 4 ][ 0 ] : 'C';

 // either deserialize() or load() the [BinaryKnapsack]Block
 Block * BKB;

 std::string iname( argv[ 1 ] );
 if( ( iname.size() > 4 ) &&
     ( iname.substr( iname.size() - 4 , 4 ) == ".nc4" ) )  // a netCDF file
  BKB = Block::deserialize( iname );
 else {                                                    // a text file
  // construct the [BinaryKnapsack]Block via the factory
  BKB = Block::new_Block( "BinaryKnapsackBlock" );
  BKB->load( iname , frmt_in );  // now load() it
  }

 // either serialize() or print() the [BinaryKnapsack]Block
 std::string oname( argv[ 3 ] );
 if( ( oname.size() > 4 ) &&
     ( oname.substr( oname.size() - 4 , 4 ) == ".nc4" ) ) {
  BKB->Block::serialize( oname , eBlockFile );
  // why the Block:: is needed completely evades me, but clang++ seems to
  // think it is
  }
 else
  BKB->print( oname , frmt_out );

 // cleanup
 delete BKB;

 // all done
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File bk2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
