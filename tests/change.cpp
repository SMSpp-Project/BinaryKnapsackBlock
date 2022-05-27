/*--------------------------------------------------------------------------*/
/*-------------------------- File change.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test for BinaryKnapsackBlockChange.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ MACROS ------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- INCLUDES -----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackBlock.h"

#include <iostream>

/*--------------------------------------------------------------------------*/
/*------------------------------- USING ------------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace std;

/*--------------------------------------------------------------------------*/
/*------------------------------- TYPES ------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv ) {
 
 // Create and load a BinaryKnapsackBlock
 auto BKB = new BinaryKnapsackBlock();
 BKB->load( 5 , 10 , { 1 , 3 , 2 , 3 , 5 } , { 3 , 2 , 5 , 2 , 5 } );

 // Create a BinaryKnapsackBlockChange from the factory
 string chgname = "BinaryKnapsackBlockChange";
 Change * CHG = Change::new_Change( chgname );
 auto BKCHG = dynamic_cast< BinaryKnapsackBlockChange * >( CHG );
 if( ! BKCHG ) {
  cerr << "Unable to cast to BinaryKnapsackBlockChange";
  exit( 1 );
 }

 BKCHG->set_type( BinaryKnapsackBlockChange::eChgCapacity );
 BKCHG->set_data( { 50 } );

 CHG->serialize( "file.nc4" );

 Change * chg;
 chg = chg->deserialize( "file.nc4" );
 auto bkchg = dynamic_cast< BinaryKnapsackBlockChange * >( chg );

 std::cout << *chg;

 delete CHG;
 delete chg;
 return( 0 );
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File change.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
