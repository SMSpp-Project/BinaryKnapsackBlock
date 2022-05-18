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
 auto CHG = Change::new_Change( chgname );
 auto BKCHG = dynamic_cast< BinaryKnapsackBlockChange * >( CHG );
 if( ! BKCHG ) {
  cerr << "Unable to cast to BinaryKnapsackBlockChange";
  exit( 1 );
 }

 // Change Capacity
 BKCHG->set_type( BinaryKnapsackBlockChange::eChgCapacity );
 BKCHG->set_data( { 50 } );
 
 cout << BKB->get_Capacity() << endl;
 auto UNDO = CHG->apply( BKB , true );
 cout << BKB->get_Capacity() << endl;
 UNDO->apply( BKB );
 cout << BKB->get_Capacity() << endl;

 // Change Sense of the objective
 BKCHG->set_type( BinaryKnapsackBlockChange::eChgSense );
 BKCHG->set_data( { 0 } );
 
 cout << BKB->get_objective_sense() << endl;
 UNDO = CHG->apply( BKB , true );
 cout << BKB->get_objective_sense() << endl;
 UNDO->apply( BKB );
 cout << BKB->get_objective_sense() << endl;

 delete CHG;
 delete BKB;

 return( 0 );
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File change.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
