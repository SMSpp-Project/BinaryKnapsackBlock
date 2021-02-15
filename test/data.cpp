#include <fstream>
#include <ctime>
#include "BinaryKnapsackBlock.h"

#define MAX_N 7
#define MAX_W 50
#define MAX_P 20

/* Generate instances for the Binary Knapsack problem.
   - txt file with the following format:
   n            Number of items
   C            Capacity of the Knapsack
   w1 ... wn    vector of the weigths
   p1 ... pn    vector of the profits               

   - netCDF file                                    */

using namespace SMSpp_di_unipi_it;

int main(){
    // initialize random seed
    srand( time( NULL ) );
    int n;
    double C;

    int ncidp;
    // create netCDF file
    nc_create( "../data/data.nc" , NC_CLOBBER , &ncidp );
    // Open and put attribute eBlockFile
    netCDF::NcFile nf("../data/data.nc", netCDF::NcFile::replace );
    nf.putAtt( "SMS++_file_type" , netCDF::NcInt() , 0 ); // eProbFile 

    std::ofstream file;
    std::ifstream ifile;

    for( int i = 0 ; i < 50 ; i++ ){
        // create and open a new txt file
        file.open("../data/data" + std::to_string(i) + ".txt" );
        
        // generate instances
        n = std::rand() % MAX_N + 3;
        C = ( MAX_W * n ) / 3;

        file << n << "\n";
        file << C << "\n";

        for( int j = 0 ; j < n ; j++ ){
            file << std::to_string( rand() %  MAX_W + 1 ) << " ";
        }

        file << "\n";

        for( int j = 0 ; j < n ; j++ ){
            file << std::to_string( rand() % MAX_P + 1 ) << " ";
        }
        file << "\n";

        
        file.close();
        
        // SERIALIZE - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ifile.open( "../data/data" + std::to_string(i) + ".txt" );
        
        // add Groups in netCDF file
        netCDF::NcGroup prob = nf.addGroup( "Prob_" + std::to_string(i) );
        prob.addGroup( "BlockConfig" );
        prob.addGroup( "BlockSolver" );
        netCDF::NcGroup group = prob.addGroup( "Block" );
        group.putAtt( "type" , "BinaryKnapsackBlock" );

        if( group.isNull() )
            throw(std::runtime_error("error creating the Block"));        

        // group.addGroup( "Solution" );

        // create a BinaryKnapsack object
        BinaryKnapsackBlock * bkb = new BinaryKnapsackBlock(); 
        
        // load content of the file
        ifile >> *bkb;
        // serialize in the netCDF file
        bkb->serialize( group );        

        // Close file and delete Knapsack
        ifile.close();
        delete bkb;
    }
    return 0;
}