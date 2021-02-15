#include "BinaryKnapsackBlock.h"
#include "CPXMILPSolver.h"
#include <fstream>
#include <random>
#include <ctime>


using namespace SMSpp_di_unipi_it;

void deserialize( BinaryKnapsackBlock * bkb , int Nprob );
void load( BinaryKnapsackBlock * bkb , int Nprob );
void printsolution( BinaryKnapsackBlock * bkb , CPXMILPSolver * sol , int status );

int main( int argc , char **argv ){

    if( argc == 1 ){
        std::cout << "Insert the number of the problem\n";
        return( 0 );
    }

    // initialize random seed
    srand( time( NULL ) );
     
    // Problem number
    int Nprob = std::stoi( argv[ 1 ] );
    std::cout << "\n\tProblem number: " << argv[ 1 ] << "\n\n";

    // create the BinaryKnapsackBlock and deserialize/load the problem
    BinaryKnapsackBlock * bkb = new BinaryKnapsackBlock();
    deserialize( bkb , Nprob );
    //load( bkb , Nprob );

    // print problem
    std::cout << *bkb;

    // SOLVER - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    CPXMILPSolver * sol = new CPXMILPSolver();
    bkb->register_Solver( sol );
    sol->set_Block( bkb );
    int status;

    sol->write_lp("model.lp");
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // MODIFICATION - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // change profit of a random item
    int item = rand() % bkb->get_NItems();
    bkb->chg_profit( bkb->get_Profit( item ) * 2 , item );
    std::cout << "New Profit of x" << item << ": " << bkb->get_Profit(item);

    // change weight of a random item
    item = rand() % bkb->get_NItems();
    bkb->chg_weight( bkb->get_Weight( item ) * 3 , item );
    std::cout << "\nNew Weight of x" << item << ": " << bkb->get_Weight(item);
    
    // change capacity - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    bkb->chg_capacity( bkb->get_Capacity() / 2 );
    std::cout << "\nNew capacity " << bkb->get_Capacity() << "\n\n";

    // print new problem
    std::cout << *bkb;

    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // fix item - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // random item
    item = rand() % bkb->get_NItems();
    
    // if its current value is 0, set it (and fix it) to 1 and viceversa
    if( bkb->get_x( item ) )
        bkb->fix_x( 0 , item );
    else
        bkb->fix_x( 1 , item );

    std::cout << "\n item x" << item << " fixed " << bkb->get_x(item) << "\n";
    
    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // unfix item - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // unfix previous fixed item
    bkb->unfix_x( item );
    std::cout << "\n item x" << item << " unfixed \n";
    
    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // change contiguous interval of weights - - - - - - - - - - - - - - - - - -
    // extreme of the interval
    int l = rand() % ( (int) bkb->get_NItems() / 2 );
    int r = l + rand() % ( bkb->get_NItems() - l ) ;
    Block::Range rng = std::make_pair( l , r );

    std::cout << "\nDouble weights of items " << l << " to " << r-1 << "\n\n";

    std::vector< double > new_weights( r - l );
    for( int i = 0 ; i < r - l ; i++ )
        new_weights[ i ] = bkb->get_Weight( i + l ) * 2 ;
    bkb->chg_weights( new_weights.begin() , rng );

    // print the problem
    std::cout << *bkb;
    
    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );
    
    // change contiguous interval of profits - - - - - - - - - - - - - - - - - -
    std::cout << "\nDouble profits of items " << l << " to " << r-1 << "\n\n";

    std::vector< double > new_profits( r - l );
    for( int i = 0 ; i < r - l ; i++ )
        new_profits[ i ] = bkb->get_Profit( i + l ) * 2 ;
    bkb->chg_profits( new_profits.begin() , std::make_pair( l , r ) );

    // print the problem
    std::cout << *bkb;
    
    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // --------------------------------------------------------------------------
    for( int i = 0 ; i < 80 ; i++ )
        std::cout << "-";
    std::cout << "\n";
    std::cout << "\n\n CHANGE ABSTRACT REPRESENTATION\n\n";

    FRowConstraint * cnst = bkb->get_Cnst();
    FRealObjective * obj = bkb->get_Obj();
    // Linear function pointer
    LinearFunction * lf;

    std::cout << "\nRanged: change coefficients of the constraint (weight)\n\n";
    
    // New weights
    std::vector<double> NWeight( r - l );
    for( int i = 0 ; i < r - l ; i++ )
        NWeight[ i ] = bkb->get_Weight( i + l ) / 3 ;
    
    // Modify abstract representation
    lf = dynamic_cast<LinearFunction*>( cnst->get_function() ); 
    lf->modify_coefficients( std::move( NWeight ) , rng );

    // print problem
    std::cout << *bkb;

    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // --------------------------------------------------------------------------
 
    std::cout << "\nSubset: change coefficients of the constraint (weight)\n";

    // vector of indexes
    Block::Subset items( bkb->get_NItems() );
    for( int i = 0 ; i < items.size() ; i++ )
        items[ i ] = i;

    // select half of the items
    Block::Subset nms( (int) bkb->get_NItems() / 2 );
    std::sample( items.begin() , items.end() , nms.begin() , 
                 nms.size() , std::mt19937{std::random_device{}()} );
    
    std::cout << "Double weights of items: ";
    for( auto i : nms )
        std::cout << i << " ";
    std::cout << "\n\n";

    std::vector< double > nw( nms.size() );
    auto nmsi = nms.begin();
    for( auto &w : nw )
        w = bkb->get_Weight( *nmsi++ ) * 2;

    lf->modify_coefficients( std::move( nw ) , std::move( nms ) );

    // print the problem
    std::cout << *bkb;

    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // --------------------------------------------------------------------------   

        std::cout << "\nRanged: change coefficients of the objective (profits)\n\n";
    
    // New profits
    std::vector<double> NProfit( r - l );
    for( int i = 0 ; i < r - l ; i++ )
        NProfit[ i ] = bkb->get_Profit( i + l ) / 3 ;
    
    // Modify abstract representation
    lf = dynamic_cast<LinearFunction*>( obj->get_function() ); 
    lf->modify_coefficients( std::move( NProfit ) , rng );

    // print problem
    std::cout << *bkb;

    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );

    // --------------------------------------------------------------------------
 
    std::cout << "\nSubset: change coefficients of the objective (profit)\n";

    // select half of the items
    Block::Subset nms1( (int) bkb->get_NItems() / 2 );
    std::sample( items.begin() , items.end() , nms1.begin() , 
                 nms1.size() , std::mt19937{std::random_device{}()} );
    
    std::cout << "Double weights of items: ";
    for( auto i : nms1 )
        std::cout << i << " ";
    std::cout << "\n\n";

    std::vector< double > np( nms1.size() );
    auto nms1i = nms1.begin();
    for( auto &p : np )
        p = bkb->get_Weight( *nms1i++ ) * 2;

    lf->modify_coefficients( std::move( np ) , std::move( nms1 ) );

    // print the problem
    std::cout << *bkb;

    // Solve again
    status = sol->compute();
    if( status )
        std::cout << "\nSolved with Cplex\n\n";
    printsolution( bkb , sol , status );


    delete sol;
    delete bkb;
    return 0;
}


void deserialize( BinaryKnapsackBlock * bkb , int Nprob ){
    // open netCDF file
    netCDF::NcFile nf("../data/data.nc", netCDF::NcFile::read );
    netCDF::NcGroup prob = nf.getGroup( "Prob_" + std::to_string(Nprob) );
    netCDF::NcGroup group = prob.getGroup( "Block" );
    if( group.isNull() )
        throw(std::runtime_error("error reading the Block"));

    bkb->deserialize( group );
    nf.close();
}

void load( BinaryKnapsackBlock * bkb , int Nprob ){
    std::ifstream file("../data/data" + std::to_string(Nprob) + ".txt");
    file >> *bkb;
    file.close();
}

void printsolution( BinaryKnapsackBlock * bkb , CPXMILPSolver * sol , int status ){

    if( status != 10 ){
        std::cout << "Not solved to optimality\n";
        for( int i = 0 ; i < 80 ; i++ )
            std::cout << "-";
        std::cout << "\n";
        return;
    }

    sol->get_var_solution();
    std::cout << "Current solution:\n";

    for( int i = 0 ; i < bkb->get_NItems() ; i++ ){
        std::cout << "x" + std::to_string(i) << "\t";
    }
    std::cout << "\n";

    for( int i = 0 ; i < bkb->get_NItems() ; i++ ){
        std::cout << bkb->get_x(i) << "\t";
    }
    std::cout << "\n\n";

    std::cout << "Lower bound: " << bkb->get_valid_lower_bound() << 
        "   "  << "Upper bound: " << bkb->get_valid_upper_bound() <<
        "   "  << "Optimal value: " << sol->get_var_value() << std::endl;

    for( int i = 0 ; i < 80 ; i++ )
        std::cout << "-";
    std::cout << "\n";

}