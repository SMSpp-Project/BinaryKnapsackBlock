/*--------------------------------------------------------------------------*/
/*------------------------------ File bkbench.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small in-process benchmark driver for a BinaryKnapsackBlock Solver over a
 * Pisinger/Jooken knapsack benchmark file (a .csv holding many instances, each
 * a block "knapPI_..._i / n N / c C / z Z / time T / <idx,profit,weight,x*> /
 * -----"). Every instance is solved with a fresh Block + Solver and the
 * optimum is cross-checked against the z recorded in the file. Running all
 * the instances in a single process avoids the heavy per-instance start-up of
 * the generic block_solver tool.
 *
 * usage: bkbench < file.csv | jooken-dir > [ <slow threshold s, default 0.5> ]
 *                [ <solver, default CoreDPBinaryKnapsackSolver> ]
 * A directory argument is treated as a Jooken benchmark root: every
 * problemInstances/<name>/test.in ("n / <id p w> lines / capacity") is solved
 * and cross-checked against the optimum recorded in <dir>/optima.csv; a
 * directory holding test.in itself is solved as a single instance (a
 * per-instance process unit, so a driver can parallelise and timeout each
 * instance independently). The solver is created through the Solver factory,
 * so any Solver that reads the BinaryKnapsackBlock physical representation
 * can be benchmarked.
 * output: one "SLOW ..." line per instance above the threshold, then a final
 *         "<class>: n=<#> ok=<#> mism=<#> tot=<s> max=<s> (<slowest>)" summary.
 *
 * \author Donato Meoli
 */
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "BinaryKnapsackBlock.h"
#include "CoreDPBinaryKnapsackSolver.h"

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

using Index = BinaryKnapsackBlock::Index;

// name (in the Solver factory) of the Solver to benchmark
static std::string solver_name = "CoreDPBinaryKnapsackSolver";

/*--------------------------------------------------------------------------*/

static std::string classname( const std::string & path )
{
 auto s = path.find_last_of( "/\\" );
 std::string f = ( s == std::string::npos ) ? path : path.substr( s + 1 );
 auto d = f.rfind( ".csv" );
 return( d == std::string::npos ? f : f.substr( 0 , d ) );
 }

/*--------------------------------------------------------------------------*/

// solve one instance with a fresh Block + Solver (load() on a reused Block
// does not re-trigger the Solver's extraction) and return the optimum

static double solve_one( Index N , double C , std::vector< double > && W ,
                         std::vector< double > && P , double & dt )
{
 auto BKB = new BinaryKnapsackBlock();
 BKB->load( N , C , std::move( W ) , std::move( P ) );
 auto slv = Solver::new_Solver( solver_name );
 slv->set_Block( BKB );

 const auto t0 = std::chrono::steady_clock::now();
 slv->compute();
 const double val = slv->get_var_value();
 const auto t1 = std::chrono::steady_clock::now();
 dt = std::chrono::duration< double >( t1 - t0 ).count();

 delete slv;
 delete BKB;
 return( val );
 }

/*--------------------------------------------------------------------------*/

// read a Jooken test.in ("n / <id p w> lines / capacity")

static bool read_jooken( const std::filesystem::path & file , Index & N ,
                         std::vector< double > & W ,
                         std::vector< double > & P , double & C )
{
 std::ifstream tin( file );
 if( ! tin )
  return( false );

 tin >> N;
 W.resize( N ); P.resize( N );
 for( Index i = 0 ; i < N ; ++i ) {
  long id; double p, w;
  tin >> id >> p >> w;
  P[ i ] = p; W[ i ] = w;
  }
 tin >> C;
 return( ( bool ) tin );
 }

/*--------------------------------------------------------------------------*/

// the optimum recorded for the named instance in <root>/optima.csv
// ("name,optimum" lines); < 0 (the benchmark's sentinel for "unknown": combo
// ran out of time / memory) also when the instance is not found

static double jooken_optimum( const std::string & root ,
                              const std::string & name )
{
 std::ifstream oin( root + "/optima.csv" );
 std::string line;
 while( std::getline( oin , line ) ) {
  const auto c = line.find( ',' );
  if( ( c != std::string::npos ) && ( line.compare( 0 , c , name ) == 0 ) )
   return( std::stod( line.substr( c + 1 ) ) );
  }
 return( -1 );
 }

/*--------------------------------------------------------------------------*/

// Jooken benchmark root: solve every problemInstances/<name>/test.in and
// cross-check against the optimum in <dir>/optima.csv

static int run_jooken( const std::string & dir , double slow )
{
 namespace fs = std::filesystem;

 long n = 0, mism = 0, skip = 0;
 double tot = 0, mx = 0;
 std::string slowest;

 std::vector< fs::path > dirs;
 for( const auto & e : fs::directory_iterator( dir + "/problemInstances" ) )
  if( e.is_directory() )
   dirs.push_back( e.path() );
 std::sort( dirs.begin() , dirs.end() );

 for( const auto & d : dirs ) {
  const std::string name = d.filename().string();
  const double z = jooken_optimum( dir , name );

  Index N = 0;
  double C = 0;
  std::vector< double > W, P;
  if( ( z < 0 ) || ( ! read_jooken( d / "test.in" , N , W , P , C ) ) ) {
   ++skip;
   continue;
   }

  double dt;
  const double val = solve_one( N , C , std::move( W ) , std::move( P ) , dt );

  ++n; tot += dt;
  if( dt > mx ) { mx = dt; slowest = name; }
  if( std::abs( val - z ) > 1e-6 * std::max( 1.0 , std::abs( z ) ) ) {
   ++mism;
   std::cout << "MISMATCH " << name << " val=" << val << " z=" << z << "\n";
   }
  if( dt > slow ) {
   std::cout << "SLOW " << name << " n=" << N << " t=" << dt << "s\n";
   std::cout.flush();
   }
  }

 std::cout << "jooken: n=" << n << " ok=" << ( n - mism ) << " mism=" << mism
           << " skip=" << skip << " tot=" << tot << "s max=" << mx
           << "s (" << slowest << ")\n";
 return( mism ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 if( argc < 2 ) {
  std::cerr << "usage: " << argv[ 0 ]
            << " <file.csv|jooken-dir> [slow_threshold_s] [solver]\n";
  return( 1 );
  }
 const std::string path = argv[ 1 ];
 const double slow = ( argc > 2 ) ? std::stod( argv[ 2 ] ) : 0.5;
 if( argc > 3 )
  solver_name = argv[ 3 ];

 if( std::filesystem::is_directory( path ) ) {
  // a single Jooken instance dir (it holds test.in), with optima.csv two
  // levels up
  if( std::filesystem::exists( std::filesystem::path( path ) / "test.in" ) ) {
   const auto p = std::filesystem::absolute( path );
   const std::string name = p.filename().string();
   const std::string root = p.parent_path().parent_path().string();

   const double z = jooken_optimum( root , name );
   if( z < 0 ) {
    std::cout << name << ": SKIP (no known optimum)\n";
    return( 0 );
    }

   Index N = 0;
   double C = 0;
   std::vector< double > W, P;
   if( ! read_jooken( p / "test.in" , N , W , P , C ) ) {
    std::cerr << "cannot open " << path << "/test.in\n";
    return( 1 );
    }

   double dt;
   const double val = solve_one( N , C , std::move( W ) , std::move( P ) ,
                                 dt );
   const bool ok = std::abs( val - z ) <= 1e-6 * std::max( 1.0 ,
                                                           std::abs( z ) );
   std::cout << name << ": " << ( ok ? "ok" : "MISMATCH" ) << " val="
             << val << " z=" << z << " t=" << dt << "s\n";
   return( ok ? 0 : 1 );
   }
  return( run_jooken( path , slow ) );
  }

 std::ifstream in( path );
 if( ! in ) { std::cerr << "cannot open " << path << "\n"; return( 1 ); }

 long n = 0, mism = 0;
 double tot = 0, mx = 0;
 std::string slowest;

 std::string line, name;
 while( std::getline( in , line ) ) {
  // start of an instance block: a name line "knapPI_..."
  if( line.compare( 0 , 6 , "knapPI" ) != 0 )
   continue;
  name = line;
  // remove a trailing carriage return, if any
  if( ! name.empty() && name.back() == '\r' ) name.pop_back();

  Index N = 0;
  double C = 0, Z = 0;
  std::vector< double > W, P;

  // header lines: n / c / z / time (order as produced by the generator)
  std::getline( in , line ); std::sscanf( line.c_str() , "n %u" , &N );
  std::getline( in , line ); std::sscanf( line.c_str() , "c %lf" , &C );
  std::getline( in , line ); std::sscanf( line.c_str() , "z %lf" , &Z );
  std::getline( in , line ); // time ...

  W.resize( N ); P.resize( N );
  for( Index i = 0 ; i < N ; ++i ) {
   std::getline( in , line );
   long idx; double p, w;
   // "idx,profit,weight,xstar"
   std::sscanf( line.c_str() , "%ld,%lf,%lf" , &idx , &p , &w );
   P[ i ] = p; W[ i ] = w;
   }
  std::getline( in , line ); // "-----"

  double dt;
  const double val = solve_one( N , C , std::move( W ) , std::move( P ) , dt );

  ++n; tot += dt;
  if( dt > mx ) { mx = dt; slowest = name; }
  if( std::abs( val - Z ) > 1e-6 * std::max( 1.0 , std::abs( Z ) ) ) {
   ++mism;
   std::cout << "MISMATCH " << name << " val=" << val << " z=" << Z << "\n";
   }
  if( dt > slow ) {
   std::cout << "SLOW " << name << " n=" << N << " t=" << dt << "s\n";
   std::cout.flush();
   }
  }

 std::cout << classname( path ) << ": n=" << n << " ok=" << ( n - mism )
           << " mism=" << mism << " tot=" << tot << "s max=" << mx
           << "s (" << slowest << ")\n";

 return( mism ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File bkbench.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
