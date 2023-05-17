/*--------------------------------------------------------------------------*/
/*----------------------- File BinaryKnapsackBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class BinaryKnapsackBlock, which
 * implements the Block concept [see Block.h] for a Knapsack problem with
 * linear objective and both integer and continuous variables.
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
/*---------------------------- DEFINITIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BinaryKnapsackBlock
 #define __BinaryKnapsackBlock 
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "Change.h"

#include "FRowConstraint.h"

#include "FRealObjective.h"

#include "LinearFunction.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class BinaryKnapsackSolution;//forward declaration of BinaryKnapsacKSolution

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASSES -----------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------- CLASS BinaryKnapsackBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Binary Knapsack problem
/** The BinaryKnapsackBlock class implements the Block concept [see Block.h]
 * for the Binary Knapsack problem.
 *
 * The data of the problem consist of two vectors W and P of the same size 
 * n = |W| = |P|, where n is the number of items. W[ i ] is the weight of
 * item i and P[ i ] is the profit of item i. Finally C is the total capacity
 * of the Knapsack.
 *
 * A formulation of the problem is:
 * \f[
 *  \max \sum_{ i = 1 }^{ n } P[ i ] X[ i ]
 * \f]
 * \f[
 *  \sum_{ i = 1 }^{ n } W[ i ] X[ i ] \leq C      (1)
 * \f]
 * \f[
 *   X[ i ] \in \{ 0 , 1 \} \quad i = 1, \dots, n  (2)
 * \f]
 * 
 * By default it is a maximization problem, but it is possible to change the 
 * sense of the objective with the method set_objective_sense() that changes 
 * the value of f_sense.
 *
 *   - f_sense = true for a maximization problem
 *
 *   - f_sense = false for a minimization problem
 * 
 * The set of items is assumed not to be changed (save for changing weights and
 * profits, and for items to be fixed or unfixed). */

class BinaryKnapsackBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

 using boolVec = std::vector< bool >;
 using c_boolVec = const boolVec;
 using boolVec_it = boolVec::iterator;
 using c_boolVec_it = boolVec::const_iterator;

 using doubleVec = std::vector< double >;
 using c_doubleVec = const doubleVec ;
 using dblVec_it = doubleVec::iterator;
 using c_dblVec_it = doubleVec::const_iterator;

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend BinaryKnapsackSolution;  ///< make BinaryKnapsackSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of BinaryKnapsackBlock, taking a pointer to the father Block
 /** Constructor of BinaryKnapsackBlock. It accepts a pointer to the father 
  * Block, which can be of any type, defaulting to nullpt so that this can 
  * also be used as the void constructor. */

 explicit BinaryKnapsackBlock( Block * father = nullptr )
  : Block( father ) , f_C( 0 ) , AR( 0 ) , f_sense( true ) , 
    f_cond_lower( - Inf< double >() ) , f_cond_upper( + Inf< double >() ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of BinaryKnapsackBlock: deletes the abstract representation

 virtual ~BinaryKnapsackBlock() { guts_of_destructor(); }

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the Binary Knapsack instance from memory, copying
 /** Loads the Binary Knapsack instance from memory, copying the data from
  * the input parameters. The parameters are:
  *
  * - n         is the number of items   
  *
  * - Capacity  is the capacity of the Knapsack
  *
  * - Weights   is the vector of weights, which must have size equal to n
  *
  * - Profits   is the vector of profits, which must have size equal to n
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * BinaryKnapsackBlock then a NBModification (the "nuclear option") is 
  * issued. */

 void load( Index n , double Capacity , c_doubleVec & Weights , 
            c_doubleVec & Profits , c_boolVec & Integrality = {} );

/*--------------------------------------------------------------------------*/
 /// loads the Binary Knapsack instance from memory, moving
 /** Loads the Binary Knapsack instance from memory, moving the data from
  * the input parameters. The parameters are:
  *
  * - n         is the number of items   
  *
  * - Capacity  is the capacity of the Knapsack
  *
  * - Weights   is the vector of weights, which must have size equal to n
  *
  * - Profits   is the vector of profits, which must have size equal to n
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * BinaryKnapsackBlock then a NBModification (the "nuclear option") is 
  * issued. */

 void load( Index n , double Capacity , doubleVec && Weights , 
            doubleVec && Profits , boolVec && Integrality = {} );

/*--------------------------------------------------------------------------*/
 /// load instance from txt file  
 /** Loads a BinaryKnapsackBlock out of std::istream. The format is the
  * following, with each element being separated by whitespaces and possibly
  * comments:
  *
  * - number of items
  *
  * - capacity of the knapsack
  *
  * - for i = 1 to n: weight of item i
  *
  * - for i = 1 to n: profit of item i
  *
  * If the stream (after having extracted whitespaces and comments) does not
  * eof() here, then
  *
  * - for i = 1 to n: integrality of item i (true if integral, false if not)
  *
  * If integrality is not specified, true is assumed for all objects.
  *
  * Since there is only one supported input format, \p frmt is ignored.
  *
  * Like load( memory ), if there is any Solver attached to this 
  * BinaryKnapsackBlock then a NBModification (the "nuclear option") is 
  * issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a BinaryKnapsackBlock. Besides what is managed by serialize() method of 
  * the base Block class, the group should contains the following:
  *
  * - the dimension "NItems" containing the number of items
  *
  * - the dimension "Capacity" containing the capacity of the Knapsack
  *
  * - the variable "Weights", of type double and indexed over the dimension
  *   "NItems"; the i-th entry of the variable is assumed to contain the 
  *   weight of the i-th item;
  *
  * - the variable "Profits", of type double and indexed over the dimension
  *   "NItems"; the i-th entry of the variable is assumed to contain the 
  *   profit of the i-th item;
  * 
  * All dimensions and variables are mandatory. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the Binary Knapsack
 /** Method that generates the abstract Variable of the Binary Knapsack. 
  *
  * There is a binary variable for each item. An std::vector< ColVariable >
  * with exactly NItems entries is generated. */

 void generate_abstract_variables( Configuration * stvv = nullptr )
  override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the Binary Knapsack.
 /** Method that generates the abstract constraint of the Binary Knapsack. 
  *
  * There is only one FRowConstraint represented as a LinearFunction whose
  * coefficients are the weights of the items. The RHS of the constraint is
  * the total capacity C of the Knapsack. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the Binary Knapsack
 /** Method that generates the objective of the Binary Knapsack.
  *
  * The objective function is represented as a LinearFunction whose 
  * coefficients are the profits of the items. The sense of the objective is
  * maximization if f_sense == true or minimization if f_sense == false. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- Methods for reading the data of the BinaryKnapsackBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the BinaryKnapsackBlock
 *  @{ */

 /// getting the current sense of the Objective

 int get_objective_sense( void ) const override final { 
  return( f_sense ? Objective::eMax : Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// getting a valid upper bound on the value of the Objective
 /** In order to compute a valid upper bound, first note:
  * 
  * If it is a maximization problem, then 
  *
  * - all the items i with negative weight W[ i ] < 0 and positive profit 
  *   P[ i ] > 0 are contained in the optimal solution.
  * 
  * - all the items i with positive weight W[ i ] > 0 and negative profit 
  *   P[ i ] < 0 are not contained in the optimal solution.
  *
  * Similarly if it is a minimization problem changing the sign of the profits
  *
  * Once all these items have been preprocessed, a valid upper bound on the 
  * optimal value of the problem is computed as a sum of the positive profits 
  * of all the remaining items (plus the sum of the profits of the selected 
  * preprocessed items).
  *
  * If the objective sense is maximization this is also a global valid upper
  * bound, otherwise it is a conditionally valid upper bound, because the 
  * problem could be empty (and it being a minimization one this would mean 
  * that its optimal value is + infinity) */

 double get_valid_upper_bound( bool conditional = false ) 
  override final {
  if( ! f_sense  )  // if the sense is minimization
   if( ( ! conditional ) && is_empty() )
    return( Inf< double >() );

  if( std::isinf( f_cond_upper ) )
   compute_conditional_bounds();

  return( f_cond_upper );
  }

/*--------------------------------------------------------------------------*/
 /// getting a lower bound on the value of the Objective
 /** In order to compute a lower bound, first note:
  * 
  * If it is a maximization problem, then 
  *
  * - all the items i with negative weight W[ i ] < 0 and positive profit 
  *   P[ i ] > 0 are contained in the optimal solution.
  * 
  * - all the items i with positive weight W[ i ] > 0 and negative profit 
  *   P[ i ] < 0 are not contained in the optimal solution.
  *
  * Similarly if it is a minimization problem changing the sign of the profits
  *
  * Once all these items have been preprocessed, a lower bound on the optimal
  * value of the problem is computed as a sum of the negative profits of all 
  * the remaining items (plus the sum of the profits of the selected 
  * preprocessed items). 
  *
  * If the objective sense is minimization this is also a global valid lower
  * bound, otherwise it is a conditionally valid lower bound, because the 
  * problem may be empty (and it being a maximization one this would mean that
  * its optimal value is - infinity). */

 double get_valid_lower_bound( bool conditional = false ) 
  override final {
  if( f_sense )  // if the sense is maximization
  if( ( ! conditional ) && is_empty() )
   return( -Inf< double >() );

  if( std::isinf( f_cond_lower ) )
   compute_conditional_bounds();

  return( f_cond_lower );
  }

/*--------------------------------------------------------------------------*/
 /// get the number of items

 Index get_NItems( void ) const { return( v_P.size() ); }

/*--------------------------------------------------------------------------*/
 /// get the capacity of the Knapsack

 double get_Capacity( void ) const { return( f_C ); }

/*--------------------------------------------------------------------------*/
 /// get the number of continuous variables

 Index get_N_Cont_Items( void ) const { return( countCont ); }

/*--------------------------------------------------------------------------*/
 /// given an index i return the weight of item i 

 double get_Weight( Index i ) const { 
  if( i >= get_NItems() )
   throw( std::invalid_argument( "BinaryKnapsackBlock: invalid item" ) );
  return( v_W[ i ] ); 
  }

/*--------------------------------------------------------------------------*/
 /// given a range rng return the vector of weigths of items in that range

 void get_Weights( dblVec_it weights, Range rng ) const { 
  
  if( rng.second >= get_NItems() )
   throw( std::invalid_argument( "Range not valid" ) );

  for( Index i = rng.first; i < rng.second ; ++i )
   *(weights++) = v_W[ i ];
  }

/*--------------------------------------------------------------------------*/
 /// given a Subset of items nms return the vector of their weigths 

 void get_Weights( dblVec_it weights , Subset && nms ) const { 
  
  for( Index i : nms ){
   if( i >= get_NItems() )
    throw( std::invalid_argument( "Invalid index in nms" ) ); 
   *(weights++) = v_W[ i ];
   }
  }

/*--------------------------------------------------------------------------*/
 /// get the vector of Weights

 c_doubleVec & get_Weights( void ) const { return( v_W ); }

/*--------------------------------------------------------------------------*/
 /// given an index i return the integrality of item i 

 bool get_Integrality( Index i ) const { 
  if( i >= get_NItems() )
   throw( std::invalid_argument( "BinaryKnapsackBlock: invalid item" ) );
  return( v_I[ i ] );
  }

/*--------------------------------------------------------------------------*/
 /// given a range rng return the integrality vector of items in that range

 void get_Integrality( boolVec_it integrality ,  Range rng ) const { 
  
  if( rng.second >= get_NItems() )
   throw( std::invalid_argument( "Range not valid" ) );


  for( Index i = rng.first; i < rng.second ; ++i )
   *(integrality++) = v_I[ i ];
  }

/*--------------------------------------------------------------------------*/
 /// given a Subset of items nms return the vector of their integrality 

 void get_Integrality( boolVec_it integrality , Subset && nms ) const { 

  for( Index i : nms ){
   
   if( i >= get_NItems() )
    throw( std::invalid_argument( "Invalid index in nms" ) );

   *(integrality++) = v_I[ i ];
   }
  }

/*--------------------------------------------------------------------------*/
 /// get the vector of Integrality 

 c_boolVec & get_Integrality( void ) const { return( v_I ); }

/*--------------------------------------------------------------------------*/
 /// given an index i return the profit of item i 

 double get_Profit( Index i ) const { 
  if( i >= get_NItems() )
   throw( std::invalid_argument( "BinaryKnapsackBlock: invalid item" ) );
  return( v_P[ i ] ); 
  }

/*--------------------------------------------------------------------------*/
 /// given a range rng return the vector of profits of items in that range

 void get_Profits( dblVec_it profits , Range rng ) const { 
  
  if( rng.second >= get_NItems() )
   throw( std::invalid_argument( "Range not valid" ) );

  for( Index i = rng.first; i < rng.second ; ++i )
   *(profits++) = v_P[ i ];
  }

/*--------------------------------------------------------------------------*/
 /// given a Subset of items nms return the vector of their profits

 void get_Profits( dblVec_it profits , Subset && nms ) const { 

  for( Index i : nms ){
   
   if( i >= get_NItems() )
    throw( std::invalid_argument( "Invalid index in nms" ) );

   *(profits++) = v_P[ i ];
   }
  }

/*--------------------------------------------------------------------------*/
 /// get the vector of Profits

 c_doubleVec & get_Profits( void ) const { return( v_P ); }

/*--------------------------------------------------------------------------*/
 /// given an index get a pointer to the corresponding variable

 ColVariable * get_Var( Index i ) const {
  if( i >= v_x.size() )
   throw( std::invalid_argument( "BinaryKnapsackBlock: invalid variable" ) );
  return( const_cast< ColVariable * >( & v_x[ i ]  ) ); 
  }

/**@} ----------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is feasible
 /** Returns true if the solution encoded in the current value of the x 
  * Variables of the BinaryKnapsackBlock is feasible. This clearly requires 
  * the Variable of the BinaryKnapsackBlock to have been defined, i.e., that 
  * generate_abstract_variables() has been called prior to this method. */

 bool is_feasible( bool useabstract = false , 
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the Binary Knapsack problem is empty.
 /** Returns true if the Binary Knapsack problem is empty. 
  * If the Capacity C of the Knapsack is positive, then x = 0 is a feasible 
  * solution and the problem is not empty. Otherwise the only way the Binary 
  * Knapsack problem could be empty is when the capacity is negative and there 
  * are not enough items with negative weights whose total weights is <= C.
  * The method is_empty() checks this condition. */

 bool is_empty( bool useabstract = false ,
                Configuration * optc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the Binary Knapsack problem is unbounded.

 bool is_unbounded( bool useabstract = false ,
                    Configuration * fsbc = nullptr ) override {
  return( false );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for R3 Blocks
 *  @{ */

 /// gets an R3 Block of BinaryKnapsackBlock currently only the copy one
 /** Gets an R3 Block of the BinaryKnapsackBlock. The list of currently 
  * supported R3 Block is:
  *
  * - r3bc == nullptr: the copy (a BinaryKnapsackBlock identical to this)
  */

 Block * get_R3_Block( Configuration *r3bc = nullptr , Block * base = nullptr , 
                       Block * father = nullptr ) override;

 /*--------------------------------------------------------------------------*/
 /// maps back the solution from a copy BinaryKnapsackBlock to the current one
 /** Maps back the solution from a copy BinaryKnapsackBlock to the current one
  */

 void map_back_solution( Block *R3B , Configuration *r3bc = nullptr ,
                         Configuration *solc = nullptr ) override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps the solution of the current BinaryKnapsackBlock to a copy.
 /** Maps the solution of the current BinaryKnapsackBlock to a copy. */

 void map_forward_solution( Block *R3B , Configuration *r3bc = nullptr ,
                            Configuration *solc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

/// returns a BinaryKnapsackSolution representing the current solution of this
/// BinaryKnapsackBlock
 /** Returns a BinaryKnapsackSolution representing the current solution status
  * of this BinaryKnapsackBlock.
  *
  * Note that BinaryKnapsackBlock may not contain the required solution,
  * if the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception, unless emptys = true, in which case the
  * BinaryKnapsackSolution object is only prepped for getting a solution, but 
  * it is not really getting one now.
  *
  * Note that, although the method clearly returns a BinaryKnapsackSolution, 
  * formally the return type is Solution *. This is because it is not possible
  * to forward declare BinaryKnapsackSolution as a derived class from 
  * Solution, nor to define BinaryKnapsackSolution before BinaryKnapsackBlock
  * because the former uses some type information declared in the latter. */ 

 Solution * get_Solution( Configuration * solc = nullptr, 
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// returns the objective value of the current solution

 RealObjective::OFValue get_objective_value( void ) {
  if( ! ( AR & HasObj ) )  // the objective is not there
   return( Inf< RealObjective::OFValue >() );
  f_obj.compute();
  return( f_obj.value() );
  }

/*--------------------------------------------------------------------------*/
 /// given an index gets the solution of the corresponding item

 double get_x( Index i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// get a contiguous interval of the solution

 void get_x( dblVec_it xSol , Range rng = INFRange ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// get the solution for an arbitrary subset of items

 void get_x( dblVec_it xSol , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// return the dual value of the knapsack constraint
 /** Returns the dual value of the knapsack constraint; this only makes sense
  * if all the variables are continuous. */

 double get_dual( void ) const { return( f_cnst.get_dual() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set the dual value of the knapsack constraint
 /** Set the dual value of the knapsack constraint; this only makes sense
  * if all the variables are continuous. */

 void set_dual( RowConstraint::c_RHSValue d_value ) {  
  f_cnst.set_dual( d_value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index, set the solution of the corresponding item

 void set_x( Index i , double value );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set a contiguous interval of the solution

 void set_x( c_dblVec_it xSol , Range rng = INFRange );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set the solution for an arbitrary subset of items

 void set_x( c_dblVec_it xSol , c_Subset & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i return true if the corresponding variable is fixed

 bool is_fixed( Index i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// get a reference to v_fxd

 const std::vector< unsigned char > & get_fxd() const { return( v_fxd ); }

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// returns true if there is any Solver "listening to this Block"
 /** Returns true if there is any Solver "listening to this 
  * BinaryKnapsackBlock", or if the BinaryKnapsackBlock has to "listen" 
  * anyway because the "abstract" representation is constructed, and therefore
  * "abstract" Modification have to be generated anyway to keep the two 
  * representations in sync.
  *
  * No, this should not be needed. In fact, if the "abstract" representation
  * is modified with the default eModBlck value of issueMod, it is issued
  * irrespectively to the value of anyone_there(); see Observer::issue_mod().
  * If the value of issueMod is anything else the  "abstract" representation
  * has been modified already and there is no point in issuing the
  * Modification.
  * Note that that Observer::issue_mod() does not check if the "abstract"
  * representation has been constructed, but this is clearly not
  * necessary, as the Modification we are speaking of are issued while
  * changing the "abstract" representation, if that has not been
  * constructed then it cannot issue Modification

 bool anyone_there( void ) const override {
  return( AR ? true : Block::anyone_there() );
  }
 */
/*--------------------------------------------------------------------------*/ 
 /** Method for handling Modification. 
  *
  * The version of BinaryKnapsackBlock has to intercept any "abstract 
  * Modification" that modifies the "abstract representation" of the 
  * BinaryKnapsackBlock, and "translate" them into both changes of the actual 
  * data structures and corresponding "physical Modification". These 
  * Modification are those for which Modification::concerns_Block() is true. 
  * Note, however, that before sending the Modification to the Solver and/or 
  * the father Block, the concerns_Block() value is set to false. This is 
  * because once it is passed through this method, the "abstract Modification" 
  * has "already done its duty" of providing the information to the 
  * BinaryKnapsackBlock, and this must not be repeated.
  * 
  * The following "abstract Modification" are handled:
  *
  * - GroupModification, that are simply unpacked into the individual
  *   sub-[Group]Modification and dealt with individually;
  *
  * - C05FunctionModRngd and C05FunctionModSbst changing coefficients coming
  *   from the LinearFunction into the Objective or from the LinearFunction 
  *   into the Constraint;
  *
  * - RowConstraintMod changing the RHS of the Constraint;
  *
  * - VariableMod fixing and un-fixing a ColVariable;
  *
  * - ObjectiveMod changing the sense of the objective function. 
  *
  * Any other Modification reaching the BinaryKnapsackBlock will lead to 
  * exception being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR PRINTING & SAVING THE BinaryKnapsackBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the BinaryKnapsackBlock
 *  @{ */

 /// print the BinaryKnapsackBlock on an ostream
 /** Print information about the BinaryKnapsackBlock with the given verbosity
  * level. */ 

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * BinaryKnapsackBlock. See BinaryKnapsackBlock::deserialize(netCDF::NcGroup) 
 * for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the BinaryKnapsack instance
 *
 * All the methods in this section have two parameters issueMod and issueAMod
 * which control if and how the, respectively, "physical Modification" and
 * "abstract Modification" corresponding to the change have to be issued, and
 * where (to which channel). The format of the parameters is that of
 * Observer::make_par(), except that the value eModBlck is ignored and
 * treated it as if it were eNoBlck [see Observer::issue_pmod()]. This is
 * because it makes no sense to issue an "abstract" Modification with
 * concerns_Block() == true, since the changes in the BinaryKnapsackBlock have 
 * surely been done already, and this is just not possible for a "physical"
 * Modification. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index and a value, fix the corresponding variable
 /** Fixes the value of the variable \p i to 0 or 1 depending on the value of
  * \p value; however, if \p is fixed already then the value of \p i is not
  * changed, even if it differs from \p value. */

 void fix_x( bool value , Index i , ModParam issueMod = eNoBlck ,
                                    ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// fix variables of a contiguos interval of indices
 /** Fixes the value of the variables whose names are in the interval
  * specified by \p rng to the values specified by \p value; that is, the
  * variable rng.first + h takes value 1 if *( value + h ) == true and 0
  * otherwise. \p value is assumed to point to an array at least as long as
  * rng.second - rng.first, except if rng.second is larger than the number
  * of variables, in which case is taken to mean "up until the end", and
  * \p value need only be correspondingly long.
  * Note that if any of the variable is fixed already then its value is not
  * changed, even if it differs from what \p value would dictate. */

 void fix_x( c_boolVec_it value , Range rng = INFRange , 
             ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// fix variables of an arbitrary subsets of indices
 /** Fixes the value of the variables whose names are in the Subset \p nms
  * to the values specified by \p value; that is, the variable nms[ h ] takes
  * value 1 if *( value + h ) == true and 0 otherwise. \p value is assumed to
  * point to an array at least as long as nms.size(). \p nms is assumed to be
  * ordered in increasing sense if ordered == true, and of course it must
  * only contain names < get_NItems().
  * Note that if any of the variable is fixed already then its value is not
  * changed, even if it differs from what \p value would dictate. */

 void fix_x( c_boolVec_it value , Subset && nms , bool ordered = false ,
             ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index, unfixes the corresponding variable

 void unfix_x( Index i , ModParam issueMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// unfix variables of a contiguos interval of indices

 void unfix_x( Range rng = INFRange , 
               ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// unfix variables of an arbitrary subsets of indices

 void unfix_x( Subset && nms , bool ordered = false ,
               ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i change the weight of item i

 void chg_weight( double NWeight , Index item , 
                  ModParam issueMod = eNoBlck ,
                  ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the weights of a contiguous interval of items

 void chg_weights( c_dblVec_it NWeight , Range rng = INFRange , 
                   ModParam issueMod = eNoBlck ,
                   ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the weights of an arbitrary subsets of items

 void chg_weights( c_dblVec_it NWeight , 
                   Subset && nms , bool ordered = false ,  
                   ModParam issueMod = eNoBlck ,
                   ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i change the boolean vector that tell which variable is 
 /// continuous and which integer (integraliy vector)

 void chg_integrality( bool NIntegrality , Index item , 
                       ModParam issueMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the integrality vector of a contiguous interval of items

 void chg_integrality( c_boolVec_it NIntegrality , Range rng = INFRange , 
                       ModParam issueMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the weights of an arbitrary subsets of items

 void chg_integrality( c_boolVec_it NIntegrality , 
                       Subset && nms , bool ordered = false ,  
                       ModParam issueMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i change the profit of item i 

 void chg_profit( double NProfit , Index item , 
                  ModParam issueMod = eNoBlck ,
                  ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the profits of a contiguous interval of items

 void chg_profits( c_dblVec_it NProfit , Range rng = INFRange , 
                   ModParam issueMod = eNoBlck ,
                   ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the profits of an arbitrary subsets of items

 void chg_profits( c_dblVec_it NProfit , 
                   Subset && nms , bool ordered = false ,  
                   ModParam issueMod = eNoBlck ,
                   ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the capacity of the Knapsack

 void chg_capacity( double NC , ModParam issueMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck ); 

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// set the sense of the objective function

 void set_objective_sense( bool sense , ModParam issueMod = eNoBlck ,
                                        ModParam issueAMod = eNoBlck ); 

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 double f_C;           ///< the Capacity of the Knapsack
 doubleVec v_W;        ///< vector of Weights          
 doubleVec v_P;        ///< vector of Profits
 boolVec v_I;          ///< vector of boolean for what variable are discrete
 Index countCont;      ///< counter for the umber of continuous variables

 unsigned char AR;     ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variable have been constructed
 static constexpr unsigned char HasObj = 2;
 ///< second bit of AR == 1 if the Objective has been constructed
 static constexpr unsigned char HasCns = 4;
 ///< third bit of AR == 1 if the Constraint has been constructed

 bool f_sense;         ///< the sense of the objective
                       /** f_sense = true for maximization
                        *  f_sense = false for minimization */

 double f_cond_lower;            ///< conditional lower bound, can be +INF
 double f_cond_upper;            ///< conditional upper bound, can be -INF

 std::vector< unsigned char > v_fxd;  ///< vector saying how the x are fixed
                                      /* < v_fxd[ i ] indicates if x_i is
                                       * fixed, with the following encoding:
                                       * 0 = not fixed , 
                                       * 1 = fixed to 0 , 
                                       * 2 = fixed to 1                     */


 std::vector< ColVariable > v_x; ///< the static binary/continuous variables
 FRowConstraint f_cnst;          ///< the static constraint 
 FRealObjective f_obj;           ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

 void guts_of_add_Modification( c_p_Mod mod , ChnlName chnl );

 void compute_conditional_bounds( void );

 int p2i_x( const Variable * var ) const {
  return( std::distance( v_x.data() ,
       static_cast< const ColVariable * >( var ) ) ); 
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h; // insert BinaryKnapsackBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }; // end( class( BinaryKnapsackBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS BinayKnapsackBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a BinaryKnapsackBlock

class BinaryKnapsackBlockMod : public Modification {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of BinaryKnapsackBlockMod
 
 enum BinaryKnapsackBlock_mod_type {
  eChgWeight = 0  ,   ///< change the item weight
  eChgProfit      ,   ///< change the item profit
  eChgIntegrality ,   ///< change the variable integrality constraint
  eChgCapacity    ,   ///< change the Knapsack capacity
  eFixX           ,   ///< fix a variable x
  eUnfixX         ,   ///< unfix a variable x
  eChgSense           ///< change the sense of the objective
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the BinaryKnapsackBlock and the type

 BinaryKnapsackBlockMod( BinaryKnapsackBlock * const fblock , int type )
                       : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~BinaryKnapsackBlockMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// returns the [BinaryKnapsack]Block to which the BinaryKnapsackBlockMod 
 /// refers

 Block * get_Block( void ) const override { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of modification

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the BinaryKnapsackBlockMod

 void print( std::ostream &output ) const override {
  output << "BinaryKnapsackBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eChgWeight ):  output << "change weight "; break;
   case( eChgProfit ):  output << "change profit "; break;
   case( eChgIntegrality ): output << "change integrality "; break;
   case( eChgCapacity ):  output << "change capacity "; break;
   case( eFixX ):  output << "fix x "; break;
   case( eUnfixX ):  output << "unfix x "; break;
   case( eChgSense ):  output << "change objective sense "; break;
   }
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 BinaryKnapsackBlock * f_Block;
  ///< pointer to the BinaryKnapsackBlock to which the modification refers

 int f_type;   ///< type of modification

/*--------------------------------------------------------------------------*/

 };  // end( class( BinaryKnapsackBlockMod ) )

/*--------------------------------------------------------------------------*/
/*----------------- CLASS BinaryKnapsackBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from BinaryKnapsackBlockMod for "ranged" modifications
/** Derived class from BinaryKnapsackBlockMod to describe "ranged" 
 * modifications to a BinaryKnapsackBlock, i.e., modifications that apply to 
 * an interval of items. */

class BinaryKnapsackBlockRngdMod : public BinaryKnapsackBlockMod
{

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the BinaryKnapsackBlock, the type, and the range

 BinaryKnapsackBlockRngdMod( BinaryKnapsackBlock * const fblock , 
           int type , Block::Range rng )
  : BinaryKnapsackBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~BinaryKnapsackBlockRngdMod() = default;   ///< destructor

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return( f_rng ); }
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the BinaryKnapsackBlockRngdMod

 void print( std::ostream &output ) const override {
  BinaryKnapsackBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;     ///< the range

/*--------------------------------------------------------------------------*/

 };  // end( class( BinaryKnapsackBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS BinaryKnapsackBlockSbstMod --------------------*/
/*--------------------------------------------------------------------------*/
/// derived from BinaryKnapsackBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 * BinaryKnapsackBlock, i.e., modifications that apply to an arbitrary subset 
 * of items. */

class BinaryKnapsackBlockSbstMod : public BinaryKnapsackBlockMod {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the BinaryKnapsackBlock, the type, and the subset
 /**< Constructor: takes the BinaryKnapsackBlock, the type, and the subset. 
  * As the the && tells, nms is "consumed" by the constructor and its resources 
  * become property of the BinaryKnapsackBlockSbstMod object. */

 BinaryKnapsackBlockSbstMod( BinaryKnapsackBlock * const fblock , 
           int type , Block::Subset && nms )
  : BinaryKnapsackBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~BinaryKnapsackBlockSbstMod() = default;  ///< destructor

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the BinaryKnapsackBlockSbstMod

 void print( std::ostream &output ) const override {
  BinaryKnapsackBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms;   ///< the subset

/*--------------------------------------------------------------------------*/

 };  // end( class( BinaryKnapsackBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS BinaryKnapsackSolution -----------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a BinaryKnapsackBlock

class BinaryKnapsackSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend BinaryKnapsackBlock;        ///< make BinaryKnapsackBlock friend

/*----------- CONSTRUCTING AND DESTRUCTING BinaryKnapsackSolution ----------*/

 explicit BinaryKnapsackSolution() {}  ///< constructor

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~BinaryKnapsackSolution() = default;  ///< destructor

/*------ METHODS DESCRIBING THE BEHAVIOR OF A BinaryKnapsackSolution ------*/

 void read( const Block * const block ) override final;

 void write( Block * const block ) override final;

 void serialize( netCDF::NcGroup & group ) const override final;

 void print( std::ostream & output );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 BinaryKnapsackSolution * scale( double factor ) const override final;

 void sum( const Solution * solution , double multiplier ) override final;

 BinaryKnapsackSolution * clone( bool empty = false ) const override final;

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/
/// vector containing the values of the current solution. It contains doubles
/// to allow to store scaled solutions of some double factor

 std::vector< double > v_x;      ///< the variables

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( BinaryKnapsackSolution ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS BinaryKnapsackBlockChange ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Change for a BinaryKnapsackBlock

class BinaryKnapsackBlockChange : public Change {

public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of BinaryKnapsackBlockChange
 
 enum BinaryKnapsackBlock_chg_type {
  eChgWeight = 0  ,   ///< change the item weight
  eChgProfit      ,   ///< change the item profit
  eChgIntegrality ,   ///< change the variable integrality constraint
  eChgCapacity    ,   ///< change the Knapsack capacity
  eFixX           ,   ///< fix a variable x
  eUnfixX         ,   ///< unfix a variable x
  eChgSense       ,   ///< change the sense of the objective
  eEmpty
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/
 
 BinaryKnapsackBlockChange( int type = eEmpty ,
                            std::vector< double > && data = {} ) : 
                            f_type( type ) , v_data( data ) {}    

 ~BinaryKnapsackBlockChange() = default;            

/*-------------------- PUBLIC METHODS OF THE CLASS -------------------------*/

 void deserialize( const netCDF::NcGroup & group ) override {
  
  // read f_type
  auto ftype = group.getAtt( "BinaryKnapsackBlock_chg_type" );
  if( ftype.isNull() )
   throw( std::invalid_argument( 
         "Error reading BinaryKnapsackBlock_chg_type" ) );
  ftype.getValues( &f_type );

  // read v_data
  netCDF::NcDim ni = group.getDim( "dim" );
  netCDF::NcVar data = group.getVar( "Data" );
  if( data.isNull() )
   v_data.clear();
  else {
   v_data.resize( ni.getSize() );
   data.getVar( v_data.data() ); 
   } 

  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void serialize( netCDF::NcGroup & group ) const override {
  
  // always call the method of the base class first
  Change::serialize( group );

  group.putAtt( "BinaryKnapsackBlock_chg_type" , netCDF::NcInt() ,  f_type );
  
  netCDF::NcDim ni = group.addDim( "dim" , v_data.size() );
  ( group.addVar( "Data" , 
                  netCDF::NcDouble() , ni ) ).putVar( v_data.data() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 Change * apply( Block * block , 
                 bool doUndo = false , 
                 ModParam issueMod = eNoBlck , 
                 ModParam issueAMod = eNoBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void load( int type , const std::vector< double > & data = {} ) {

  if( type == eEmpty ) {
   f_type = type;         // update type
   v_data = {};           // discard data (if provided)
   return;
   }

  if( ( type == eChgSense ) || ( type == eChgCapacity ) ) {
   
   // check that data contains at least 1 entry (new sense or new capacity)
   if( ! data.size() )
    throw( std::invalid_argument( "Data not provided" ) );
   
   // update
   f_type = type;
   v_data = data;
   }
  else
   throw( std::invalid_argument( "type not supported" ) );  
  
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of change

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the BinaryKnapsackBlockChange

 void print( std::ostream &output ) const override {
  output << "BinaryKnapsackBlockChange[" << this << "]: ";
  switch( f_type ) {
   case( eChgWeight ):  output << "change weight "; break;
   case( eChgProfit ):  output << "change profit "; break;
   case( eChgIntegrality ): output << "change integrality "; break;
   case( eChgCapacity ):  output << "change capacity "; break;
   case( eFixX ):  output << "fix x "; break;
   case( eUnfixX ):  output << "unfix x "; break;
   case( eChgSense ):  output << "change objective sense "; break;
   case( eEmpty ): output << "Empty Change "; break;
   }
  
  output << std::endl << "New data: ";
  for( auto x : v_data )
   output << x << " ";
  output << std::endl;
  
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 int f_type;                      ///< type of change

 std::vector< double > v_data;    ///< data to change

private:

 SMSpp_insert_in_factory_h;

}; // end( class( BinaryKnapsackBlockChange ) )

/*--------------------------------------------------------------------------*/
/*----------------- CLASS BinaryKnapsackBlockRngdChange --------------------*/
/*--------------------------------------------------------------------------*/

class BinaryKnapsackBlockRngdChange : public BinaryKnapsackBlockChange {

public:

 static constexpr auto INFRange = BinaryKnapsackBlock::INFRange;

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/
 
 BinaryKnapsackBlockRngdChange( int type = eEmpty ,
                                std::vector< double > && data = {} ,
                                Block::Range rng = INFRange ) : 
                                BinaryKnapsackBlockChange( type , 
                                std::move( data ) ) , f_rng( rng ) {}  

 ~BinaryKnapsackBlockRngdChange() = default;            

/*-------------------- PUBLIC METHODS OF THE CLASS -------------------------*/

 void deserialize( const netCDF::NcGroup & group ) override {
  // read f_rng
  netCDF::NcDim ni = group.getDim( "RangeDim" );
  netCDF::NcVar rng = group.getVar( "Range" );
  rng.getVar( &f_rng );
  BinaryKnapsackBlockChange::deserialize( group ); 
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void serialize( netCDF::NcGroup & group ) const override {
  BinaryKnapsackBlockChange::serialize( group );
  netCDF::NcDim ni = group.addDim( "RangeDim" , 2 );
  ( group.addVar( "Range" , netCDF::NcInt() , ni ) ).putVar( &f_rng );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void load( int type , const std::vector< double > & data = {} , 
            Block::Range rng = INFRange ) {

  if( type == eEmpty ) {
   f_type = type;       // update type
   v_data = {};         // discard data (if provided)
   f_rng = INFRange;    // discard range (if provided)
   return;
   }
  
  if( ( type == eChgWeight ) || ( type == eChgProfit ) || 
      ( type == eChgIntegrality ) || ( type == eFixX ) || 
      ( type == eUnfixX ) ) {
   
   // check data size and range 
   if( data.size() != rng.second - rng.first )
    throw( std::invalid_argument( "data.size() != range size" ) );

   // update all the data
   f_type = type;
   v_data = data;
   f_rng = rng;
   }
  else
   throw( std::invalid_argument( "type not supported" ) );

  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 Change * apply( Block * block , 
                 bool doUndo = false , 
                 ModParam issueMod = eNoBlck , 
                 ModParam issueAMod = eNoBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to f_rng

 Block::Range rng( void ) const { return( f_rng ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the BinaryKnapsackBlockRngdChange

 void print( std::ostream &output ) const override {
  BinaryKnapsackBlockChange::print( output );
  output << "Range: ";
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;                      ///< type of change

private:

 SMSpp_insert_in_factory_h;

}; // end( class( BinaryKnapsackBlockRngdChange ) )

/*--------------------------------------------------------------------------*/
/*----------------- CLASS BinaryKnapsackBlockSbstChange --------------------*/
/*--------------------------------------------------------------------------*/

class BinaryKnapsackBlockSbstChange : public BinaryKnapsackBlockChange {

public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/
 
 BinaryKnapsackBlockSbstChange( int type = eEmpty ,
                                std::vector< double > && data = {} ,
                                Block::Subset && nms = {} ) : 
                                BinaryKnapsackBlockChange( type , 
                                std::move( data ) ) , 
                                v_nms( std::move( nms ) ) {}

 ~BinaryKnapsackBlockSbstChange() = default;            

/*-------------------- PUBLIC METHODS OF THE CLASS -------------------------*/

 void deserialize( const netCDF::NcGroup & group ) override {
  // read v_nms
  netCDF::NcDim ni = group.getDim( "SubsetDim" );
  netCDF::NcVar nms = group.getVar( "Subset" );
  nms.getVar( &v_nms ); 
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void serialize( netCDF::NcGroup & group ) const override {
  BinaryKnapsackBlockChange::serialize( group );
  netCDF::NcDim ni = group.addDim( "SubsetDim" , v_nms.size() );
  ( group.addVar( "Subset" , netCDF::NcInt() , ni ) ).putVar( &v_nms );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 void load( int type , const std::vector< double > & data = {} , 
            Block::Subset && nms = {} ) {

  if( type == eEmpty ) {
   f_type = type;       // update type
   v_data = {};         // discard data (if provided)
   v_nms = {};          // discard nms (if provided)
   return;
   }
  
  if( ( type == eChgWeight ) || ( type == eChgProfit ) || 
      ( type == eChgIntegrality ) || ( type == eFixX ) || 
      ( type == eUnfixX ) ) {
   
   // check data size and subset size 
   if( data.size() != nms.size() )
    throw( std::invalid_argument( "data.size() != nms.size()" ) );

   // update all the data
   f_type = type;
   v_data = data;
   v_nms = std::move( nms );
   }
  else
   throw( std::invalid_argument( "type not supported" ) );

  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 Change * apply( Block * block , 
                 bool doUndo = false , 
                 ModParam issueMod = eNoBlck , 
                 ModParam issueAMod = eNoBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to f_nms

 Block::c_Subset & nms( void ) const { return( v_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the BinaryKnapsackBlockSbstChange

 void print( std::ostream &output ) const override {
  BinaryKnapsackBlockChange::print( output );
  output << "Subset: ";
  output << "(# " << v_nms.size() << ")" << std::endl;
  for( auto i : v_nms )
   output << i << " ";
  output << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset v_nms;                      ///< subset of items to change

private:

 SMSpp_insert_in_factory_h;

}; // end( class( BinaryKnapsackBlockSbstChange ) )

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif /* BinaryKnapsackBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File BinaryKnapsackBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
