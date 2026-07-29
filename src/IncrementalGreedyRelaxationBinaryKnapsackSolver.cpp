/*--------------------------------------------------------------------------*/
/*--------------------- File IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class IncrementalGreedyRelaxationBinaryKnapsackSolver,
 *  which
 * implements the RelaxationSolver concept [see RelaxationSolver.h] and the
 * Solver concept [see Solver.h] for solving the continuous relaxation of a
 * Knapsack problem as represented by a BinaryKnapsackBlock.
 *
 * \author 	Filippo Magi
 * 	   		Dipartimento di Informatica \n
 * 			Universita' di Pisa \n
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * Copyright &copy by Filippo Magi, Federica Di Pasquale, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "IncrementalGreedyRelaxationBinaryKnapsackSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using idx_type = ThinComputeInterface::idx_type;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register IncrementalGreedyRelaxationBinaryKnapsackSolver  and GreedyHeuristicSolver to the factory

SMSpp_insert_in_factory_cpp_1(IncrementalGreedyChangeBinaryKnapsackSolver);
SMSpp_insert_in_factory_cpp_1(IncrementalGreedyRelaxationBinaryKnapsackSolver);

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF IncrementalGreedyRelaxationBinaryKnapsackSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::set_Block(Block *block)
{

	if (block == f_Block) // nothing to do
		return;

	Solver::set_Block(block); // attach to the new Block

	load(); // load Binary Knapsack instance
}

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int IncrementalGreedyChangeBinaryKnapsackSolver::compute(bool changedvars)
{
	lock();

	process_outstanding_Modification();

	if (changedData)
	{
		obj = -Inf<double>();
		initializeVariables();
		C = presolveCapacity;
	}
	changedData = false;
	// preprocessing

	// C,P and startingSearchIndex updated during apply

	// presolve capacity consider
	if (presolveCapacity < 0)
	{
		unlock();
		return (kInfeasible);
	}
	// solve continuous knapsack

	// f_ci = sortedVar.back(); // index of the critical item
	f_ciVal = 1; // variable value of the critical item

	reachTheEnd = true;
	for (Index j = startingSearchIndex; j < f_N; ++j)
	{
		Index i = sortedVar[j];
		if (skip[i])
			continue;
		if (C - v_W[i] < 0)
		{
			// goneBack = varToSorted[f_ci] > j;
			f_ci = i;
			reachTheEnd = false;
			f_ciVal = C / v_W[i];
			startingSearchIndex = j;
			break;
		}
		C -= v_W[i];
		P += v_P[i];
	}

	f_ci = reachTheEnd ? sortedVar.back() : f_ci;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	obj = P + (f_ciVal == 0 || f_ciVal == 1 ? 0 : v_P[f_ci] * f_ciVal); // update objective value

	// std::cout << "Obj " << obj << " true lb: " << get_true_lb() << std::endl;

Return_OK:

	unlock(); // unlock the mutex

	return (kOK);
}

/*--------------------------------------------------------------------------*/

// in case we have multiple items, we consider in last position the profit of the father solution and, in the second to last position the capacity of the father solution, to be used in the branching change
// gives the value for wthich f_i is fixed to, the residual capacity after the fix and the profit of that node considering the fix
std::vector<Change *> IncrementalGreedyRelaxationBinaryKnapsackSolver::branch()
{

	// branch on the critical item
	if (presolveCapacity - v_W[f_ci] < 0)
	{
		std::vector<Change *> branches(1);
		std::vector<double> zero = {0.0, C, P};
		branches[0] = new BinaryKnapsackBlockRngdChange(
			BinaryKnapsackBlockChange::eFixX,
			std::move(zero),
			std::make_pair(f_ci, f_ci + 1));
		return branches;
	}
	std::vector<Change *> branches(2);
	std::vector<double> zero;
	std::vector<double> one;

	zero = {0.0, C, P};
	one = {1.0, C - v_W[f_ci], P + v_P[f_ci]};

	branches[0] = new BinaryKnapsackBlockRngdChange(
		BinaryKnapsackBlockChange::eFixX,
		std::move(one),
		std::make_pair(f_ci, f_ci + 1));

	branches[1] = new BinaryKnapsackBlockRngdChange(
		BinaryKnapsackBlockChange::eFixX,
		std::move(zero),
		std::make_pair(f_ci, f_ci + 1));

	return branches;
}

/*--------------------------------------------------------------------------*/
Change *IncrementalGreedyChangeBinaryKnapsackSolver::apply(Change *chg, bool doUndo)
{
	auto CHG = dynamic_cast<BinaryKnapsackBlockChange *>(chg);

	if (!CHG)
		throw(std::invalid_argument("Change must be a BinaryKnapsackBlockChange"));

	if (CHG->type() != BinaryKnapsackBlockChange::eFixX && CHG->type() != BinaryKnapsackBlockChange::eUnfixX)
		return CHG->apply(f_Block, doUndo);
	// subset change (from external fixing)
	else if (BinaryKnapsackBlockSbstChange *change = dynamic_cast<BinaryKnapsackBlockSbstChange *>(CHG))
	{
		Block::Subset subset = change->nms();
		// TODO valutare se è giusto usare la move
		std::vector<double> data = std::move(change->data());
		Change *undo = nullptr;
		// branching case
		if (subset.size() == data.size() - 2)
		{
			// in last position there is the profit in case we branch
			P = data.back();
			data.pop_back();
			// in second to last position there is the residual capacity in case we branch
			C = data.back();
			data.pop_back();
		}
		// create the undo change
		if (doUndo)
		{
			// don't insert information about C and P, since there will be a fix that will give that values before compute
			// std::vector<double> inversed_data(subset.size());
			// TODO capire perché faccio inversed data
			/* 			for (Index i = 0; i < subset.size(); ++i)
						{
							inversed_data[i] = (1 - data[i]);
						} */
			undo = new BinaryKnapsackBlockSbstChange(
				(change->type() == BinaryKnapsackBlockChange::eFixX) ? BinaryKnapsackBlockChange::eUnfixX : BinaryKnapsackBlockChange::eFixX,
				std::move(std::vector<double>(data)),
				std::move(Block::Subset(subset)));
		}
		// unfixing
		if (change->type() == BinaryKnapsackBlockChange::eUnfixX)
		{
			for (auto i : subset)
			{
				skip[i] = false;
				presolveCapacity += data[i] * v_W[i];
				startingSearchIndex = startingSearchIndex > varToSorted[i] ? varToSorted[i] : startingSearchIndex;
			}
			return undo;
		}
		// fixing
		else if (change->type() == BinaryKnapsackBlockChange::eFixX)
		{
			Index idx = 0;
			for (auto i : subset)
			{
				f_ci = i;
				v_x[i] = complemented[i] ? 1 - (data[idx] == 1) : (data[idx] == 1);
				presolveCapacity -= data[idx] * v_W[i];
				skip[i] = true;
				if (data[idx] == 1 && varToSorted[f_ci] == startingSearchIndex)
					for (int j = startingSearchIndex; j >= 0; j--)
					{
						if (C >= 0)
						{
							startingSearchIndex = j;
							break;
						}
						Index z = sortedVar[j];
						if (skip[z])
							continue;
						P -= v_P[z];
						C += v_W[z];
					}
				// TODO capire se e come è necessario trattare il caso fisso a 1  AND f_ci != startingSearchIndex (quindi ho avuto unfix)
				else if (data[idx] == 1 && startingSearchIndex > varToSorted[f_ci])
				{
					for (Index j = varToSorted[f_ci]; j < startingSearchIndex; ++j)
					{
						Index z = sortedVar[j];
						if (skip[z])
							continue;
						P -= v_P[z];
						C += v_W[z];
					}
					P -= data[startingSearchIndex] * v_P[i];
					C += data[startingSearchIndex] * v_W[i];
					startingSearchIndex = varToSorted[i];
				}
				else if (data[idx] == 0)
				{
					if (startingSearchIndex > varToSorted[i])
					{

						for (Index j = varToSorted[i]; j < startingSearchIndex; ++j)
						{
							Index z = sortedVar[j];
							if (skip[z])
								continue;
							P -= v_P[z];
							C += v_W[z];
						}
						P -= data[startingSearchIndex] * v_P[i];
						C += data[startingSearchIndex] * v_W[i];
						startingSearchIndex = varToSorted[i];
					}
				}
				else
					throw(std::invalid_argument("Data for fixing in BinaryKnapsackBlockChange must be binary"));
				++idx;
			}
			return undo;
		}
		else
		{
			throw(std::invalid_argument("Change type not supported"));
		}
	}
	// ranged change with double (from external fixing)
	else if (BinaryKnapsackBlockRngdChange *change =
				 dynamic_cast<BinaryKnapsackBlockRngdChange *>(CHG))
	{
		Block::Range rng = change->rng();
		std::vector<double> data = change->data();
		Change *undo = nullptr;
		// branching case
		if (rng.second - rng.first == data.size() - 2)
		{

			// in last position there is the profit in case we branch
			P = data.back();
			data.pop_back();
			// in second to last position there is the residual capacity in case we branch
			C = data.back();
			data.pop_back();
		}
		// create the undo change
		if (doUndo)
		{
			std::vector<double> inversed_data(rng.second - rng.first);
			undo = new BinaryKnapsackBlockRngdChange(
				(change->type() == BinaryKnapsackBlockChange::eFixX) ? BinaryKnapsackBlockChange::eUnfixX : BinaryKnapsackBlockChange::eFixX,
				std::move(std::vector<double>(data)),
				std::move(Block::Range(rng)));
		}

		// unfixing
		if (change->type() == BinaryKnapsackBlockChange::eUnfixX)
		{
			for (Index i = rng.first; i < rng.second; ++i)
			{
				skip[i] = false;
				presolveCapacity += data[i - rng.first] * v_W[i];
				// update with minimum starting index
				startingSearchIndex = startingSearchIndex > varToSorted[i] ? varToSorted[i] : startingSearchIndex;
			}
			return undo;
		}
		// fixing
		else if (change->type() == BinaryKnapsackBlockChange::eFixX)
		{
			std::vector<double> data = change->data();
			// obtain data for change
			for (Index i = rng.first; i < rng.second; ++i)
			{
				v_x[i] = complemented[i] ? 1 - (data[i - rng.first] == 1) : (data[i - rng.first] == 1);
				presolveCapacity -= data[i - rng.first] * v_W[i];
				skip[i] = true;
				f_ci = i;

				//   fixing to 1
				// Main idea: if I fix to 1 after and I don't have a enough capacity, I come back to search the first element that permit to have a feasible region, and then I restart from that element my new search

				//   TODO capire se è possibile gone back in questo caso (probabilmente no)
				if (data[i - rng.first] == 1 && /* varToSorted[i] == startingSearchIndex  &&*/ C < 0)
				{
					for (int j = varToSorted[i]; j >= 0; --j)
					{
						Index z = sortedVar[j];
						if (skip[z])
							continue;
						P -= v_P[z];
						C += v_W[z];
						if (C >= 0)
						{
							startingSearchIndex = j;
							break;
						}
					}
				}
				// TODO capire se e come è necessario trattare il caso fisso a 1  AND f_ci != startingSearchIndex (quindi ho avuto unfix)
				// Main idea: in this case I fix after at least one unfix, then I have to update correctly my node
				else if (data[i - rng.first] == 1 /* && startingSearchIndex > varToSorted[i] */)
				{
					if (startingSearchIndex > varToSorted[i])
					{
						/* 						for (Index j = varToSorted[i]; j <= startingSearchIndex; ++j)
												{
													Index z = sortedVar[j];
													if (skip[z])
														continue;
													P -= v_P[z];
													C += v_W[z];
												} */
						startingSearchIndex = varToSorted[i];
					}
					else
						for (Index j = startingSearchIndex; j <= varToSorted[i]; ++j)
						{
							Index z = sortedVar[j];
							if (skip[z])
								continue;
							P -= v_P[z];
							C += v_W[z];
						}
				}
				// per capire cosa succede per davvore, dovrebbe venire coperto alla fine
				/* 				else if (data[i - rng.first] == 1 && startingSearchIndex < varToSorted[i])
								{
									for (Index j = startingSearchIndex; j <= varToSorted[i]; ++j)
									{
										Index z = sortedVar[j];
										if (skip[z])
											continue;
										P -= v_P[z];
										C += v_W[z];
									}
								} */

				//  fixing to zero
				else if (data[i - rng.first] == 0)
				{
					// I move back
					if (startingSearchIndex > varToSorted[i])
					{
						// TODO undestand why I change the value of P and C accordly (actually commented)
						/* 						for (Index j = varToSorted[i]; j <= startingSearchIndex; ++j)
												{
													Index z = sortedVar[j];
													if (skip[z])
														continue;
													P -= v_P[z];
													C += v_W[z];
												} */
						// P -= data[startingSearchIndex] * v_P[i];
						// C += data[startingSearchIndex] * v_W[i];
						startingSearchIndex = varToSorted[i];
					}
					// main idea, i have to start from SSI, I update the value C and P from the position of the fix to SSI
					else if (startingSearchIndex < varToSorted[i])
						for (Index j = startingSearchIndex; j <= varToSorted[i]; ++j)
						{
							Index z = sortedVar[j];
							if (skip[z])
								continue;
							P -= v_P[z];
							C += v_W[z];
						}
				}

				else
					throw(std::invalid_argument("Data for fixing in BinaryKnapsackBlockChange must be binary"));
			}
			return undo;
		}
		else
		{
			throw(std::invalid_argument("Change type not supported"));
		}
	}
	else
	{
		throw(std::invalid_argument("Change must be a BinaryKnapsackBlockChange"));
	}
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::get_var_solution(Configuration *solc)
{

	auto BKB = static_cast<BinaryKnapsackBlock *>(f_Block);

	BKB->set_x(v_x.begin()); // write the solution in BinaryKnapsackBlock
	// modify the critical item
	BKB->set_x(f_ci, (v_W[f_ci] < 0 ? 1 : 0));

} // end( IncrementalGreedyChangeBinaryKnapsackSolver::get_var_solution )

Solution *IncrementalGreedyChangeBinaryKnapsackSolver::get_Solution(Configuration *solc)
{
	std::vector<double> sol_x(v_x);
	sol_x[f_ci] = sol_x[f_ci] == 1 ? 1 : (v_W[f_ci] < 0 ? 1 : 0);
	BinaryKnapsackSolution *sol = new BinaryKnapsackSolution(std::move(sol_x));
	return sol;
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::add_Modification(sp_Mod &mod)
{

	if (f_no_Mod)
		return;

	// try to acquire lock, spin on failure
	while (f_mod_lock.test_and_set(std::memory_order_acquire))
		;

	// if NBModification, reload BinaryKnapsack instance and clear modifications
	if (const auto tmod = std::dynamic_pointer_cast<NBModification>(mod))
	{
		load();
		v_mod.clear();
	}
	else
		v_mod.push_back(mod);

	// release lock
	f_mod_lock.clear(std::memory_order_release);

} // end( add_Modification() )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::initializeVariables()
{
	// order sortedVar according to profit/weight ratio
	C = 0;
	P = 0;
	startingSearchIndex = 0;
	sortedVar.resize(f_N);
	varToSorted.resize(f_N);
	// complemented.resize(f_N);
	std::iota(sortedVar.begin(), sortedVar.end(), 0);
	// sort indeces in order of profit/weight
	sort(sortedVar.begin(), sortedVar.end(), [&](auto a, auto b)
		 { return (v_P[a] / v_W[a] > v_P[b] / v_W[b]); });
	for (int i = 0; i < f_N; ++i)
	{
		varToSorted[sortedVar[i]] = i;
	}
	f_ci = sortedVar.front();
	v_x.assign(f_N, 0);
	presolveCapacity = f_C;
	// lastFix.first = -1;
	skip.assign(f_N, false);
	for (int i = 0; i < f_N; ++i)
	{

		if (v_fxd[i] == 1)
		{
			v_x[i] = 0;
			skip[i] = true;
			continue;
		}

		// if the variable is fixed to 1 by the user
		if (v_fxd[i] == 2)
		{
			v_x[i] = 1;
			skip[i] = true;
			// controllare che abbia senso
			P += v_P[i];
			presolveCapacity -= v_W[i];
			continue;
		}

		// if the item has positive weight and negative profit
		if (v_W[i] >= 0 && v_P[i] <= 0)
		{
			v_x[i] = 0;
			// v_fxd[i] = 1;
			skip[i] = true;
		}

		// if the item has negative weight and positive profit
		else if (v_W[i] <= 0 && v_P[i] >= 0)
		{
			v_x[i] = 1;
			// v_fxd[i] = 2;
			skip[i] = true;
			presolveCapacity -= v_W[i]; // update weight
			P += v_P[i];				// update total Profit
		}
		else if (v_W[i] < 0 && v_P[i] < 0)
		{
			// v_x[i] = 1;
			presolveCapacity -= v_W[i]; // update weight
			P += v_P[i];				// update total Profit
			v_W[i] = -v_W[i];
			v_P[i] = -v_P[i];
			complemented[i] = true;
			// skip[i] = false;
		}
		else if (complemented[i])
		{
			// since it have already positive values, we have to consider the opposite operation
			presolveCapacity += v_W[i]; // update weight
			P -= v_P[i];				// update total Profit
		}
	}
	for (int i = 0; i < f_N; ++i)
	{
		if (!skip[i] && v_W[i] > presolveCapacity)
		{
			skip[i] = true;
			v_x[i] = complemented[i] ? 1 : 0;
		}
	}
}

void IncrementalGreedyChangeBinaryKnapsackSolver::load()
{

	if (!f_Block)
	{
		f_N = 0;
		v_P.clear();
		v_W.clear();
		v_x.clear();
		v_fxd.clear();
		sortedVar.clear();
		varToSorted.clear();
		// lastFix = std::make_pair(-1, false);
		complemented.clear();
		skip.clear();
		v_x.clear();
		presolveCapacity = 0;
		C = 0;
		P = 0;
		startingSearchIndex = 0;
		return;
	}

	auto BKB = dynamic_cast<BinaryKnapsackBlock *>(f_Block);
	if (!BKB)
		throw(std::invalid_argument("Block must be a BinaryKnapsackBlock"));

	// (try to) lock the BinaryKnapsackBlock
	bool owned = BKB->is_owned_by(f_id);
	if ((!owned) && (!BKB->read_lock()))
		throw(std::runtime_error("Unable to lock the Block"));

	// load Binary Knapsack instance - - - - - - - - - - - - - - - - - - - - -

	// get scalar data: sense of the objective, number of items and capacity
	f_sense = (BKB->get_objective_sense() == Objective::eMax);
	f_N = BKB->get_NItems();
	f_C = BKB->get_Capacity();

	// resize all
	v_P.resize(f_N);
	v_W.resize(f_N);
	v_fxd.resize(f_N);
	complemented.resize(f_N);
	// get references to profits, weights and integrality values
	const auto &P = BKB->get_Profits();
	const auto &W = BKB->get_Weights();
	const auto &FXD = BKB->get_fxd();

	for (Index i = 0; i < f_N; ++i)
	{
		v_P[i] = f_sense ? P[i] : -P[i];
		v_W[i] = W[i];
		v_fxd[i] = FXD[i];
		complemented[i] = false;
		if (v_W[i] < 0 && v_P[i] < 0)
		{
			v_W[i] = -v_W[i];
			v_P[i] = -v_P[i];
			complemented[i] = true;
		}
	}
	changedData = true;
	if (!owned)
		BKB->read_unlock();

} // end( IncrementalGreedyChangeBinaryKnapsackSolver::load() )

/*--------------------------------------------------------------------------*/

void IncrementalGreedyChangeBinaryKnapsackSolver::process_outstanding_Modification()
{

	// copy v_mod in a temporary list of modifications - - - - - - - - - - - - -

	Lst_sp_Mod v_mod_tmp; // temporary list of modifications

	// try to acquire lock, spin on failure
	while (f_mod_lock.test_and_set(std::memory_order_acquire))
		;

	for (auto mod : v_mod)
		v_mod_tmp.push_back(mod); // copy v_mod in v_mod_tmp

	v_mod.clear();

	f_mod_lock.clear(std::memory_order_release); // release lock

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	auto BKB = static_cast<BinaryKnapsackBlock *>(f_Block);

	// Any changes in Profits must be processed only AFTER checking the changes
	// on the sense of the Objective. Hence v_mod_tmp is scanned twice, checking
	// modifications on Objective sense (and also Capacity) first

	auto mod = v_mod_tmp.begin();

	while (mod != v_mod_tmp.end())
	{
		changedData = true;

		// BinaryKnapsackBlockMod - - - - - - - - - - - - - - - - - - - - - - - - -
		if (const auto tmod = dynamic_cast<BinaryKnapsackBlockMod *>(mod->get()))
		{
			switch (tmod->type())
			{

			case (BinaryKnapsackBlockMod::eChgCapacity):
			{
				f_C = BKB->get_Capacity();
				mod = v_mod_tmp.erase(mod);
				break;
			}

			case (BinaryKnapsackBlockMod::eChgSense):
				f_sense = (BKB->get_objective_sense() == Objective::eMax);
				for (Index i = 0; i < f_N; ++i)
				{
					if (complemented[i])
					{
						complemented[i] = false;
						v_W[i] = -v_W[i];
					}
					else
						v_P[i] = -v_P[i];
				}
				mod = v_mod_tmp.erase(mod);
				break;

			default:
				mod++;
			}
		}
		else
			mod = v_mod_tmp.erase(mod); // it is not a physical modification

	} // end( while() )

	for (auto mod : v_mod_tmp)
	{
		changedData = true;
		// BinaryKnapsackBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - -
		if (const auto tmod = dynamic_cast<BinaryKnapsackBlockRngdMod *>(
				mod.get()))
		{
			switch (tmod->type())
			{

			case (BinaryKnapsackBlockMod::eChgProfit):
				for (Index i = tmod->rng().first; i < tmod->rng().second; i++)
				{
					double oldP = v_P[i];
					v_P[i] = f_sense ? BKB->get_Profit(i) : -BKB->get_Profit(i);
					if (complemented[i] && oldP * v_P[i] >= 0)
					{
						complemented[i] = false;
						v_W[i] = -v_W[i];
					}
					else if (complemented[i])
						v_P[i] = -v_P[i];
				}
				break;

			case (BinaryKnapsackBlockMod::eFixX):
				for (Index i = tmod->rng().first; i < tmod->rng().second; i++)
				{
					// check if the variable is fixed to 0 or to 1
					if (BKB->is_fixed(i) && std::abs(BKB->get_x(i)) < 1e-6)
						v_fxd[i] = 1;
					else if (BKB->is_fixed(i) && std::abs(BKB->get_x(i) - 1) < 1e-6)
						v_fxd[i] = 2;
					if (complemented[i])
					{
						v_P[i] = -v_P[i];
						v_W[i] = -v_W[i];
						complemented[i] = false;
					}
				}
				break;

			case (BinaryKnapsackBlockMod::eUnfixX):
				for (Index i = tmod->rng().first; i < tmod->rng().second; i++)
					v_fxd[i] = 0;
				break;

			case (BinaryKnapsackBlockMod::eChgWeight):
				for (Index i = tmod->rng().first; i < tmod->rng().second; i++)
				{
					double oldW = v_W[i];
					v_W[i] = BKB->get_Weight(i);
					if (complemented[i] && oldW * v_W[i] >= 0)
					{
						complemented[i] = false;
						v_P[i] = -v_P[i];
					}
					else if (complemented[i])
						v_W[i] = -v_W[i];
				}
				break;
			}
		}

		// BinaryKnapsackBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - -
		if (const auto tmod = dynamic_cast<BinaryKnapsackBlockSbstMod *>(
				mod.get()))
		{
			changedData = true;

			switch (tmod->type())
			{

			case (BinaryKnapsackBlockMod::eChgProfit):
				for (auto i : tmod->nms())
				{
					double oldP = v_P[i];
					v_P[i] = f_sense ? BKB->get_Profit(i) : -BKB->get_Profit(i);
					if (complemented[i] && oldP * v_P[i] >= 0)
					{
						complemented[i] = false;
						v_W[i] = -v_W[i];
					}
					else if (complemented[i])
						v_P[i] = -v_P[i];
				}
				break;

			case (BinaryKnapsackBlockMod::eFixX):
				for (auto i : tmod->nms())
				{
					// check if the variable is fixed to 0 or to 1
					if (BKB->is_fixed(i) && std::abs(BKB->get_x(i)) < 1e-6)
						v_fxd[i] = 1;
					else if (BKB->is_fixed(i) && std::abs(BKB->get_x(i) - 1) < 1e-6)
						v_fxd[i] = 2;
					if (complemented[i])
					{
						v_P[i] = -v_P[i];
						v_W[i] = -v_W[i];
						complemented[i] = false;
					}
				}
				break;

			case (BinaryKnapsackBlockMod::eUnfixX):
				for (auto i : tmod->nms())
					v_fxd[i] = 0;

				break;

			case (BinaryKnapsackBlockMod::eChgWeight):
				for (auto i : tmod->nms())
				{
					double oldW = v_W[i];
					v_W[i] = BKB->get_Weight(i);
					if (complemented[i] && oldW * v_W[i] >= 0)
					{
						complemented[i] = false;
						v_P[i] = -v_P[i];
					}
					else if (complemented[i])
						v_W[i] = -v_W[i];
				}
				break;
			}
		}
	} // end( for(  ) )

	v_mod_tmp.clear(); // clear the temporary list of Modification

} // end( process_outstanding_Modification() )

/*--------------------------------------------------------------------------*/
/*----------------- End File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
