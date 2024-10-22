#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <list>
#include "CbcEventHandler.hpp"
#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "CoinModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiSolverInterface.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"
#include "CoinTime.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveGuided.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcHeuristic.hpp"
#include "CoinError.hpp"
#include "OsiCuts.hpp"
#include "CglCutGenerator.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"
#include "ClpSimplex.hpp"
#include "ClpPresolve.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinBuild.hpp"
#include "CbcBranchDynamic.hpp"
#include "CbcBranchDecision.hpp"
#include "CbcBranchDefaultDecision.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglSimpleRounding.hpp"
#include "CbcHeuristicDINS.hpp"
#include "CbcHeuristicDive.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDW.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
#include "CbcHeuristicRandRound.hpp"
#include "CbcHeuristicRENS.hpp"
#include "CbcHeuristicVND.hpp"
#include "CglGMI.hpp"
#include "CglRedSplit2.hpp"
#include "CglResidualCapacity.hpp"
#include "CglZeroHalf.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcCompare.hpp"
#include "CbcBranchActual.hpp"
#include "ClpSolve.hpp"
#include "CoinWarmStart.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
//#include <memory>



//=================================================================================================================
using namespace std;
int b;
bool isLog = true;
bool doWhsScaling = true;
bool whsSafe = false;
string jCbc_version = "2.9.10";

int Y2_simple(CbcModel* model, OsiClpSolverInterface* solver);
//==================================================================================================================	
class MyEventHandler3 : public CbcEventHandler {

public:
	/**@name Overrides */
	//@{
	virtual CbcAction event(CbcEvent whichEvent);
	//@}

	/**@name Constructors, destructor etc*/
	//@{
	/** Default constructor. */
	MyEventHandler3();
	/// Constructor with pointer to model (redundant as setEventHandler does)
	MyEventHandler3(CbcModel* model);
	/** Destructor */
	virtual ~MyEventHandler3();
	/** The copy constructor. */
	MyEventHandler3(const MyEventHandler3& rhs);
	/// Assignment
	MyEventHandler3& operator=(const MyEventHandler3& rhs);
	/// Clone
	virtual CbcEventHandler* clone() const;
	//@}


protected:
	// data goes here
};
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3()
	: CbcEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3(const MyEventHandler3& rhs)
	: CbcEventHandler(rhs)
{
}

// Constructor with pointer to model
MyEventHandler3::MyEventHandler3(CbcModel* model)
	: CbcEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyEventHandler3::~MyEventHandler3()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyEventHandler3&
MyEventHandler3::operator=(const MyEventHandler3& rhs)
{
	if (this != &rhs) {
		CbcEventHandler::operator=(rhs);
	}
	return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CbcEventHandler* MyEventHandler3::clone() const
{
	return new MyEventHandler3(*this);
}

CbcEventHandler::CbcAction
MyEventHandler3::event(CbcEvent whichEvent)
{
	if (b == 1) {
		return stop; // say finished
	}
	else {
		return noAction; // carry on

	}

}

//=============================================================================================================
// CoinModel build objects
int addCol(CoinModel* build, double collb, double colub, double obj, const char* name, bool isInt) {	
	
	try { 
		build->addCol(0, NULL, NULL, collb, colub, obj, name, isInt); 
	}
	catch (CoinError e) { e.print(); return 91; }
	catch (...) { cout << "\nUnknown exception caught in addCol!\n"; return 92; }
	return 0;
}

int addRow(CoinModel* build, int numberInRow, long int index[], double values[], double rowlb, double rowup, const char* name) {

	try { 
		build->addRow(numberInRow, (const int*)index, values, rowlb, rowup, name); 
	}
	catch (CoinError e) { e.print(); return 91; }
	catch (...) { cout << "\nUnknown exception caught in addRow!\n"; return 92; }
	return 0;
}


// Add rows at once to OsiClpSolverInterface
int addRows(OsiClpSolverInterface* solver, CoinModel* build) { 

	try { 
		solver->loadFromCoinModel(*build); 
	}
	catch (CoinError e) { e.print(); return 91; }
	catch (...) { cout << "\nUnknown exception caught in addRows!\n"; return 92; }
	return 0;

}

//set Integer variables in OsiClpSolverInterface
void setInteger(OsiClpSolverInterface* solver, int i) { solver->setInteger(i); }

// Assign OsiClpSolverInterface to CbcModel
void assignSolver(CbcModel* model, OsiClpSolverInterface* solver) {
	OsiSolverInterface* solver1 = solver;
	model->assignSolver(solver1);
	delete[] solver1;
}

// Reading/Writing Lp/Mps with OsiClpSolverInterface
void readMps(OsiClpSolverInterface* solver, const char* name) { solver->readMps(name, ""); }
int readLp(OsiClpSolverInterface* solver, const char* name) {
		
		try{
			solver->setLogLevel(0); solver->readLp(name);
		}
		catch(CoinError e){e.print();return 91;}
		catch (...) {cout << "\nUnknown exception caught in readLp!\n";return 92;}
		return 0;
}
void writeMps(OsiClpSolverInterface* solver, const char* name) { solver->writeMps(name); }
void writeLp(OsiClpSolverInterface* solver, const char* name) { solver->writeLp(name); }




// CbcModel functions
void branchAndBound(CbcModel* model) { model->branchAndBound(); }
void setLogLevel(CbcModel* model, int i) { model->setLogLevel(i); }


//=========================================================================================================
//solvers
void solve(CbcModel* model, OsiClpSolverInterface* solver, int logLevel = 0) {
	CbcMain0(*model);
	model->setLogLevel(logLevel);
	const char* argv2[] = { "driver3","-strong","5","-heuristicsOnOff","off","-presolve","off","-cutsOnOff","off","-primalS","-preprocess","off","-tune","100000000","-passc","1","-feas","off","-rins","off","-solve","-quit" };
	CbcMain1(22, argv2, *model);
}

void solve_1(CbcModel* model, OsiClpSolverInterface* solver, int logLevel = 0) {
	// Set strategy - below is == CbcStrategyDefault()
	OsiSolverInterface* solver2 = solver;
	CbcModel model1(*solver2);

	/*model->assignSolver(solver2);
	model->solver()->setIntParam(OsiNameDiscipline,1);
  CbcStrategyDefault strategy(true,5,5);
  model->setStrategy(strategy);
  // Do complete search
  model->setLogLevel(logLevel);
  model->branchAndBound();*/

	model1.branchAndBound();
	model->gutsOfCopy(model1, 2);
	model->setProblemStatus(model1.status());
	model->setSecondaryStatus(model1.secondaryStatus());


}

double round(double a, int n) {
	int temp = pow(10., n);
	double temp1 = a * temp;
	double temp2 = temp1 - floor(temp1);
	if (temp2 <= 0.5)
		return floor(temp1) / temp;
	else
		return ceil(temp1) / temp;

}

void cutoff(CbcModel* model, int n) {
	const CoinPackedMatrix* B = model->solver()->getMatrixByCol();
	int nnz = B->getNumElements();
	const double* el = B->getElements();
	double* el_new = new double[nnz];
	for (int i = 0; i < nnz; i++) {
		el_new[i] = round(el[i], 11);
	}
	CoinPackedMatrix* A = new CoinPackedMatrix(true, B->getMinorDim(), B->getMajorDim(), nnz, el_new, B->getIndices(), B->getVectorStarts(), B->getVectorLengths());
	model->solver()->replaceMatrix(*A);
}


int solve_2(CbcModel* model, OsiClpSolverInterface* solver, int logLevel = 0) {

	try {
	solver->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
	solver->setHintParam(OsiDoScale, false, OsiHintTry); 
	solver->setLogLevel(logLevel);


	//int n = solver->getNumCols();
	//cout << "-----------NumVars = " << n << " --------------\n"; 
	//std::string * Names = new std::string[n];
	 //double time = CoinCpuTime();

	//int n = model->solver()->getNumCols();
	//for (int i = 0; i < n; i++) {
	//	if (model->solver()->getColName(i) == "shrtg_64_xa") {
	//		cout << "++++++++++++++++++++;" << "\n";
	//		solver->setColBounds(i, 0, 0);

	//	}
	//}
	OsiClpSolverInterface solver1 = OsiClpSolverInterface(*solver);


	OsiSolverInterface* solver2 = &solver1;
	CglPreProcess* process = new CglPreProcess();
	process->messageHandler()->setLogLevel(logLevel);
	solver2 = process->preProcess(solver1, false, 1);
	//model.assignSolver(solver2);
	if (solver2) {
		//cout << " preprocess says feasible! " << "\n";
		solver2->resolve();
	}
	else {
		//model most probably infeasible but double check!
		//cout << " preprocess says infeasible! Let's make sure by solving without preprocessing! " << "\n";
		solver->setHintParam(OsiDoScale, true, OsiHintTry);

		CbcModel model1(*solver);

		model1.setIntegerTolerance(model->getIntegerTolerance());
		//model1.branchAndBound();
		model1.setLogLevel(0);
		callCbc("-loglevel 0 -solv ", model1);
		if (model1.secondaryStatus() == 7) {
			model->setProblemStatus(0);
			model->setSecondaryStatus(7);
			delete process, solver2, model1;
			return 1;
		}

		else if (model1.isProvenInfeasible()) {
			//model->gutsOfCopy(model1, 2);
			model->setProblemStatus(0);
			model->setSecondaryStatus(1);
			//cout << " \nmodel infeasible\n";
			delete process, solver2, model1;
			return 1;
		}
		else if (model1.isProvenOptimal()) {


			const double* sol = model1.bestSolution();
			double Obj = model1.getObjValue();
			int nCols = model1.getNumCols();
			model->gutsOfCopy(model1, 2);
			model->setBestSolution(sol, nCols, Obj, false);

			int violationStatus = Y2_simple(model, solver);
			if (violationStatus) {
				model->setProblemStatus(0);
				model->setSecondaryStatus(9000 + violationStatus);
				delete process, solver2, model1, sol;
				return 9;
			}

			model->setProblemStatus(0);
			model->setSecondaryStatus(0);
			delete process, solver2, model1, sol;
			//cout << "=====================" <<model->isProvenInfeasible() << " =======" << model->isProvenOptimal() <<endl;
			return 0;

		}

	}


	//solver2->initialSolve();
	//solver2->resolve();

	CbcModel model1(*solver2);

	if (true) {
		model1.solver()->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
		model1.solver()->setHintParam(OsiDoScale, false, OsiHintTry);
		model1.setLogLevel(logLevel);
		model1.messageHandler()->setLogLevel(logLevel);
		model1.solver()->messageHandler()->setLogLevel(logLevel);

		model1.setIntegerTolerance(model->getIntegerTolerance());

		model1.initialSolve();


		//================================================================================================================
		/*
		CglProbing generator1;
		generator1.setUsingObjective(true);
		generator1.setMaxPass(1);
		generator1.setMaxPassRoot(5);
		// Number of unsatisfied variables to look at
		generator1.setMaxProbe(10);
		generator1.setMaxProbeRoot(1000);
		// How far to follow the consequences
		generator1.setMaxLook(50);
		generator1.setMaxLookRoot(500);
		// Only look at rows with fewer than this number of elements
		generator1.setMaxElements(200);
		generator1.setRowCuts(3);

		CglGomory generator2;
		// try larger limit
		generator2.setLimit(300);

		CglKnapsackCover generator3;

		CglRedSplit generator4;
		// try larger limit
		generator4.setLimit(200);

		CglClique generator5;
		generator5.setStarCliqueReport(false);
		generator5.setRowCliqueReport(false);

		CglMixedIntegerRounding2 mixedGen;
		CglFlowCover flowGen;

		CglGMI cut1;
		CglMixedIntegerRounding2 cut2;
		CglOddHole cut3;
		CglSimpleRounding cut4;
		CglResidualCapacity cut5;
		CglTwomir cut6;
		CglZeroHalf cut7;


		// Add in generators
		// Experiment with -1 and -99 etc
		model1.addCutGenerator(&generator1,-1,"Probing");
		model1.addCutGenerator(&generator2,-1,"Gomory");
		model1.addCutGenerator(&generator3,-1,"Knapsack");
		model1.addCutGenerator(&generator4,-1,"RedSplit");
		model1.addCutGenerator(&generator5,-1,"Clique");
		model1.addCutGenerator(&flowGen,-1,"FlowCover");
		model1.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");
		model1.addCutGenerator(&cut1,-1,"GMI");
		model1.addCutGenerator(&cut2,-1,"MixedIntegerRounding2");
		model1.addCutGenerator(&cut3,-1,"OddHole");
		model1.addCutGenerator(&cut4,-1,"SimpleRounding");
		model1.addCutGenerator(&cut5,-1,"ResidualCapacity");
		model1.addCutGenerator(&cut6,-1,"Twomir");
		model1.addCutGenerator(&cut7,-1,"ZeroHalf");

	   */
	   // Uncommenting this should switch off all CBC messages
	   // model1.messagesPointer()->setDetailMessages(10,10000,NULL);
	   // Allow rounding heuristic

	   //CbcRounding heuristic1(model1);
	   //CbcHeuristicLocal heuristic2(model1);

		 //CbcHeuristicDiveLineSearch heuristic3(model1);
		 //CbcHeuristicDivePseudoCost heuristic4(model1);
		 //CbcHeuristicDiveVectorLength heuristic5(model1);
		 //CbcHeuristicDW heuristic6(model1);

		 //CbcHeuristicLocal heuristic8(model1);
		 //CbcHeuristicPivotAndFix heuristic9(model1);
		 //CbcHeuristicRENS heuristic10(model1);
		 //CbcHeuristicVND heuristic11(model1);



	  //model1.addHeuristic(&heuristic1);
	   //model1.addHeuristic(&heuristic2);
		 //model1.addHeuristic(&heuristic3);
		 //model1.addHeuristic(&heuristic4);
		 //model1.addHeuristic(&heuristic5);
		 //model1.addHeuristic(&heuristic6);

		 //model1.addHeuristic(&heuristic8);
		 //model1.addHeuristic(&heuristic9);
		 //model1.addHeuristic(&heuristic10);
		 //model1.addHeuristic(&heuristic11);

	   // Do initial solve to continuous
	   //model1.initialSolve();

	   // Could tune more
		double objValue = model1.solver()->getObjSense() * model1.solver()->getObjValue();
		double minimumDropA = CoinMin(1.0, fabs(objValue) * 1.0e-3 + 1.0e-4);
		double minimumDrop = fabs(objValue) * 1.0e-4 + 1.0e-4;
		//printf("min drop %g (A %g)\n",minimumDrop,minimumDropA);
		model1.setMinimumDrop(minimumDrop);



		model1.setMaximumCutPassesAtRoot(50); // use minimum drop
		model1.setMaximumCutPasses(1000);
		// Switch off strong branching if wanted
		// model1.setNumberStrong(0);
		// Do more strong branching if small


		model1.setNumberStrong(8);

		model1.setNumberBeforeTrust(5);

		model1.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);



		model1.messageHandler()->setLogLevel(logLevel);
		model1.solver()->messageHandler()->setLogLevel(logLevel);


		// Default strategy will leave cut generators as they exist already
		// so cutsOnlyAtRoot (1) ignored
		// numberStrong (2) is 5 (default)
		// numberBeforeTrust (3) is 5 (default is 0)
		// printLevel (4) defaults (0)
		//CbcStrategyDefault strategy(5);

		//  strategy.setupPreProcessing(2);
		//model1.setStrategy(strategy);
		model1.setLogLevel(logLevel);

		model1.findCliques(true, 10, 20);
		model1.setTypePresolve(2);
		//model1.setSpecialOptions(4);
		model1.setMoreSpecialOptions(16777216);
		model1.setMoreSpecialOptions(32768);
		model1.setMoreSpecialOptions(1024);
		model1.setMoreSpecialOptions2(1);
		model1.setMoreSpecialOptions2(128);
		//model1.setMoreSpecialOptions2(8192); 
		//model1.setThreadMode(2);
		 //model1.setNumberThreads(4);
		  //================================================================================================================
		  //CbcTreeLocal localTree(&model,model.solver()->getColSolution(),10,0,50,10000,2000);
		  //model.passInTreeHandler(localTree);

	}
	model1.solver()->setDblParam(OsiPrimalTolerance, 1e-09);
	model1.branchAndBound();
	//cout << "1;;" << model1.status() << ";;" << model1.secondaryStatus()  << ";;" << model1.isProvenOptimal() << "\n";
	OsiSolverInterface* final;
	process->postProcess(*model1.solver());
	final = &solver1;//now a solution for the original model has been made and copied to the final
	//I am using final to make sure everythink is correct. If it cause memory overflow, use solver instead of final in the following lines
	// to get the solution and objective.

	//model1 = *process->originalModel();
	//cout << "2;;" << model1.getNumCols() << ";;" << model->getNumCols() << "\n";
	//cout << "3;;" << model1.getObjValue() << ";" << model->getObjValue() << ";" << solver->getObjValue() << ";" << model1.solver()->getObjValue() << ";" << model->solver()->getObjValue() << "\n";
	//model->gutsOfCopy(model1, 2);
	//cout << "4;;" << model1.getNumCols() << ";;" << model->getNumCols() << "\n";
	//cout << "5;;" << model1.getObjValue() << ";" << model->getObjValue() << ";" << solver->getObjValue() << ";" << model1.solver()->getObjValue() << ";" << model->solver()->getObjValue();
	//model->setProblemStatus(model1.status());
	//model->setSecondaryStatus(model1.secondaryStatus());



	//CbcModel *temp = model;
	//temp = model;
	//model1.assignSolver(solver3);
	//model = model1.clone(true);
	//model = model1;


	//model1.gutsOfDestructor();




	if (model1.isProvenOptimal()) {

		const double* sol = final->getColSolution();
		const double Obj = model1.getObjValue();
		int nCols = final->getNumCols();
		model->setBestSolution(sol, nCols, Obj, false);
		model->setObjValue(Obj);

		int violationStatus = Y2_simple(model, solver);
		if (violationStatus) {
			model->setProblemStatus(0);
			model->setSecondaryStatus(9000 + violationStatus);
			delete process, solver2, model1, sol;
			return 9;
		}

		model->setProblemStatus(0);
		model->setSecondaryStatus(0);
		delete process, solver2, model1, sol;
		return 0;

	}

	else if (model1.isProvenInfeasible()) {

		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
		delete process, solver2, model1;
		return 1;

	}
	else {
		model->setProblemStatus(-1);
		model->setSecondaryStatus(-1);
		delete process, solver2, model1;
		return 2;
	}




	}
	catch (CoinError e) {
		if (isLog) { e.print(); }
		return 91;
	}
	catch (...) {
		if (isLog) { cout << "\nUnknown exception caught in Solve_2\n"; }
		return 92;
	}
}



int solve_3(CbcModel* model, OsiClpSolverInterface* solver, int logLevel = 0, double presolve_tolerance = 1e-07) {

	model->setLogLevel(logLevel);
	solver->setLogLevel(logLevel);
	ClpSimplex* simplex = solver->getModelPtr();

	ClpPresolve pinfo;
	double time_2 = CoinCpuTime();
	ClpSimplex* simplex2 = pinfo.presolvedModel(*simplex, presolve_tolerance, true, 3, false, false, 0, simplex->integerInformation());
	double time_3 = CoinCpuTime() - time_2;
	//cout << "Presolve= "<< time_3 << "\n";


	bool status = false;
	if (simplex2) {
		//cout << "there is simplex2\n";
		simplex2->setLogLevel(logLevel);
		status = simplex2->isProvenPrimalInfeasible();

	}
	if (!simplex2 || status) {

		//cout << "Presolve says infeasible, let's make sure!\n";
		CbcModel model1(*solver);

		model1.branchAndBound();

		if (model1.secondaryStatus() == 7) {
			model->setProblemStatus(0);
			model->setSecondaryStatus(7);

			return 1;
		}

		if (model1.isProvenInfeasible()) {
			//model->gutsOfCopy(model1, 2);
			model->setProblemStatus(0);
			model->setSecondaryStatus(1);
			//cout << " \nmodel infeasible\n";

			return 1;
		}
		if (model1.isProvenOptimal()) {


			const double* sol = model1.bestSolution();
			double Obj = model1.getObjValue();
			int nCols = model1.getNumCols();
			//model->gutsOfCopy(model1, 2);
			model->setBestSolution(sol, nCols, Obj, false);


			model->setProblemStatus(0);
			model->setSecondaryStatus(0);

			int violationStatus = Y2_simple(model, solver);
			if (violationStatus) {
				model->setSecondaryStatus(9000 + violationStatus);
				return 9;
			}


			//cout << "=====================" <<model->isProvenInfeasible() << " =======" << model->isProvenOptimal() <<endl;
			return 0;
		}
	}




	OsiClpSolverInterface solver2(simplex2);
	solver2.messageHandler()->setLogLevel(logLevel);

	solver2.setHintParam(OsiDoScale, false, OsiHintTry);
	CbcModel model1(solver2);
	model1.setLogLevel(logLevel);
	model1.solver()->messageHandler()->setLogLevel(logLevel);

	model1.setIntegerTolerance(model->getIntegerTolerance());

	double time = CoinCpuTime();
	//model1.setLogLevel(logLevel);
	CbcMain0(model1);

	model1.setLogLevel(logLevel);
	const char* argv2[] = { "driver3","-preprocess","off","-rins","off" ,"-feas","off","-solve","-quit" };
	CbcMain1(9, argv2, model1);
	double time1 = CoinCpuTime() - time;


	if (model1.isInitialSolveProvenPrimalInfeasible()) {
		//cout << "\nStopped due to numerical instabilities.\n";
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
		delete simplex2;
		return 1;
	}


	if (model1.solver()->isProvenPrimalInfeasible() == 1) {
		//cout << "\nModel Infeasible\n";
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
		delete simplex2;


		return 1;

	}

	OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (model1.solver());
	assert(clpSolver);

	clpSolver->setLogLevel(logLevel);
	ClpSimplex* clp = clpSolver->getModelPtr();

	clp->setLogLevel(logLevel);
	*simplex2 = *clp;
	double time4 = CoinCpuTime();
	simplex2->checkSolution();

	pinfo.postsolve(true);
	//cout << "postsolve" <<  CoinCpuTime() - time4;

	if (true) {
		const int* original = pinfo.originalColumns();
		double* lower2 = simplex2->columnLower();
		double* upper2 = simplex2->columnUpper();
		const char* info2 = simplex2->integerInformation();
		double* lower = simplex->columnLower();
		double* upper = simplex->columnUpper();
		int i;

		for (i = 0; i < simplex2->numberColumns(); i++) {
			if (info2[i]) {
				int iSeq = original[i];
				upper[iSeq] = upper2[i];
				lower[iSeq] = lower2[i];
			}
		}
	}

	simplex->initialSolve();
	//model->gutsOfCopy(model1, 2);
	simplex->checkSolution(1);

	model->setProblemStatus(model1.status());
	model->setSecondaryStatus(model1.secondaryStatus());


	if (model1.isProvenOptimal()) {
		const double* sol = simplex->getColSolution();

		double Obj = model1.getObjValue();
		int nCols = simplex->getNumCols();
		//model->gutsOfCopy(model1, 2);

		model->setBestSolution(sol, nCols, Obj, true);


		if (model->bestSolution()) {

			model->setProblemStatus(0);
			model->setSecondaryStatus(0);
			int violationStatus = Y2_simple(model, solver);
			if (violationStatus) {
				model->setSecondaryStatus(9000 + violationStatus);
				return 9;
			}
			return 0;
		}
		else {
			model->setProblemStatus(-1);
			model->setSecondaryStatus(-1);
			return 2;
		}


	}

	else if (model1.isProvenInfeasible()) {
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
	}


	delete simplex2;


	//} catch(CoinError e) {
	//e.print();
	//}

}


//============================================================================================================================

// OsiClpSolverInterface functions
void initialSolve(OsiClpSolverInterface* solver) { solver->initialSolve(); }
void setLogLevel(OsiClpSolverInterface* solver, int i) { solver->setLogLevel(i); }
int isInteger(OsiClpSolverInterface* solver, int i) { return solver->isInteger(i); }
int isBinary(OsiClpSolverInterface* solver, int i) { return solver->isBinary(i); }

// All these methods can also be done using OsiClpSolverInterface, with CbcModel it will probably take longer time.
const double* getSol(CbcModel* model) { return model->solver()->getColSolution(); }
int isInteger(CbcModel* model, int i) { return model->isInteger(i); }

void addRows(CbcModel* model, CoinModel* build) {
	OsiClpSolverInterface* solver1;
	solver1 = new OsiClpSolverInterface;
	OsiSolverInterface* solver = solver1;
	solver->loadFromCoinModel(*build);
	model->assignSolver(solver);
	delete[] solver, solver1;
}

void readMps(CbcModel* model, const char* name) { model->solver()->readMps(name, ""); }
void readLp(CbcModel* model, const char* name) { model->solver()->readLp(name); }
void writeMps(CbcModel* model, const char* name) { model->solver()->writeMps(name); }
void writeLp(CbcModel* model, const char* name) { model->solver()->writeLp(name); }
void setInteger(CbcModel* model, int i) { model->solver()->setInteger(i); }

//Names
std::string getRowName(OsiClpSolverInterface* solver, int i) { return solver->getRowName(i); }
std::string getColName(OsiClpSolverInterface* solver, int i) { return solver->getColName(i); }
std::string getRowName(CbcModel* model, int i) { return model->solver()->getRowName(i); }
std::string getColName(CbcModel* model, int i) { return model->solver()->getColName(i); }
std::string getVersion() {
	std::stringstream buffer;
	buffer << jCbc_version << " (" << __DATE__ << ")";
	return buffer.str();
}
void setWhsScaling(bool solvWhs_scaling ) { 
	doWhsScaling = solvWhs_scaling;
}
void setWhsSafe(bool solvWhs_safe) {
	whsSafe = solvWhs_safe;
}
void setModelName(OsiClpSolverInterface* solver, std::string name) { solver->setStrParam(OsiProbName, name); }
std::string  getModelName(OsiClpSolverInterface* solver) {
	std::string temp;
	solver->getStrParam(OsiProbName, temp);
	return temp;
}

//problem status
int status(CbcModel* model) { return model->status(); }//returns 0 if finished (which includes the case 
													//when the algorithm is finished because it has been proved infeasible), 
													//1 if stopped by user, and 2 if difficulties arose.
int isProvenOptimal(CbcModel* model) { return model->isProvenOptimal(); }
int isProvenInfeasible(CbcModel* model) { return model->isProvenInfeasible(); }
int secondaryStatus(CbcModel* model) { return model->secondaryStatus(); }
/*  cbc secondary status of problem
		-1 unset (status_ will also be -1)
	0 search completed with solution
	1 linear relaxation not feasible (or worse than cutoff)
	2 stopped on gap
	3 stopped on nodes
	4 stopped on time
	5 stopped on user event
	6 stopped on solutions
	7 linear relaxation unbounded*/

	//problem stats
int getNumRows(OsiClpSolverInterface* solver) { return solver->getNumRows(); }
int getNumCols(OsiClpSolverInterface* solver) { return solver->getNumCols(); }
int getNumRows(CbcModel* model) { return model->getNumRows(); }
int getNumCols(CbcModel* model) { return model->getNumCols(); }

int  numberIntegers(CbcModel* model) { return model->numberIntegers(); }

//solution
const double* getColSolution(OsiClpSolverInterface* solver) { return solver->getColSolution(); }
const double* getRowPrice(OsiClpSolverInterface* solver) { return solver->getRowPrice(); }
const double* getRowActivity(OsiClpSolverInterface* solver) { return solver->getRowActivity(); }
const double* getReducedCost(OsiClpSolverInterface* solver) { return solver->getReducedCost(); }

const double* getColSolution(CbcModel* model) { return model->bestSolution(); }
const double* getRowPrice(CbcModel* model) { return model->getRowPrice(); }
const double* getRowActivity(CbcModel* model) { return model->getRowActivity(); }
const double* getReducedCost(CbcModel* model) { return model->getReducedCost(); }

const double getObjValue(CbcModel* model) { return model->getObjValue(); }
const double getObjValue(OsiClpSolverInterface* solver) { return solver->getObjValue(); }

//time
double getCoinCpuTime() { return CoinCpuTime(); }

//CoinPackedMatrix
void loadProblem(OsiClpSolverInterface* solver, CoinPackedMatrix* byRow, double  columnLower[], double columnUpper[],
	double objective[], double rowLower[], double rowUpper[]) {

	solver->loadProblem(*byRow, columnLower, columnUpper, objective, rowLower, rowUpper);

}

CoinPackedMatrix* createMatrix(int numberColumns, int numberRows, int numberElements,
	double elementByRow[], long int column[], long int rowStart[]) {

	CoinPackedMatrix* M = new CoinPackedMatrix(false, numberColumns, numberRows,
		numberElements, elementByRow, (const int*)column, (const int*)rowStart, NULL);
	return M;
}

void setRowName(OsiClpSolverInterface* solver, int i, const char* name) {
	solver->setRowName(i, name);
}


//tolerances
void setPrimalTolerance(CbcModel* model, double a) {
	model->solver()->setDblParam(OsiPrimalTolerance, a);
}

void setDualTolerance(CbcModel* model, double a) {
	model->solver()->setDblParam(OsiDualTolerance, a);
}

void setIntegerTolerance(CbcModel* model, double a) {
	model->setIntegerTolerance(a);
}

void test(std::pair<std::string, int> a) {
	//cout << a.first << " + " << a.second;
}

void test1(std::vector<std::pair<std::string, int> > a) {
	//cout << a.at(0).first << " + " << a.at(1).second;
}

const char* bool_cast(const bool b) {
	return b ? "true" : "false";
}
int callCbcJ(std::string a, CbcModel *model, OsiClpSolverInterface* solver);
int solve_whs(CbcModel* model, OsiClpSolverInterface* solver, std::string names[], int values[], int intvars, int logLevel = 0, double presolve_tolerance = 1e-07) {

	try {

	std::vector< std::pair< std::string, double > > A;

	model->setLogLevel(logLevel);
	solver->setLogLevel(logLevel);
	if (!doWhsScaling) { solver->setHintParam(OsiDoScale, false, OsiHintTry); }
	//solver->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
	ClpSimplex* simplex = solver->getModelPtr();


	int numCols = solver->getNumCols();
	ClpPresolve pinfo;
	ClpSimplex* simplex2 = pinfo.presolvedModel(*simplex, presolve_tolerance, true, 3, false, false, 0, simplex->integerInformation());
	bool status_isInfeasible = false;
	if (simplex2) {
		//cout << "there is simplex2\n";
		simplex2->setLogLevel(logLevel);
		status_isInfeasible = simplex2->isProvenPrimalInfeasible();

	}
	//cout << "926 " << bool_cast(status)<<"\n";

	if (!simplex2 || status_isInfeasible) {

		if (whsSafe) {
			if (simplex2) { delete simplex2; }
			simplex = NULL;
			int returnV = -99;
			returnV = callCbcJ("-log 0 -primalT 1e-9 -integerT 1e-9 -solve", model, solver);
			return returnV;
		}


		//cout << "Presolve says infeasible, let's make sure!\n";
		CbcModel model1(*solver);

		model1.branchAndBound();

		if (model1.secondaryStatus() == 7) {

			if (simplex2) { delete simplex2; }
			simplex = NULL;
			model->setProblemStatus(0);
			model->setSecondaryStatus(7);
			return 1;
		}

		else if (model1.isProvenInfeasible()) {
			//model->gutsOfCopy(model1, 2);

			//cout << " \nmodel infeasible\n";
			if (simplex2) { delete simplex2; }
			simplex = NULL;
			model->setProblemStatus(0);
			model->setSecondaryStatus(1);
			return 1;
		}
		else if (model1.isProvenOptimal()) {


			const double* sol = model1.bestSolution();
			double Obj = model1.getObjValue();
			int nCols = model1.getNumCols();
			model->gutsOfCopy(model1, 2);
			model->setBestSolution(sol, nCols, Obj, false);

			int violationStatus = Y2_simple(model, solver);
			if (violationStatus) {

				if (simplex2) { delete simplex2; }
				simplex = NULL;
				model->setProblemStatus(0);
				model->setSecondaryStatus(9000 + violationStatus);
				return 9;
			}


			//cout << "=====================" <<model->isProvenInfeasible() << " =======" << model->isProvenOptimal() <<endl;
			if (simplex2) { delete simplex2; }
			simplex = NULL;
			model->setProblemStatus(0);
			model->setSecondaryStatus(0);
			return 0;
		}
		else {

			if (simplex2) { delete simplex2; }
			simplex = NULL;
			model->setProblemStatus(-1);
			model->setSecondaryStatus(-1);
			return 2;
		}

	}
	//cout << "Presolve says feasible, continue.\n";
	bool whs = false;
	const double* sol0 = NULL;

	ClpSimplex simplex_cpy(*simplex2);
	OsiClpSolverInterface solver1(&simplex_cpy);


	solver1.setLogLevel(logLevel);

	//int priority[solver1.getNumCols()];
	//int* priority = new int[solver1.getNumCols()];

	//for (int i = 0; i < solver1.getNumIntegers(); i++)
	//	priority[i] = 0;

	if (intvars > 0) {

		const double* UB = solver1.getColUpper();
		const double* LB = solver1.getColLower();

		for (int i = 0; i < solver1.getNumCols(); i++) {
			//solver1.setObjCoeff(i,0);
			if (solver1.isInteger(i)) {
				for (int k = 0; k < intvars; k++) {
					if (solver1.getColName(i) == names[k]) {
						if (UB[i] > LB[i]) {
							solver1.setColLower(i, values[k]);
							solver1.setColUpper(i, values[k]);

							//if (values[k] == 1)
							//	priority[i] = -1;
							//else
							//	priority[i] = -1;
						}
					}
				}
			}
		}

		solver1.initialSolve();



		if (!solver1.isProvenPrimalInfeasible() && !solver1.isProvenDualInfeasible()) {
			sol0 = solver1.getColSolution();
			whs = true;
			//delete[] priority;
			//CbcModel model_fake(solver1);
			//A = model_fake.getMIPStart();
		}
		else {
			//cout << "\\warm start solution infeasible\\";
			//delete[] priority;
			delete simplex2; 
			simplex = NULL;
			model->setProblemStatus(-1);
			model->setSecondaryStatus(-1);
			return 2;
		}

	}

	//const double *a = simplex2->getObjCoefficients();
	//int b[simplex2->getNumCols()]; 
	//for (int i = 0;i<simplex2->getNumCols();i++)
	//b[i]=i;
	//simplex2->addRow(simplex2->getNumCols(),b,a,-COIN_DBL_MAX,solver1.getObjValue());

	OsiClpSolverInterface solver2(simplex2);
	solver2.messageHandler()->setLogLevel(logLevel);
	if (!doWhsScaling) { solver2.setHintParam(OsiDoScale, false, OsiHintTry); }
	//solver2.setHintParam(OsiDoDualInInitial, false, OsiHintTry);



	CbcModel model1(solver2);
	model1.setLogLevel(logLevel);
	model1.solver()->messageHandler()->setLogLevel(logLevel);
	model1.setLogLevel(logLevel);

	model1.setIntegerTolerance(model->getIntegerTolerance());
	//model1.initialSolve();


	/*
	for (int i=0;i<model1.solver()->getNumCols();i++){

	if (model1.solver()->isInteger(i)){
	//cout << model1.solver()->getColName(i) << " = " << sol_tmp[i]<<endl;
	if (model1.solver()->getColName(i) =="int_hands" || model1.solver()->getColName(i) =="c_sac093intflow" ||model1.solver()->getColName(i) =="int_sac_below" ||model1.solver()->getColName(i) =="int_whs_stor_chg" ||model1.solver()->getColName(i) =="int_ibu_uwfe" )
	priority[i]=0;
	else if (model1.solver()->getColName(i) =="int_sweir"  || model1.solver()->getColName(i) =="int_fweir1"  || model1.solver()->getColName(i) =="int_fweir2"  || model1.solver()->getColName(i) =="int_tweir")
	priority[i]=0;
	else
	priority[i]=1;

	}
	}

	*/


	if (whs) {

		vector<pair<string, double> > A;
		for (int i = 0; i < solver1.getNumCols(); i++) {
			A.push_back(make_pair(solver1.getColName(i), sol0[i]));
		}

		//model1.setHotstartSolution(sol0, priority);
		model1.setMIPStart(A);
		//cout << "Successfully set the hot start\n";
		//model1.setBestSolution(sol0, solver1.getNumCols(), solver1.getObjValue(),true);
		//model1.setBestObjectiveValue(solver1.getObjValue());


	}

	//model1.setSpecialOptions(4);
	//model1.setMoreSpecialOptions(1024);
	//model1.setThreadMode(2);
	//model1.setNumberThreads(4);
	//model1.passInPriorities(priority, false);

	double time = CoinCpuTime();

	//CbcMain0(model1);
	//model1.solver()->writeMps("babak.mps");
	//model1.setLogLevel(logLevel);
	//const char * argv2[]={"-heuristicsOnOff off","-presolve","off","-cutsOnOff","off","-primalS","-preprocess","off","-tune","100000000","-passc","1","-feas","off","-rins","off","-solve","-quit"};
	//CbcMain1(18,argv2,model1);

	//stringstream convert;
	//convert << logLevel;
	//string Result = convert.str();
	//string log = "-log " + Result;

	CbcMain0(model1);
	model1.setLogLevel(logLevel);
	const char* argv2[] = { "driver3","-preprocess","off","-passC","1","-rins","off" ,"-feas","off","-tune","1000000","-solve","-quit" };
	CbcMain1(*(&argv2 + 1) - argv2, argv2, model1);
	//cout << "1073." << bool_cast(model1.isProvenOptimal()) << "\n";

	if (model1.isInitialSolveProvenPrimalInfeasible()) {
		//cout << "\nStopped due to numerical instabilities.\n";

		delete simplex2; 			
		simplex = NULL;
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
		return 1;
	}


	else if (model1.solver()->isProvenPrimalInfeasible()) {
		//cout << "\nModel Infeasible\n";

		delete simplex2;
		simplex = NULL;
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);

		return 1;

	}

	double time1 = CoinCpuTime() - time;


	OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (model1.solver());
	assert(clpSolver);

	clpSolver->setLogLevel(logLevel);
	ClpSimplex* clp = clpSolver->getModelPtr();

	clp->setLogLevel(logLevel);
	*simplex2 = *clp;
	simplex2->checkSolution();
	pinfo.postsolve(true);

	//simplex2->checkSolution();
	//model->gutsOfCopy(model1, 2);
	model->setProblemStatus(model1.status());
	model->setSecondaryStatus(model1.secondaryStatus());


	if (model1.isProvenOptimal()) {

		const double* sol = simplex->getColSolution();

		double Obj = model1.getObjValue();
		int nCols = simplex->getNumCols();
		//model->gutsOfCopy(model1, 2);

		model->setBestSolution(sol, nCols, Obj, true);


		if (model->bestSolution()) {

			int violationStatus = Y2_simple(model, solver);
			if (violationStatus) {

				delete simplex2; simplex = NULL;
				model->setProblemStatus(0);
				model->setSecondaryStatus(9000 + violationStatus);
				return 9;
			}

			delete simplex2; simplex = NULL;
			model->setProblemStatus(0);
			model->setSecondaryStatus(0);
			return 0;
		}
		else {

			delete simplex2; simplex = NULL;
			model->setProblemStatus(-1);
			model->setSecondaryStatus(-1);
			return 2;
		}


	}
	//else if (model1.isProvenInfeasible()) {
	else {
		delete simplex2; simplex = NULL;
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
		return 1;
	}


	} catch(CoinError e) {
		if (isLog) { e.print(); }
		return 91;
	} catch(...) {
		if (isLog) {cout << "\nUnknown exception caught in Solve_whs\n";}
		return 92;
	}

}



int solve_4(CbcModel* model, OsiClpSolverInterface* solver, int logLevel) {


	solver->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
	solver->setLogLevel(logLevel);
	OsiClpSolverInterface* temp = new OsiClpSolverInterface(*solver);

	OsiSolverInterface* solver2 = solver;
	CglPreProcess* process = new CglPreProcess;
	process->messageHandler()->setLogLevel(logLevel);
	solver2 = process->preProcess(*solver, false, 2);
	//model.assignSolver(solver2);
	if (solver2)
		cout << "";
	else {
		//model most probably infeasible but double check!

		OsiSolverInterface* temp1 = temp;
		CbcModel model1(*temp1);
		model1.branchAndBound();

		if (model1.isProvenInfeasible()) {
			//model->gutsOfCopy(model1, 2);
			model->setProblemStatus(0);
			model->setSecondaryStatus(1);
			//cout << " \nmodel infeasible\n";
			delete process, solver2, model1;
			return 0;
		}
		if (model1.isProvenOptimal()) {
			*solver = OsiClpSolverInterface(*temp);
			const double* sol = model1.getColSolution();
			double Obj = model1.getObjValue();
			int nCols = model1.getNumCols();
			model->setBestSolution(sol, nCols, Obj, true);
			model->gutsOfCopy(model1, 2);


			model->setProblemStatus(0);
			model->setSecondaryStatus(0);
			delete process, solver2, model1, sol;
			//cout << "=====================" <<model->isProvenInfeasible() << " =======" << model->isProvenOptimal() <<endl;
			return 1;
		}
	}


	CbcModel model1(*solver2);
	model1.setLogLevel(logLevel);
	model1.messageHandler()->setLogLevel(logLevel);
	model1.solver()->messageHandler()->setLogLevel(logLevel);

	model1.initialSolve();



	model1.setTypePresolve(2);
	model1.setSpecialOptions(4);
	model1.setMoreSpecialOptions(16777216);
	model1.setMoreSpecialOptions(32768);
	model1.setMoreSpecialOptions(1024);
	model1.setMoreSpecialOptions2(1);
	model1.setMoreSpecialOptions2(128);
	model1.setThreadMode(2);


	CbcMain0(model1);
	model->setLogLevel(logLevel);
	const char* argv2[] = { "driver3","-strong","5","-trust","0","-heuristicsOnOff","off","-presolve","off","-cutsOnOff","off","-primalS","-preprocess","off","-tune","100000000","-passc","1","-feas","off","-rins","off","-solve","-quit" };
	CbcMain1(24, argv2, model1);

	OsiSolverInterface* solver3;

	process->postProcess(*model1.solver());

	solver3 = solver;


	model->setProblemStatus(model1.status());
	model->setSecondaryStatus(model1.secondaryStatus());

	model->gutsOfCopy(model1, 2);




	/*
   if (model1.solver()->isProvenOptimal()){
	   model->setProblemStatus(0);
	   model->setSecondaryStatus(0);
   }
	if(!model1.isProvenInfeasible()){
			model->setProblemStatus(0);
			 model->setSecondaryStatus(1);
	   }*/
	delete process, solver3, solver2, model1;
	return 2;



}

int getNumIntegers(OsiClpSolverInterface* solver) {
	return solver->getNumIntegers();
}

CoinWarmStart* getWarmStart(CbcModel* model) {
	return model->solver()->getWarmStart();
}

void writeMps1(CbcModel* model, const char* name, int 	formatType = 0, int numberAcross = 2, double objSense = 0.0) {
	model->solver()->writeMpsNative(name, NULL, NULL, formatType, numberAcross, objSense, 0);
}

void writeLp1(CbcModel* model, const char* name, double epsilon = 1e-5, int decimals = 5) {
	model->solver()->writeLp(name, "lp", epsilon, 10, decimals, 0.0, true);
}

double minCoeff(CbcModel* model) {
	double* minimumNegative = new double;
	double* maximumNegative = new double;
	double* minimumPositive = new double;
	double* maximumPositive = new double;
	model->solver()->statistics(*minimumNegative, *maximumNegative,
		*minimumPositive, *maximumPositive,
		0);
	return min(fabs(*minimumNegative), fabs(*minimumPositive));
	delete minimumNegative;
	delete maximumNegative;
	delete minimumPositive;
	delete maximumPositive;
}

int callCbcJ(std::string a, CbcModel* model, OsiClpSolverInterface* solver) {
	try {
	callCbc(a, *model);

	if (model->secondaryStatus() == 7) {

		// one more try with branch and bound
		model->branchAndBound();
		//model->setProblemStatus(0);
		//model->setSecondaryStatus(7);
	}

	if (model->isProvenOptimal()) {
		model->setProblemStatus(0);
		model->setSecondaryStatus(0);

		int violationStatus = Y2_simple(model, solver);
		if (violationStatus) {
			model->setSecondaryStatus(9000 + violationStatus);
			return 9;
		}
		return 0;
	}
	else if (model->isProvenInfeasible() == 1) {
		model->setProblemStatus(0);
		model->setSecondaryStatus(1);
		return 1;
	}
	else {
		model->setProblemStatus(-1);
		model->setSecondaryStatus(-1);
		return 2;
	}

	}
	catch (CoinError e) {
		if (isLog) { e.print(); }
		return 91;
	}
	catch (...) {
		if (isLog) { cout << "\nUnknown exception caught in callCbc\n"; }
		return 92;
	}
	
}

double maxCoeff(CbcModel* model) {
	double* minimumNegative = new double;
	double* maximumNegative = new double;
	double* minimumPositive = new double;
	double* maximumPositive = new double;
	model->solver()->statistics(*minimumNegative, *maximumNegative,
		*minimumPositive, *maximumPositive,
		0);
	return max(fabs(*minimumNegative), fabs(*minimumPositive));


}

double minRHS(CbcModel* model) {
	double* minRHS = new double;
	model->solver()->getMinRHS(*minRHS);
	return *minRHS;

}

int solve_unified(CbcModel* model, OsiClpSolverInterface* solver, std::string names[] = NULL, int values[] = NULL, int intvars = 0, int logLevel = 0) {

	OsiClpSolverInterface* solver1 = new OsiClpSolverInterface(*solver);
	CbcModel* model1 = new CbcModel(*solver1);
	model1->setIntegerTolerance(1e-9);
	model1->solver()->setDblParam(OsiPrimalTolerance, 1e-9);
	solver1->setDblParam(OsiPrimalTolerance, 1e-9);
	int f = 0;
	if (intvars == 0) {
		//delete solver1;
		delete model1;
		goto try_solve_2;
	}
	else {
		f = solve_whs(model1, solver1, names, values, intvars, logLevel);
		if (model1->status() == 0 && model1->secondaryStatus() == 1) {
			//cout <<  model1.isProvenInfeasible() << "  " <<  model1.isProvenOptimal() << "\n";
			//cout << "whs inf, now trying solve_2\n";
			//delete solver1;
			delete model1;
			goto try_solve_2;
		}
		else if (f == 3) {
			//cout << "\nprovided solution for whs inf, now trying solve_2\n";
			//delete solver1;
			delete model1;
			goto try_solve_2;
		}
		else {
			*solver = *solver1;
			model->gutsOfCopy(*model1, 2);
			model->setProblemStatus(0);
			model->setSecondaryStatus(0);
			//delete solver1;
			delete model1;
			return 6;
		}
	}

try_solve_2:

	OsiClpSolverInterface* solver2 = new OsiClpSolverInterface(*solver);
	CbcModel* model2 = new CbcModel(*solver2);
	model2->setIntegerTolerance(1e-9);
	model2->solver()->setDblParam(OsiPrimalTolerance, 1e-8);
	solver2->setDblParam(OsiPrimalTolerance, 1e-8);

	solve_2(model2, solver2, logLevel);
	if (model2->status() == 0 && model2->secondaryStatus() == 1) {
		//cout <<  model2.isProvenInfeasible() << "  " <<  model2.isProvenOptimal() << "\n";
		//cout << "solve_2 inf, now trying solve_2 with tol 1e-4\n";
		//delete solver2;
		delete model2;
		goto try_solve_2_tol;

	}
	else {
		*solver = *solver2;
		model->gutsOfCopy(*model2, 2);
		model->setProblemStatus(0);
		model->setSecondaryStatus(0);
		//delete solver2;
		delete model2;
		return 5;
	}
try_solve_2_tol:

	OsiClpSolverInterface* solver3 = new OsiClpSolverInterface(*solver);
	CbcModel* model3 = new CbcModel(*solver3);
	assignSolver(model3, solver3);
	//model3->setIntegerTolerance(1e-9);

	model3->solver()->setDblParam(OsiPrimalTolerance, 1e-4);
	solver3->setDblParam(OsiPrimalTolerance, 1e-4);
	solve_2(model3, solver3, logLevel);

	if (model3->status() == 0 && model3->secondaryStatus() == 1) {
		//cout << "Model infeasible!\n";
		//cout <<  model3.isProvenInfeasible() << "  " <<  model3.isProvenOptimal() << "\n";
		//delete solver3;
		delete model3;
		return 1;
	}
	else {
		cout << "Optimal Solution found!\n";
		*solver = *solver3;
		model->gutsOfCopy(*model3, 2);
		//cout << model->getObjValue() << "   " << solver->getObjValue();
		model->setProblemStatus(0);
		model->setSecondaryStatus(0);
		//delete solver3;
		delete model3;
		return 4;
	}
}

void iis(OsiClpSolverInterface* solver) {

	OsiClpSolverInterface* solver0 = new OsiClpSolverInterface(*solver);
	int n = solver0->getNumRows();
	int m = 0;
	string temp;
	double time = CoinCpuTime();
	for (int i = 0; i < n; i++) {
		temp = solver->getRowName(i);
		OsiClpSolverInterface* solver1 = new OsiClpSolverInterface(*solver0);
		solver1->setLogLevel(0);
		const int* j = &m;
		solver1->deleteRows(1, j);
		solver1->initialSolve();
		if (solver1->isProvenPrimalInfeasible()) {
			solver0->deleteRows(1, j);
		}
		else {
			m++;
			//cout << temp << "\n";
		}
		delete j;
		delete solver1;
	}
	double time1 = CoinCpuTime() - time;
	solver0->writeLp("iis");
	cout << "\nIIS is written to iis.lp. The size is " << solver0->getNumRows() << ". It took " << time1 / 60. << " minutes.\n\nList of constraints causing infeasibility:\n\n";
	for (int j = 0; j < solver0->getNumRows(); j++)
		cout << solver0->getRowName(j) + "\n";
	delete solver0;

}

double max_var_violation = 0;
double max_int_violation = 0;
double max_row_violation = 0;

std::string* var_viol_names;
double* var_viol_bound;
int var_viol_count = 0;

std::string* row_viol_names;
int* row_viol_dir;
int row_viol_count = 0;

std::string* int_var_viol_names;
int* int_var_viol_round;
int int_var_viol_count = 0;

double get_max_var_violation() {
	return max_var_violation;
}

double get_max_int_violation() {
	return max_int_violation;
}

double get_max_row_violation() {
	return max_row_violation;
}

std::string* get_var_viol_names() {
	return var_viol_names;
}

double* get_var_viol_bound() {
	return var_viol_bound;
}

int get_var_viol_count() {
	return var_viol_count;
}

std::string* get_row_viol_names() {
	return row_viol_names;
}

int* get_row_viol_dir() {
	return row_viol_dir;
}

int get_row_viol_count() {
	return row_viol_count;
}

std::string* get_int_var_viol_names() {
	return int_var_viol_names;
}

int* get_int_var_viol_round() {
	return int_var_viol_round;
}

int get_int_var_viol_count() {
	return int_var_viol_count;
}

void clean() {

	max_int_violation = 0;
	max_var_violation = 0;
	max_row_violation = 0;


	if (var_viol_count > 0){ 
		delete[] var_viol_names; delete[] var_viol_bound;
	}

	if (row_viol_count > 0) {
		delete[] row_viol_names; delete[] row_viol_dir;
	}

	if (int_var_viol_count > 0) {
		delete[] int_var_viol_names; delete[] int_var_viol_round;
	}
	
	int_var_viol_count = 0;
	var_viol_count = 0;
	row_viol_count = 0;
}



//  0: no violations
// 100: integer violation
//  10: var violation
//   1: row violation

int Y2(CbcModel* model, OsiClpSolverInterface* solver, int write_to_file = 0, const char* address = NULL, int terminal_output = 1, double threshold = 0, double threshold_int = 0) {
	
	//cout << " Y2 called \n";

	clean();

	//OsiSolverInterface* solver = solver_;

	int return_val = 0;
	int n = solver->getNumCols();
	double* solution_new = model->bestSolution();

	std::list<string> var_viol_names1 = {};
	std::list<double> var_viol_bound1 = {};
	std::list<string> row_viol_names1 = {};
	std::list<int> row_viol_dir1 = {};
	std::list<string> int_var_viol_names1 = {};
	std::list<int> int_var_viol_round1 = {};

	ofstream myfile;
	if (write_to_file) myfile.open(address);

	int numberColumns = solver->getNumCols();
	int numberRows = solver->getNumRows();

	const double* columnLower = solver->getColLower(); // Column Bounds
	const double* columnUpper = solver->getColUpper();
	const double* rowLower = solver->getRowLower();
	const double* rowUpper = solver->getRowUpper();
	const double* objective = solver->getObjCoefficients(); // In code we also use min/max
	//const double integerTolerance = model->getIntegerTolerance();
	//double primalTolerance;
	//model->solver()->getDblParam(OsiPrimalTolerance, primalTolerance);

	const CoinPackedMatrix* matrix_ = model->getMatrixByCol();
	const double* element = matrix_->getElements();
	const int* row = matrix_->getIndices();
	const CoinBigIndex* columnStart = matrix_->getVectorStarts();
	const int* columnLength = matrix_->getVectorLengths();

	double* rowActivity = new double[numberRows];
	for (int i = 0; i < numberRows; i++) {
		rowActivity[i] = 0;
	}
	double solutionValue = 0;
	for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
		CoinBigIndex j;
		double value = solution_new[iColumn];

		if (value) {
			double cost = objective[iColumn];

			solutionValue += value * cost;

			for (j = columnStart[iColumn];
				j < columnStart[iColumn] + columnLength[iColumn]; j++) {
				int iRow = row[j];
				rowActivity[iRow] += value * element[j];
			}
		}
	}



	cout.precision(ceil(max(-log10(threshold), -log10(threshold_int))) + 2);
	if (write_to_file) myfile.precision(ceil(max(-log10(threshold), -log10(threshold_int))) + 2);
	int row_violation = 0;

	for (int iRow = 0; iRow < numberRows; iRow++) {
		if (rowActivity[iRow] < rowLower[iRow] - threshold || rowActivity[iRow] > rowUpper[iRow] + threshold) {
			row_violation = 1;
			if (rowActivity[iRow] < rowLower[iRow] - threshold) {
				max_row_violation = max(max_row_violation, std::abs(rowActivity[iRow] - rowLower[iRow]));
			}
			else {
				max_row_violation = max(max_row_violation, std::abs(rowActivity[iRow] - rowUpper[iRow]));
			}

			row_viol_names1.push_back(solver->getRowName(iRow));
			row_viol_dir1.push_back((rowActivity[iRow] < rowLower[iRow] - threshold) ? -1 : 1);
			row_viol_count++;
			if (terminal_output) cout << "Row " << iRow << "(" << solver->getRowName(iRow) << ") with activity " << rowActivity[iRow] << " violates one of LB (" << rowLower[iRow] << ") or UB (" << rowUpper[iRow] << ")\n";
			if (write_to_file) myfile << "Row " << iRow << "(" << solver->getRowName(iRow) << ") with activity " << rowActivity[iRow] << " violates one of LB (" << rowLower[iRow] << ") or UB (" << rowUpper[iRow] << ")\n";
		}
	}

	int col_violation = 0;
	int int_violation = 0;
	int rounded = -999;

	for (int iCol = 0; iCol < numberColumns; iCol++) {
		//if (solver->getColName(iCol) == "ir_64_xa_gp") {
		//	cout << "!!! " << solver->getColName(iCol) << "\n";
		//}
		if (solver->isInteger(iCol)) {
			rounded = round(solution_new[iCol]);
			if ( std::abs(rounded - solution_new[iCol]) > threshold_int) {
				int_violation = 1;
				max_int_violation = max(max_int_violation, std::abs(round(solution_new[iCol]) - solution_new[iCol]));
				int_var_viol_names1.push_back(solver->getColName(iCol));
				int_var_viol_round1.push_back(rounded);
				int_var_viol_count++;
				if (terminal_output) cout << std::setprecision(12) << "Integer Variable" << iCol << " (" << solver->getColName(iCol) << ") with value " << solution_new[iCol] << " violates threshold_int " << threshold_int << "\n";
				if (write_to_file) myfile << "Integer Variable" << iCol << " (" << solver->getColName(iCol) << ") with value " << solution_new[iCol] << " violates threshold_int " << threshold_int << "\n";
			}
			else {
				solution_new[iCol] = rounded;
			}
		}
		else { // non integer var
			
			if (   solution_new[iCol] < columnLower[iCol] - threshold 
				|| solution_new[iCol] > columnUpper[iCol] + threshold) {

				col_violation = 1;
				if (solution_new[iCol] < columnLower[iCol]) {
					max_var_violation = max(max_var_violation, std::abs(solution_new[iCol] - columnLower[iCol]));
				}
				else if (solution_new[iCol] > columnUpper[iCol]) {
					max_var_violation = max(max_var_violation, std::abs(solution_new[iCol] - columnUpper[iCol]));
				}
				var_viol_names1.push_back(solver->getColName(iCol));
				var_viol_bound1.push_back((solution_new[iCol] < columnLower[iCol]) ? columnLower[iCol] : columnUpper[iCol]);
				var_viol_count++;
				if (terminal_output) cout << "Variable " << iCol << " (" << solver->getColName(iCol) << ") with value " << solution_new[iCol] << " violates one of LB (" << columnLower[iCol] << ") or UB(" << columnUpper[iCol] << ")\n";
				if (write_to_file) myfile << "Variable " << iCol << " (" << solver->getColName(iCol) << ") with value " << solution_new[iCol] << " violates one of LB (" << columnLower[iCol] << ") or UB(" << columnUpper[iCol] << ")\n";

			}
			else if (solution_new[iCol] < columnLower[iCol]) {
				//cout << "Variable: " << iCol << " adjusted to " << columnLower[iCol] << "\n";
				solution_new[iCol] = columnLower[iCol];

			}
			else if (solution_new[iCol] > columnUpper[iCol]) {

				solution_new[iCol] = columnUpper[iCol];

			}
		}
	}

	//cout << row_feasible << "   " << col_feasible << "   " << int_feasible << "   " << primalTolerance << "         " << integerTolerance  << "     " << max_viol << "\n";

	//bool all_pass = row_feasible && col_feasible && int_feasible;
	//bool threshold_ok = max_row_violation <= threshold && max_var_violation <= threshold && max_int_violation <= int_threshold;
	return_val = 100 * int_violation + 10 * col_violation + row_violation;

	if (return_val == 0) {
		if (terminal_output) cout << "max violation < threshold!\n";
		if (write_to_file) myfile << "max violation < threshold!\n";
		if (write_to_file) myfile.close();
		return return_val;
	}
	else {

		if (col_violation) {
			var_viol_names = new string[var_viol_count];
			var_viol_bound = new double[var_viol_count];
			int i = 0;
			for (std::list<string>::iterator it = var_viol_names1.begin(); it != var_viol_names1.end(); it++) {
				var_viol_names[i] = *it;
				i++;
			}
			i = 0;
			for (std::list<double>::iterator it = var_viol_bound1.begin(); it != var_viol_bound1.end(); it++) {
				var_viol_bound[i] = *it;
				i++;
			}
		}

		if (row_violation) {
			row_viol_names = new string[row_viol_count];
			row_viol_dir = new int[row_viol_count];
			int i = 0;
			for (std::list<string>::iterator it = row_viol_names1.begin(); it != row_viol_names1.end(); it++) {
				row_viol_names[i] = *it;
				i++;
			}
			i = 0;
			for (std::list<int>::iterator it = row_viol_dir1.begin(); it != row_viol_dir1.end(); it++) {
				row_viol_dir[i] = *it;
				i++;
			}
		}

		if (int_violation) {
			int_var_viol_names = new string[int_var_viol_count];
			int_var_viol_round = new int[int_var_viol_count];
			int i = 0;
			for (std::list<string>::iterator it = int_var_viol_names1.begin(); it != int_var_viol_names1.end(); it++) {
				int_var_viol_names[i] = *it;
				i++;
			}
			i = 0;
			for (std::list<int>::iterator it = int_var_viol_round1.begin(); it != int_var_viol_round1.end(); it++) {
				int_var_viol_round[i] = *it;
				i++;
			}
		}

		if (write_to_file) {
			myfile << "Max int violation: " << max_int_violation << " Max var violation:" << max_var_violation << " Max row violation:" << max_row_violation << "\n";
			myfile.close();

		}

	}

	return return_val;
}

int Y2_simple(CbcModel* model, OsiClpSolverInterface* solver) {

	const double integerTolerance = model->getIntegerTolerance();
	double primalTolerance;
	model->solver()->getDblParam(OsiPrimalTolerance, primalTolerance);

	return Y2(model, solver, 0, NULL, 0, primalTolerance * 100, integerTolerance * 10);
	//return Y2(model, solver, 0, NULL, 0, primalTolerance, integerTolerance);
}