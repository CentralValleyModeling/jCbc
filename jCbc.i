%module jCbc
%{
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
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





/* Put header files here or function declarations like below */
extern int addCol(CoinModel *build, double collb, double colub, double obj,const char *name, bool isInt);
extern void addRows(CbcModel *model, CoinModel *build);
extern int addRow(CoinModel *build, int numberInRow, long int index [], double values [], double rowlb, double rowup, const char *name);
extern void branchAndBound(CbcModel * model);
extern void readMps(CbcModel *model, const char *name);
extern void readLp(CbcModel *model, const char *name);
extern void writeMps(CbcModel *model, const char *name);
extern void writeLp(CbcModel *model, const char *name);
extern const double * getSol(CbcModel *model);
extern void setInteger(CbcModel *model, int i);
extern int addRows(OsiClpSolverInterface *solver, CoinModel *build);
extern void assignSolver(CbcModel *model, OsiClpSolverInterface *solver);
extern void readMps(OsiClpSolverInterface *solver, const char *name);
extern int readLp(OsiClpSolverInterface *solver, const char *name);
extern void setInteger(OsiClpSolverInterface *solver, int i);
extern void initialSolve(OsiClpSolverInterface *solver);
extern void setLogLevel(OsiClpSolverInterface *solver, int i);
extern void setLogLevel(CbcModel * model, int i);
extern void solve(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0);
extern void solve_1(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0);
extern int solve_2(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0);
extern int solve_3(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0,double presolve_tolerance=1e-07);
extern void writeMps(OsiClpSolverInterface *solver, const char *name);
extern void writeLp(OsiClpSolverInterface *solver, const char *name);
extern int isInteger(OsiClpSolverInterface *solver, int i);
extern int isInteger(CbcModel * model, int i);
extern std::string getRowName(OsiClpSolverInterface *solver, int i);
extern std::string getColName(OsiClpSolverInterface *solver, int i);
extern std::string getRowName(CbcModel *model, int i);
extern std::string getColName(CbcModel *model, int i);
extern std::string getVersion();
extern void setWhsScaling(bool solvWhs_scaling);
extern void setWhsSafe(bool solvWhs_safe);
extern int getNumRows(OsiClpSolverInterface *solver);
extern int getNumCols(OsiClpSolverInterface *solver);
extern int getNumRows(CbcModel *model);
extern int getNumCols(CbcModel *model);
extern const double * getColSolution(OsiClpSolverInterface *solver);
extern const double * getRowPrice(OsiClpSolverInterface *solver);
extern const double * getRowActivity(OsiClpSolverInterface *solver);
extern const double * getReducedCost(OsiClpSolverInterface *solver);
extern const double * getColSolution(CbcModel *model);
extern const double * getRowPrice(CbcModel *model);
extern const double * getRowActivity(CbcModel *model);
extern const double * getReducedCost(CbcModel *model);
extern int  numberIntegers(CbcModel *model);
extern const double getObjValue(CbcModel *model);
extern int status(CbcModel *model);
extern int isProvenOptimal(CbcModel *model);
extern int isProvenInfeasible(CbcModel *model);
extern void setModelName(OsiClpSolverInterface *solver, std::string name);
extern std::string getModelName(OsiClpSolverInterface *solver);
extern int isBinary(OsiClpSolverInterface *solver, int i);
extern int secondaryStatus(CbcModel *model);
extern double getCoinCpuTime();
extern const double getObjValue(OsiClpSolverInterface *solver);
extern void loadProblem(OsiClpSolverInterface *solver, CoinPackedMatrix *byRow, double columnLower[], double columnUpper[], 
			double objective[], double rowLower[], double rowUpper[]);
extern CoinPackedMatrix * createMatrix( int numberColumns, int numberRows, int numberElements ,
								double elementByRow[], long int column[], long int rowStart[]);

extern void setRowName(OsiClpSolverInterface *solver, int i , const char *name);
extern void setPrimalTolerance(CbcModel *model, double a);
extern void setDualTolerance(CbcModel *model, double a);
extern void setIntegerTolerance(CbcModel *model, double a);
extern void test(std::pair<std::string, int> a);
extern void test1(std::vector<std::pair<std::string, int> > a);
extern int solve_whs(CbcModel *model, OsiClpSolverInterface *solver, std::string names[], int values[],int intvars,  int logLevel= 0,double presolve_tolerance=1e-07);
extern int getNumIntegers(OsiClpSolverInterface *solver);
extern CoinWarmStart * getWarmStart(CbcModel *model);
extern void writeLp1(CbcModel *model, const char *name, double epsilon=1e-5, int decimals=5 );
extern void writeMps1(CbcModel *model, const char *name, int 	formatType = 0, int numberAcross = 2,double objSense = 0.0);
extern int solve_4(CbcModel *model, OsiClpSolverInterface *solver, int logLevel);
extern int callCbcJ(std::string a,CbcModel *model, OsiClpSolverInterface *solver);
extern double minCoeff(CbcModel *model);
extern double maxCoeff(CbcModel *model);
//extern int callCbc(std::string a,CbcModel *model);
extern double minRHS(CbcModel *model);
extern void cutoff(CbcModel *model, int n);
extern int solve_unified(CbcModel *model, OsiClpSolverInterface *solver, std::string names[] = NULL, int values[] = NULL, int intvars = 0, int logLevel= 0 );
extern void iis(OsiClpSolverInterface *solver);
extern std::string* get_var_viol_names() ;
extern double* get_var_viol_bound() ;
extern int get_var_viol_count() ;
extern std::string* get_row_viol_names() ;
extern int* get_row_viol_dir();
extern int get_row_viol_count() ;
extern std::string* get_int_var_viol_names() ;
extern int* get_int_var_viol_round();
extern int get_int_var_viol_count() ;
extern void clean();
//extern double Y(CbcModel* model, OsiClpSolverInterface* solver, double solution[], std::string names[], int n, int write_to_file = 0, int terminal_output = 1, double threshold = 0);

extern int Y2(CbcModel* model, OsiClpSolverInterface* solver, int write_to_file = 0, const char* address = NULL, int terminal_output = 1, double threshold = 0, double int_threshold = 0);
extern int Y2_simple(CbcModel* model, OsiClpSolverInterface* solver);
extern double get_max_var_violation();
extern double get_max_int_violation();
extern double get_max_row_violation();
//extern int alt_opt(CbcModel *model, OsiClpSolverInterface *solver, int priorities[] = NULL);
//extern double* get_sol(int i); 
//extern void clean_sol();
%}

%include cpointer.i
%include "typemaps.i"
#include <string>
%include carrays.i
%include <std_pair.i>
%include <std_vector.i>
%include <std_string.i>


%array_functions( double, jarray_double);
%array_functions( int, jarray_int);
%array_functions(std::string, jarray_string);
%pointer_functions(CbcModel, jCbcModel);
%pointer_functions(CoinModel, jCoinModel);
%pointer_functions(OsiClpSolverInterface, jOsiClpSolverInterface);
%pointer_functions(CoinPackedMatrix, jCoinPackedMatrix);


%typemap(jtype) double values[] "double[]"
%typemap(jstype) double values[] "double[]"
%typemap(javain) double values[] "$javainput"
%typemap(jni) double values[] "jdoubleArray"
%typemap(in) double values[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double values[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) long int index[] "int[]"
%typemap(jstype) long int index[] "int[]"
%typemap(javain) long int index[] "$javainput"
%typemap(jni) long int index[] "jintArray"
%typemap(in) long int index[] {
  jboolean isCopy;
  $1 = JCALL2(GetIntArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) long int index[] {
  JCALL3(ReleaseIntArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) double elementByRow[] "double[]"
%typemap(jstype) double elementByRow[] "double[]"
%typemap(javain) double elementByRow[] "$javainput"
%typemap(jni) double elementByRow[] "jdoubleArray"
%typemap(in) double elementByRow[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double elementByRow[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) long int column[] "int[]"
%typemap(jstype) long int column[] "int[]"
%typemap(javain) long int column[] "$javainput"
%typemap(jni) long int column[] "jintArray"
%typemap(in) long int column[] {
  jboolean isCopy;
  $1 = JCALL2(GetIntArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) long int column[] {
  JCALL3(ReleaseIntArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) long int rowStart[] "int[]"
%typemap(jstype) long int rowStart[] "int[]"
%typemap(javain) long int rowStart[] "$javainput"
%typemap(jni) long int rowStart[] "jintArray"
%typemap(in) long int rowStart[] {
  jboolean isCopy;
  $1 = JCALL2(GetIntArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) long int rowStart[] {
  JCALL3(ReleaseIntArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) double objective[] "double[]"
%typemap(jstype) double objective[] "double[]"
%typemap(javain) double objective[] "$javainput"
%typemap(jni) double objective[] "jdoubleArray"
%typemap(in) double objective[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double objective[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}


%typemap(jtype) double columnUpper[] "double[]"
%typemap(jstype) double columnUpper[] "double[]"
%typemap(javain) double columnUpper[] "$javainput"
%typemap(jni) double columnUpper[] "jdoubleArray"
%typemap(in) double columnUpper[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double columnUpper[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}


%typemap(jtype) double columnLower[] "double[]"
%typemap(jstype) double columnLower[] "double[]"
%typemap(javain) double columnLower[] "$javainput"
%typemap(jni) double columnLower[] "jdoubleArray"
%typemap(in) double columnLower[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double columnLower[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) double rowLower[] "double[]"
%typemap(jstype) double rowLower[] "double[]"
%typemap(javain) double rowLower[] "$javainput"
%typemap(jni) double rowLower[] "jdoubleArray"
%typemap(in) double rowLower[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double rowLower[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}

%typemap(jtype) double rowUpper[] "double[]"
%typemap(jstype) double rowUpper[] "double[]"
%typemap(javain) double rowUpper[] "$javainput"
%typemap(jni) double rowUpper[] "jdoubleArray"
%typemap(in) double rowUpper[] {
  jboolean isCopy;
  $1 = JCALL2(GetDoubleArrayElements, jenv, $input, &isCopy);
}
%typemap(freearg) double rowUpper[] {
  JCALL3(ReleaseDoubleArrayElements, jenv, $input, $1, 0);
}



extern int addCol(CoinModel *build, double collb, double colub, double obj,const char *name, bool isInt);
extern void addRows(CbcModel *model, CoinModel *build);
extern int addRow(CoinModel *build, int numberInRow, long int index [], double values [], double rowlb, double rowup, const char *name);
extern void branchAndBound(CbcModel * model);
extern void readMps(CbcModel *model, const char *name);
extern void readLp(CbcModel *model, const char *name);
extern void writeMps(CbcModel *model, const char *name);
extern void writeLp(CbcModel *model, const char *name);
extern const double * getSol(CbcModel *model);
extern void setInteger(CbcModel *model, int i);
extern int addRows(OsiClpSolverInterface *solver, CoinModel *build);
extern void assignSolver(CbcModel *model, OsiClpSolverInterface *solver);
extern void readMps(OsiClpSolverInterface *solver, const char *name);
extern int readLp(OsiClpSolverInterface *solver, const char *name);
extern void setInteger(OsiClpSolverInterface *solver, int i);
extern void initialSolve(OsiClpSolverInterface *solver);
extern void setLogLevel(OsiClpSolverInterface *solver, int i);
extern void setLogLevel(CbcModel * model, int i);
extern void writeMps(OsiClpSolverInterface *solver, const char *name);
extern void writeLp(OsiClpSolverInterface *solver, const char *name);
extern void solve(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0);
extern void solve_1(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0);
extern int solve_2(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0);
extern int solve_3(CbcModel *model, OsiClpSolverInterface *solver, int logLevel= 0,double presolve_tolerance=1e-07);
extern int isInteger(OsiClpSolverInterface *solver, int i);
extern int isInteger(CbcModel * model, int i);
extern std::string getRowName(OsiClpSolverInterface *solver, int i);
extern std::string getColName(OsiClpSolverInterface *solver, int i);
extern std::string getRowName(CbcModel *model, int i);
extern std::string getColName(CbcModel *model, int i);
extern std::string getVersion();
extern void setWhsScaling(bool solvWhs_scaling);
extern void setWhsSafe(bool solvWhs_safe);
extern int getNumRows(OsiClpSolverInterface *solver);
extern int getNumCols(OsiClpSolverInterface *solver);
extern int getNumRows(CbcModel *model);
extern int getNumCols(CbcModel *model);
extern const double * getColSolution(OsiClpSolverInterface *solver);
extern const double * getRowPrice(OsiClpSolverInterface *solver);
extern const double * getRowActivity(OsiClpSolverInterface *solver);
extern const double * getReducedCost(OsiClpSolverInterface *solver);
extern const double * getColSolution(CbcModel *model);
extern const double * getRowPrice(CbcModel *model);
extern const double * getRowActivity(CbcModel *model);
extern const double * getReducedCost(CbcModel *model);
extern int  numberIntegers(CbcModel *model);
extern const double getObjValue(CbcModel *model);
extern int status(CbcModel *model);
extern int isProvenOptimal(CbcModel *model);
extern int isProvenInfeasible(CbcModel *model);
extern void setModelName(OsiClpSolverInterface *solver, std::string name);
extern std::string getModelName(OsiClpSolverInterface *solver);
extern int isBinary(OsiClpSolverInterface *solver, int i);
extern int secondaryStatus(CbcModel *model);
extern double getCoinCpuTime();
extern const double getObjValue(OsiClpSolverInterface *solver);
extern void loadProblem(OsiClpSolverInterface *solver, CoinPackedMatrix *byRow, double columnLower[], double columnUpper[], 
			double objective[], double rowLower[], double rowUpper[]);
extern CoinPackedMatrix * createMatrix( int numberColumns, int numberRows, int numberElements , double elementByRow[], long int column[], long int rowStart[]);
extern void setRowName(OsiClpSolverInterface *solver, int i , const char *name);
extern void setPrimalTolerance(CbcModel *model, double a);
extern void setDualTolerance(CbcModel *model, double a);
extern void setIntegerTolerance(CbcModel *model, double a);
extern void test(std::pair<std::string, int> a);
extern void test1(std::vector<std::pair<std::string, int> > a);
extern int solve_whs(CbcModel *model, OsiClpSolverInterface *solver, std::string names[], int values[],int intvars, int logLevel= 0,double presolve_tolerance=1e-07);
extern int getNumIntegers(OsiClpSolverInterface *solver);
extern CoinWarmStart * getWarmStart(CbcModel *model);
extern void writeLp1(CbcModel *model, const char *name, double epsilon=1e-5, int decimals=5 );
extern void writeMps1(CbcModel *model, const char *name, int 	formatType = 0, int numberAcross = 2,double objSense = 0.0);
extern int solve_4(CbcModel *model, OsiClpSolverInterface *solver, int logLevel);
extern int callCbcJ(std::string a,CbcModel *model, OsiClpSolverInterface *solver);
extern double minCoeff(CbcModel *model);
extern double maxCoeff(CbcModel *model);
//extern int callCbc(std::string a,CbcModel *model);
extern double minRHS(CbcModel *model);
extern void cutoff(CbcModel *model, int n);
extern int solve_unified(CbcModel *model, OsiClpSolverInterface *solver, std::string names[] = NULL, int values[] = NULL, int intvars = 0, int logLevel= 0 );
extern void iis(OsiClpSolverInterface *solver);
extern std::string* get_var_viol_names() ;
extern double* get_var_viol_bound() ;
extern int get_var_viol_count() ;
extern std::string* get_row_viol_names() ;
extern int* get_row_viol_dir();
extern int get_row_viol_count() ;
extern std::string* get_int_var_viol_names() ;
extern int* get_int_var_viol_round();
extern int get_int_var_viol_count() ;
extern void clean();
//extern double Y(CbcModel* model, OsiClpSolverInterface* solver, double solution[], std::string names[], int n, int write_to_file = 0, int terminal_output = 1, double threshold = 0);

extern int Y2(CbcModel* model, OsiClpSolverInterface* solver, int write_to_file = 0, const char* address = NULL, int terminal_output = 1, double threshold = 0, double int_threshold = 0);
extern int Y2_simple(CbcModel* model, OsiClpSolverInterface* solver);
extern double get_max_var_violation();
extern double get_max_int_violation();
extern double get_max_row_violation();
//extern int alt_opt(CbcModel *model, OsiClpSolverInterface *solver, int priorities[] = NULL);
//extern double* get_sol(int i); 
//extern void clean_sol();


