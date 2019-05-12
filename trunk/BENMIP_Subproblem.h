#pragma once
using namespace std;
#include <vector>

#include "UtilParameters.h"
#include "DecompParam.h"
#include "UtilMacrosDecomp.h"
#include "DecompAlgoD.h"
class BENMIP_Subproblem 
{
public:

	BENMIP_Subproblem()
	{
	}


	//BENMIP_Subproblem::~BENMIP_Subproblem()
	//{
	//}



	vector<double*> getDualRaysOsi(int maxNumRays);
	bool isDualRayInfProof(const double* dualRay, const CoinPackedMatrix* rowMatrix, const double* colLB, const double* colUB, const double* rowRhs, ostream* os);
	bool isDualRayInfProofCpx(const double* dualRay, const CoinPackedMatrix* rowMatrix, const double* colLB, const double* colUB, const double* rowRhs, ostream* os);
	void printBasisInfo(OsiSolverInterface* si, ostream* os);
	vector<double*> getDualRaysCpx(int maxNumRays);
	vector<double*> getDualRays(int maxNumRays);


	/**
	* Store the name of the class (for logging/debugging) - "who am I?"
	*/
	std::string m_classTag;

	/**
	* Parameters
	*/
	DecompParam m_param;
	UtilParameters* m_utilParam;

	/**
	* Stream for log file (default to stdout).
	*/
	std::ostream* m_osLog;
	

	/**
	* Solver interface(s) for subproblems (P').
	*/
	//vector<OsiSolverInterface*> m_subprobSI;

	/**
	* Solver interface(s) for master problem (Q'').
	*   CPM: holds model core (and optionally relaxed)
	*        in original space
	*   PC : holds model core in reformulated space
	*/
	OsiSolverInterface* m_subproblemOSI;

	/**
	* Pointer to current active DECOMP application.
	*/
	DecompApp* m_app;
	int NumBlocks;
	double* m_objective;
	/**
	*  Model data: the core model (A'')
	*/
	DecompConstraintSet*  modelSubproblem;

	/** Original constraint matrix for the instance */

	const CoinPackedMatrix* m_matrix;
	CoinPackedMatrix*  varBounds;
	std::vector<double> varLBs, varUBs;
	std::vector<string>  varBounsNames;

	BENMIP_Subproblem(UtilParameters& utilParam) :
		m_classTag("D-APP"),
		m_osLog(&std::cout),
		NumBlocks(0),
		m_utilParam(&utilParam),
		m_objective(NULL),
		//m_modelSubproblem(utilParam),
		m_matrix(NULL)
		//m_modelC(NULL),
		//m_threadIndex(0)
	{
		//---
		//--- get application parameters
		//---
		m_param.getSettings(utilParam);

		if (m_param.LogLevel >= 1) {
			m_param.dumpSettings();
		}

		//startupLog();
	};

	/**
	* Destructor.
	*/
	virtual ~BENMIP_Subproblem() {
		UTIL_DELARR(m_objective);
		//UtilDeleteMapPtr(m_modelR);
		//UTIL_DELPTR(m_modelC);
	};
};

