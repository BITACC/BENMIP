#include "BENMIP_Subproblem.h"



vector<double*> BENMIP_Subproblem::getDualRays(int maxNumRays)
{
	if (m_param.DecompLPSolver == "CPLEX"){
		return(getDualRaysCpx(maxNumRays));
	}
	else if (m_param.DecompLPSolver == "Clp" ||
		m_param.DecompLPSolver == "Gurobi"){
		return(getDualRaysOsi(maxNumRays));
	}
	else{
		throw UtilException("Unknown solver selected.",
			"getDualRays", "DecompAlgo");
	}
}


//===========================================================================//
vector<double*> BENMIP_Subproblem::getDualRaysCpx(int maxNumRays)
{
#ifdef DIP_HAS_CPX
	bool useMultiRay = true;
	if (useMultiRay){
		OsiCpxSolverInterface* siCpx
			= dynamic_cast<OsiCpxSolverInterface*>(m_subproblemOSI);
		const int m = m_subproblemOSI->getNumRows();
		const int n = m_subproblemOSI->getNumCols();
		const double* rowRhs = m_subproblemOSI->getRightHandSide();
		const char*    rowSense = m_subproblemOSI->getRowSense();
		int r, b, c;
		vector<double*> rays;
		//Ax + Is = b
		// ax     <= b
		// ax + s  = b, s >= 0
		// ax     >= b
		// ax + s  = b, s <= 0
		UTIL_DEBUG(m_param.LogDebugLevel, 5,

			for (r = 0; r < m; r++) {
				(*m_osLog) << "Row r: " << r << " sense: " << rowSense[r]
					<< " rhs: " << rowRhs[r] << endl;
			}
		);
		m_subproblemOSI->enableSimplexInterface(false);
		double* tabRhs = new double[m];
		int*     basics = new int[m];
		double* yb = new double[m];
		double* bInvRow = new double[m];
		double* bInvARow = new double[n];
		//STOP ============================================
		//tabRhs and yb do NOT match up.... is this an issue?
		//have to hand adjust or use tabRhs since proof is based on B-1
		//which matches up with bhead - what to do in the case of CLP?
		//but, we are multiplying this by A'' later on which is based on
		//original variable space, not the one adjusted by simplex - so if
		//we return the dual ray directly from B-1 then do B-1A by hand -
		//do we have a problem?
		//need to add a check that B-1A matches my dualray.A calculation
		//in generate vars... it might be ok and yb not ok, because the
		//adjustments in simplex might only be related to rhs...
		//i don't think Osi returns tabRhs... that should be changed
		CPXgetbhead(siCpx->getEnvironmentPtr(),
			siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
			basics, tabRhs);
		//as a sanity check print out the basis status next to the yb vs tabRhs
		//calculation.... let's see why and where things don't match up...
		//yb, where y is a row of B-1 (note, can get from bhead?)
		UTIL_DEBUG(m_param.LogDebugLevel, 6,
			(*m_osLog) << "\nB-1:";

		for (r = 0; r < m; r++) {
			yb[r] = 0.0;
			m_subproblemOSI->getBInvRow(r, bInvRow);
			(*m_osLog) << "\nB-1Row r: " << r << ": " << endl;

			for (b = 0; b < m; b++) {
				yb[r] += bInvRow[b] * rowRhs[b];
				(*m_osLog) << setw(6) << "bind: "
					<< setw(4) << basics[b]
					<< setw(12) << bInvRow[b]
					<< " ["
					<< setw(12) << rowRhs[b]
					<< "] "
					<< setw(8) << " +=: "
					<< setw(12) << bInvRow[b] * rowRhs[b]
					<< setw(8) << " yb: "
					<< setw(12) << yb[r]
					<< setw(8) << " tabRhs: "
					<< setw(12) << tabRhs[r] << endl;
			}

			if (!UtilIsZero(yb[r] - tabRhs[r])) {
				(*m_osLog) << " DIFF is " << yb[r] - tabRhs[r] << endl;
			}

			assert(UtilIsZero(yb[r] - tabRhs[r], 1.0e-4));
		}
		);

		for (r = 0; r < m; r++) {
			yb[r] = 0.0;
			m_subproblemOSI->getBInvRow(r, bInvRow);

			for (b = 0; b < m; b++) {
				yb[r] += bInvRow[b] * rowRhs[b];//(B-1)_r.b
			}

			if (!UtilIsZero(yb[r] - tabRhs[r])) {
				(*m_osLog) << " DIFF is " << yb[r] - tabRhs[r] << endl;
				(*m_osLog) << "\nB-1Row r: " << r << ": basics[r]=" << basics[r]
					<< endl;
				yb[r] = 0.0;

				for (b = 0; b < m; b++) {
					if (UtilIsZero(bInvRow[b])) {
						continue;
					}

					yb[r] += bInvRow[b] * rowRhs[b];
					(*m_osLog) << setw(6) << "bind: "
						<< setw(4) << basics[b]
						<< setw(12) << bInvRow[b]
						<< " ["
						<< setw(12) << rowRhs[b];

					if (basics[b] < 0) { //== -rowIndex-1
						(*m_osLog) << " sense = " << rowSense[-(basics[b] + 1)];
					}

					(*m_osLog) << "] "
						<< setw(8) << " +=: "
						<< setw(12) << bInvRow[b] * rowRhs[b]
						<< setw(8) << " yb: "
						<< setw(12) << yb[r]
						<< setw(8) << " tabRhs: "
						<< setw(12) << tabRhs[r] << endl;
				}
			}

			//assert(UtilIsZero(yb[r] - tabRhs[r], 1.0e-4));
		}

		for (r = 0; r < m; r++) {
			if (UtilIsZero(tabRhs[r])) {
				continue;
			}

			//all pos case? if yb < 0 (then we want to minimize B-1Ax, x in P')
			//all neg case? if yb > 0 (then we want to maximize B-1Ax, x in P')
			UTIL_DEBUG(m_param.LogDebugLevel, 6,
				(*m_osLog) << "\nB-1A:";
			);

			if (tabRhs[r] > 0) {  //instead of yb
				//Ted also checks that it is a slack var here - why?
				bool allneg = true;
				m_subproblemOSI->getBInvARow(r, bInvARow);
				UTIL_DEBUG(m_param.LogDebugLevel, 6,
					(*m_osLog) << "\nB-1ARow r: " << r << ": ";
				);
				allneg = true;

				for (c = 0; c < n; c++) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << bInvARow[c] << " ";
					);

					if (bInvARow[c] >= DecompEpsilon) {
						allneg = false;
						break;
					}
				}

				if (allneg) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << " ---> allneg";
					);
					double* dualRay = new double[m];
					m_subproblemOSI->getBInvRow(r, dualRay);
					transform(dualRay, dualRay + m, dualRay, negate<double>());
					rays.push_back(dualRay);
				}
			}
			else {
				bool allpos = true;
				m_subproblemOSI->getBInvARow(r, bInvARow);
				UTIL_DEBUG(m_param.LogDebugLevel, 6,
					(*m_osLog) << "\nB-1ARow r: " << r << ": ";
				);
				allpos = true;

				for (c = 0; c < n; c++) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << bInvARow[c] << " ";
					);

					if (bInvARow[c] <= -DecompEpsilon) {
						allpos = false;
						break;
					}
				}

				if (allpos) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << " ---> allpos";
					);
					double* dualRay = new double[m];
					m_subproblemOSI->getBInvRow(r, dualRay);
					rays.push_back(dualRay);
				}
			}
		}

		UTIL_DELARR(tabRhs);
		UTIL_DELARR(basics);
		UTIL_DELARR(yb);
		UTIL_DELARR(bInvRow);
		UTIL_DELARR(bInvARow);
		m_subproblemOSI->disableSimplexInterface();
		printf("rays.size = %d\n", static_cast<int>(rays.size()));

		if (rays.size() <= 0) {
			printf("NO RAYS using standard lookup - try dualfarkas\n");
			double   proof_p;
			double* dualRay = new double[m];
			CPXdualfarkas(siCpx->getEnvironmentPtr(),
				siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
				dualRay, &proof_p);
			(*m_osLog) << "After dual farkas proof_p = " << proof_p << "\n";
			transform(dualRay, dualRay + m, dualRay, negate<double>());

			for (int i = 0; i < m; i++) {
				printf("dualRay[%d]: %g\n", i, dualRay[i]);
			}

			rays.push_back(dualRay);
		}

		//NOTE: you will have dup rays here - need to filter out...
		printf("rays.size = %d", static_cast<int>(rays.size()));

		for (size_t i = 0; i < rays.size(); i++) {
			bool isProof = isDualRayInfProof(rays[i],
				m_subproblemOSI->getMatrixByRow(),
				m_subproblemOSI->getColLower(),
				m_subproblemOSI->getColUpper(),
				m_subproblemOSI->getRightHandSide(),
				NULL);

			if (!isProof) {
				isDualRayInfProof(rays[i],
					m_subproblemOSI->getMatrixByRow(),
					m_subproblemOSI->getColLower(),
					m_subproblemOSI->getColUpper(),
					m_subproblemOSI->getRightHandSide(),
					m_osLog);
			}

			//assert(isProof);
		}

		assert(rays.size() > 0);
		return rays;
	}
	else{//useMultiRay == false
		//TEST THIS
		OsiCpxSolverInterface* siCpx
			= dynamic_cast<OsiCpxSolverInterface*>(m_subproblemOSI);
		const int m = m_subproblemOSI->getNumRows();
		const int n = m_subproblemOSI->getNumCols();
		double proof_p;
		bool   isProof;
		vector<double*> rays;
		double* ray = new double[m];
		int err
			= CPXdualfarkas(siCpx->getEnvironmentPtr(),
			siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
			ray, &proof_p);//proof_p

		if (err) {
			cerr << "CPXdualfarkas returns err " << err << endl;
			abort();
		}

		cout << "After dual farkas proof_p = " << proof_p << "\n";
		//We have to flip because in this context we want to max B-1Ax, x in P'
		double* pneg = new double[m];
		transform(ray, ray + m, pneg, negate<double>());
		rays.push_back(pneg);
#if 1
		UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
			bool isProof = isDualRayInfProof(rays[0],
			m_subproblemOSI->getMatrixByRow(),
			m_subproblemOSI->getColLower(),
			m_subproblemOSI->getColUpper(),
			m_subproblemOSI->getRightHandSide(),
			NULL);
		printf("isProof = %d\n", isProof);
		printBasisInfo(m_subproblemOSI, m_osLog);
		fflush(stdout);

		if (!isProof) {
			isDualRayInfProof(ray,
				m_subproblemOSI->getMatrixByRow(),
				m_subproblemOSI->getColLower(),
				m_subproblemOSI->getColUpper(),
				m_subproblemOSI->getRightHandSide(),
				m_osLog);
			printBasisInfo(m_subproblemOSI, m_osLog);
			fflush(stdout);
		}
		);
		assert(isDualRayInfProof(ray,
			m_subproblemOSI->getMatrixByRow(),
			m_subproblemOSI->getColLower(),
			m_subproblemOSI->getColUpper(),
			m_subproblemOSI->getRightHandSide(),
			NULL));
#endif
		return rays;
	}
#else
	throw UtilException("CPLEX function called when CPLEX is not available",
		"getDualRaysCpx", "DecompAlgo");
#endif
}

//===========================================================================//
//STOP - try this...
vector<double*> BENMIP_Subproblem::getDualRaysOsi(int maxNumRays)
{
	if (m_param.UseMultiRay){
		const int m = m_subproblemOSI->getNumRows();
		const int n = m_subproblemOSI->getNumCols();
		const double* rowRhs = m_subproblemOSI->getRightHandSide();
		const char*    rowSense = m_subproblemOSI->getRowSense();
		int i, r, b, c;
		vector<double*> rays;
		UtilPrintFuncBegin(m_osLog, m_classTag,
			"getDualRays()", m_param.LogDebugLevel, 2);
		UTIL_DEBUG(m_param.LogDebugLevel, 5,

			for (r = 0; r < m; r++) {
				(*m_osLog) << "Row r: " << r << " sense: " << rowSense[r]
					<< " rhs: " << rowRhs[r] << endl;
			}
		);
		m_subproblemOSI->enableSimplexInterface(false);
		//with simplex interface, this is slightly different...
		const double* primSolution = m_subproblemOSI->getColSolution();
		const double* rowAct = m_subproblemOSI->getRowActivity(); //==slacks?
		double* tabRhs = new double[m]; //osi_clp does not give this?
		//B-1b just equals x, but what if art column then is slack var
		int*     basics = new int[m];
		double* yb = new double[m];
		double* bInvRow = new double[m];
		double* bInvARow = new double[n];
		m_subproblemOSI->getBasics(basics);

		for (r = 0; r < m; r++) {
			i = basics[r];

			if (i < n) {
				tabRhs[r] = primSolution[i]; //should == B-1b
				//printf("tabRhs[c:%d]: %g\n", i, tabRhs[r]);
			}
			else {
				//this really should be slack vars...
				//assuming clp does Ax-Is = b, s = ax-b ??? nope...
				//tabRhs[r] = rowAct[i - n] - rowRhs[i - n];
				tabRhs[r] = rowRhs[i - n] - rowAct[i - n];
				//printf("tabRhs[r:%d]: %g [act: %g rhs: %g sense: %c]\n",
				//	i-n, tabRhs[r], rowAct[i-n], rowRhs[i-n], rowSense[i-n]);
			}
		}

		//as a sanity check print out the basis status next to the yb vs tabRhs
		//calculation.... let's see why and where things don't match up...
		//yb, where y is a row of B-1 (note, can get from bhead?)
		//B-1b is tab rhs, is this equivalent to x for struct columns?
		UTIL_DEBUG(m_param.LogDebugLevel, 6,
			(*m_osLog) << "\nB-1:";

		for (r = 0; r < m; r++) {
			if (UtilIsZero(tabRhs[r])) {
				continue;
			}

			yb[r] = 0.0;
			m_subproblemOSI->getBInvRow(r, bInvRow);
			(*m_osLog) << "\nB-1Row r: " << r << ": " << endl;

			for (b = 0; b < m; b++) {
				yb[r] += bInvRow[b] * rowRhs[b];
				(*m_osLog) << setw(6) << "bind: "
					<< setw(4) << basics[b]
					<< setw(12) << bInvRow[b]
					<< " ["
					<< setw(12) << rowRhs[b]
					<< "] "
					<< setw(8) << " +=: "
					<< setw(12) << bInvRow[b] * rowRhs[b]
					<< setw(8) << " yb: "
					<< setw(12) << yb[r]
					<< setw(8) << " tabRhs: "
					<< setw(12) << tabRhs[r]
					<< endl;
			}

			if (!UtilIsZero(yb[r] - tabRhs[r])) {
				(*m_osLog) << " DIFF is " << yb[r] - tabRhs[r] << endl;
			}

			assert(UtilIsZero(yb[r] - tabRhs[r], 1.0e-4));
		}
		);

		for (r = 0; r < m; r++) {
			if (UtilIsZero(tabRhs[r])) {
				continue;
			}

			//all pos case? if yb < 0 (then we want to minimize B-1Ax, x in P')
			//all neg case? if yb > 0 (then we want to maximize B-1Ax, x in P')
			if (tabRhs[r] > 0) { //instead of yb
				//Ted also checks that it is a slack var here - why?
				bool allneg = true;
				//not getting back slacks part here... need?
				m_subproblemOSI->getBInvARow(r, bInvARow);
				UTIL_DEBUG(m_param.LogDebugLevel, 6,
					(*m_osLog) << "B-1ARow r: " << r << ": ";
				);
				allneg = true;

				for (c = 0; c < n; c++) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << bInvARow[c] << " ";
					);

					if (bInvARow[c] >= DecompEpsilon) {
						allneg = false;
						break;
					}
				}

				if (allneg) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << " ---> allneg";
					);
					double* dualRay = new double[m];
					m_subproblemOSI->getBInvRow(r, dualRay);
					transform(dualRay, dualRay + m, dualRay, negate<double>());
					rays.push_back(dualRay);
				}
			}
			else {
				bool allpos = true;
				m_subproblemOSI->getBInvARow(r, bInvARow);
				UTIL_DEBUG(m_param.LogDebugLevel, 6,
					(*m_osLog) << "B-1ARow r: " << r << ": ";
				);
				allpos = true;

				for (c = 0; c < n; c++) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << bInvARow[c] << " ";
					);

					if (bInvARow[c] <= -DecompEpsilon) {
						allpos = false;
						break;
					}
				}

				if (allpos) {
					UTIL_DEBUG(m_param.LogDebugLevel, 6,
						(*m_osLog) << " ---> allpos";
					);
					double* dualRay = new double[m];
					m_subproblemOSI->getBInvRow(r, dualRay);
					rays.push_back(dualRay);
				}
			}

			UTIL_DEBUG(m_param.LogDebugLevel, 6,
				(*m_osLog) << endl;
			);
		}

		UTIL_DELARR(basics);
		UTIL_DELARR(yb);
		UTIL_DELARR(bInvRow);
		UTIL_DELARR(bInvARow);
		m_subproblemOSI->disableSimplexInterface();
		/*
		if(rays.size() <= 0){
		double   proof_p;
		double * dualRay = new double[m];
		CPXdualfarkas(siCpx->getEnvironmentPtr(),
		siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
		dualRay, &proof_p);
		(*m_osLog) << "After dual farkas proof_p = " << proof_p << "\n";
		transform(dualRay, dualRay + m, dualRay, negate<double>());
		for(int i = 0; i < m; i++){
		printf("dualRay[%d]: %g\n", i, dualRay[i]);
		}
		rays.push_back(dualRay);
		}
		*/
		//NOTE: you will have dup rays here - need to filter out...
		UTIL_DEBUG(m_param.LogDebugLevel, 5,
			(*m_osLog) << "Number of Rays = " << rays.size() << endl;
		);

		for (int i = 0; i < (int)rays.size(); i++) {
			bool isProof = isDualRayInfProof(rays[i],
				m_subproblemOSI->getMatrixByRow(),
				m_subproblemOSI->getColLower(),
				m_subproblemOSI->getColUpper(),
				m_subproblemOSI->getRightHandSide(),
				NULL);

			if (!isProof) {
				isDualRayInfProof(rays[i],
					m_subproblemOSI->getMatrixByRow(),
					m_subproblemOSI->getColLower(),
					m_subproblemOSI->getColUpper(),
					m_subproblemOSI->getRightHandSide(),
					m_osLog);
			}

			assert(isProof);
		}

		assert(rays.size() > 0);
		UTIL_DELARR(tabRhs);
		UtilPrintFuncEnd(m_osLog, m_classTag,
			"getDualRays()", m_param.LogDebugLevel, 2);
		return rays;
	}
	else{//m_param.UseMultiRay == false

		UtilPrintFuncBegin(m_osLog, m_classTag,
			"getDualRays()", m_param.LogDebugLevel, 2);
		vector<double*> raysT = m_subproblemOSI->getDualRays(maxNumRays);
		const double* rayT = raysT[0];
		assert(rayT);
		//stop
		//what is yb, that will tell me if i want to opt over uA or -uA
		//y^T b
		int   i;
		const CoinPackedMatrix* rowMatrix = m_subproblemOSI->getMatrixByRow();
		const double*            rowRhs = m_subproblemOSI->getRightHandSide();
		const int                m = rowMatrix->getNumRows();
		double yb = 0.0;

		for (i = 0; i < m; i++) {
			yb += rayT[i] * rowRhs[i]; //safe to use rowRhs? or flips in tab going on
		}

		(*m_osLog) << " yb = " << yb << endl;
		//need tabRhs if doing this way?
		//see Clp/examples/decompose.cpp
		//   he flips the infeasibility ray (always...)
		//---    yA >= 0, yb < 0, or  --> find a yAs <= 0 (min)
		//---    yA <= 0, yb > 0 ??   --> find a yAs >= 0 (max <--> -min)
		vector<double*> rays;

		if (yb > 0) {
			double* pneg = new double[m];
			transform(rayT, rayT + m, pneg, negate<double>());
			rays.push_back(pneg);
		}
		else {
			rays.push_back(raysT[0]);
		}

		for (int i = 0; i < m; i++)
			cout << rays[0][i] << endl;
		cout << endl;

#if 1
		//UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
		//	const double* ray = rays[0];
		//assert(ray);
		//bool isProof = isDualRayInfProof(ray,
		//	m_subproblemOSI->getMatrixByRow(),
		//	m_subproblemOSI->getColLower(),
		//	m_subproblemOSI->getColUpper(),
		//	m_subproblemOSI->getRightHandSide(),
		//	NULL);
		//printf("isProof = %d\n", isProof);
		//fflush(stdout);

		//if (!isProof) {
		//	isDualRayInfProof(ray,
		//		m_subproblemOSI->getMatrixByRow(),
		//		m_subproblemOSI->getColLower(),
		//		m_subproblemOSI->getColUpper(),
		//		m_subproblemOSI->getRightHandSide(),
		//		m_osLog);
		//	printBasisInfo(m_subproblemOSI, m_osLog);
		//	fflush(stdout);
		//}
		//assert(isDualRayInfProof(ray,
		//	m_subproblemOSI->getMatrixByRow(),
		//	m_subproblemOSI->getColLower(),
		//	m_subproblemOSI->getColUpper(),
		//	m_subproblemOSI->getRightHandSide(),
		//	NULL));
		//);;
#endif
		UtilPrintFuncEnd(m_osLog, m_classTag,
			"getDualRays()", m_param.LogDebugLevel, 2);
		return rays;
	}
}

//===========================================================================//

//===========================================================================//
bool BENMIP_Subproblem::isDualRayInfProof(const double*            dualRay,
	const CoinPackedMatrix* rowMatrix,
	const double*            colLB,
	const double*            colUB,
	const double*            rowRhs,
	ostream*                 os)
{
	//---
	//--- Does dualRay provide a proof according to Farkas Lemma?
	//---    yA >= 0, yb < 0, or
	//---    yA <= 0, yb > 0 ??
	//---
	int      i;
	double   yb;
	bool     isProof = true;
	bool     ybPos = true;
	double* yA = 0;
	const int m = rowMatrix->getNumRows();
	const int n = rowMatrix->getNumCols();
	//y^T b
	yb = 0.0;

	for (i = 0; i < m; i++) {
		yb += dualRay[i] * rowRhs[i];

		if (os) {
			(*os) << "i : " << i << " dualRay = " << dualRay[i]
				<< " rowRhs = " << rowRhs[i] << " yb = " << yb << endl;
		}
	}

	//TODO: tol
	if (yb > 1.0e-10) {
		ybPos = true;
	}
	else if (yb < -1.0e-10) {
		ybPos = false;
	}
	else {
		return isProof;
	}

	yA = new double[n];
	rowMatrix->transposeTimes(dualRay, yA);     //y^T A

	for (i = 0; i < n; i++) {
		if (os) {
			(*os) << "yA[" << i << "]:\t" << yA[i];
		}

		//TODO: tol 1.0e-6 is too tight?
		if ((ybPos && (yA[i] >  1.0e-2)) ||
			(!ybPos && (yA[i] < -1.0e-2))) {
			if (os) {
				(*os) << " -->isProof (false)" << endl;
			}

			isProof = false;
		}
		else if (os) {
			(*os) << endl;
		}
	}

	UTIL_DELARR(yA);
#if 0

	//sanity check
	if (!isProof)
		isProof
		= isDualRayInfProofCpx(dualRay, rowMatrix, colLB, colUB, rowRhs, os);

#endif
	return isProof;
}



//===========================================================================//
bool BENMIP_Subproblem::isDualRayInfProofCpx(const double*            dualRay,
	const CoinPackedMatrix* rowMatrix,
	const double*            colLB,
	const double*            colUB,
	const double*            rowRhs,
	ostream*                 os)
{
	//---
	//--- Assume:
	//---   Ax     >= b
	//---   y^T Ax >= y^T b, y >= 0 (for >=)
	//---
	//--- Let z[j] = u[j], if y^T A[j] > 0
	//---          = l[j], if y^T A[j] < 0
	//---          = arbitrary, otherwise
	//---
	//--- Then, WHY?
	//---   y^T b - y^T A z > 0 ==> contradiction
	//---
	//---  proof_p = y^T b - y^T A z > 0
	//---
	//---  So, we want to maximize y^T A x to break the proof.
	//---
	int       i, j;
	double   yb, yAz;
	double* yA = 0;
	double*   z = 0;
	const int m = rowMatrix->getNumRows();
	const int n = rowMatrix->getNumCols();
	//TODO: check for out-of-mem conditions?
	yA = new double[n];
	UtilFillN(yA, n, 0.0);
	double* yA2 = new double[n];
	rowMatrix->transposeTimes(dualRay, yA2);     //y^T A

	for (i = 0; i < m; i++) {
		double yA_i = 0;
		CoinShallowPackedVector pv = rowMatrix->getVector(i);
		const int*     indI = pv.getIndices();
		const double* elsI = pv.getElements();
		const int      lenI = pv.getNumElements();

		for (int j = 0; j < lenI; j++) {
			yA_i += dualRay[indI[j]] * elsI[j];
			printf("i: %d, j: %d, indIj: %d, elsIj: %g ray: %g yA_i: %g\n",
				i, j, indI[j], elsI[j], dualRay[indI[j]], yA_i);
		}

		yA[i] = yA_i;

		if (!UtilIsZero(yA[i] - yA2[i])) {
			printf(" ---> yA: %g, yA2: %g\n", yA[i], yA2[i]);
		}

		fflush(stdout);
		CoinAssert(UtilIsZero(yA[i] - yA2[i]));
	}

	z = new double[n];

	for (j = 0; j < n; j++) {
		if (yA[j] >= 0) {
			z[j] = CoinMin(1.0e20, colUB[j]);
		}
		else {
			z[j] = colLB[j];
		}
	}

	//y^T b
	yb = 0.0;

	for (i = 0; i < m; i++) {
		yb += dualRay[i] * rowRhs[i];

		if (os)
			(*os) << "\ni : " << i << " dualRay = " << dualRay[i]
			<< " rowRhs = " << rowRhs[i] << " yb = " << yb;
	}

	//y^T A z
	yAz = 0.0;

	for (j = 0; j < n; j++) {
		yAz += yA[j] * z[j];

		if (os)
			(*os) << "\nj : " << j << " yA = " << yA[j]
			<< " z = " << z[j] << " yAz = " << yAz;
	}

	if (os) {
		(*os) << "\nyb - yAz = " << yb - yAz << endl;
	}

	UTIL_DELARR(yA);
	UTIL_DELARR(z);

	//TODO: tol
	if (yb - yAz > 1.0e-3) {
		return true;
	}
	else {
		return false;
	}
}

//===========================================================================//


//===========================================================================//
void BENMIP_Subproblem::printBasisInfo(OsiSolverInterface* si,
	ostream*              os)
{
	int      b, r, c;
	int*     basics = 0;
	int*     rstat = 0;
	int*     cstat = 0;
	double* bInvRow = 0;
	double* bInvARow = 0;
	const int n = si->getNumCols();
	const int m = si->getNumRows();
	char type[4] = { 'F', 'B', 'U', 'L' };
	//TODO: have to check sense?
	const double* rowRhs = si->getRightHandSide();
	basics = new int[m];
	bInvRow = new double[m];
	bInvARow = new double[n];
	rstat = new int[m];
	cstat = new int[n];
	si->enableSimplexInterface(false);
	si->getBasics(basics);
	(*os) << "\n\nBasics: ";

	for (b = 0; b < m; b++) {
		(*os) << basics[b] << " ";
	}

	si->getBasisStatus(cstat, rstat);
	(*os) << "\ncstat: ";

	for (c = 0; c < n; c++) {
		(*os) << type[cstat[c]];
	}

	(*os) << "\n";
	(*os) << "rstat: ";

	for (r = 0; r < m; r++) {
		(*os) << type[rstat[r]];
	}

	(*os) << "\n";
	//yb, where y is a row of B-1
	double yb = 0.0;
	(*os) << "\nB-1:";

	for (r = 0; r < m; r++) {
		yb = 0.0;
		si->getBInvRow(r, bInvRow);
		(*os) << "\nB-1Row r: " << r << ": ";

		for (b = 0; b < m; b++) {
			(*os) << bInvRow[b] << " ";
			//rowRhs is just orig row rhs? or change based on who is basic?
			yb += bInvRow[b] * rowRhs[b];
		}

		(*os) << " ---> yb: " << yb;
	}

	//all pos case? if yb < 0
	//all neg case? if yb > 0
	//  what if yb=0?
	(*os) << "\nB-1A:";
	bool allpos = true;
	bool allneg = true;

	for (r = 0; r < m; r++) {
		si->getBInvARow(r, bInvARow);
		(*os) << "\nB-1ARow r: " << r << ": ";
		allpos = true;
		allneg = true;

		for (c = 0; c < n; c++) {
			(*os) << bInvARow[c] << " ";

			if (bInvARow[c] < 0) {
				allpos = false;
			}

			if (bInvARow[c] > 0) {
				allneg = false;
			}
		}

		if (allpos) {
			(*os) << " ---> allpos";
		}

		if (allneg) {
			(*os) << " ---> allneg";
		}
	}

	UTIL_DELARR(basics);
	UTIL_DELARR(bInvRow);
	UTIL_DELARR(bInvARow);
	UTIL_DELARR(rstat);
	UTIL_DELARR(cstat);
	si->disableSimplexInterface();
	//if you do this and want dual ray back, you need to resolve
	si->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
	si->resolve();
	si->setHintParam(OsiDoPresolveInResolve, true, OsiHintDo);
}
