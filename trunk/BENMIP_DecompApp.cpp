 //===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2015, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//

#include "Decomp.h"
#include "DecompAlgo.h"
#include "BENMIP_DecompApp.h"

//===========================================================================//
void BENMIP_DecompApp::initializeApp(UtilParameters & utilParam)  {
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---   
   m_appParam.getSettings(utilParam);   
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings();



   //---
   //--- read BENMIP instance (mps format)
   //---
   string fileName;
   if (m_appParam.DataDir != "") {
	   fileName = m_appParam.DataDir + UtilDirSlash() + m_appParam.Instance;
   } else {
      fileName = m_appParam.Instance;
   }

   m_lpIO.messageHandler()->setLogLevel(m_param.LogLpLevel);
   m_lpIO.readMps(fileName.c_str());
   int rstatus = 1;
   if(rstatus < 0){
      cerr << "Error: Filename = " << fileName << " failed to open." << endl;
      throw UtilException("I/O Error.", "initalizeApp", "BENMIP_DecompApp");
   }
   if(m_appParam.LogLevel >= 2)
      (*m_osLog) << "Objective Offset = " 
                 << UtilDblToStr(m_lpIO.objectiveOffset()) << endl;

   //---
   //--- read the master variable file
   //---
   
   //readMasterVars();


   ////---
   ////--- set best known lb/ub
   ////---
   //double offset = m_mpsIO.objectiveOffset();
   //setBestKnownLB(m_appParam.BestKnownLB + offset);
   //setBestKnownUB(m_appParam.BestKnownUB + offset);

   ////---
   ////--- read block file
   ////---
   //readBlockFile();

   //---
   //--- create models
   //---
   createModels();

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initializeApp()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void BENMIP_DecompApp::createModels(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);
   

   UTIL_DEBUG(m_param.LogDebugLevel, 2,
	   std::cout << "extracting statistics from the loaded MPS file" << endl;
   );
   //---
   //--- extract elements in the model
   //---
   const int      nbRows = m_lpIO.getNumRows();
   const int      m_nbColsOrigModel = m_lpIO.getNumCols();
   const CoinPackedMatrix * M = m_lpIO.getMatrixByRow();
   const int              * ind = M->getIndices();
   const int              * beg = M->getVectorStarts();
   const int              * len = M->getVectorLengths();
   const double           * elem = M->getElements();
   const int              * indR = NULL;
   m_nbMasterRows = 0;
   m_nbMasterCols = 0;
   m_nbMasterColsPlusEta = 0;
   m_nbSubproblemRows = 0;
   m_nbSubproblemCols = 0;
   const double * rowLB = m_lpIO.getRowLower();
   const double * rowUB = m_lpIO.getRowUpper();
   const double * colLBOrigModel = m_lpIO.getColLower();
   const double * colUBOrigModel = m_lpIO.getColUpper();
   const double * m_objCoefsOrigModel = m_lpIO.getObjCoefficients();

   m_rowNamesOrigModel = std::vector<string>();
   m_colNamesOrigModel = std::vector<string>();
   m_colsMarkerOrigModel = std::vector<int>(m_nbColsOrigModel, -2);
   m_rowsMarkerOrigModel = std::vector<int>(nbRows, -3);
   
   // for all master variables supplied by user,
   // make a map between the original index and 
   // the index in the master problem 

   std::vector<double> _UB, _LB;
   _LB.assign(colLBOrigModel, colLBOrigModel + m_nbColsOrigModel);
   _UB.assign(colUBOrigModel, colUBOrigModel + m_nbColsOrigModel);

   m_masterVarNames;
   readColumns(m_nbColsOrigModel, m_colsMarkerOrigModel, m_objCoefsOrigModel, m_colNamesOrigModel);




   UTIL_DEBUG(m_param.LogDebugLevel, 2,
	   std::cout << "Stats:" << endl;
   std::cout << "m_modelMaster.getNumInts(): " << m_modelMaster.getNumInts() << endl;
   std::cout << "m_modelMaster.getNumCols(): " << m_modelMaster.getNumCols() << endl;
   std::cout << "m_modelMaster.getNumRows(): " << m_modelMaster.getNumRows() << endl;
   std::cout << "----------------------------------------------------------" << endl;
   //std::cout << "m_modelSubproblem.getNumInts(): " << m_SP.modelSubproblem->getNumInts() << endl;
   //std::cout << "m_modelSubproblem.getNumCols(): " << m_SP.modelSubproblem->getNumCols() << endl;
   //std::cout << "m_modelSubproblem.getNumRows(): " << m_SP.modelSubproblem->getNumRows() << endl;
   std::cout << "----------------------------------------------------------" << endl;
   std::cout << "m_modelRelax.getNumInts(): " << m_modelRelax.getNumInts() << endl;
   std::cout << "m_modelRelax.getNumCols(): " << m_modelRelax.getNumCols() << endl;
   std::cout << "m_modelRelax.getNumRows(): " << m_modelRelax.getNumRows() << endl;
   std::cout << "----------------------------------------------------------" << endl;
   );

   initializeModels();

   m_SP.varBounds = new CoinPackedMatrix(false, 0.0, 0.0);

   cout << "partitioning variables of the original model. " << endl;
   for (int j = 0; j < m_nbColsOrigModel; j++)
   {
	   cout << "\r" << j << " of " << m_nbColsOrigModel;

	   m_masterOldIndexToNewIndex;
	   int _idx1 = getMPVarNewIndex(j);
	   if (_idx1 > -100){
		   assert(_idx1 > -1);
		   m_modelMaster.colLB.push_back(colLBOrigModel[j]);
		   m_modelMaster.colUB.push_back(colUBOrigModel[j]);
		   m_modelMaster.colNames.push_back(m_colNamesOrigModel[j]);
		   m_modelMaster.integerVars.push_back(_idx1);
		   m_modelMaster.integerMark.push_back('I');
	   }
	   else
	   {
		   m_subproblemOldIndexToNewIndex;
		   int idx = getSPVarNewIndex(j);
		   assert(idx > -1);

		   //*
		   // if lower bound is not 0 and the upper bound is not DecompInf
		   // add them as constraint
		   //*
		   //m_SP.modelSubproblem->colLB.push_back(colLBOrigModel[j]);
		   //m_SP.modelSubproblem->colUB.push_back(colUBOrigModel[j]);
		   bool hasUB = false;
		   bool hasLB = false;
		   if (fabs(colLBOrigModel[j]) > DecompEpsilon)
		   {
			   // add the constraint
			   std::vector<int> indexes{ getSPVarNewIndex(j) };
			   std::vector<double> coeffs{ 1 };
			   m_SP.varBounds->appendRow(indexes.size(), &(indexes[0]), &(coeffs[0]));
			   m_SP.varLBs.push_back(colLBOrigModel[j]);
			   m_SP.varUBs.push_back(DecompInf);
			   cout << colLBOrigModel[j] << "\n";
			   //m_SP.modelSubproblem->rowLB.push_back(colLBOrigModel[j]);
			   //m_SP.modelSubproblem->rowUB.push_back(DecompInf);
			   string name = "LB_" + m_colNamesOrigModel[j];
			   m_SP.varBounsNames.push_back(name);
			   hasLB = true;
			   //m_masterVariablesInMixedConstraints.push_back(BENMIP_MixedConstraints(++m_nbSubproblemRows));
		   }
		   if (fabs(colUBOrigModel[j]) < DecompInf - DecompEpsilon)
		   {
			   // add the constraint
			   std::vector<int> indexes{ getSPVarNewIndex(j) };
			   std::vector<double> coeffs{ 1 };
			   m_SP.varBounds->appendRow(indexes.size(), &(indexes[0]), &(coeffs[0]));
			   m_SP.varLBs.push_back(-DecompInf);
			   m_SP.varUBs.push_back(colUBOrigModel[j]);
			   cout << colUBOrigModel[j] << "\n";
			   //m_SP.modelSubproblem->rowLB.push_back(-DecompInf);
			   //m_SP.modelSubproblem->rowUB.push_back(colUBOrigModel[j]);
			   string name = "UB_" + m_colNamesOrigModel[j];
			   m_SP.varBounsNames.push_back(name);
			   hasUB = true;
			   //m_masterVariablesInMixedConstraints.push_back(BENMIP_MixedConstraints(++m_nbSubproblemRows));
		   }
		   {
			   if (hasUB && !hasLB){
				   m_SP.modelSubproblem->colLB.push_back(0);//colLBOrigModel[j]
				   m_SP.modelSubproblem->colUB.push_back(DecompInf);//colUBOrigModel[j]
			   }
			   else if (hasLB && !hasUB){
				   m_SP.modelSubproblem->colLB.push_back(-DecompInf);//colLBOrigModel[j]
				   m_SP.modelSubproblem->colUB.push_back(DecompInf);//colUBOrigModel[j]
			   }
			   else if (!hasLB && !hasUB){
				   m_SP.modelSubproblem->colLB.push_back(0);//colLBOrigModel[j]
				   m_SP.modelSubproblem->colUB.push_back(DecompInf);//colUBOrigModel[j]
			   }
			   else if (hasLB && hasUB){
				   m_SP.modelSubproblem->colLB.push_back(-DecompInf);//colLBOrigModel[j]
				   m_SP.modelSubproblem->colUB.push_back(DecompInf);//colUBOrigModel[j]
			   }
			   m_SP.modelSubproblem->colNames.push_back(m_colNamesOrigModel[j]);
			   m_SP.modelSubproblem->integerMark.push_back('C');
		   }
	   }
   }
   cout << endl;
   if (m_SP.modelSubproblem->colNames.size() < 1){
	   std::cerr << "BENMIP has not yet implemented IP decomposition" << endl;
	   exit(1);
   }



   m_modelMaster.colUB.push_back(DecompInf);
   m_modelMaster.colLB.push_back(-DecompInf);
   m_modelMaster.colNames.push_back("eta");
   m_masterVarNames.push_back("eta");

   m_SP.modelSubproblem->colUB.push_back(0);
   m_SP.modelSubproblem->colLB.push_back(0);
   m_SP.modelSubproblem->colNames.push_back("delta");




   //m_modelRandCore.M->submatrixOf(*m_mpsIO.getMatrixByRow(),
   // nRowsCore, rowsCore);

   m_indexOfIntegersVarsInMixedConstraints = std::vector<std::vector<int> >(nbRows);

   readRows(nbRows, m_rowNamesOrigModel, M, beg, len, ind, m_rowsMarkerOrigModel, m_colNamesOrigModel, rowLB, rowUB);

   m_SP.m_osLog = m_osLog;

   m_modelMaster.rowLB.push_back(0);//-1000
   m_modelMaster.rowUB.push_back(DecompInf);
   std::vector<int > inds(1, m_nbMasterColsPlusEta - 1);
   std::vector<double > elms(1, 1);
   m_modelMaster.M->appendRow(1, &inds[0], &elms[0]);






   UTIL_DEBUG(m_param.LogDebugLevel, 1,
	   std::cout << "Stats:" << endl;
	   std::cout << "m_modelMaster.getNumInts(): " << m_modelMaster.getNumInts() << endl;
	   std::cout << "m_modelMaster.getNumCols(): " << m_modelMaster.getNumCols() << endl;
	   std::cout << "m_modelMaster.getNumRows(): " << m_modelMaster.getNumRows() << endl;
	   std::cout << "----------------------------------------------------------" << endl;
	   //std::cout << "m_modelSubproblem.getNumInts(): " << m_SP.modelSubproblem->getNumInts() << endl;
	   std::cout << "m_modelSubproblem.getNumCols(): " << m_SP.modelSubproblem->getNumCols() << endl;
	   std::cout << "m_modelSubproblem.getNumRows(): " << m_SP.modelSubproblem->getNumRows() << endl;
	   std::cout << "----------------------------------------------------------" << endl;
	   std::cout << "m_modelRelax.getNumInts(): " << m_modelRelax.getNumInts() << endl;
	   std::cout << "m_modelRelax.getNumCols(): " << m_modelRelax.getNumCols() << endl;
	   std::cout << "m_modelRelax.getNumRows(): " << m_modelRelax.getNumRows() << endl;
	   std::cout << "----------------------------------------------------------" << endl;
   );


   // create modelRelax composed of eta>=0
   m_modelRelax.colLB = m_modelMaster.colLB;
   m_modelRelax.colUB = m_modelMaster.colUB;
   m_modelRelax.colNames = m_modelMaster.colNames;
   m_modelRelax.integerVars = m_modelMaster.integerVars;
   m_modelRelax.integerMark = m_modelMaster.integerMark;

   m_modelRelax.rowLB.push_back(0);
   m_modelRelax.rowUB.push_back(DecompInf);
   m_modelRelax.M->appendRow(1, &inds[0], &elms[0]);



   int nCols = m_SP.modelSubproblem->getNumCols();



   //====================================================

   // add the bound constraints
   for (int i = 0; i < m_SP.varLBs.size(); i++)
   {
	   //int n = m_SP.varBounds->getVector(i).getNumElements();
	   //const int * ids = m_SP.varBounds->getVector(i).getIndices();
	   //const double * els = m_SP.varBounds->getVector(i).getElements();
	   //std::vector<int> idss;
	   //std::vector<double> elss;
	   //idss.assign(ids, ids + n);
	   //elss.assign(els, els + n);
	   //cout << i << "\t" << idss[0] << "\t" << elss[0] << endl;

	   m_SP.modelSubproblem->M->appendRow(
		   m_SP.varBounds->getVector(i).getNumElements(),
		   m_SP.varBounds->getVector(i).getIndices(),
		   m_SP.varBounds->getVector(i).getElements()
		   );
   }
   m_SP.modelSubproblem->rowLB.insert(m_SP.modelSubproblem->rowLB.end(), m_SP.varLBs.begin(), m_SP.varLBs.end());
   m_SP.modelSubproblem->rowUB.insert(m_SP.modelSubproblem->rowUB.end(), m_SP.varUBs.begin(), m_SP.varUBs.end());
   m_SP.modelSubproblem->rowNames.insert(m_SP.modelSubproblem->rowNames.end(), m_SP.varBounsNames.begin(), m_SP.varBounsNames.end());
   
   nCols = m_SP.modelSubproblem->getNumCols();


	   //====================================================


   m_SP.modelSubproblem;
   m_modelRelax;
   m_modelMaster;
  

   UTIL_DEBUG(m_param.LogDebugLevel, 2,
	   std::cout << "Stats:" << endl;
   std::cout << "m_modelMaster.getNumInts(): " << m_modelMaster.getNumInts() << endl;
   std::cout << "m_modelMaster.getNumCols(): " << m_modelMaster.getNumCols() << endl;
   std::cout << "m_modelMaster.getNumRows(): " << m_modelMaster.getNumRows() << endl;
   std::cout << "----------------------------------------------------------" << endl;
   //std::cout << "m_modelSubproblem.getNumInts(): " << m_SP.modelSubproblem->getNumInts() << endl;
   std::cout << "m_modelSubproblem.getNumCols(): " << m_SP.modelSubproblem->getNumCols() << endl;
   std::cout << "m_modelSubproblem.getNumRows(): " << m_SP.modelSubproblem->getNumRows() << endl;
   std::cout << "----------------------------------------------------------" << endl;
   std::cout << "m_modelRelax.getNumInts(): " << m_modelRelax.getNumInts() << endl;
   std::cout << "m_modelRelax.getNumCols(): " << m_modelRelax.getNumCols() << endl;
   std::cout << "m_modelRelax.getNumRows(): " << m_modelRelax.getNumRows() << endl;
   std::cout << "----------------------------------------------------------" << endl;
   );




   //---
   //--- Construct the objective function of Master.
   //---
   m_objective = new double[m_nbMasterColsPlusEta];
   if (!m_objective)
	   throw UtilExceptionMemory("createModels", "BENMIP_DecompApp");
  
   //memcpy(m_objective, 
   //       m_mpsIO.getObjCoefficients(), nCols * sizeof(double));
   if (m_appParam.ObjectiveSense == 1){
	   for (int i = 0; i < m_nbMasterCols ; i++)
	   {
		   m_objective[i] = m_objCoefsOrigModel[getMPVarOldIndex(i)];
		   UTIL_DEBUG(m_param.LogDebugLevel, 5,
			   std::cout << m_masterVarNames[i] << "\t" << m_objective[i] << endl;
		   )
		   
	   }

	   m_objective[m_nbMasterColsPlusEta - 1] = 1;
   }
  else
	   if (m_appParam.ObjectiveSense == -1){
		   for (int i = 0; i < m_nbMasterColsPlusEta; i++)
			   m_objective[i] *= -1;
	   }
   // for eta
   setModelObjective(&m_objective[0], m_nbMasterColsPlusEta);


   //for (int i = 0; i < m_nbMasterColsPlusEta; i++)
	  // cout << m_objective[i] << endl;

   //---
   //--- Construct the objective function of subproblem.
   //---
   m_SP.m_objective = new double[m_nbSubproblemColsPlusDummy];
   //m_objective_SP = new double[m_nbSubproblemCols];
   if (!m_SP.m_objective)
	   throw UtilExceptionMemory("createModels", "BENMIP_DecompApp");

 /*  for (int j = 0; j < m_nbColsOrigModel; j++)
   {
	   cout << m_objCoefsOrigModel[j] << endl;
   }*/

   //memcpy(m_objective, 
   //       m_mpsIO.getObjCoefficients(), nCols * sizeof(double));
   if (m_appParam.ObjectiveSense == 1){
	   for (int i = 0; i < m_nbSubproblemCols; i++)
	   {
		   m_SP.m_objective[i] = m_objCoefsOrigModel[getSPVarOldIndex(i)];

		   UTIL_DEBUG(m_param.LogDebugLevel, 2,
			   std::cout << m_subproblemVarNames[i] << "\t" << m_objCoefsOrigModel[getSPVarOldIndex(i)] << endl;
		   )
	   }
	   m_SP.m_objective[m_nbSubproblemCols] = 10;
	   //UTIL_DEBUG(m_param.LogDebugLevel, 5,
	//	   std::cout << m_masterVarNames[getSPVarOldIndex(-)] << "\t" << m_objCoefsOrigModel[getSPVarOldIndex(i)] << endl;
	   //)
   }
   else
	   std::cout << "error: obj sense" << endl;
   // for eta



   //---
   //--- Construct the core matrix.
   //---
   //DecompConstraintSet * modelCore = &m_modelMaster;
   //createModelPart(modelCore, nRowsCore, rowsCore);
      
   setModelCore(&m_modelMaster, "MasterProblem");
   //setModelRelax(&m_modelRelax, "ModelRelaxOnlyEta");
   //m_models.push_back(model);
   //setModelRelax(&m_modelSubproblem, "SubProblem");

   
   //---
   //--- save a pointer so we can delete it later
   //---
   //m_modelC = modelCore;
   m_MPOSI = createOsiProblem(&m_modelMaster, m_objective, m_MPOSI, "MP");
   m_subprobOSI = createOsiProblem(m_SP.modelSubproblem, m_SP.m_objective, m_subprobOSI, "SP");

   //for (int i = 0; i < m_nbSubproblemCols; i++)
	  // cout << m_SP.m_objective[i] << endl;


	UTIL_DEBUG(m_param.LogDebugLevel, 2,
		std::cout << "Stats:" << endl;
	std::cout << "m_modelMaster.getNumInts(): " << m_modelMaster.getNumInts() << endl;
	std::cout << "m_modelMaster.getNumCols(): " << m_modelMaster.getNumCols() << endl;
	std::cout << "m_modelMaster.getNumRows(): " << m_modelMaster.getNumRows() << endl;
	std::cout << "----------------------------------------------------------" << endl;
	//std::cout << "m_modelSubproblem.getNumInts(): " << m_SP.modelSubproblem->getNumInts() << endl;
	std::cout << "m_modelSubproblem.getNumCols(): " << m_SP.modelSubproblem->getNumCols() << endl;
	std::cout << "m_modelSubproblem.getNumRows(): " << m_SP.modelSubproblem->getNumRows() << endl;
	std::cout << "----------------------------------------------------------" << endl;
	std::cout << "m_modelRelax.getNumInts(): " << m_modelRelax.getNumInts() << endl;
	std::cout << "m_modelRelax.getNumCols(): " << m_modelRelax.getNumCols() << endl;
	std::cout << "m_modelRelax.getNumRows(): " << m_modelRelax.getNumRows() << endl;
	std::cout << "----------------------------------------------------------" << endl;
	);






	m_param.DecompLPSolver = "CPLEX";
	m_param.DecompIPSolver = "CPLEX";
	m_param.DoInteriorPoint = true;
//#######################

	//writeRelaxedModel();


   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModels()", m_appParam.LogLevel, 2);   
   //exit(1);
}


bool BENMIP_DecompApp::solveSubproblem(const double * x, DecompCutList    & _newCuts)
{

	std::vector<double> xxx(x, x + m_nbMasterColsPlusEta);

	int nbRows = m_SP.modelSubproblem->getNumRows();

	std::vector<double> UBs = m_SP.modelSubproblem->rowUB;
	std::vector<double> LBs = m_SP.modelSubproblem->rowLB;

	// get all the master solutions
	m_SP.modelSubproblem->boundsToSenses();

	m_SP.m_subproblemOSI = m_subprobOSI;
	m_SP.m_subproblemOSI = createOsiProblem(m_SP.modelSubproblem, m_SP.m_objective, m_SP.m_subproblemOSI, "SP");


	double oldUB = 0;
	for (int i = 0; i < nbRows - m_SP.varLBs.size(); i++)
	{
		int nbMPVarsInRow = m_masterVariablesInMixedConstraints[i].getNbMPVars(i);
		if (nbMPVarsInRow > 0)
		{
			//cout << m_modelSubproblem.rowSense[i] << endl;
			switch (m_SP.modelSubproblem->rowSense[i])
			{
			case 'L':
				oldUB = m_SP.modelSubproblem->rowUB[i];
				for (int j = 0; j < nbMPVarsInRow; j++)
				{
					oldUB -= x[m_masterVariablesInMixedConstraints[i].m_masterVarIndexes[j]] *
						m_masterVariablesInMixedConstraints[i].m_masterVarCoeffs[j];
				}
				m_SP.modelSubproblem->rowUB[i] = oldUB;
				break;
			case 'G':
				oldUB = m_SP.modelSubproblem->rowLB[i];
				for (int j = 0; j < nbMPVarsInRow; j++)
				{
					oldUB -= x[m_masterVariablesInMixedConstraints[i].m_masterVarIndexes[j]] *
						m_masterVariablesInMixedConstraints[i].m_masterVarCoeffs[j];
				}
				m_SP.modelSubproblem->rowUB[i] = oldUB;
				break;
			case 'E':
				oldUB = m_SP.modelSubproblem->rowUB[i];
				for (int j = 0; j < nbMPVarsInRow; j++)
				{
					oldUB -= x[m_masterVariablesInMixedConstraints[i].m_masterVarIndexes[j]] *
						m_masterVariablesInMixedConstraints[i].m_masterVarCoeffs[j];
				}
				m_SP.modelSubproblem->rowUB[i] = oldUB;
				m_SP.modelSubproblem->rowLB[i] = oldUB;
				break;
			default:
				break;
			}

		}


	}

	m_SP.m_subproblemOSI = m_subprobOSI;
	m_SP.m_subproblemOSI = createOsiProblem(m_SP.modelSubproblem, m_SP.m_objective, m_SP.m_subproblemOSI, "SP");

	
	//set the dummy variable bounds to 0
	m_SP.modelSubproblem->colUB[m_SP.modelSubproblem->colUB.size() - 1] = 0;
	m_SP.m_subproblemOSI->setColUpper(&m_SP.modelSubproblem->colUB[0]);

	m_SP.modelSubproblem->colLB[m_SP.modelSubproblem->colLB.size() - 1] = 0;
	m_SP.m_subproblemOSI->setColLower(&m_SP.modelSubproblem->colLB[0]);


	// update SP RHS using the master index vector


	// solve 
	double primal_time = 0.0;
	double subp_time = 0.0;
	double start_time = 0.0;
	double end_time = 0.0;

	// Solve the (relaxation of the) problem
	start_time = CoinCpuTime();
	m_SP.m_subproblemOSI->setHintParam(OsiDoPresolveInInitial, false);
	m_SP.m_subproblemOSI->initialSolve();

	cout<< m_SP.m_subproblemOSI->getObjValue() << endl;

	

	end_time = CoinCpuTime();
	primal_time += end_time - start_time;

	cutType cuttype = OPTIMALITY;

	

	std::vector<double>  _dualSolOrRay;
	if (m_SP.m_subproblemOSI->isProvenOptimal())//->optimalBasisIsAvailable())//
	{
		cuttype = OPTIMALITY;
		m_SPStatus = OPTIMAL;
		const double * duals = m_SP.m_subproblemOSI->getRowPrice();
		_dualSolOrRay = std::vector<double>(duals, duals + m_SP.modelSubproblem->getNumRows());

		m_SPObjectiveValue = m_SP.m_subproblemOSI->getObjValue();
		cout << "m_SPObjectiveValue: " << m_SPObjectiveValue << endl;
		//getchar();
	}
	else
	{
		////set the dummy variable bounds to 0
		//m_SP.modelSubproblem->colUB[m_SP.modelSubproblem->colUB.size() - 1] = DecompInf;
		//m_SP.m_subproblemOSI->setColUpper(&m_SP.modelSubproblem->colUB[0]);

		//m_SP.modelSubproblem->colLB[m_SP.modelSubproblem->colLB.size() - 1] = -DecompInf;
		//m_SP.m_subproblemOSI->setColLower(&m_SP.modelSubproblem->colLB[0]);

		m_SP.m_subproblemOSI->setHintParam(OsiDoPresolveInInitial, false);
		m_SP.m_subproblemOSI->initialSolve();

		m_SPStatus = INFEASIBILE;
		cuttype = FEASIBILITY;

		//const double * duals = m_SP.m_subproblemOSI->getRowPrice();
		//_dualSolOrRay = std::vector<double>(duals, duals + m_SP.modelSubproblem->getNumRows());



		m_SPObjectiveValue = DecompInf;
		//m_modelRelax.M->appendRow(1, &inds[0], &elms[0]);
		//std::vector<double*> rays = m_subproblem.m_subproblemOSI->getDualRays(1);
		std::vector<double*> rays = m_SP.m_subproblemOSI->getDualRays(1);  //m_SP.getDualRays(1);
		_dualSolOrRay = std::vector<double>(rays[0], rays[0] + m_SP.m_subproblemOSI->getNumRows());
		if (strcmp(m_param.DecompLPSolver.c_str(), "CPLEX" ) == 0)
			std::transform(&_dualSolOrRay[0], &_dualSolOrRay[0] + _dualSolOrRay.size(), &_dualSolOrRay[0], std::negate<double>());
	

		//std::vector<double>::iterator min_elem = std::min_element(std::begin(_dualSolOrRay), std::end(_dualSolOrRay));
		//std::vector<double>::iterator max_elem = std::max_element(std::begin(_dualSolOrRay), std::end(_dualSolOrRay));
		//double val = fabs(*min_elem) > fabs(*max_elem) ? fabs(*min_elem) : fabs(*max_elem);
		//int divisor = log10(val)+1;

		//for (int i = 0; i < nbRows; i++)
		//{
		//	_dualSolOrRay[i] /= pow(10, divisor);// scale(rays[0][i], *min_elem, *max_elem, 0, 1);
		//}
	}




	// verify if the duals are correct



	//double objTest = 0;

	//for (int i = 0; i < m_masterVariablesInMixedConstraints.size(); i++)
	//{
	//	//_dualSolOrRay[i] = -_dualSolOrRay[i];
	//	cout << "u(" << i << ")= " << _dualSolOrRay[i] << endl;
	//	switch (m_modelSubproblem.rowSense[i])
	//	{
	//	case 'L':
	//		objTest += m_modelSubproblem.rowUB[i] * _dualSolOrRay[i];
	//		break;
	//	case 'G':
	//		objTest += m_modelSubproblem.rowLB[i] * _dualSolOrRay[i];
	//		break;
	//	case 'E':
	//		objTest += m_modelSubproblem.rowUB[i] * _dualSolOrRay[i];
	//		break;
	//	}
	//}
	//std::cout << "objTest: " << objTest << endl;


	// reset the bounds
	m_SP.modelSubproblem->rowUB = UBs;
	m_SP.modelSubproblem->rowLB = LBs;

	// generate cut
	std::vector<int > inds(m_nbMasterColsPlusEta, 0);
	std::vector<double > elms(m_nbMasterColsPlusEta, 0);

	double constantTerm = 0;
	for (int i = 0; i < m_masterVariablesInMixedConstraints.size(); i++)
	{
		//inds[i] = i;
		int nbMPVarsInRow = m_masterVariablesInMixedConstraints[i].getNbMPVars(i);
		switch (m_SP.modelSubproblem->rowSense[i])
		{

		case 'L':
			constantTerm += m_SP.modelSubproblem->rowUB[i] * _dualSolOrRay[i];
			if (nbMPVarsInRow > 0)
				for (int j = 0; j < nbMPVarsInRow; j++){
					elms[m_masterVariablesInMixedConstraints[i].m_masterVarIndexes[j]] -=
						m_masterVariablesInMixedConstraints[i].m_masterVarCoeffs[j] * _dualSolOrRay[i];
				}
			break;

		case 'G':
			constantTerm += m_SP.modelSubproblem->rowLB[i] * _dualSolOrRay[i];
			if (nbMPVarsInRow > 0)
				for (int j = 0; j < nbMPVarsInRow; j++){
					elms[m_masterVariablesInMixedConstraints[i].m_masterVarIndexes[j]] -=
						m_masterVariablesInMixedConstraints[i].m_masterVarCoeffs[j] * _dualSolOrRay[i];
				}
			break;
		case 'E':
			constantTerm += m_SP.modelSubproblem->rowUB[i] * _dualSolOrRay[i];
			if (nbMPVarsInRow > 0)
				for (int j = 0; j < nbMPVarsInRow; j++){
					elms[m_masterVariablesInMixedConstraints[i].m_masterVarIndexes[j]] -=
						m_masterVariablesInMixedConstraints[i].m_masterVarCoeffs[j] * _dualSolOrRay[i];
				}
			break;
		}

	}
	for (int i = m_masterVariablesInMixedConstraints.size(); i < nbRows; i++)
		constantTerm += m_SP.modelSubproblem->rowUB[i] * _dualSolOrRay[i];



	if (cuttype == OPTIMALITY)
		elms[m_nbMasterColsPlusEta - 1] = -1;

	
	//cout << constantTerm << endl;

	BENMIP_BendersCut*  bendersCut;
	if (cuttype == OPTIMALITY)
		bendersCut = new BENMIP_BendersCut(inds, elms, constantTerm, m_masterVarNames, OPTIMALITY);
	else
		bendersCut = new BENMIP_BendersCut(inds, elms, constantTerm, m_masterVarNames, FEASIBILITY);
	bendersCut->setBounds(-DecompInf, -constantTerm);


	UTIL_DEBUG(m_param.LogDebugLevel, 3,
		bendersCut->print();
	);



	_newCuts.push_back(bendersCut);
	//m_cutPool.push_back(bendersCut);

	return false;

	
}

bool BENMIP_DecompApp::isInteger(const double* x)
{
	bool is_integer = true;
	std::vector<double> xxx(x, x + m_nbMasterColsPlusEta);
	for (int i = 0; i < xxx.size(); i++)
	{
		if (fabs(xxx[i] - ceil(xxx[i])) > DecompEpsilon || fabs(xxx[i] - floor(xxx[i])) > DecompEpsilon)
			is_integer &= false;
		if (is_integer == false)
			return false;
	}

	return is_integer;
}

////===========================================================================//
//int BENMIP_DecompApp::generateInitVars(DecompVarList & initVars){	
//   UtilPrintFuncBegin(m_osLog, m_classTag,
//		      "generateInitVars()", m_appParam.LogLevel, 2);
//
//   readInitSolutionFile(initVars);
//   
//   UtilPrintFuncEnd(m_osLog, m_classTag,
//                    "generateInitVars()", m_appParam.LogLevel, 2);  
//   return static_cast<int>(initVars.size());
//}
//
//===========================================================================//

//===========================================================================//
int BENMIP_DecompApp::generateCuts(const double              * x,
	DecompCutList             & newCuts){
	static int enter = 0;
	static int last_node ;
	//if (isInteger(x) == false)
	//	return 0;
	DecompAlgo *algo = getDecompAlgo();

	//if (last_node > 0 && last_node == algo->getCurrentNode()->getIndex())
	//	return 0;


	UtilPrintFuncBegin(m_osLog, m_classTag,
		"generateCuts()", m_appParam.LogLevel, 1);

	UTIL_DEBUG(m_param.LogDebugLevel, 3,
		printSolution(x);
	);




	bool isFeas = solveSubproblem(x, newCuts);
	double LB = 0;
	double UB = 0;
	//if (algo->getObjBestBoundLB() > m_globalLB)
	//	LB = algo->getObjBestBoundLB();
	//if (algo->getObjBestBoundLB() + getSubProblemObj() - x[m_nbMasterColsPlusEta - 1] < m_globalUB)
	//UB = algo->getObjBestBoundLB() + getSubProblemObj() - x[m_nbMasterColsPlusEta - 1];

	//fbounds << LB << "\t" << UB << endl;
	if (isInteger(x) == true)
		UB = algo->getObjBestBoundLB() + getSubProblemObj() - x[m_nbMasterColsPlusEta - 1];

	enter++;
	if (LB>0 && enter == 1)
		m_modelMaster.rowLB[m_modelMaster.rowLB.size()-1] = 0;

	cout << "++++++++++++generateCuts+++++++++++++" << endl;
	cout << "UB: " << setw(6) << setprecision(6) << UB << endl;
	cout << "LB: " << setw(6) << setprecision(6) << m_globalLB << endl;
	cout << "gLB: " << setw(6) << setprecision(6) << algo->getObjBestBoundLB() << endl;
	cout << "+++++++++++++++++++++++++++++++++++++" << endl;

	//delete algo;
	//if ((m_globalUB - m_globalLB) < DecompEpsilon && isInteger(x) == true)
	//	return 0;
	m_nbBendersCuts++;
	// if the number of added cuts at the root node is more than ALPHA and 
	// and the current solution is fractional then do not add cut

	int nodeNumber = algo->getCurrentNode()->getIndex();
	last_node = nodeNumber;
	cout << "Node: " << nodeNumber << endl;
	if (nodeNumber == 0 && m_nbBendersCuts > 100)
		return 0;
	

	UtilPrintFuncEnd(m_osLog, m_classTag,
		"generateCuts()", m_appParam.LogLevel, 1);
	return 1;
}
//===========================================================================//
bool BENMIP_DecompApp::APPisUserFeasible(const double * x,
	const int      n_cols,
	const double   tolZero){

	m_nbBendersCuts++;



	UtilPrintFuncBegin(m_osLog, m_classTag,
		"APPisUserFeasible()", m_appParam.LogLevel, 1);

	UTIL_DEBUG(m_param.LogDebugLevel, 3,
		printSolution(x);
	);


	DecompCutList            newCuts;
	bool isFeas = solveSubproblem(x, newCuts);

	//TODO: feasibility cuts?
	//(*m_osLog) << "Not Feasible: disconnected, n_comp : " << n_comp << endl;

	DecompAlgo *algo = getDecompAlgo();

	if (algo->getObjBestBoundLB() > m_globalLB)
		m_globalLB = algo->getObjBestBoundLB();
	//if (algo->getObjBestBoundLB() + getSubProblemObj() - x[m_nbMasterColsPlusEta - 1] < m_globalUB)
	m_globalUB = algo->getObjBestBoundLB() + getSubProblemObj() - x[m_nbMasterColsPlusEta - 1];
	
	int nodeNumber = algo->getCurrentNode()->getIndex();
	cout << "Node: " << nodeNumber << endl;

	cout << "++++++++++APPisUserFeasible++++++++++" << endl;
	cout << "UB: " << setw(6) << setprecision(6) << m_globalUB << endl;
	cout << "LB: " << setw(6) << setprecision(6) << m_globalLB << endl;
	cout << "gLB: " << setw(6) << setprecision(6) << algo->getObjBestBoundLB() << endl;
	cout << "+++++++++++++++++++++++++++++++++++++" << endl;


	if (fabs(m_globalUB - algo->getObjBestBoundLB()) < DecompEpsilon && (m_globalUB < DecompInf && algo->getObjBestBoundLB() > -DecompInf))
		return true;





	UtilPrintFuncEnd(m_osLog, m_classTag,
		"APPisUserFeasible()", m_appParam.LogLevel, 1);

	return isFeas;
}

OsiSolverInterface* BENMIP_DecompApp::createOsiProblem( DecompConstraintSet* _model, const double * _objective, OsiSolverInterface* _OSI, char* _modelName, bool _exportLP/*=true*/ )
{
	//TODO: design question, we are assuming that master solver is
	//  an LP solver and relaxed solver is an IP - it really should
	//  be a generic object and an LP or IP solver is just one option
	//  for a solver
	_OSI = NULL;
	DecompConstraintSet* model = _model;// &m_modelMaster;

	if (!model || !model->M) {
		//---
		//--- if not using built-in solver, make sure user has
		//---   provided a solver function
		//--- TODO: how?
		//---
		//const DecompApp * app = getDecompApp();
		return NULL;
	}

	UtilPrintFuncBegin(m_osLog, m_classTag,
		"createOsiProblem()", m_param.LogDebugLevel, 2);
	int nInts = model->getNumInts();
	int nCols = model->getNumCols();
	int nRows = model->getNumRows();
	/*
	if (nInts) {
	subprobSI = new OsiIpSolverInterface();
	} else {
	subprobSI = new OsiLpSolverInterface();
	}
	*/
//#if defined(DIP_HAS_CLP)
//	_OSI = new OsiClpSolverInterface();
//#elif defined(DIP_HAS_CPX)
//	_OSI = new OsiCpxSolverInterface();
//#elif defined(DIP_HAS_SYMPHONY)
//	m_subprobOSI = new OsiSymSolverInterface();
//#elif defined(DIP_HAS_GRB)
//	m_subprobOSI = new OsiGrbSolverInterface();
//#endif

if ( strcmp( (m_param.DecompLPSolver).c_str() , "Clp") == 0)
	_OSI = new OsiClpSolverInterface();
else if (strcmp((m_param.DecompLPSolver).c_str(), "CPLEX") == 0)
	_OSI = new OsiCpxSolverInterface();
//else if (m_param.DecompLPSolver == "Sym")
//	m_subprobOSI = new OsiSymSolverInterface();
//else if (m_param.DecompLPSolver == "Grb")
//	m_subprobOSI = new OsiGrbSolverInterface();




	assert(_OSI);
	_OSI->messageHandler()->setLogLevel(m_param.LogLpLevel);
	//TODO: use assign vs load? just pass pointers?
	_OSI->loadProblem(*model->getMatrix(),
		model->getColLB(),
		model->getColUB(),
		//NULL, //null objective
		_objective,
		model->getRowLB(),
		model->getRowUB());

	if (nInts > 0) {
		_OSI->setInteger(model->getIntegerVars(), nInts);
//#if defined(__DECOMP_IP_CPX__) && defined (__DECOMP_LP_CPX__)
//		OsiCpxSolverInterface* osiCpx
//			= dynamic_cast<OsiCpxSolverInterface*>(_OSI);
//		osiCpx->switchToMIP();
//#endif
}
	//m_param.DecompIPSolver = "Clp";
	//m_param.DecompLPSolver = "Cbc";
	////---
	//--- set column and row names (if they exist)
	//---
	string           objName = "objective";
	vector<string>& colNames = model->colNames;
	vector<string>& rowNames = model->rowNames;
	_OSI->setIntParam(OsiNameDiscipline, 2);	//1=Lazy, 2=Full

	if (colNames.size()) {
		_OSI->setColNames(colNames, 0, nCols, 0);
	}

	if (rowNames.size()) {
		_OSI->setRowNames(rowNames, 0, nRows, 0);
	}

	_OSI->setObjName(objName);
	UTIL_DEBUG(m_param.LogDebugLevel, 5,
		int i;

	for (i = 0; i < nCols; i++) {
		(*m_osLog) << "User column name (" << i << ") = "
			<< colNames[i] << endl;
	}
	for (i = 0; i < nCols; i++) {
		(*m_osLog) << "OSI  column name (" << i << ") = "
			<< _OSI->getColName(i) << endl;
	}
	);

	if (_exportLP)
		_OSI->writeLp(_modelName);
	//---
	//--- set subproblem pointer
	//---
	//model_SP.setOsi(subprobSI);

	//double primal_time = 0.0;
	//double subp_time = 0.0;
	//double start_time = 0.0;
	//double end_time = 0.0;

	//// Solve the (relaxation of the) problem
	//start_time = CoinCpuTime();
	//subprobSI->initialSolve();
	//end_time = CoinCpuTime();
	//primal_time += end_time - start_time;

	UtilPrintFuncEnd(m_osLog, m_classTag,
		"createOsiSubProblem()", m_param.LogDebugLevel, 2);

	return _OSI;
}


//===========================================================================//
//DecompSolverStatus
//BENMIP_DecompApp::solveRelaxed(const int          whichBlock,
//const double     * redCostX,
//const double       convexDual,
//DecompVarList    & varList){
//
//	UtilPrintFuncBegin(m_osLog, m_classTag,
//		"solveRelaxed()", m_appParam.LogLevel, 2);
//
//	DecompSolverStatus   solverStatus = DecompSolStatNoSolution;
//	//TODO: BranchEnforceInSubProb option?
//
//	
//
//	UtilPrintFuncEnd(m_osLog, m_classTag,
//		"solveRelaxed()", m_appParam.LogLevel, 2);
//	return solverStatus;
//}

/*
#if 0
//===========================================================================//
DecompSolverStatus 
BENMIP_DecompApp::solveRelaxedNest(const int          whichBlock,
				      const double     * redCostX,
				      const double       convexDual,
				      DecompVarList    & varList){
   
   //---
   //--- solve full model heuristically  as IP
   //---   if get incumbent, break them out into approriate blocks
   //---   and return those partial columns
   //---


   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxedNest()", m_appParam.LogLevel, 2);
   
   /////STOP
   //---
   //--- this allows user direct access to access methods in 
   //---   algorithm interface (in case they want to use any 
   //---   of its data)
   //---
   //--- for example, if the user wants to enforce the branching
   //---   decisions in the oracle
   //--- TODO: can this be done using mcknap solver?
   //---
   //const DecompAlgo * decompAlgo = getDecompAlgo();
   //const double     * colLBNode  = decompAlgo->getColLBNode();
   //const double     * colUBNode  = decompAlgo->getColUBNode();

   //---
   //--- in the case where oracle=MCKP0, we have a specialized solver
   //--- in the case where oracle=MDKP,  we do not have a specialized solver
   //---   so, we just return with no solution and the framework will 
   //---   attempt to use the built-in MILP solver
   //---
   DecompSolverStatus solverStatus = DecompSolStatNoSolution;
   if(m_appParam.ModelNameRelax == "MCKP0"){
      vector<int>           solInd;
      vector<double>        solEls;
      double                varRedCost    = 0.0;
      double                varOrigCost   = 0.0;
      MMKP_MCKnap         * mcknapK       = m_mcknap[whichBlock];
      //TODO: check status return codes here
      mcknapK->solveMCKnap(redCostX, m_objective,
                           solInd, solEls, varRedCost, varOrigCost);
      assert(static_cast<int>(solInd.size()) == m_instance.getNGroupRows());
      
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 printf("PUSH var with k = %d RC = %g origCost = %g\n", 
                        whichBlock, varRedCost - convexDual, varOrigCost);
                 );
      
      //the user should not have to calculate orig cost too
      //  the framework can do this... in fact the framework shoudl
      //  calculate red-cost too... but the user might want to check this stuff
      
      //this is way too confusing for user to remember they need -alpha!
      //  let framework do that - also setting the block id - framework!
      DecompVar * var = new DecompVar(solInd, solEls, 
                                      varRedCost - convexDual, varOrigCost);
      var->setBlockId(whichBlock);
      varList.push_back(var);      
      solverStatus = DecompSolStatOptimal;
   }
           
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solveRelaxedNest()", m_appParam.LogLevel, 2);
   return solverStatus;
}
#endif
*/
