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

#ifndef BENMIP_DECOMPAPP_INCLUDED
#define BENMIP_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"
#include "BENMIP_Param.h"
#include "BENMIP_Utilities.h"
#include "BENMIP_MixedConstraints.h"
#include "BENMIP_BendersCut.h"
#include "BENMIP_Subproblem.h"
//===========================================================================//
#include "CoinMpsIO.hpp"
#include "VarVecs.h"
//===========================================================================//

/*!
 * \class BENMIP_DecompApp
 * A DecompApp to illustrate a basic usage of Decomp.
 * 
 * \see
 * DecompApp
 */
//===========================================================================//
class BENMIP_DecompApp : public DecompApp{
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;
   
   /** MPS object for reading BENMIP instances. */
   CoinMpsIO m_lpIO;
   
   /** Application specific parameters. */
   BENMIP_Param m_appParam;

   /** The model objective coefficients (original space). */
   //double * m_objective;
   double * m_objective;
   //double * m_objective_SP;
   
   /** The model constraint systems used for different algos. */
   DecompConstraintSet *          m_modelC;
   map<int, DecompConstraintSet*> m_modelR;

   /** Definition of blocks (by rows). */
   map<int, vector<int> > m_blocks;

   /**master variable file name**/
   vector<string>	m_masterVarNames;
   vector<int>		m_masterVarIndexes;
   vector<string>	m_subproblemVarNames;
   vector<int>		m_subproblemVarIndexes;
   std::vector<double>	m_masterObjCoeffs;
   std::vector<double>	m_subprobemObjCoeffs;

   /**master and subproblem models**/
   DecompConstraintSet m_modelMaster;
   //DecompConstraintSet m_modelSubproblem;
   DecompConstraintSet m_modelRelax;

   std::vector<std::pair<int, int> > m_masterOldIndexToNewIndex;
   std::vector<std::pair<int, int> > m_subproblemOldIndexToNewIndex;
   int				m_nbMPVars;
   int				m_nbSPVars;
   std::vector<BENMIP_MixedConstraints > m_masterVariablesInMixedConstraints;
   std::vector<double>		m_masterRowLBs, m_masterRowUBs, m_subproblemRowLBs, m_subproblemRowUBs;

   std::vector<string>		m_rowNamesOrigModel, m_colNamesOrigModel;
   std::vector<int> m_colsMarkerOrigModel;
   std::vector<int> m_rowsMarkerOrigModel;

   BENMIP_Subproblem	m_SP;

private:

   /* @name Inherited (from virtual) methods. */

   /** Generate init columns. */
   //virtual int generateInitVars(DecompVarList & initVars);

   int generateCuts(const double * x, DecompCutList & newCuts);

   void printSolution(const double * x)
   {
	   cout << "==================================================" << endl;
	   for (int i = 0; i < m_nbMasterColsPlusEta; i++)
	   {
		   cout << m_masterVarNames[i] << " : " << x[i] << endl;
	   }
	   cout << "==================================================" << endl;
   }

   bool APPisUserFeasible(const double * x, const int n_cols, const double tolZero);
   OsiSolverInterface* createOsiProblem(DecompConstraintSet* _model, const double * _objective, OsiSolverInterface* _OSI, char* _modelName, bool _exportLP=true);

   /** Solve the relaxed problem. */
   //virtual DecompSolverStatus solveRelaxed(const int          whichBlock,
	  // const double     * redCostX,
	  // const double       convexDual,
	  // DecompVarList    & varList);

   /** @name Helper functions (private). */

   /** Initialize application. */
   void initializeApp(UtilParameters & utilParam);

   std::vector<std::vector<int> > m_indexOfIntegersVarsInMixedConstraints;

   void readMasterVars()
   {
	   string       MasterVarsFile = m_appParam.MasterVarsFile;
	   ifstream varNamefile(MasterVarsFile.c_str());
	   if (varNamefile.is_open())
	   {
		   (*m_osLog) << "Reading master variable file" << endl;
		   string line;
		   while (!varNamefile.eof())
		   {
			   getline(varNamefile, line);
			   m_masterVarNames.push_back(trim(line, " "));
		   }
	   }
	   else
	   {
		   cerr << "The Master Variable file is not accessible" << endl;
	   }
   }
   OsiSolverInterface*   m_subprobOSI;
   OsiSolverInterface*   m_relaxOSI;
   OsiSolverInterface*   m_MPOSI;
   /** Create model parts. */
   void                  createModels();

   int					m_nbBendersCuts;

   void initializeModels()
   {
	   //---
	   //--- Construct the master and sub problem matrices.
	   //---
	   m_modelMaster.M = new CoinPackedMatrix(false, 0.0, 0.0);
	   if (!m_modelMaster.M)
		   throw UtilExceptionMemory("createModels", "MILP_DecompApp");
	   //m_modelMaster.reserve(m_nbMasterRows, m_nbMasterColsPlusEta);
	   m_modelMaster.reserve(0, m_nbMasterColsPlusEta);

	   //---
	   //--- Construct the subproblem matrix.
	   //---
	   m_SP.modelSubproblem = new DecompConstraintSet();
	   m_SP.modelSubproblem->M = new CoinPackedMatrix(false, 0.0, 0.0);
	   if (!m_SP.modelSubproblem->M)
		   throw UtilExceptionMemory("createModels", "MILP_DecompApp");
	   //m_modelSubproblem.reserve(m_nbSubproblemRows, m_nbSubproblemCols);
	   m_SP.modelSubproblem->reserve(0, m_nbSubproblemColsPlusDummy);


	   m_modelRelax.M = new CoinPackedMatrix(false, 0.0, 0.0);
	   if (!m_modelRelax.M)
		   throw UtilExceptionMemory("createModels", "MILP_DecompApp");
	   m_modelRelax.reserve(1, m_nbMasterColsPlusEta);
   }

   void readRows(const int nbRows, std::vector<string> &rowNames, const CoinPackedMatrix * M, const int * beg, const int * len, const int * ind, std::vector<int> &rowsMarker, std::vector<string> colNames, const double * rowLB, const double * rowUB)
   {
	   cout << "partitioning rows of the original model. " << endl;

	   for (int i = 0; i < nbRows; i++)
	   {
		   cout << "\r" << i << " of " << nbRows ;
		   rowNames.push_back(m_lpIO.rowName(i));

		   CoinShallowPackedVector row = M->getVector(i);
		   const double *vecelem = row.getElements();
		   const int *vecind = row.getIndices();
		   const int vecsize = row.getNumElements();


		   //---
		   //--- which columns are present in this  rows
		   //---
		   int startIndex = beg[i];
		   int endIndex = beg[i] + len[i];

		   bool onlyInteger = true;
		   bool hasInteger = false;
		   for (int k = 0; k < len[i]; k++)
		   {
			   if ( isIn(m_masterVarIndexes, ind[startIndex + k]) ){
				   onlyInteger &= false;
			   }
			   else
			   {
				   hasInteger |= true;
				   m_indexOfIntegersVarsInMixedConstraints[i].push_back(ind[startIndex + k]);
			   }
		   }

		   if (onlyInteger == true)
		   {
			   rowsMarker[i] = UNIQUELY_INTEGER;  //-1 corresponds to the Master problem
			   m_nbMasterRows++;

			   //m_modelMaster.M->appendRow(vecsize, vecind, vecelem);
			   int *vecNewInd = (int*)row.getIndices();
			   double *vecNewElms = (double*)row.getElements();
			   //std::vector<int> inds;
			   //inds.assign(vecNewInd, vecNewInd + row.getNumElements());
			   //std::vector<double> elms;
			   //elms.assign(vecNewElms, vecNewElms + row.getNumElements());

			   //for (int kk = 0; kk < vecsize; kk++)
				  // cout << vecind[kk] << "\t" << vecelem[kk] << endl;


			   for (int k = 0; k < vecsize; k++)
			   {
				   int idx = getMPVarNewIndex(vecind[k]);
				   assert(idx > -1);
				   vecNewInd[k] = idx;

				   UTIL_DEBUG(m_param.LogDebugLevel, 5,
					   std::cout << colNames[vecind[k]] << "\t";
				   );

			   }

			   UTIL_DEBUG(m_param.LogDebugLevel, 5,
				   std::cout << endl;
			   );
			   m_modelMaster.M->appendRow(vecsize, vecNewInd, vecelem);
			   m_masterRowLBs.push_back(rowLB[i]);
			   m_masterRowUBs.push_back(rowUB[i]);
			   m_modelMaster.rowLB.push_back(rowLB[i]);
			   m_modelMaster.rowUB.push_back(rowUB[i]);
			   m_modelMaster.rowNames.push_back(rowNames[i]);
		   }
		   else
		   {
			   if (hasInteger && onlyInteger == false)
			   {
				   // for all the mixed constraints, we create new constrains with 
				   // mapped indexes and do not insert any of the MP variables but we
				   // store all the indexes of MP variables in the mixed constraints 
				   // w.r.t the new indexes. 
				   // ATTN: new constraint do not contain MP variables.

				   rowsMarker[i] = MIXED; // index of subproblem
				   m_nbSubproblemRows++;
				   m_masterVariablesInMixedConstraints.push_back(BENMIP_MixedConstraints(m_nbSubproblemRows - 1));

				   std::vector<int> indexes;
				   std::vector<double> coeffs;
				   for (int k = 0; k < vecsize; k++)
				   {
					   if (m_colsMarkerOrigModel[vecind[k]] == CONTINUOUS)
					   {
						   int idx = getSPVarNewIndex(vecind[k]);
						   assert(idx > -1);
						   UTIL_DEBUG(m_param.LogDebugLevel, 5,
							   std::cout << m_subproblemVarNames[idx] << "\t" << vecelem[k] << endl;
						   );
						   indexes.push_back(idx);
						   coeffs.push_back(vecelem[k]);
					   }
					   else
					   {
						   int idx = getMPVarNewIndex(vecind[k]);
						   assert(idx > -1);
						   m_masterVariablesInMixedConstraints[m_nbSubproblemRows - 1].m_masterVarIndexes.push_back(idx);
						   m_masterVariablesInMixedConstraints[m_nbSubproblemRows - 1].m_masterVarCoeffs.push_back(vecelem[k]);
						   m_masterVariablesInMixedConstraints[m_nbSubproblemRows - 1].m_masterVarName.push_back(colNames[vecind[k]]);

						   UTIL_DEBUG(m_param.LogDebugLevel, 5,
							   std::cout << m_masterVarNames[idx] << "\t" << vecelem[k] << endl;
						   )
					   }

				   }
				   UTIL_DEBUG(m_param.LogDebugLevel, 5,
					   std::cout << "\n";
				   )

					// add dummy variable
				   int idx = getSPVarNewIndex(-1);
				   indexes.push_back(idx);
				   coeffs.push_back(-1);

					   // the new row to be added is not of the same size as before,
					   // 
				   m_SP.modelSubproblem->M->appendRow(indexes.size(), &(indexes[0]), &(coeffs[0]));
				   m_masterVariablesInMixedConstraints[m_nbSubproblemRows - 1].m_LB = rowLB[i];
				   m_masterVariablesInMixedConstraints[m_nbSubproblemRows - 1].m_UB = rowUB[i];

				   m_subproblemRowLBs.push_back(rowLB[i]);
				   m_subproblemRowUBs.push_back(rowUB[i]);
				   m_SP.modelSubproblem->rowLB.push_back(rowLB[i]);
				   m_SP.modelSubproblem->rowUB.push_back(rowUB[i]);
				   m_SP.modelSubproblem->rowNames.push_back(rowNames[i]);

			   }
			   else
			   {

				   rowsMarker[i] = UNIQUELY_CONTINUOUS;
				   m_nbSubproblemRows++;
				   m_masterVariablesInMixedConstraints.push_back(BENMIP_MixedConstraints(m_nbSubproblemRows - 1));

				   std::vector<int> indexes;
				   std::vector<double> coeffs;
				   int *vecNewInd = (int*)row.getIndices();

				   for (int i = 0; i < vecsize; i++)
				   {
					   int idx = getSPVarNewIndex(vecind[i]);
					   vecNewInd[i] = idx;
					   indexes.push_back(idx);
					   coeffs.push_back(vecelem[i]);
				   }

				   //m_SP.modelSubproblem->M->appendRow(vecsize, vecNewInd, vecelem);
				   // add dummy variable
				   int idx = getSPVarNewIndex(-1);
				   indexes.push_back(idx);
				   coeffs.push_back(-1);
				   m_SP.modelSubproblem->M->appendRow(indexes.size(), &(indexes[0]), &(coeffs[0]));

				   m_subproblemRowLBs.push_back(rowLB[i]);
				   m_subproblemRowUBs.push_back(rowUB[i]);

				   m_SP.modelSubproblem->rowLB.push_back(rowLB[i]);
				   m_SP.modelSubproblem->rowUB.push_back(rowUB[i]);
				   m_SP.modelSubproblem->rowNames.push_back(rowNames[i]);
			   }

		   }

	   }
   }
   std::vector<VarVecs> m_MPVars, m_allVars;
   void readColumns(const int nbCols, std::vector<int> &colsMarker, const double * objCoefs, std::vector<string> &colNames)
   {
	   m_nbMPVars = -1;
	   m_nbSPVars = -1;
	   cout << "reading columns in the original model " << endl;
	   for (int j = 0; j < nbCols; j++)
	   {
		   cout << "\r" << j << " of " << nbCols;
		   VarVecs var;
		   string colName = trim(m_lpIO.columnName(j));
		   int idx = m_lpIO.columnIndex(colName.c_str());
		   
		   //-----------------------------
		   if (m_lpIO.isInteger(j)) var.setInteger(true);


		   if (m_lpIO.isInteger(j))
			   m_masterVarNames.push_back(colName);
		   else
			   m_subproblemVarNames.push_back(colName);



		   if (std::find(m_masterVarNames.begin(), m_masterVarNames.end(), colName) != m_masterVarNames.end())
		   {

			   m_masterVarIndexes.push_back(idx);
			   colsMarker[j] = INTEGER;
			   m_masterOldIndexToNewIndex.push_back(std::make_pair(idx, ++m_nbMPVars));
			   m_masterObjCoeffs.push_back(objCoefs[idx]);

		   }
		   else
		   {
			   m_subproblemVarIndexes.push_back(idx);
			   colsMarker[j] = CONTINUOUS;
			   m_subproblemOldIndexToNewIndex.push_back(std::make_pair(idx, ++m_nbSPVars));
			   m_subprobemObjCoeffs.push_back(objCoefs[idx]);
		   }
		   colNames.push_back(colName);
	   }
	   //	eta variable

	   m_nbSubproblemCols = nbCols - m_masterVarIndexes.size();
	   m_nbMasterCols = m_masterVarIndexes.size();
	   m_nbMasterColsPlusEta = m_masterVarIndexes.size() + 1;

	   colNames.push_back("eta");
	   m_masterOldIndexToNewIndex.push_back(std::make_pair(-1, ++m_nbMPVars));
	   m_masterObjCoeffs.push_back(1);

	   // subproblem dummy
	   m_nbSubproblemColsPlusDummy = m_nbSubproblemCols + 1;
	   m_subproblemVarIndexes.push_back(-1);
	   m_subproblemVarNames.push_back("delta");
	   m_subproblemOldIndexToNewIndex.push_back(std::make_pair(-1, ++m_nbSPVars));
   }

   DecompCutList             m_cutPool;

  

   int getMPVarNewIndex(const int i)
   {
	   int idx = -100;
	   std::vector<std::pair<int, int> >::iterator mapIter = m_masterOldIndexToNewIndex.begin();
	   while (mapIter != m_masterOldIndexToNewIndex.end()){
		   if (mapIter->first == i)
		   {
			   idx = mapIter->second;
			   break;
		   }
		   mapIter++;
	   }
	   return idx;
   }
   int getSPVarNewIndex(const int  i)
   {
	   int idx = -100;
	   std::vector<std::pair<int, int> >::iterator mapIter = m_subproblemOldIndexToNewIndex.begin();
	   while (mapIter != m_subproblemOldIndexToNewIndex.end()){
		   if (mapIter->first == i)
		   {
			   idx = mapIter->second;
			   break;
		   }
		   mapIter++;

	   }
	   return idx;
   }
   int getMPVarOldIndex(const int i)
   {
	   int idx = -1;
	   std::vector<std::pair<int, int> >::iterator mapIter = m_masterOldIndexToNewIndex.begin();
	   while (mapIter != m_masterOldIndexToNewIndex.end()){
		   if (mapIter->second == i)
		   {
			   idx = mapIter->first;
			   break;
		   }
		   mapIter++;
	   }
	   return idx;
   }
   int getSPVarOldIndex(const int  i)
   {
	   int idx = -1;
	   std::vector<std::pair<int, int> >::iterator mapIter = m_subproblemOldIndexToNewIndex.begin();
	   while (mapIter != m_subproblemOldIndexToNewIndex.end()){
		   if (mapIter->second == i)
		   {
			   idx = mapIter->first;
			   break;
		   }
		   mapIter++;
	   }
	   return idx;
   }
   DecompConstraintSet * createModelPart(const int   nRowsPart,
                                         const int * rowsPart);

   Status   m_SPStatus;
   Status	getSubProblemStatus(){ return m_SPStatus; }

   double	m_SPObjectiveValue;
   double	getSubProblemObj(){ return m_SPObjectiveValue; }

   //void createModelPart(DecompConstraintSet * model,
			//const int             nRowsPart,
			//const int           * rowsPart);
   //void createModelPartSparse(DecompConstraintSet * model,
			//      const int             nRowsPart,
			//      const int           * rowsPart);   
   //void                  createModelMasterOnlys(vector<int> & masterOnlyCols);
   //void                  readInitSolutionFile(DecompVarList & initVars);

   double					m_globalLB;
   double					m_globalUB;

   int						m_nbMasterRows ;
   int						m_nbMasterCols  ;
   int						m_nbMasterColsPlusEta  ;
   int						m_nbSubproblemRows  ;
   int						m_nbSubproblemCols  ;
   int						m_nbSubproblemColsPlusDummy;

   /** Read block file. */
   void readBlockFile();
   ofstream fbounds;
   ///** Find the active columns for some block. */
   //void findActiveColumns(const vector<int> & rowsPart,
   //                       set<int>          & activeColsSet);
   
public:
   /** User access methods. */
   
   /** Get instance name. */
   const string getInstanceName(){
      return m_appParam.Instance;
   }
  
public:
   /** @name Constructor and Destructor */
   BENMIP_DecompApp(UtilParameters & utilParam) : 
      DecompApp  (utilParam),
      m_classTag ("BENMIP-APP"),
      m_objective(NULL),
	  m_globalLB(-DecompInf),
	  m_globalUB(DecompInf)
   {
      initializeApp(utilParam); //can there be a default?
	  m_nbBendersCuts = 0;
	  fbounds.open("bounds.txt");
   }
   
   virtual ~BENMIP_DecompApp() {
      //UTIL_DELARR(m_objective);
      UTIL_DELPTR(m_modelC);
      UtilDeleteMapPtr(m_modelR);
   }
   bool  solveSubproblem(const double * x, DecompCutList    & _newCuts);
   bool isInteger(const double* x);
};

#endif
