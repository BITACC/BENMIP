// $Id: TSP_SubtourCut.hpp,v 1.9 2004/08/10 03:43:15 magh Exp $

/*-------------------------------------------------------------------------
  Author: Matthew Galati (magh@lehigh.edu)

  (c) Copyright 2004 Lehigh University. All Rights Reserved.

  This software is licensed under the Common Public License. Please see 
  accompanying file for terms.
  ---------------------------------------------------------------------------*/

#ifndef TSP_SUBTOUR_CUT_INCLUDED
#define TSP_BENDERS_CUT_INCLUDED

/*----------------------------------------------------------------------
  TSP_SubtourCut: TSP subtour elimination constraint
  
  (i)  ACROSS: sum{e in delta(S)} x_e >= 2
  (ii) SIDE  : sum{e in E(S)}     x_e <= |S| - 1
  ---------------------------------------------------------------------- */

#include "Decomp.h"
#include "DecompCut.h"
#include "UtilMacros.h"
#include "BENMIP_Utilities.h"
#include <vector>
using namespace std;

/*---------------------------------------------------------------------------*/
class BENMIP_BendersCut : public DecompCut {
public:
   enum storageType {VECTOR, BITSET, BOTH};

private:
   //vector<int>       m_S;
   //vector<bool>      m_inS;
   cutType           m_type;
   storageType       m_storage;
   std::vector<int>		m_idx;
   std::vector<double>		m_coeffs	;
   std::vector<string>		m_names;
   double					m_constTerm;

public:
   //these (pure virutal) methods are inherited from DecompCut
   virtual void expandCutToRow(CoinPackedVector * row);
   virtual void setBounds(double _lb = -DecompInf, double _ub = DecompInf);

   //these (virutal) methods are inherited from DecompCut
   virtual void print(ostream * os = &cout) const;
   virtual bool isSame(const DecompCut * cut) const;

public:
   void init();
   void setCutType();
   void create_bitset();
   void create_vector();

public:
	BENMIP_BendersCut( vector<int> & inds, const vector<double> & elems,  double constTerm, const vector<string> & names,
		  const cutType        type = OPTIMALITY){
      m_type    = type;
      setBounds(-DecompInf, -constTerm);
      init();
	  m_names = names;
	  m_coeffs = elems;
	  m_idx = inds;
   };
  
	//BENMIP_BendersCut(const vector<int> & inds, const vector<double> & elems, const vector<string> & names,
	//	  const cutType        type){
 //     m_type    = type;
 //     setBounds();
 //     init();
 //  };

	BENMIP_BendersCut(const vector<int> & inds, const vector<double> & elems, const vector<string> & names ){
      setCutType();
      setBounds();
      init();
   };

   virtual ~BENMIP_BendersCut() {}
};



   
  


#endif
