// $Id: TSP_SubtourCut.cpp,v 1.6 2004/03/03 01:00:19 magh Exp $

/*-------------------------------------------------------------------------
  Author: Matthew Galati (magh@lehigh.edu)

  (c) Copyright 2004 Lehigh University. All Rights Reserved.

  This software is licensed under the Common Public License. Please see 
  accompanying file for terms.
  ---------------------------------------------------------------------------*/

#include "BENMIP_BendersCut.h"

#include "Decomp.h"
#include "UtilMacros.h"

/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::init(){
   //setCutType();
   //switch(m_storage){
   //case VECTOR:
   //   create_bitset();
   //   break;
   //case BITSET:
   //   create_vector();
   //   break;
   //case BOTH:
   //   break;
   //default:
   //   //throw exception
   //   assert(0);
   //   return;
   //}
}

/*-------------------------------------------------------------------------*/
bool BENMIP_BendersCut::isSame(const DecompCut * cut) const{
   const BENMIP_BendersCut * sec_cut = dynamic_cast<const BENMIP_BendersCut*>(cut);
   if(!sec_cut)
      return false;

   if(m_type != sec_cut->m_type)
      return false;
   
   return m_idx == sec_cut->m_idx;
   return false;
}

/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::setCutType(){
   //what does concorde do?
   //sec_type = getSize() <= (n_vertices - 1)/2 ? SIDE : ACROSS;
   m_type = OPTIMALITY;
}

//another function that probably belongs in a utility class    
/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::create_bitset(){
   //create bitset from vector
   //for(int i = 0; i < m_nverts; i++)
   //   m_inS.push_back(false);//??
   //for(vector<int>::iterator it = m_S.begin(); it != m_S.end(); it++)
   //   m_inS[*it] = true;
   //m_storage = m_storage == VECTOR ? BOTH : BITSET;
}
  
//another function that probably belongs in a utility class    
/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::create_vector(){
  // //create vector from bistet
  // //m_S.reserve(m_inS.count());//is this worth it? or is count costly?
  // for(unsigned int u = 0; u < m_inS.size(); u++)
  //    if(m_inS[u]) 
	 //m_S.push_back(u);
  // m_storage = m_storage == BITSET ? BOTH : VECTOR;
}

/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::setBounds(double _lb , double _ub )
{
   switch(m_type){
   case OPTIMALITY:
	   setLowerBound(_lb);
	   setUpperBound(_ub);
      break;
   case FEASIBILITY:
	   setLowerBound(_lb);
	   setUpperBound(_ub);
      break;
   default:
      assert(0);
   }   
}

/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::expandCutToRow(CoinPackedVector * row){
   vector<int>    indices;
   vector<double> elements;
  
   for (int i = 0; i < m_idx.size(); i++)
   {
	   if (fabs(m_coeffs[i]) > DecompEpsilon)
	   {
		   indices.push_back(i);
		   elements.push_back(m_coeffs[i]);
	   }
   }

     
   row->setVector(static_cast<int>(indices.size()), 
		  &indices[0], &elements[0], false);
}

/*-------------------------------------------------------------------------*/
void BENMIP_BendersCut::print(ostream * os) const{
   double lb, ub;

   (*os) << "==================================================" << endl; 

   switch(m_type){
   case OPTIMALITY:      
	   (*os) << "Optimality Cut: ";
      break;
   case FEASIBILITY:
      (*os) << "Feasibility Cut: ";
      break;
   }

   /*  switch(m_storage){
	 case VECTOR:
	 {
	 (*os) << "S: ";
	 for(vector<int>::const_iterator it = m_S.begin();
	 it != m_S.end(); it++)
	 (*os) << *it << " ";
	 }
	 break;
	 case BITSET:
	 case BOTH:
	 {
	 (*os) << "S: ";
	 for(int i = 0; i < m_nverts; i++)
	 if(m_inS[i])
	 (*os) << i << " ";
	 }
	 break;
	 default:
	 cerr << "ERROR in print - BAD cut storage_type" << endl;
	 abort();
	 }*/

	   
  
   
   lb = getLowerBound();
   ub = getUpperBound();
   if(lb > -DecompInf){
	   (*os) << lb << " <=";
   }
   else{
      (*os) << " -INF <=";
   }

   for (int i = 0; i < m_names.size(); i++)
	   (*os) << m_coeffs[i] << " " << m_names[i] << " + ";

   if(ub < DecompInf){
      (*os) << "<= " << ub;
   }
   else{
      (*os) << "<= INF";
   }
   (*os) << endl;

   (*os) << "==================================================" << endl;

}
