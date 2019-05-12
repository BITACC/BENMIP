#pragma once
#include "BaseAdvancedVector.h"
using namespace std;
#include <string>

class VarVecs :
	public BaseAdvancedVector
{
public:
	VarVecs();
	~VarVecs();
	int getOriginalIndex(){ return index_in_original_model; }
	void setOriginalIndex(int _idx){ index_in_original_model = _idx; }
	int getSPIndex(){ return index_in_sp; }
	void setIndex(int _idx){ index_in_sp = _idx; }
	void setName(string colName);
	void setInteger(bool param1);
	void setObjCoef(const double objCoefs);
private:
	int index_in_sp;
	int index_in_original_model;
	std::string m_colName;
	bool m_isInteger;
	double m_objCoefs;
};

