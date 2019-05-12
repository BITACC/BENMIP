#pragma once
#include <vector>
#include <string>
using namespace std;
class BENMIP_MixedConstraints
{
public:
	BENMIP_MixedConstraints();
	~BENMIP_MixedConstraints();
	BENMIP_MixedConstraints(int idx) { m_masterRowIndex = idx; }
	std::vector<int> m_masterVarIndexes;
	std::vector<double> m_masterVarCoeffs;
	std::vector<string> m_masterVarName;
	double m_LB;
	double m_UB;
	int    m_masterRowIndex;

	int getNbMPVars(int i) { return m_masterVarIndexes.size(); };
	

};

