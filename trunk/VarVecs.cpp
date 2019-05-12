#include "VarVecs.h"



VarVecs::VarVecs()
{
}


VarVecs::~VarVecs()
{
}

void VarVecs::setName(string colName)
{
	m_colName = colName;
}

void VarVecs::setInteger(bool param1)
{
	m_isInteger = param1;
}

void VarVecs::setObjCoef(const double objCoefs)
{
	m_objCoefs = objCoefs;
}
