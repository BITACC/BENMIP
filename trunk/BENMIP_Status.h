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

#ifndef GAP_STATUS_INCLUDED
#define GAP_STATUS_INCLUDED

//===========================================================================//
#include <string>
//===========================================================================//
//---
//--- return codes
//---
//===========================================================================//
enum BendersStatus {
   BendersStatusOk = 0,
   BendersStatusError,      //general error
   BendersStatusFileIO,     //error in i/o
   BendersStatusOutOfMemory //out of memory
};
const std::string BendersStatusStr[4] = {
   "BendersStatusOk",
   "BendersStatusError",
   "BendersStatusFileIO",
   "BendersStatusOutOfMemory"
};
#endif
