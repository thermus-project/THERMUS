/************************************************************************
 * Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 *                                                                      *
 * Author: The THERMUS Project (A Thermal Model Package for ROOT).      *
 * Contributors (UCT-IPHC) are mentioned in the code where appropriate. *
 *                                                                      *
 * Permission to use, copy, modify and distribute this software and its *
 * documentation strictly for non-commercial purposes is hereby granted *
 * without fee, provided that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice *
 * appear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is        *
 * provided "as is" without express or implied warranty.                *
 ************************************************************************/

// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Base class of parameter set objects. All derived classes must define  
// Double_t GetRadius(), since this is used by this class to calculate 
// the volume.
//

#include <TTMParameterSet.h>

ClassImp(TTMParameterSet)

TTMParameterSet::TTMParameterSet()
{
  fPar = (TTMParameter *) 0;
  fConstraintInfo = "Parameters unconstrained";
}
