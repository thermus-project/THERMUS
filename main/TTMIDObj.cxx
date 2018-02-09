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
// ID object for storage in ROOT container. 
//  

#include <TTMIDObj.h>

ClassImp(TTMIDObj)

//__________________________________________________________________________
TTMIDObj::TTMIDObj()
{
  SetID(0);                    
}

//__________________________________________________________________________
TTMIDObj::TTMIDObj(Int_t x)
{
  SetID(x);
}

//__________________________________________________________________________
void TTMIDObj::SetID(Int_t x)
{
  // Sets ID  
  //
	
  fID = x;
}

//__________________________________________________________________________
TTMIDObj& TTMIDObj::operator=(const TTMIDObj& obj)
{
  if (this == &obj) return *this;

  SetID(obj.GetID());
  return *this;
}
