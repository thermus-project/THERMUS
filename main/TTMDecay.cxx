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
// Stores parent & daughter id's as well as the branching ratio as a 
// fraction 
//

#include <TTMDecay.h>

ClassImp(TTMDecay)

//__________________________________________________________________________
TTMDecay::TTMDecay() 
{
  fParentID = 0;
  fDaughterID = 0;
  fBRatio = 0.;
}

//__________________________________________________________________________
TTMDecay::TTMDecay(Int_t parent, Int_t daughter, Double_t fraction)
{
  fParentID = parent;
  fDaughterID = daughter;
  fBRatio = fraction;
}

//__________________________________________________________________________
void TTMDecay::List()
{
  cout<<" BRatio: "<<fBRatio<<" Parent: "<<fParentID<<" Daughter: "
      <<fDaughterID<<endl;
}

//__________________________________________________________________________
TTMDecay& TTMDecay::operator=(const TTMDecay& obj)
{
  if (this == &obj) return *this;

  fParentID = obj.GetParentID();
  fDaughterID = obj.GetDaughterID();
  fBRatio = obj.GetBRatio();
  return *this;
}
