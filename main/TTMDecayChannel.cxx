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
// Stores daughter id's as well as the branching ratio as a 
// fraction 
//

#include <TTMDecayChannel.h>
#include <TTMIDObj.h>

ClassImp(TTMDecayChannel)

//__________________________________________________________________________
TTMDecayChannel::TTMDecayChannel() 
{
  fBRatio = 0.;
  fDaughters = (TList *)0;
}

//__________________________________________________________________________
TTMDecayChannel::TTMDecayChannel(Double_t fraction, TList *list)
{
  fBRatio = fraction;
  fDaughters = list;
}

//__________________________________________________________________________
void TTMDecayChannel::List()
{
  cout<<"\t\t BRatio: "<<fBRatio;

  if(fDaughters){
  TIter next(fDaughters);
  cout<<"\t\t"<<"Daughters: "<<"\t";
  TTMIDObj *did;
  while((did = (TTMIDObj *)next())){
     cout<<did->GetID()<<"\t";
  }
  }
  cout<<endl;
}

//__________________________________________________________________________
TTMDecayChannel& TTMDecayChannel::operator=(const TTMDecayChannel& obj)
{
  if (this == &obj) return *this;

  fBRatio = obj.GetBRatio();
  fDaughters = obj.GetDaughterList();
  return *this;
}
