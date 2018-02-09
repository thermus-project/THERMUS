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
// Canonical thermal fit class
//  

#include <TThermalFitCanBSQ.h>

ClassImp(TTMThermalFitCanBSQ)

//__________________________________________________________________________
TTMThermalFitCanBSQ::TTMThermalFitCanBSQ():TTMThermalFit()
{
  fDescriptor = "BSQCanonical";
  fParm = (TTMParameterSetCanBSQ *) 0;
  fWidth = true;
  fModel = (TTMThermalModelCanBSQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitCanBSQ::TTMThermalFitCanBSQ(TTMParticleSet *set, TTMParameterSetCanBSQ *par, char *file):TTMThermalFit() 
{
  fDescriptor = "BSQCanonical";
  fParm = par;
  fPartSet = set;
  fWidth = true;
  InputExpYields(file);            
  fModel = (TTMThermalModelCanBSQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelCanBSQ* TTMThermalFitCanBSQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelCanBSQ(set,fParm,fWidth);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitCanBSQ::~TTMThermalFitCanBSQ()
{
  if(fModel){
    delete fModel;
  }
}

