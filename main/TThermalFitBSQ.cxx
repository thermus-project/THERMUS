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
// BSQ thermal fit class
//  

#include <TThermalFitBSQ.h>

ClassImp(TTMThermalFitBSQ)

//__________________________________________________________________________
TTMThermalFitBSQ::TTMThermalFitBSQ():TTMThermalFit()
{
  fDescriptor = "GCanonical";
  fParm = (TTMParameterSetBSQ *) 0;
  fQStats = true;
  fWidth = true;
  fExclVol = false;
  fModel = (TTMThermalModelBSQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitBSQ::TTMThermalFitBSQ(TTMParticleSet *set, TTMParameterSetBSQ *par, const char *file):TTMThermalFit() 
{
  fDescriptor = "GCanonical";
  fParm = par;
  fPartSet = set;
  fQStats = true;
  fWidth = true;
  fExclVol = false;
  InputExpYields(file);            
  fModel = (TTMThermalModelBSQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelBSQ* TTMThermalFitBSQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelBSQ(set,fParm,fQStats,fWidth);
  fModel->SetExcludedVolume(fExclVol);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitBSQ::~TTMThermalFitBSQ()
{
  if(fModel){
    delete fModel;
  }
}

