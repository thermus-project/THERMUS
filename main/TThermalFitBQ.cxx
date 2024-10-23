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
// Strangeness-canonical thermal fit class
//  

#include <TThermalFitBQ.h>

ClassImp(TTMThermalFitBQ)
	
//__________________________________________________________________________
TTMThermalFitBQ::TTMThermalFitBQ():TTMThermalFit()
{
  fDescriptor = "SCanonical";
  fParm = (TTMParameterSetBQ *) 0;
  fWidth = true;
  fNonStrangeQStats = true;
  fModel = (TTMThermalModelBQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitBQ::TTMThermalFitBQ(TTMParticleSet *set, TTMParameterSetBQ *par, char *file):TTMThermalFit() 
{
  fDescriptor = "SCanonical";
  fParm = par;
  fPartSet = set;
  fWidth = true;
  fNonStrangeQStats = true;
  InputExpYields(file);            
  fModel = (TTMThermalModelBQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelBQ* TTMThermalFitBQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelBQ(set,fParm,fWidth);
  fModel->SetNonStrangeQStats(fNonStrangeQStats);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitBQ::~TTMThermalFitBQ()
{
  if(fModel){
    delete fModel;
  }
}
