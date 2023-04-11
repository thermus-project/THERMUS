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

//Authors: Natasha Sharma and Jean Cleymans 2 February 2021 (modification of the corresponding
//BQ file, TThermalFitBQ.cxx,  by Spencer Wheaton 14 July 2004)
//__________________________________________________________________________
// Baryon-canonical thermal fit class

#include <TThermalFitSQ.h>

ClassImp(TTMThermalFitSQ)
	
//__________________________________________________________________________
TTMThermalFitSQ::TTMThermalFitSQ():TTMThermalFit()
{
  fDescriptor = "SCanonical";
  fParm = (TTMParameterSetSQ *) 0;
  fWidth = true;
  fNonBaryonQStats = true;
  fModel = (TTMThermalModelSQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitSQ::TTMThermalFitSQ(TTMParticleSet *set, TTMParameterSetSQ *par, char *file):TTMThermalFit() 
{
  fDescriptor = "SCanonical";
  fParm = par;
  fPartSet = set;
  fWidth = true;
  fNonBaryonQStats = true;
  InputExpYields(file);            
  fModel = (TTMThermalModelSQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelSQ* TTMThermalFitSQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelSQ(set,fParm,fWidth);
  fModel->SetNonBaryonQStats(fNonBaryonQStats);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitSQ::~TTMThermalFitSQ()
{
  if(fModel){
    delete fModel;
  }
}
