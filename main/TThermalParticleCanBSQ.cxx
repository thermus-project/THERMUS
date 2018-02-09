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
// Thermal Particle class based on the full canonical approach (i.e. 
// strangeness, baryon number and charge treated canonically).
// The canonical particle and energy densities and pressure in the Boltzmann
// approximation differ from those of the complete grand-canonical approach 
// by a multiplicative factor calculated from the densities of all particles 
// in the fireball. Thus at this level, this factor cannot be calculated 
// but rather has to be specified (fCorrFactor). Only when a particle set 
// is specified can it be determined.
//  

#include <TThermalParticleCanBSQ.h>

ClassImp(TTMThermalParticleCanBSQ)

//__________________________________________________________________________
TTMThermalParticleCanBSQ::TTMThermalParticleCanBSQ():TTMThermalParticle()
{
  fParameters = (TTMParameterSetCanBSQ *) 0;
  fCorrFactor = 1.;
}

//__________________________________________________________________________
TTMThermalParticleCanBSQ::TTMThermalParticleCanBSQ(TTMParticle *part,
                                               TTMParameterSetCanBSQ *parm, 
                                               Double_t correction):TTMThermalParticle()
{
  fParticle = part;
  fParameters = parm;
  fCorrFactor = correction;
}
  
//__________________________________________________________________________
TTMThermalParticleCanBSQ::TTMThermalParticleCanBSQ(TTMThermalParticleCanBSQ& obj)
{
  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
}

//__________________________________________________________________________
void TTMThermalParticleCanBSQ::UpdateMembers()
{	
  fDeg = fParticle->GetDeg();
  fM = fParticle->GetMass();
  fT = fParameters->GetT();

  fMu = 0.;
   
  Double_t SContent = fParticle->GetSContent();
  Double_t gammas = fParameters->GetGammas();

  if (SContent != 0.0) {
    fG = pow(gammas, SContent);
  } else {
    fG = 1.0;
  }
}

//__________________________________________________________________________
TTMThermalParticleCanBSQ& TTMThermalParticleCanBSQ::operator=(TTMThermalParticleCanBSQ& obj)
{
  if (this == &obj) return *this;

  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
  return *this;
}
