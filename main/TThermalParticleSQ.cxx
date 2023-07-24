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
//BQ file, TThermalParticleBQ.cxx,  by Spencer Wheaton 14 July 2004)
//__________________________________________________________________________
// Thermal Particle class based on the baryon-canonical approach (i.e. 
// baryon canonical, strangeness and charge grand-canonical).
// The baryon-canonical particle and energy densities and pressure in the 
// Boltzmann approximation differ from those of the complete grand-
// canonical approach by a multiplicative factor calculated from the 
// densities of all particles in the fireball. Thus at this level, this 
// factor cannot be calculated, but rather has to be specified (fCorrFactor). 
// Only when a particle set is specified can it be determined.
//  

#include <TThermalParticleSQ.h>

ClassImp(TTMThermalParticleSQ)

//__________________________________________________________________________
TTMThermalParticleSQ::TTMThermalParticleSQ():TTMThermalParticle()
{
  fParameters = (TTMParameterSetSQ *) 0;
  fCorrFactor = 1.;
}

//__________________________________________________________________________
TTMThermalParticleSQ::TTMThermalParticleSQ(TTMParticle *part,
                                       TTMParameterSetSQ *parm, Double_t correction)
  :TTMThermalParticle()
{
  fParticle = part;
  fParameters = parm;
  fCorrFactor = correction;
}

//__________________________________________________________________________
TTMThermalParticleSQ::TTMThermalParticleSQ(TTMThermalParticleSQ& obj)
{
  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
} 

//__________________________________________________________________________
void TTMThermalParticleSQ::UpdateMembers()
{	
  fDeg = fParticle->GetDeg();
  fM = fParticle->GetMass();
  fT = fParameters->GetT();

  Double_t S = fParticle->GetS();
  Double_t Q = fParticle->GetQ();

  Double_t muS = fParameters->GetMuS();
  Double_t muQ = fParameters->GetMuQ();

  fMu = S * muS + Q * muQ;

  Double_t BContent = fParticle->GetBContent();
  Double_t gammas = fParameters->GetGammas();

  if (BContent != 0.0) {
    fG = pow(gammas, BContent);
  } else {
    fG = 1.0;
  }
}

//__________________________________________________________________________
TTMThermalParticleSQ& TTMThermalParticleSQ::operator=(TTMThermalParticleSQ& obj)
{
  if (this == &obj) return *this;

  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
  return *this;
}
