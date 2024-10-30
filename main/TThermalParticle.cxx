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
// Base class for thermal particles 
//

#include <TThermalParticle.h>
#include <TF2.h>
#include <TF1.h>
#include <FncsThermalModel.h>
#include <TMath.h>

ClassImp(TTMThermalParticle)

const  Double_t TTMThermalParticle::eps = 1e-5;   // required accuracy of integrals

//__________________________________________________________________________
TTMThermalParticle::TTMThermalParticle()
{
  fParticle = (TTMParticle*) 0;
}

//__________________________________________________________________________
Double_t TTMThermalParticle::DensityBoltzmannNoWidth()
{
  // Primordial Particle density assuming no width and Boltzmann stats
  // 
	
  UpdateMembers();
   
  return fCorrFactor / (2. * pow(TMath::Pi(), 2)) * fG * fDeg * pow(fM, 2) * 
    fT * TMath::BesselK(2,fM / fT) * exp(fMu / fT) / pow(0.197, 3.);
}

//__________________________________________________________________________
Double_t TTMThermalParticle::DensityBoltzmannWidth()
{
  // Primordial Particle density assuming finite width and Boltzmann stats
  // 

  UpdateMembers();
  Double_t lDensity = 0;
  Double_t width = fParticle->GetWidth();
  
  if(width != 0){

    Double_t threshold = fParticle->GetThreshold();

    Double_t a = TMath::Max(fM - 2.*width,threshold);

    TF1 *fn = new TF1("n Boltzmann Width",FcnDensBoltzmannWidth,0.,(fM + 3.*width)/fT,5);
    fn->SetParameters(fMu/fT,fM/fT,fG,fDeg,width/fT);
  
    Double_t n = IntegrateLegendre40(fn,a/fT,(fM + 2.*width)/fT);

    TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
    fnorm->SetParameters(fM/fT,width/fT);

    Double_t norm = IntegrateLegendre40(fnorm,a/fT,(fM + 2.*width)/fT);

    lDensity = fCorrFactor * pow(fT,3.) * n / norm / pow(0.197, 3.);
    delete fn;
    delete fnorm;

  }else{

    lDensity = DensityBoltzmannNoWidth();

  }
  return lDensity;
}

//__________________________________________________________________________
Double_t TTMThermalParticle::EnergyBoltzmannNoWidth()
{
  // Primordial Energy density contribution assuming no width and 
  // Boltzmann stats
  //
	
  UpdateMembers();
   
  return fCorrFactor / (2. * pow(TMath::Pi(), 2)) * fG * fDeg * pow(fM, 3)
    * fT * exp(fMu / fT) * (TMath::BesselK(1, fM / fT) + 3. * fT / fM * 
                            TMath::BesselK(2, fM / fT)) / pow(0.197, 3.);
}

//__________________________________________________________________________
Double_t TTMThermalParticle::EnergyBoltzmannWidth()
{
  // Primordial Energy density contribution assuming finite width and 
  // Boltzmann stats
  //

  UpdateMembers();	
  Double_t lEnergy = 0;
  Double_t width = fParticle->GetWidth();
  
  if(width != 0){

    Double_t threshold = fParticle->GetThreshold();

    Double_t a = TMath::Max(fM - 2.*width,threshold);

    TF1 *fe = new TF1("e Boltzmann Width",FcnEnergyBoltzmannWidth,0.,(fM + 3.*width)/fT,5);
    fe->SetParameters(fMu/fT,fM/fT,fG,fDeg,width/fT);
  
    Double_t e = IntegrateLegendre40(fe,a/fT,(fM + 2.*width)/fT);

    TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
    fnorm->SetParameters(fM/fT,width/fT);

    Double_t norm = IntegrateLegendre40(fnorm,a/fT,(fM + 2.*width)/fT);

    lEnergy = fCorrFactor * pow(fT,4.) * e / norm / pow(0.197, 3.);
    delete fe;
    delete fnorm;

  }else{

    lEnergy = EnergyBoltzmannNoWidth();

  }
  return lEnergy;
}

//__________________________________________________________________________
Double_t TTMThermalParticle::PressureBoltzmannNoWidth()
{
  // Pressure contribution assuming no width and Boltzmann stats
  //
	
  UpdateMembers();
   
  Double_t dens = DensityBoltzmannNoWidth();

  return dens * fT ;
}

//__________________________________________________________________________
Double_t TTMThermalParticle::PressureBoltzmannWidth()
{
  // Pressure contribution assuming finite width and Boltzmann stats
  //
	
  UpdateMembers();

  Double_t dens = DensityBoltzmannWidth();

  return dens * fT;
}
