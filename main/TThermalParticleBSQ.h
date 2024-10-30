/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 26 April 2010 //

#include <TROOT.h>
#include <TObject.h>
#include <math.h>

#ifndef Thermal_TTMThermalParticle
#include <TThermalParticle.h>
#endif

#ifndef Thermal_TTMParameterSetBSQ
#include <TTMParameterSetBSQ.h>
#endif

#ifndef Thermal_TTMThermalParticleBSQ
#define Thermal_TTMThermalParticleBSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalParticleBSQ                                       		//
//                                                                  	//
// A thermal particle in the complete grand-canonical formalism        	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalParticleBSQ:public TTMThermalParticle {

 private:

  TTMParameterSetBSQ *fParameters;    // pointer to thermal parameters
   
  void UpdateMembers();
   
 public:

  TTMThermalParticleBSQ();
  TTMThermalParticleBSQ(TTMParticle *part, TTMParameterSetBSQ *parm);
  TTMThermalParticleBSQ(TTMThermalParticleBSQ& obj);

  TTMParameterSetBSQ* GetParameters() {return fParameters;}
  void SetParameters(TTMParameterSetBSQ *x) {fParameters = x;}
   
  Double_t DensityQStatNoWidth();
  Double_t DensityQStatWidth();
   
  Double_t EnergyQStatNoWidth();
  Double_t EnergyQStatWidth();
  
  Double_t EntropyBoltzmannNoWidth();
  Double_t EntropyBoltzmannWidth();
  Double_t EntropyQStatNoWidth();
  Double_t EntropyQStatWidth();
   
  Double_t PressureQStatNoWidth();
  Double_t PressureQStatWidth();

  Bool_t ParametersAllowed();	      // checks for Bose-Einstein Condensation

  TTMThermalParticleBSQ& operator=(TTMThermalParticleBSQ& obj);

  ClassDef(TTMThermalParticleBSQ,1) // The complete grand-canonical approach

};

#endif
