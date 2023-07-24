/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

//Authors: Natasha Sharma and Jean Cleymans 2 February 2021 (modification of the corresponding
//BQ file, TThermalParticleBQ.h,  by Spencer Wheaton 14 July 2004)
#include <TROOT.h>
#include <TObject.h>
#include <math.h>

#ifndef Thermal_TTMThermalParticle
#include <TThermalParticle.h>
#endif

#ifndef Thermal_TTMParameterSetSQ
#include <TTMParameterSetSQ.h>
#endif

#ifndef Thermal_TTMThermalParticleSQ
#define Thermal_TTMThermalParticleSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalParticleSQ                                        		//
//                                                                  	//
// A thermal particle in the formalism treating baryon canonically	//
// and charge and strangeness number grand-canonically. The Boltzmann 	//
// approximation is used since in this case the correction to the 	//
// grand-canonical results amounts to a simple multiplicative factor.	//
// The entropy cannot be split into the sum of particle entropies, so	//
// no entropy calculation is made here.				       	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalParticleSQ:public TTMThermalParticle {

 private:

  TTMParameterSetSQ *fParameters;	//pointer to thermal parameters
   
  void UpdateMembers();

 public:

  TTMThermalParticleSQ();
  TTMThermalParticleSQ(TTMParticle *part, TTMParameterSetSQ *parm,
                     Double_t correction);
  TTMThermalParticleSQ(TTMThermalParticleSQ& obj);

  TTMParameterSetSQ* GetParameters() {return fParameters;}
  void SetParameters(TTMParameterSetSQ *x) {fParameters = x;}
  void SetCorrFactor(Double_t x) {fCorrFactor = x;}
  
  TTMThermalParticleSQ& operator=(TTMThermalParticleSQ& obj);
 
  ClassDef(TTMThermalParticleSQ,1) // The baryon-canonical approach

};

#endif
