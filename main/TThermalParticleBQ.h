/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <math.h>

#ifndef Thermal_TTMThermalParticle
#include <TThermalParticle.h>
#endif

#ifndef Thermal_TTMParameterSetBQ
#include <TTMParameterSetBQ.h>
#endif

#ifndef Thermal_TTMThermalParticleBQ
#define Thermal_TTMThermalParticleBQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalParticleBQ                                        		//
//                                                                  	//
// A thermal particle in the formalism treating strangeness canonically	//
// and charge and baryon number grand-canonically. The Boltzmann 	//
// approximation is used since in this case the correction to the 	//
// grand-canonical results amounts to a simple multiplicative factor.	//
// The entropy cannot be split into the sum of particle entropies, so	//
// no entropy calculation is made here.				       	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalParticleBQ:public TTMThermalParticle {

 private:

  TTMParameterSetBQ *fParameters;	//pointer to thermal parameters
   
  void UpdateMembers();

 public:

  TTMThermalParticleBQ();
  TTMThermalParticleBQ(TTMParticle *part, TTMParameterSetBQ *parm,
                     Double_t correction);
  TTMThermalParticleBQ(TTMThermalParticleBQ& obj);

  TTMParameterSetBQ* GetParameters() {return fParameters;}
  void SetParameters(TTMParameterSetBQ *x) {fParameters = x;}
  void SetCorrFactor(Double_t x) {fCorrFactor = x;}
  
  TTMThermalParticleBQ& operator=(TTMThermalParticleBQ& obj);
 
  ClassDef(TTMThermalParticleBQ,1) // The strangeness-canonical approach

};

#endif
