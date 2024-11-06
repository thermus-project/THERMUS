/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <math.h>

#ifndef Thermal_TTMThermalParticle
#include <TThermalParticle.h>
#endif

#ifndef Thermal_TTMParameterSetCanBSQ
#include <TTMParameterSetCanBSQ.h>
#endif

#ifndef Thermal_TTMThermalParticleCanBSQ
#define Thermal_TTMThermalParticleCanBSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalParticleCanBSQ                                        	//
//                                                                  	//
// A thermal particle in the formalism treating strangeness, charge 	//
// and baryon number canonically. The Boltzmann approximation is used	//
// since in this case the correction to the grand-canonical results 	//
// amounts to a simple multiplicative factor. The entropy cannot be 	//
// split into the sum of particle entropies, so no entropy calculation	//
// is made here.						       	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalParticleCanBSQ:public TTMThermalParticle {

 private:

  TTMParameterSetCanBSQ *fParameters;	//pointer to thermal parameters
   
  void UpdateMembers();

 public:

  TTMThermalParticleCanBSQ();
  TTMThermalParticleCanBSQ(TTMParticle *part, TTMParameterSetCanBSQ *parm,
                         Double_t correction);
  TTMThermalParticleCanBSQ(TTMThermalParticleCanBSQ& obj);

  TTMParameterSetCanBSQ* GetParameters() {return fParameters;}
  void SetParameters(TTMParameterSetCanBSQ *x) {fParameters = x;}
  void SetCorrFactor(Double_t x) {fCorrFactor = x;}
  
  TTMThermalParticleCanBSQ& operator=(TTMThermalParticleCanBSQ& obj);
 
  ClassDef(TTMThermalParticleCanBSQ,1) // The canonical BSQ approach

};

#endif
