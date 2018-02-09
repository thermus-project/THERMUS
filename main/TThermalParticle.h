/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>

#ifndef Thermal_TTMParticle
#include <TTMParticle.h>
#endif

#ifndef Thermal_TTMParameterSet
#include <TTMParameterSet.h>
#endif

#ifndef Thermal_TTMThermalParticle
#define Thermal_TTMThermalParticle

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalParticle                                          		//
//                                                                  	//
// Base class 	 				     	            	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalParticle:public TObject {

 protected:

  TTMParticle *fParticle; 	// pointer to particle object
   
  Double_t fCorrFactor;	// w.r.t. grand-canonical results (in the
                        // Boltzmann approximation there is a simple
                        // multiplicative correction factor for e, n 
                        // and P)
   
  Double_t fG;			// heavy flavour suppression 
  Double_t fDeg;		// copy of particle's degeneracy
  Double_t fM;			// copy of particle's mass
  Double_t fT;			// copy of temperature
  Double_t fMu;			// copy of particle's chemical potential

  const static Double_t eps;	// required accuracy of integrals

 public:

  TTMThermalParticle();

  virtual TTMParameterSet* GetParameters() = 0;
  virtual void UpdateMembers() = 0;
   
  Double_t DensityBoltzmannNoWidth();
  Double_t DensityBoltzmannWidth();

  Double_t EnergyBoltzmannNoWidth();
  Double_t EnergyBoltzmannWidth();

  Double_t PressureBoltzmannNoWidth();
  Double_t PressureBoltzmannWidth();

  Double_t GetCorrFactor() const {return fCorrFactor;}
  TTMParticle* GetParticle() const {return fParticle;}
  TString GetPartName() const {return fParticle->GetPartName();}
  Int_t GetID() const {return fParticle->GetID();}
  TList* GetDecays() const {return fParticle->GetDecaySummary();};

  void SetParticle(TTMParticle *x) {fParticle = x;}
   
  ClassDef(TTMThermalParticle,1) // Base Class 

};

#endif
