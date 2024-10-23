/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 26 April 2010 //

#include <TROOT.h>
#include <TObject.h>

#ifndef Thermal_TTMThermalModel
#include <TThermalModel.h>
#endif

#ifndef Thermal_TTMThermalParticleBSQ
#include <TThermalParticleBSQ.h>
#endif

#ifndef Thermal_TTMThermalModelBSQ
#define Thermal_TTMThermalModelBSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalModelBSQ                                  			//
//                                                                  	//
// Grand Canonical Thermal Model Class					//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalModelBSQ:public TTMThermalModel {

 private:

  TTMParameterSetBSQ *fParm;	 // pointer to thermal parameters
   
  Bool_t fQStats;		 // true if Quantum Stats required 

  Bool_t fExclVolCorrection;	 // true if excluded volume is to be taken into account 
  Double_t fExclVolPressure;     // excluded volume pressure responsible for shift in mu's
  Double_t fExclVolDenominator;  // denominator correction in n, e and s 

  void CalcExclVolPressure();

  TTMThermalParticleBSQ* ShiftedParticle(TTMParticle* CopyPart, 
                                         TTMParameterSetBSQ* CopyParm, Double_t P);


 public:

  Double_t ExclVolShiftedPressure(Double_t x);
  Int_t PrimPartDens();
  TTMThermalModelBSQ();
  TTMThermalModelBSQ(TTMParticleSet *particles, 
                     TTMParameterSetBSQ *parameters, Bool_t qstats = true, 
                     Bool_t width = true);

  friend void BSQfuncS(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncQ(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQ(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQC(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQCb(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncBSQ(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncQQ(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQEN(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQPercolation(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSPercolation(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSEN(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQST3(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSST3(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSBDens(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQBDens(Int_t n, Float_t x[], Float_t f[]);
  friend void BSQfuncSQNetBDens(Int_t n, Float_t x[], Float_t f[]);

  friend Double_t FindExclVolPressure(TTMThermalModelBSQ *model, Double_t limit);
  friend Float_t ExclVolPressureFunc(Float_t x);

  Int_t ConstrainEoverN(Double_t EoverN);
  Int_t ConstrainSoverT3(Double_t SoverT3);
  Int_t ConstrainTotalBaryonDensity(Double_t nb);
  Int_t ConstrainNetBaryonDensity(Double_t nb);
  Int_t ConstrainBSQ(Double_t B, Double_t S, Double_t Q);
  Int_t ConstrainQ(Double_t Q);
  Int_t ConstrainPercolation();

  void SetQStats(Bool_t x) {fQStats = x;}
  Bool_t GetQStats() const {return fQStats;}

  void SetExcludedVolume(Bool_t x) {fExclVolCorrection = x;}
  Bool_t GetExcludedVolume() const {return fExclVolCorrection;}

  Double_t GetExclVolPressure() const {return fExclVolPressure;}
  Double_t GetExclVolDenominator() const {return fExclVolDenominator;}

  void SetParameterSet(TTMParameterSetBSQ *parameters) {fParm = parameters;}
   
  TTMParameterSetBSQ *GetParameterSet() {return fParm;} 
   
  Int_t GenerateParticleDens();
  void GenerateEnergyDens();
  void GenerateEntropyDens();
  void GeneratePressure();

  void ListInfo();

  ClassDef(TTMThermalModelBSQ,1) // Grand-canonical thermal model class

};

#endif
