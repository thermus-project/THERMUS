/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>

#ifndef Thermal_TTMThermalModel
#include <TThermalModel.h>
#endif

#ifndef Thermal_TTMThermalParticleBQ
#include <TThermalParticleBQ.h>
#endif

#ifndef Thermal_TTMThermalModelBQ
#define Thermal_TTMThermalModelBQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalModelBQ							//
//            	                                 			//
// Strangeness Canonical Thermal Model Class				// 
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalModelBQ:public TTMThermalModel {

 private:

  TTMParameterSetBQ *fParm;	// pointer to thermal parameters
  
  Double_t flnZtot;	// log(total partition function)
  Double_t flnZ0;	// log(non-strange part of partition function)
  Double_t fExactMuS;   // exact strangeness chemical potential
  Double_t fCorrP1; 	// canonical correction factor for S=+1 particles
  Double_t fCorrP2; 	// canonical correction factor for S=+2 particles
  Double_t fCorrP3;	// canonical correction factor for S=+3 particles
  Double_t fCorrM1; 	// canonical correction factor for S=-1 particles
  Double_t fCorrM2; 	// canonical correction factor for S=-2 particles
  Double_t fCorrM3;	// canonical correction factor for S=-3 particles

  Bool_t fNonStrangeQStats;	// true if S=0 hadrons are to be treated 
				// with quantum statistics


  void Term(Double_t *x, Double_t *y, Int_t m, Int_t n, Double_t *t);

 public:

  Int_t PrimPartDens();

  TTMThermalModelBQ(TTMParticleSet *particles, TTMParameterSetBQ *parameters,
                  Bool_t width = true);
  TTMThermalModelBQ();

  friend void BQfuncQ(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncEN(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncQEN(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncQPercolation(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncST3(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncQST3(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncBDens(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncQBDens(Int_t n, Float_t x[], Float_t f[]);
  friend void BQfuncQNetBDens(Int_t n, Float_t x[], Float_t f[]);

  void SetParameterSet(TTMParameterSetBQ* parm){fParm = parm;}
  TTMParameterSetBQ* GetParameterSet() {return fParm;} 
   
  Int_t GenerateParticleDens();
  void GenerateEnergyDens();
  void GenerateEntropyDens();
  void GeneratePressure();

  Int_t ConstrainEoverN(Double_t eovern);
  Int_t ConstrainSoverT3(Double_t SoverT3);
  Int_t ConstrainTotalBaryonDensity(Double_t nb);
  Int_t ConstrainNetBaryonDensity(Double_t nb);
  Int_t ConstrainPercolation();

  Double_t GetlnZtot() const {return flnZtot;}
  Double_t GetlnZ0() const {return flnZ0;}
  Double_t GetExactMuS() const {return fExactMuS;}
   
  Double_t GetCorrP1() const {return fCorrP1;}
  Double_t GetCorrP2() const {return fCorrP2;}
  Double_t GetCorrP3() const {return fCorrP3;}
  Double_t GetCorrM1() const {return fCorrM1;}
  Double_t GetCorrM2() const {return fCorrM2;}
  Double_t GetCorrM3() const {return fCorrM3;}

  Bool_t GetNonStrangeQStats() const {return fNonStrangeQStats;}
  void SetNonStrangeQStats(Bool_t x){fNonStrangeQStats = x;}
  
  void ListInfo();

  ClassDef(TTMThermalModelBQ,1) // Strangeness-canonical thermal model class

};

#endif
