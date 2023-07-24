/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

//Authors: Natasha Sharma and Jean Cleymans 2 February 2021 (modification of the corresponding
//BQ file, TThermalModelBQ.h,  by Spencer Wheaton 14 July 2004)

#include <TROOT.h>
#include <TObject.h>

#ifndef Thermal_TTMThermalModel
#include <TThermalModel.h>
#endif

#ifndef Thermal_TTMThermalParticleSQ
#include <TThermalParticleSQ.h>
#endif

#ifndef Thermal_TTMThermalModelSQ
#define Thermal_TTMThermalModelSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalModelSQ							//
//            	                                 			//
// Baryon Canonical Thermal Model Class				// 
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalModelSQ:public TTMThermalModel {

 private:

  TTMParameterSetSQ *fParm;	// pointer to thermal parameters
  
  Double_t flnZtot;	// log(total partition function)
  Double_t flnZ0;	// log(non-baryonic part of partition function)
  Double_t fExactMuB;   // exact baryon chemical potential
  Double_t fCorrP1; 	// canonical correction factor for B=+1 particles
  Double_t fCorrP2; 	// canonical correction factor for B=+2 particles
  Double_t fCorrP3;	// canonical correction factor for B=+3 particles
  Double_t fCorrM1; 	// canonical correction factor for B=-1 particles
  Double_t fCorrM2; 	// canonical correction factor for B=-2 particles
  Double_t fCorrM3;	// canonical correction factor for B=-3 particles

  Bool_t fNonBaryonQStats;	// true if B=0 hadrons are to be treated 
				// with quantum statistics


  void Term(Double_t *x, Double_t *y, Int_t m, Int_t n, Double_t *t);

 public:

  Int_t PrimPartDens();

  TTMThermalModelSQ(TTMParticleSet *particles, TTMParameterSetSQ *parameters,
                  Bool_t width = true);
  TTMThermalModelSQ();

  friend void SQfuncQ(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncEN(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncQEN(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncQPercolation(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncST3(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncQST3(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncBDens(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncQBDens(Int_t n, Float_t x[], Float_t f[]);
  friend void SQfuncQNetBDens(Int_t n, Float_t x[], Float_t f[]);

  void SetParameterSet(TTMParameterSetSQ* parm){fParm = parm;}
  TTMParameterSetSQ* GetParameterSet() {return fParm;} 
   
  Int_t GenerateParticleDens();
  void GenerateEnergyDens();
  void GenerateEntropyDens();
  void GeneratePressure();

  Int_t ConstrainEoverN(Double_t eovern);
  Int_t ConstrainSoverT3(Double_t SoverT3);
  Int_t ConstrainTotalBaryonDensity(Double_t nb);
  Int_t ConstrainNetStrangeDensity(Double_t nb);  // Int_t ConstrainNetBaryonDensity(Double_t nb);
  Int_t ConstrainPercolation();

  Double_t GetlnZtot() const {return flnZtot;}
  Double_t GetlnZ0() const {return flnZ0;}
  Double_t GetExactMuB() const {return fExactMuB;}
   
  Double_t GetCorrP1() const {return fCorrP1;}
  Double_t GetCorrP2() const {return fCorrP2;}
  Double_t GetCorrP3() const {return fCorrP3;}
  Double_t GetCorrM1() const {return fCorrM1;}
  Double_t GetCorrM2() const {return fCorrM2;}
  Double_t GetCorrM3() const {return fCorrM3;}

  Bool_t GetNonBaryonQStats() const {return fNonBaryonQStats;}
  void SetNonBaryonQStats(Bool_t x){fNonBaryonQStats = x;}
  
  void ListInfo();

  ClassDef(TTMThermalModelSQ,1) // Baryon-canonical thermal model class

};

#endif
