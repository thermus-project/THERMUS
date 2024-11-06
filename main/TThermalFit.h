/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TMinuit.h>
#include <fstream>
#include <iostream>

#ifndef Thermal_TTMParticleSet
#include <TTMParticleSet.h>
#endif

#ifndef Thermal_TTMParameterSet
#include <TTMParameterSet.h>
#endif

#ifndef Thermal_TTMThermalModel
#include <TThermalModel.h>
#endif

#ifndef Thermal_TTMYield
#include <TTMYield.h>
#endif

#ifndef Thermal_TTMThermalFit
#define Thermal_TTMThermalFit

using namespace std;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalFit							//
//            	                                 			//
// Base class for thermal fit objects. 					// 
//									//
//////////////////////////////////////////////////////////////////////////

class TTMThermalFit:public TObject {

 protected:

  TTMParticleSet *fPartSet;	// pointer to BASE Particle Set
  TList *fYields; 		// Container for yields of interest 
  Double_t fChiSquare;		// chi-squared 
  Double_t fQuadDev;		// quadratic deviation
  TMinuit *fMinuit;		//! pointer to TMinuit obj (not written to file)

  TString fDescriptor;		// string describing fit

  TString GetTMName(Int_t id1,Int_t id2=0,TString descr="");

  virtual TTMThermalModel* GenerateThermalModel(TTMParticleSet *set) = 0;

 public:

  TTMThermalFit();  
   
  ~TTMThermalFit();
   
  friend void Minuit_fcn(Int_t & npar, Double_t * gin, Double_t & f,
                         Double_t * par, Int_t iflag);
  
  virtual TTMParameterSet* GetParameterSet() = 0;

  void InputExpYields(const char *file);
  void FitData(Int_t flag = 0);
  void GenerateYields();

  TList* GetYields(){return fYields;}	
  Double_t GetChiSquare(){return fChiSquare;}
  Double_t GetQuadDev(){return fQuadDev;}
  TMinuit* GetMinuit(){return fMinuit;}

  TString GetDescriptor() const {return fDescriptor;}

  void SetMinuit(TMinuit *x){fMinuit = x;}
  void AddYield(TTMYield *yield);
  void RemoveYield(Int_t id1,Int_t id2=0,TString descr = "");

  TTMYield* GetYield(Int_t id1,Int_t id2=0,TString descr = "");
  void SetParticleSet(TTMParticleSet* set){fPartSet = set;}

  TTMParticleSet* GetParticleSet() {return fPartSet;}
  void ListYields();
  void ListMinuitInfo();

  ClassDef(TTMThermalFit,1) // Base class for thermal fit objects

};

#endif
