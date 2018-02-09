/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>

using namespace std;

#ifndef Thermal_TTMDensObj
#define Thermal_TTMDensObj

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMDensObj                                          			//
//                                                                  	//
// Object containing densities for storage in container class. 		// 
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMDensObj:public TNamed {

 private:

  Int_t fID;                   // particle ID
  Double_t fEnergy;            // primordial energy density
  Double_t fEntropy;           // primordial entropy density
  Double_t fPressure;          // primordial pressure
  Double_t fDensity;           // primordial particle density
  Double_t fDecayDensity;      // decay density contribution

 public:

  TTMDensObj();
  TTMDensObj(Int_t x);
  ~TTMDensObj() { } 
   
  void List();

  void SetID(Int_t x);

  void SetPrimEnergy(Double_t x) {fEnergy = x;}
  void SetPrimEntropy(Double_t x) {fEntropy = x;}
  void SetPrimPressure(Double_t x) {fPressure = x;}
  void SetPrimDensity(Double_t x) {fDensity = x;}
  void SetDecayDensity(Double_t x) {fDecayDensity = x;}

  Int_t GetID() const {return fID;}
  Double_t GetPrimEnergy() const {return fEnergy;}
  Double_t GetPrimEntropy() const {return fEntropy;}
  Double_t GetPrimPressure() const {return fPressure;}
  Double_t GetPrimDensity() const {return fDensity;}
  Double_t GetDecayDensity() const {return fDecayDensity;}
  Double_t GetFinalDensity() const {return fDensity + fDecayDensity;}
   
  TTMDensObj& operator=(const TTMDensObj& obj);

  ClassDef(TTMDensObj,1) // Density Object for storage in ROOT container class

};

#endif
