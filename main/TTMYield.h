/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>

using namespace std;

#ifndef Thermal_TTMParticleSet
#include <TTMParticleSet.h>
#endif

#ifndef Thermal_TTMYield
#define Thermal_TTMYield

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMYield                                               		//
//  									//
// Object for storage of experimental yields and model predictions 	//
// together with the decay chain relevant to the measurement in the 	//
// form of a particle set					        //
//									//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMYield:public TNamed {

 private:

   Int_t fID1;                  // id for yield; numerator id for ratio
   Int_t fID2;                  // denominator id for ratio   
   
   Bool_t fFit;			// true if yield is to be fitted (predicted if false)
   
   TTMParticleSet *fSet1;	// particle set relevant to yield (numerator) 
   TTMParticleSet *fSet2;	// particle set relevant to denominator
   
   Double_t fExpValue;          // experimental value
   Double_t fExpError;          // experimental error
   Double_t fModelValue;        // model value
   Double_t fModelError;        // model error

 public:

   TTMYield();
   TTMYield(TString name, Double_t exp_val, Double_t exp_err,
             Int_t id1, Int_t id2 = 0, Bool_t fit = true);

   void SetTMName(TString x) {fName = x;}

   void Fit() {fFit = true;}
   void Predict() {fFit = false;}

   void SetID(Int_t x, Int_t y) 
   {
      fID1 = x;
      fID2 = y;
   }
   
   void SetPartSet(TTMParticleSet *x, TTMParticleSet *y = (TTMParticleSet *) 0) 
   {
      fSet1 = x;
      fSet2 = y;
   }
   
   TTMParticleSet* GetPartSet1() const {return fSet1;}
   TTMParticleSet* GetPartSet2() const {return fSet2;}
   
   Bool_t GetFit() const {return fFit;}
   
   void SetExpValue(Double_t x) {fExpValue = x;}
   void SetExpError(Double_t x) {fExpError = x;}
   void SetModelValue(Double_t x) {fModelValue = x;}
   void SetModelError(Double_t x) {fModelError = x;}

   TString GetTMName() {return fName;}
   Int_t GetID1() const {return fID1;}
   Int_t GetID2() const {return fID2;}
   Double_t GetExpValue() const {return fExpValue;}
   Double_t GetExpError() const {return fExpError;}
   Double_t GetModelValue() const {return fModelValue;}
   Double_t GetModelError() const {return fModelError;}

   Double_t GetStdDev() const
   {
      return (fModelValue - fExpValue) / fExpError;
   }

   Double_t GetQuadDev() const
   {
      return (fModelValue - fExpValue) / fModelValue;
   }

   void List();

   TTMYield operator=(TTMYield obj);

   ClassDef(TTMYield, 1) // Object for storage of model and experimental yields

};

#endif
