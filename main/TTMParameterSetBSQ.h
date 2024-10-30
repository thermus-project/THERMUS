/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 26 April 2010 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>

#ifndef Thermal_TTMParameter
#include <TTMParameter.h>
#endif

#ifndef Thermal_TTMParameterSet
#include <TTMParameterSet.h>
#endif

#ifndef Thermal_TTMParameterSetBSQ
#define Thermal_TTMParameterSetBSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMParameterSetBSQ    						//
// 									//
// The parameter set used when Baryon Content, Strangeness, Charge, 	//
// Charm & Beauty are to be treated grand-canonically:			//   
//     									//
//     									//
//     	 fParArray[0]:T 		Temperature			//
//     	 fParArray[1]:muB 		Baryon Chemical Potential	//
//     	 fParArray[2]:muS 		Strangeness Chemical Potential	//
//     	 fParArray[3]:muQ 		Charge Chemical Potential	//
//     	 fParArray[4]:gammas 		Strangeness Suppression Factor	//
//     	 fParArray[5]:radius		Fireball Radius			//
//     	 fParArray[6]:muC 		Charm Chemical Potential	//
//       fParArray[7]:gammac		Charm Suppression Factor	//
//	 fParArray[8]:mub		Beauty Chemical Potential	//
//	 fParArray[9]:gammab		Beauty Suppression Factor	//
//     	 								//
//   The initial strangeness density may be used to constrain muS.	//
//   The ratio Baryon/(2*Charge) may be used to constrain muQ.		//
//   The initial charm density may be used to constrain muC.		//
//   The initial beauty density may be used to constrain mub.		//
//     		 							//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMParameterSetBSQ:public TTMParameterSet {

 private:

  TTMParameter fParArray[10];

  Bool_t fMuSConstrain;        // true if muS must be constrained
  Bool_t fMuQConstrain;        // true if muQ must be constrained
  Bool_t fMuCConstrain;	       // true if muC must be constrained
  Bool_t fMubConstrain;	       // true if mub must be constrained

  Double_t fSDens;	       // the initial strangeness density
  Double_t fCDens;	       // the initial charm density 	
  Double_t fbDens;	       // the initial beauty density 	
  Double_t fB2Q;	       // the initial B/2Q ratio

 public:

  TTMParameterSetBSQ(Double_t temp, Double_t mub,
                     Double_t mus, Double_t muq, Double_t gs,
                     Double_t r = 0., Double_t muc = 0., Double_t gc = 1., Double_t mubeauty = 0., Double_t gb = 1., 
		     Double_t b2q = 0.,
                     Double_t s = 0., Double_t c = 0., Double_t beauty = 0., Double_t temp_error = 0.,
                     Double_t mub_error = 0., Double_t mus_error = 0.,
                     Double_t muq_error = 0., Double_t gs_error = 0.,
                     Double_t r_error = 0., Double_t muc_error = 0., Double_t gc_error = 0., 
                     Double_t mubeauty_error = 0., Double_t gb_error = 0.);

  TTMParameterSetBSQ();
  ~TTMParameterSetBSQ() { } 
   
  Bool_t GetMuSConstrain() const {return fMuSConstrain;}
  Bool_t GetMuQConstrain() const {return fMuQConstrain;}
  Bool_t GetMuCConstrain() const {return fMuCConstrain;}
  Bool_t GetMubConstrain() const {return fMubConstrain;}

  Double_t GetSDens() const {return fSDens;}
  Double_t GetCDens() const {return fCDens;}
  Double_t GetbDens() const {return fbDens;}
  Double_t GetB2Q() const {return fB2Q;}

  Double_t GetT() const {return fParArray[0].GetValue();}
  TTMParameter* GetTPar() {return &fParArray[0];}
  Double_t GetMuB() const {return fParArray[1].GetValue();}
  TTMParameter* GetMuBPar() {return &fParArray[1];}
  Double_t GetMuS() const {return fParArray[2].GetValue();}
  TTMParameter* GetMuSPar() {return &fParArray[2];}
  Double_t GetMuQ() const {return fParArray[3].GetValue();}
  TTMParameter* GetMuQPar() {return &fParArray[3];}
  Double_t GetGammas() const {return fParArray[4].GetValue();}
  TTMParameter* GetGammasPar() {return &fParArray[4];}
  Double_t GetRadius() const {return fParArray[5].GetValue();}
  TTMParameter* GetRadiusPar() {return &fParArray[5];}
  Double_t GetMuC() const {return fParArray[6].GetValue();}
  TTMParameter* GetMuCPar() {return &fParArray[6];}
  Double_t GetGammac() const {return fParArray[7].GetValue();}
  TTMParameter* GetGammacPar() {return &fParArray[7];}
  Double_t GetMub() const {return fParArray[8].GetValue();}
  TTMParameter* GetMubPar() {return &fParArray[8];}
  Double_t GetGammab() const {return fParArray[9].GetValue();}
  TTMParameter* GetGammabPar() {return &fParArray[9];}


  void SetT(Double_t x) {fParArray[0].SetValue(x);}
  void SetMuB(Double_t x) {fParArray[1].SetValue(x);}
  void SetMuS(Double_t x) {fParArray[2].SetValue(x);}
  void SetMuQ(Double_t x) {fParArray[3].SetValue(x);}
  void SetGammas(Double_t x) {fParArray[4].SetValue(x);}
  void SetRadius(Double_t x) {fParArray[5].SetValue(x);}
  void SetMuC(Double_t x) {fParArray[6].SetValue(x);}
  void SetGammac(Double_t x) {fParArray[7].SetValue(x);}
  void SetMub(Double_t x) {fParArray[8].SetValue(x);}
  void SetGammab(Double_t x) {fParArray[9].SetValue(x);}

  void SetSDens(Double_t x) {fSDens = x;}
  void SetCDens(Double_t x) {fCDens = x;}
  void SetbDens(Double_t x) {fbDens = x;}
  void SetB2Q(Double_t x) {fB2Q = x;}

  void ConstrainMuS(Double_t x = 0.);	
  void ConstrainMuQ(Double_t x);
  void ConstrainMuC(Double_t x = 0.);
  void ConstrainMub(Double_t x = 0.);

  void FitT(Double_t start, Double_t min = 0.050,
            Double_t max = 0.180, Double_t step = 0.001);
  void FitMuB(Double_t start, Double_t min = 0.000,
              Double_t max = 0.500, Double_t step = 0.001);
  void FitMuS(Double_t start, Double_t min = 0.000,
              Double_t max = 0.500, Double_t step = 0.001);
  void FitMuQ(Double_t start, Double_t min = -0.200,
              Double_t max = 0.200, Double_t step = 0.001);
  void FitGammas(Double_t start, Double_t min = 0.,
                 Double_t max = 1.5, Double_t step = 0.001);
  void FitRadius(Double_t start, Double_t min = 0.,
                 Double_t max = 15., Double_t step = 0.01);
  void FitMuC(Double_t start, Double_t min = -0.200,
                 Double_t max = 0.200, Double_t step = 0.001);
  void FitGammac(Double_t start, Double_t min = 0.,
                 Double_t max = 1.5, Double_t step = 0.001);
  void FitMub(Double_t start, Double_t min = -0.200,
                 Double_t max = 0.200, Double_t step = 0.001);
  void FitGammab(Double_t start, Double_t min = 0.,
                 Double_t max = 1.5, Double_t step = 0.001);

  void FixT(Double_t value, Double_t error = 0.);
  void FixMuB(Double_t value, Double_t error = 0.);
  void FixMuS(Double_t value, Double_t error = 0.);
  void FixMuQ(Double_t value, Double_t error = 0.);
  void FixGammas(Double_t value, Double_t error = 0.);
  void FixRadius(Double_t value, Double_t error = 0.);
  void FixMuC(Double_t value, Double_t error = 0.);
  void FixGammac(Double_t value, Double_t error = 0.);
  void FixMub(Double_t value, Double_t error = 0.);
  void FixGammab(Double_t value, Double_t error = 0.);

  void List();

  TTMParameterSetBSQ& operator=(const TTMParameterSetBSQ& obj);

  ClassDef(TTMParameterSetBSQ, 1) // Grand-canonical parameter set

};

#endif
