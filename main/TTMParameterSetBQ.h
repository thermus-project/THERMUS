/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>

#ifndef Thermal_TTMParameter
#include <TTMParameter.h>
#endif

#ifndef Thermal_TTMParameterSet
#include <TTMParameterSet.h>
#endif

#ifndef Thermal_TTMParameterSetBQ
#define Thermal_TTMParameterSetBQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMParameterSetBQ    						//
//									//
// The parameter set used when Baryon# and Charge are to be treated 	//
// grand-canonically, and strangeness canonically 			//   
//     									//
//     									//
//     	 fParArray[0]:T 		Temperature			//
//     	 fParArray[1]:muB 		Baryon Chemical Potential	//
//     	 fParArray[2]:muQ 		Charge Chemical Potential	//
//     	 fParArray[3]:gammas 		Strangeness Suppression Factor	//
//	 fParArray[4]:Can. radius	Correlation radius		//
//     	 fParArray[5]:radius		Fireball Radius			//
//     	 								//
//   The ratio Baryon/(2*Charge) may be used to constrain muQ.		//
//     		 							//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////


class TTMParameterSetBQ:public TTMParameterSet {

 private:

  TTMParameter fParArray[6];

  Bool_t fMuQConstrain;         // true if muQ must be constrained
  Bool_t fCorrRConstrain;	// true if the correlation and fireball 
   				// radii are the same

  Double_t fS;	        	// the required strangeness inside the correlation volume	
  Double_t fB2Q;		// the initial B/2Q ratio

 public:

  TTMParameterSetBQ(Double_t temp, Double_t mub,
                    Double_t muq, Double_t gs, Double_t can_r,
                    Double_t r = 0., Double_t b2q = 0., Double_t S = 0.,
                    Double_t temp_error = 0., Double_t mub_error = 0.,
                    Double_t muq_error = 0., Double_t gs_error = 0.,
                    Double_t can_r_error = 0., Double_t r_error = 0.);

  TTMParameterSetBQ();
  ~TTMParameterSetBQ() { }
   
  Bool_t GetMuQConstrain() const {return fMuQConstrain;}
  Bool_t GetCorrRConstrain() const {return fCorrRConstrain;}

  Double_t GetS() const {return fS;}
  Double_t GetB2Q() const {return fB2Q;}
   
  Double_t GetT() const {return fParArray[0].GetValue();}
  TTMParameter* GetTPar() {return &fParArray[0];}
  Double_t GetMuB() const {return fParArray[1].GetValue();}
  TTMParameter* GetMuBPar() {return &fParArray[1];}
  Double_t GetMuQ() const {return fParArray[2].GetValue();}
  TTMParameter* GetMuQPar() {return &fParArray[2];}
  Double_t GetGammas() const {return fParArray[3].GetValue();}
  TTMParameter* GetGammasPar() {return &fParArray[3];}
  Double_t GetCanRadius() const {return fParArray[4].GetValue();}
  TTMParameter* GetCanRadiusPar() {return &fParArray[4];}
  Double_t GetRadius() const {return fParArray[5].GetValue();}
  TTMParameter* GetRadiusPar() {return &fParArray[5];}

  void SetT(Double_t x) {fParArray[0].SetValue(x);}
  void SetMuB(Double_t x) {fParArray[1].SetValue(x);}
  void SetMuQ(Double_t x) {fParArray[2].SetValue(x);}
  void SetGammas(Double_t x) {fParArray[3].SetValue(x);}
  void SetCanRadius(Double_t x) {fParArray[4].SetValue(x);}
  void SetRadius(Double_t x) {fParArray[5].SetValue(x);}

  void SetS(Double_t x) {fS = x;}
  void SetB2Q(Double_t x) {fB2Q = x;}

  void ConstrainMuQ(Double_t x);
  void ConserveSGlobally();
   
  void FitT(Double_t start, Double_t min = 0.050,
            Double_t max = 0.180, Double_t step = 0.001);
  void FitMuB(Double_t start, Double_t min = 0.000,
              Double_t max = 0.500, Double_t step = 0.001);
  void FitMuQ(Double_t start, Double_t min = -0.200,
              Double_t max = 0.200, Double_t step = 0.001);
  void FitGammas(Double_t start, Double_t min = 0.,
                 Double_t max = 1.5, Double_t step = 0.001);
  void FitCanRadius(Double_t start, Double_t min = 0.1,
                    Double_t max = 15., Double_t step = 0.01);
  void FitRadius(Double_t start, Double_t min = 0.,
                 Double_t max = 15., Double_t step = 0.01);

  void FixT(Double_t value, Double_t error = 0.);
  void FixMuB(Double_t value, Double_t error = 0.);
  void FixMuQ(Double_t value, Double_t error = 0.);
  void FixGammas(Double_t value, Double_t error = 0.);
  void FixCanRadius(Double_t value, Double_t error = 0.);
  void FixRadius(Double_t value, Double_t error = 0.);

  void List();

  TTMParameterSetBQ& operator=(const TTMParameterSetBQ& obj);

  ClassDef(TTMParameterSetBQ, 1) // Canonical strangeness parameter set

};

#endif
