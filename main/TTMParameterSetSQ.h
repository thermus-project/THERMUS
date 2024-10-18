/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

//Authors: Natasha Sharma and Jean Cleymans 2 February 2021 (modification of the corresponding
//BQ file, TTMParameterSetBQ.h,  by Spencer Wheaton 14 July 2004)

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>

#ifndef Thermal_TTMParameter
#include <TTMParameter.h>
#endif

#ifndef Thermal_TTMParameterSet
#include <TTMParameterSet.h>
#endif

#ifndef Thermal_TTMParameterSetSQ
#define Thermal_TTMParameterSetSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMParameterSetSQ    						//
//									//
// The parameter set used when Strangeness# and Charge are to be treated 	//
// grand-canonically, and baryon canonically 			//   
//     									//
//     									//
//     	 fParArray[0]:T 		Temperature			//
//     	 fParArray[1]:muS 		Strangeness Chemical Potential	//
//     	 fParArray[2]:muQ 		Charge Chemical Potential	//
//     	 fParArray[3]:gammas 		Strangeness Suppression Factor	//
//	 fParArray[4]:Can. radius	Correlation radius		//
//     	 fParArray[5]:radius		Fireball Radius			//
//     	 								//
//                                                       		//
//     		 							//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////


class TTMParameterSetSQ:public TTMParameterSet {

 private:

  TTMParameter fParArray[6];

//  
  Bool_t fCorrRConstrain;	// true if the correlation and fireball 
   				// radii are the same

  Double_t fB;	        	// the required baryon inside the correlation volume	
//  

 public:

  TTMParameterSetSQ(Double_t temp, Double_t mus,
                    Double_t muq, Double_t gs, Double_t can_r,
            //       
                    Double_t r = 0., Double_t B = 0.,
                    Double_t temp_error = 0., Double_t mus_error = 0.,
                    Double_t muq_error = 0., Double_t gs_error = 0.,
                    Double_t can_r_error = 0., Double_t r_error = 0.);

  TTMParameterSetSQ();
  ~TTMParameterSetSQ() { }
   
// 
  Bool_t GetCorrRConstrain() const {return fCorrRConstrain;}

  Double_t GetB() const {return fB;}
//  
   
  Double_t GetT() const {return fParArray[0].GetValue();}
  TTMParameter* GetTPar() {return &fParArray[0];}
  Double_t GetMuS() const {return fParArray[1].GetValue();}
  TTMParameter* GetMuSPar() {return &fParArray[1];}
  Double_t GetMuQ() const {return fParArray[2].GetValue();}
  TTMParameter* GetMuQPar() {return &fParArray[2];}
  Double_t GetGammas() const {return fParArray[3].GetValue();}
  TTMParameter* GetGammasPar() {return &fParArray[3];}
  Double_t GetCanRadius() const {return fParArray[4].GetValue();}
  TTMParameter* GetCanRadiusPar() {return &fParArray[4];}
  Double_t GetRadius() const {return fParArray[5].GetValue();}
  TTMParameter* GetRadiusPar() {return &fParArray[5];}

  void SetT(Double_t x) {fParArray[0].SetValue(x);}
  void SetMuS(Double_t x) {fParArray[1].SetValue(x);}
  void SetMuQ(Double_t x) {fParArray[2].SetValue(x);}
  void SetGammas(Double_t x) {fParArray[3].SetValue(x);}
  void SetCanRadius(Double_t x) {fParArray[4].SetValue(x);}
  void SetRadius(Double_t x) {fParArray[5].SetValue(x);}

  void SetB(Double_t x) {fB = x;}
//  

  void ConstrainMuQ(Double_t x);
  void ConserveBGlobally();
   
  void FitT(Double_t start, Double_t min = 0.050,
            Double_t max = 0.180, Double_t step = 0.001);
  void FitMuS(Double_t start, Double_t min = 0.000,
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
  void FixMuS(Double_t value, Double_t error = 0.);
  void FixMuQ(Double_t value, Double_t error = 0.);
  void FixGammas(Double_t value, Double_t error = 0.);
  void FixCanRadius(Double_t value, Double_t error = 0.);
  void FixRadius(Double_t value, Double_t error = 0.);

  void List();

  TTMParameterSetSQ& operator=(const TTMParameterSetSQ& obj);

  ClassDef(TTMParameterSetSQ, 1) // Canonical baryon parameter set

};

#endif
