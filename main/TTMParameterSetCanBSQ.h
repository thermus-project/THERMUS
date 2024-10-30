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

#ifndef Thermal_TTMParameterSetCanBSQ
#define Thermal_TTMParameterSetCanBSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMParameterSetCanBSQ    						//
//									//
// The parameter set used when Baryon#, Strangeness and Charge are	//
// to be treated canonically:						//   
//     									//
//     									//
//     	 fParArray[0]:T 	Temperature				//
//     	 fParArray[1]:B 	Total Baryon Number			//
//       fParArray[2]:S         Total Strangeness			//
//     	 fParArray[3]:Q 	Total Charge 				//
//     	 fParArray[4]:gammas 	Strangeness Suppression Factor		//
//     	 fParArray[5]:radius	Fireball Radius				//
//     	 								//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////


class TTMParameterSetCanBSQ:public TTMParameterSet {

 public:

  TTMParameter fParArray[6];

  TTMParameterSetCanBSQ(Double_t temp, Int_t b, Double_t s,
                        Double_t q, Double_t gs, Double_t r = 0.,
                        Double_t temp_error = 0., Double_t b_error = 0.,
                        Double_t s_error = 0., Double_t q_error = 0.,
                        Double_t gs_error = 0., Double_t r_error = 0.);

  TTMParameterSetCanBSQ();
  ~TTMParameterSetCanBSQ() { }
   
  Double_t GetT() const {return fParArray[0].GetValue();}
  TTMParameter* GetTPar() {return &fParArray[0];}
  Int_t GetB() const {return (Int_t)fParArray[1].GetValue();}
  TTMParameter* GetBPar() {return &fParArray[1];}
  Double_t GetS() const {return fParArray[2].GetValue();}
  TTMParameter* GetSPar() {return &fParArray[2];}
  Double_t GetQ() const {return fParArray[3].GetValue();}
  TTMParameter* GetQPar() {return &fParArray[3];}
  Double_t GetGammas() const {return fParArray[4].GetValue();}
  TTMParameter* GetGammasPar() {return &fParArray[4];}
  Double_t GetRadius() const {return fParArray[5].GetValue();}
  TTMParameter* GetRadiusPar() {return &fParArray[5];}

  void SetT(Double_t x) {fParArray[0].SetValue(x);}
  void SetB(Int_t x) {fParArray[1].SetValue(x);}
  void SetS(Double_t x) {fParArray[2].SetValue(x);}
  void SetQ(Double_t x) {fParArray[3].SetValue(x);}
  void SetGammas(Double_t x) {fParArray[4].SetValue(x);}
  void SetRadius(Double_t x) {fParArray[5].SetValue(x);}

  void FitT(Double_t start, Double_t min = 0.050,
            Double_t max = 0.180, Double_t step = 0.001);
  void FitB(Int_t start, Int_t min = 0,
            Int_t max = 250, Int_t step = 1);
  void FitS(Double_t start, Double_t min = 0,
            Double_t max = 250, Double_t step = 1);
  void FitQ(Double_t start, Double_t min = 0,
            Double_t max = 250, Double_t step = 1);
  void FitGammas(Double_t start, Double_t min = 0.,
                 Double_t max = 1.5, Double_t step = 0.001);
  void FitRadius(Double_t start, Double_t min = 0.,
                 Double_t max = 15., Double_t step = 0.01);

  void FixT(Double_t value, Double_t error = 0.);
  void FixB(Int_t value, Double_t error = 0.);
  void FixS(Double_t value, Double_t error = 0.);
  void FixQ(Double_t value, Double_t error = 0.);
  void FixGammas(Double_t value, Double_t error = 0.);
  void FixRadius(Double_t value, Double_t error = 0.);

  void List();

  TTMParameterSetCanBSQ& operator=(const TTMParameterSetCanBSQ& obj);

  ClassDef(TTMParameterSetCanBSQ, 1) // Canonical BSQ parameter set

};

#endif
