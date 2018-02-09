/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <fstream>
#include <iostream>

using namespace std;

#ifndef Thermal_TTMDecay
#define Thermal_TTMDecay

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMDecay                                                  		//
//                                                                  	//
// Stores parent & daughter id's and Branching ratio in an object for	// 
// storage in a ROOT container class				 	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMDecay:public TObject {

 private:

  Int_t fParentID;             // Parent ID
  Int_t fDaughterID;           // Daughter ID
  Double_t fBRatio;            // Branching ratio (fraction i.e. NOT %!)

 public:

  TTMDecay(); 
  TTMDecay(Int_t parent, Int_t daughter, Double_t fraction);
  ~TTMDecay() { }
   
  Int_t GetParentID() const {return fParentID;}
  Int_t GetDaughterID() const {return fDaughterID;}
  Double_t GetBRatio() const {return fBRatio;}

  void SetParentID(Int_t a) {fParentID = a;}
  void SetDaughterID(Int_t b) {fDaughterID = b;}
  void SetBRatio(Double_t x) {fBRatio = x;}

  void List();

  TTMDecay& operator=(const TTMDecay& obj);
 
  ClassDef(TTMDecay,1) // Decay object for storage in a ROOT container class 

};

#endif

    
