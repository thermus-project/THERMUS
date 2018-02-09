/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 26 April 2010 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <TList.h>
#include <fstream>
#include <iostream>

using namespace std;

TString Int_2_String(Int_t x);
class TTMDecay;
class TTMDecayChannel;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMParticle                                                  	//
//                                                                  	//
// Object containing particle properties relevant to thermal model.	//
// Inherits from TNamed to allow TTMParticle objects to be stored in 	//
// ROOT hash tables. The particle's Monte Carlo ID is used to set	// 
// fName.		 				 		//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

#ifndef Thermal_TTMParticle
#define Thermal_TTMParticle

class TTMParticle: public TNamed{

 friend class TTMParticleSet;

 private:

  TString fPartName;           // particle name (not fName!!!!!)
  Int_t fID;                   // Monte Carlo ID
  Bool_t fStable;              // true if stable
  Int_t fDeg;                  // spin degeneracy 
  Int_t fStat;                 // 0:Boltzmann +1:Fermi-Dirac -1:Bose-Einstein 
  Int_t fB;                    // baryon number 
  Int_t fS;                    // strangeness
  Int_t fQ;                    // charge
  Int_t fCharm;		       // charm
  Int_t fBeauty;	       // beauty
  Int_t fTop;		       // top
  Double_t fMass;              // mass 
  Double_t fWidth;             // width
  Double_t fThreshold;         // threshold 

  Double_t fRadius;	       // hard sphere excluded volume radius

  Double_t fThresholdCalc;     // calculated threshold	
  Bool_t fThresholdFlag;       // true if threshold has been calculated

  Double_t fSContent;          // |S_i|= #s + #s-bar quarks
  Double_t fCContent;	       // |C_i|= #c + #c-bar quarks
  Double_t fbContent;	       // |b_i|= #b + #b-bar quarks
  Double_t fTContent;	       // |T_i|= #t + #t-bar quarks

  TList* fDecaySummary;        // pointer to summary list of decays
  TList* fDecayChannels;       // pointer to list of decay channels

  void SetDecaySummary(TList* x);
  void SetThresholdCalc(Double_t x) {fThresholdCalc = x;}

 public:

  TTMParticle();
  TTMParticle(const TTMParticle &obj);
  ~TTMParticle();

  void SetPartName(TString x) {fPartName = x;} 
  void SetID(Int_t x);
  void SetStable() {fStable = true;}
  void SetUnstable() {fStable = false;}
  void SetDeg(Int_t x) {fDeg = x;}
  void SetStat(Int_t x) {fStat = x;}
  void SetB(Int_t x) {fB = x;}
  void SetS(Int_t x) {fS = x;}
  void SetQ(Int_t x) {fQ = x;}
  void SetCharm(Int_t x) {fCharm = x;}
  void SetBeauty(Int_t x) {fBeauty = x;}
  void SetTop(Int_t x) {fTop = x;}

  void SetMass(Double_t x) {fMass = x;}
  void SetWidth(Double_t x) {fWidth = x;}

  void SetThreshold(Double_t x) {fThreshold = x;}

  void SetRadius(Double_t x) {fRadius = x;}

  void SetThresholdFlag(Bool_t x) {fThresholdFlag = x;}

  void SetSContent(Double_t x) {fSContent = x;}
  void SetCContent(Double_t x) {fCContent = x;}
  void SetbContent(Double_t x) {fbContent = x;}
  void SetTContent(Double_t x) {fTContent = x;}

  void SetDecayChannels(const char *file, Bool_t ScaleBRatios = false);
  void SetDecayChannels(TList* x);

  void SetDecayChannelEfficiency(Int_t channel, Double_t eff);
  
  TString GetPartName() const {return fPartName;}

  Int_t GetID() const {return fID;}
  Bool_t GetStable() const {return fStable;}
  Int_t GetDeg() const {return fDeg;}
  Int_t GetStat() const {return fStat;}
  Int_t GetB() const {return fB;}
  Int_t GetS() const {return fS;}
  Int_t GetQ() const {return fQ;}
  Int_t GetCharm() const {return fCharm;}
  Int_t GetBeauty() const {return fBeauty;}
  Int_t GetTop() const {return fTop;}

  Double_t GetMass() const {return fMass;}
  Double_t GetWidth() const {return fWidth;}

  Double_t GetThreshold() const {return fThreshold;}

  Double_t GetRadius() const {return fRadius;}

  Double_t GetThresholdCalc() const {return fThresholdCalc;}

  Bool_t GetThresholdFlag() const {return fThresholdFlag;}

  Double_t GetSContent() const {return fSContent;}
  Double_t GetCContent() const {return fCContent;}
  Double_t GetbContent() const {return fbContent;}
  Double_t GetTContent() const {return fTContent;}

  TList* GetDecaySummary() const {return fDecaySummary;};
  TList* GetDecayChannels() const {return fDecayChannels;};

  void UpdateDecaySummary();

  TTMDecay* GetDecay(Int_t daughter_id);
  TTMDecayChannel* GetDecayChannel(Int_t channel);	 

  void List();

  TTMParticle& operator=(const TTMParticle& obj);

  ClassDef(TTMParticle,1) // Particle object 

};

#endif
