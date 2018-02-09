/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TH2F.h>

#ifndef Thermal_TTMThermalModel
#include <TThermalModel.h>
#endif

#ifndef Thermal_TTMParameterSetCanBSQ
#include <TTMParameterSetCanBSQ.h>
#endif

#ifndef Thermal_TTMThermalParticleCanBSQ
#include <TThermalParticleCanBSQ.h>
#endif

#ifndef Thermal_TTMThermalModelCanBSQ
#define Thermal_TTMThermalModelCanBSQ

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalModelCanBSQ						//
//            	                                 			//
// Canonical Thermal Model Class in the Boltzmann approximation		// 
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMThermalModelCanBSQ:public TTMThermalModel {

 private:

  TTMParameterSetCanBSQ *fParm;	// pointer to thermal parameters
  
  Double_t fMuB,fMuS,fMuQ;  		// associated chemical potentials

  Double_t flnZtot;	     // log(total canonical partition function)  

  Double_t fCorrpip;        // correction factor for "pi+ like" particles
  Double_t fCorrpim;        // correction factor for "pi- like" particles
  Double_t fCorrkm;         // correction factor for "K- like" particles
  Double_t fCorrkp;         // correction factor for "K+ like" particles
  Double_t fCorrk0;         // correction factor for "K0 like" particles
  Double_t fCorrak0;        // correction factor for "a-K0 like" particles
  Double_t fCorrproton;     // correction factor for "p like" particles
  Double_t fCorraproton;    // correction factor for "a-p like" particles
  Double_t fCorrneutron;    // correction factor for "n like" particles
  Double_t fCorraneutron;   // correction factor for "a-n like" particles
  Double_t fCorrlambda;     // correction factor for "Lambda like" particles
  Double_t fCorralambda;    // correction factor for "a-La like" particles
  Double_t fCorrsigmap;     // correction factor for "Sigma+ like" particles
  Double_t fCorrasigmap;    // correction factor for "a-Sigma+ like" particles
  Double_t fCorrsigmam;     // correction factor for "Sigma- like" particles
  Double_t fCorrasigmam;    // correction factor for "a-Sigma- like" particles
  Double_t fCorrdeltam;     // correction factor for "Delta- like" particles
  Double_t fCorradeltam;    // correction factor for "a-Delta- like" particles
  Double_t fCorrdeltapp;    // correction factor for "Delta++ like" particles
  Double_t fCorradeltapp;   // correction factor for "a-Delta++ like" particles
  Double_t fCorrksim;       // correction factor for "Ksi- like" particles
  Double_t fCorraksim;      // correction factor for "a-Ksi- like" particles
  Double_t fCorrksi0;       // correction factor for "Ksi0 like" particles
  Double_t fCorraksi0;      // correction factor for "a-Ksi0 like" particles
  Double_t fCorromega;      // correction factor for "Omega like" particles
  Double_t fCorraomega;     // correction factor for "a-Omega like" particles

  Double_t fTol;
  Double_t fRadius;
  Double_t fStepScale;
  Int_t fMin;
  Int_t fOrder;

  Int_t PrimPartDens();

  Double_t par[17];

  void UpdatePartitionFcnParameters();

 public:

  void PopulateZHistograms(TH2F *hZ[27]);

  TTMThermalModelCanBSQ(TTMParticleSet *particles, 
                        TTMParameterSetCanBSQ *parameters, Bool_t width = true);
  TTMThermalModelCanBSQ();

  Bool_t fListCorrFactorOutput;
  Bool_t fTF2IntCorrFactorCalc;

  void SetParameterSet(TTMParameterSetCanBSQ* parm){fParm = parm;}
  TTMParameterSetCanBSQ* GetParameterSet() {return fParm;} 
   
  Int_t GenerateParticleDens();
  void GenerateEnergyDens();
  void GenerateEntropyDens();
  void GeneratePressure();

  Bool_t GetListCorrFactorOutput() const {return fListCorrFactorOutput;}
  Bool_t GetTF2IntCorrFactorCalc() const {return fTF2IntCorrFactorCalc;}  

  Double_t GetCorrFactor(TTMParticle *part) const;

  Double_t GetlnZtot() const {return flnZtot;}
  Double_t GetMuB() const {return fMuB;}
  Double_t GetMuS() const {return fMuS;}
  Double_t GetMuQ() const {return fMuQ;}

  Double_t GetCorrpip() const {return fCorrpip;}        
  Double_t GetCorrpim() const {return fCorrpim;}          
  Double_t GetCorrkm() const {return fCorrkm;}           
  Double_t GetCorrkp() const {return fCorrkp;}           
  Double_t GetCorrk0() const {return fCorrk0;}           
  Double_t GetCorrak0() const {return fCorrak0;}          
  Double_t GetCorrproton() const {return fCorrproton;}       
  Double_t GetCorraproton() const {return fCorraproton;}      
  Double_t GetCorrneutron() const {return fCorrneutron;}      
  Double_t GetCorraneutron() const {return fCorraneutron;}     
  Double_t GetCorrlambda() const {return fCorrlambda;}       
  Double_t GetCorralambda() const {return fCorralambda;}      
  Double_t GetCorrsigmap() const {return fCorrsigmap;}       
  Double_t GetCorrasigmap() const {return fCorrasigmap;}      
  Double_t GetCorrsigmam() const {return fCorrsigmam;}       
  Double_t GetCorrasigmam() const {return fCorrasigmam;}      
  Double_t GetCorrdeltam() const {return fCorrdeltam;}       
  Double_t GetCorradeltam() const {return fCorradeltam;}      
  Double_t GetCorrdeltapp() const {return fCorrdeltapp;}      
  Double_t GetCorradeltapp() const {return fCorradeltapp;}     
  Double_t GetCorrksim() const {return fCorrksim;}         
  Double_t GetCorraksim() const {return fCorraksim;}        
  Double_t GetCorrksi0() const {return fCorrksi0;}         
  Double_t GetCorraksi0() const {return fCorraksi0;}        
  Double_t GetCorromega() const {return fCorromega;}        
  Double_t GetCorraomega() const {return fCorraomega;}    

  void SetListCorrFactorOutput(Bool_t x){fListCorrFactorOutput = x;}
  void SetTF2IntCorrFactorCalc(Bool_t x){fTF2IntCorrFactorCalc = x;}

  void SetTolerance(Double_t x){fTol = x;}
  void SetIntRadius(Double_t x){fRadius = x;}
  void SetStepScale(Double_t x){fStepScale = x;}
  void SetMinSteps(Int_t x){fMin = x;}

  void SetOrder(Int_t x){

    if( x==4 || x==8 || x==16 || x==20 || x==32){ 
      fOrder = x;
    }else {
      fOrder=8; 
      cout<<"Allowed values: 4, 8, 16, 20, 32. Setting to 8"<<endl;
    }

  }

  Double_t GetTolerance() const {return fTol;}
  Double_t GetIntRadius() const {return fRadius;}
  Double_t GetStepScale() const {return fStepScale;}
  Int_t GetMinSteps() const {return fMin;}
  Int_t GetOrder() const {return fOrder;}

  void ListInfo();

  ClassDef(TTMThermalModelCanBSQ,1) // Canonical thermal model class

    };

#endif
