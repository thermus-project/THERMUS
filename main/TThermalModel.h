/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 26 April 2010 //

#include <TROOT.h>
#include <TObject.h>
#include <THashTable.h>
#include <TString.h>
#include <TList.h>
#include <fstream>
#include <iostream>

#ifndef Thermal_TTMParameterSet
#include <TTMParameterSet.h>
#endif

#ifndef Thermal_TTMParticleSet
#include <TTMParticleSet.h>
#endif

#ifndef Thermal_TTMDecay
#include <TTMDecay.h>
#endif

#ifndef Thermal_TTMDensObj
#include "TTMDensObj.h"
#endif

#ifndef Thermal_TTMThermalModel
#define Thermal_TTMThermalModel

using namespace std;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalModel							//
//            	                                 			//
// Base class for thermal model objects. Thermal densities are stored	// 
// as TTMDensObj objects in a hash table. All derived classes must 	//
// define functions to calculate the primordial energy, entropy and 	//
// particle densities and pressure.					//  	
//									//
//////////////////////////////////////////////////////////////////////////

class TTMThermalModel:public TObject {

 protected:

   TTMParticleSet *fPartSet;	// pointer to particle set
   Bool_t fWidth;		// True if width is to be taken into account 
   THashTable *fDensTable;	// pointer to density hash table
  
   TString fDescriptor;		// string describing model
 
   Double_t fEnergy;		// total energy density in model
   Double_t fEntropy;		// total entropy density in model
   Double_t fPressure;		// total pressure in model
   Double_t fWroblewski;	// Wroblewski factor

   Double_t fStrange;		// total strangeness density in model
   Double_t fBaryon;		// total baryon density in model
   Double_t fCharge;		// total charge density in model
   Double_t fCharm;		// total charm density in model
   Double_t fBeauty;		// total beauty density in model
   Double_t fDensity;		// total particle density in model

   Double_t fSplus;
   Double_t fSminus;

   Double_t fBplus;
   Double_t fBminus;

   Double_t fQplus;
   Double_t fQminus;

   Double_t fCplus;
   Double_t fCminus;

   Double_t fbplus;
   Double_t fbminus;

   void CalcWroblewski();

   virtual Int_t PrimPartDens() = 0;

 public:

   TTMThermalModel() 
   {
      fDensTable = new THashTable();
      fDensTable->SetOwner(kTRUE);
   }  
   
   ~TTMThermalModel();
   
   virtual TTMParameterSet* GetParameterSet() = 0;
   
   void GenerateDecayPartDens();
   void GenerateDecayPartDens(Int_t id);
   void ListDecayContributions(Int_t d_id);
   void ListDecayContribution(Int_t p_id,Int_t d_id);
 
   virtual Int_t GenerateParticleDens() = 0;
   virtual void GenerateEnergyDens() = 0;
   virtual void GenerateEntropyDens() = 0;
   virtual void GeneratePressure() = 0;
   virtual void ListInfo() = 0;
  
   void SetParticleSet(TTMParticleSet* set){fPartSet = set;}
   void SetWidth(Bool_t x){fWidth = x;}

   Bool_t GetWidth() const {return fWidth;}
   TTMParticleSet* GetParticleSet() const {return fPartSet;}
   THashTable* GetParticleTable() const {return fPartSet->GetParticleTable();}   
   TTMDensObj* GetDensities(Int_t ID);
   THashTable* GetDensityTable() const {return fDensTable;}

   TString GetDescriptor() const {return fDescriptor;}

   Double_t GetEnergy() const {return fEnergy;}
   Double_t GetEntropy() const {return fEntropy;}
   Double_t GetPressure() const {return fPressure;}
   Double_t GetDensity() const {return fDensity;}
   Double_t GetBaryon() const {return fBaryon;}
   Double_t GetCharge() const {return fCharge;}
   Double_t GetStrange() const {return fStrange;}
   Double_t GetCharm() const {return fCharm;}
   Double_t GetBeauty() const {return fBeauty;}
   Double_t GetWroblewski() const {return fWroblewski;}
   
   Double_t GetSplus() const {return fSplus;}
   Double_t GetSminus() const {return fSminus;}

   Double_t GetBplus() const {return fBplus;}
   Double_t GetBminus() const {return fBminus;}

   Double_t GetQplus() const {return fQplus;}
   Double_t GetQminus() const {return fQminus;}

   Double_t GetCplus() const {return fCplus;}
   Double_t GetCminus() const {return fCminus;}

   Double_t Getbplus() const {return fbplus;}
   Double_t Getbminus() const {return fbminus;}

   void ListStableParticles();
   void ListStableDensities();
   
   ClassDef(TTMThermalModel,1) // Base class for thermal model objects

};

#endif
