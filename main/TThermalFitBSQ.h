/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <fstream>
#include <iostream>

#ifndef Thermal_TTMThermalFit
#include <TThermalFit.h>
#endif

#ifndef Thermal_TTMParameterSetBSQ
#include <TTMParameterSetBSQ.h>
#endif

#ifndef Thermal_TTMThermalModelBSQ
#include <TThermalModelBSQ.h>
#endif

#ifndef Thermal_TTMThermalFitBSQ
#define Thermal_TTMThermalFitBSQ

using namespace std;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalFitBSQ							//
//            	                                 			//
// Fit class applying the grand canonical formalism 			// 
//									//
//////////////////////////////////////////////////////////////////////////

class TTMThermalFitBSQ:public TTMThermalFit {

 protected:
  
  TTMParameterSetBSQ *fParm;	
  Bool_t fQStats;
  Bool_t fWidth;
  Bool_t fExclVol;
  TTMThermalModelBSQ *fModel;	 	

  TTMThermalModelBSQ* GenerateThermalModel(TTMParticleSet *set);

 public:

  TTMThermalFitBSQ();
  TTMThermalFitBSQ(TTMParticleSet *set, TTMParameterSetBSQ *par, const char *file);
   
  ~TTMThermalFitBSQ();

  TTMParameterSetBSQ* GetParameterSet(){return fParm;}

  void SetParameterSet(TTMParameterSetBSQ *par){fParm = par;}
  void SetQStats(Bool_t x){fQStats = x;}
  void SetWidth(Bool_t x){fWidth = x;}
  void SetExclVol(Bool_t x){fExclVol = x;}

  Bool_t GetQStats(){return fQStats;}
  Bool_t GetWidth(){return fWidth;}
  Bool_t GetExclVol(){return fExclVol;}

  ClassDef(TTMThermalFitBSQ,1) // Grand-Canonical thermal fit class

};

#endif
