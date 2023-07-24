/************************************************************************
 * Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 *                                                                      *
 * Author: The THERMUS Project (A Thermal Model Package for ROOT).      *
 * Contributors (UCT-IPHC) are mentioned in the code where appropriate. *
 *                                                                      *
 * Permission to use, copy, modify and distribute this software and its *
 * documentation strictly for non-commercial purposes is hereby granted *
 * without fee, provided that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice *
 * appear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is        *
 * provided "as is" without express or implied warranty.                *
 ************************************************************************/

//Authors: Natasha Sharma and Jean Cleymans 2 February 2021 (modification of the corresponding
//BQ file, TThermalModelBQ.cxx,  by Spencer Wheaton 14 July 2004)
//__________________________________________________________________________
// Baryon Canonical thermal model class.

#include <TThermalModelSQ.h>
#include <TThermalParticleBSQ.h>
#include <FncsConstrain.h>

ClassImp(TTMThermalModelSQ)

//__________________________________________________________________________
TTMThermalModelSQ::TTMThermalModelSQ(TTMParticleSet *particles,
                                       TTMParameterSetSQ *parameters, 
                                       Bool_t width)
{
  fDescriptor = "BCanonical";
  fPartSet = particles;
  fParm = parameters;
  fWidth = width;
  fDensTable->Delete();
  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = 0.;
  fNonBaryonQStats = true;
}

//__________________________________________________________________________
TTMThermalModelSQ::TTMThermalModelSQ()
{
  fDescriptor = "BCanonical";
  fPartSet = (TTMParticleSet *) 0;
  fParm = (TTMParameterSetSQ *) 0;
  fWidth = true;
  fDensTable->Delete();
  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = 0.;
  fNonBaryonQStats = true;
}

//__________________________________________________________________________
void TTMThermalModelSQ::Term(Double_t *x, Double_t *y, Int_t m, Int_t n, Double_t *t)
{
  // Private fcn which calculates normalised terms for the summations used 
  // to calculate the correction factors  
  //
  // 	t[0] : term for full partition function sum
  // 	t[1] : term for +1 sum
  // 	t[2] : term for +2 sum
  // 	t[3] : term for +3 sum
  // 	t[4] : term for -1 sum
  // 	t[5] : term for -2 sum
  // 	t[6] : term for -3 sum
  //

  Double_t B = fParm->GetB();
 
  // if b and/or c are zero, then all are zero //

  if(TMath::BesselI(TMath::Abs(n),x[1])==0. || TMath::BesselI(TMath::Abs(m),x[2])==0.){

    t[0] = 0.;
    t[1] = 0.;
    t[2] = 0.;
    t[3] = 0.;
    t[4] = 0.;
    t[5] = 0.;
    t[6] = 0.;

  }else{

    Double_t b = TMath::Log10(TMath::BesselI(TMath::Abs(n),x[1])/TMath::BesselI(0,x[1]));
    Double_t c = TMath::Log10(TMath::BesselI(TMath::Abs(m),x[2])/TMath::BesselI(0,x[2]));
    Double_t d = m*TMath::Log10(y[2]/pow(y[0],3));
    Double_t e = n*TMath::Log10(y[1]/pow(y[0],2));

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B),x[0]) == 0.){
      t[0] = 0;
    }else{
      Double_t a0 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B),x[0])/TMath::BesselI(0,x[0]));
      t[0] = pow(10.,a0 + b + c + d + e + B * TMath::Log10(y[0]));
    }

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B+1),x[0]) == 0.){
      t[1] = 0.;
    }else{
      Double_t ap1 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B+1),x[0])/TMath::BesselI(0,x[0]));
      t[1] = pow(10.,ap1 + b + c + d + e + (B - 1.)*TMath::Log10(y[0]));
    }

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B+2),x[0]) == 0.){
      t[2] = 0.;
    }else{
      Double_t ap2 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B+2),x[0])/TMath::BesselI(0,x[0]));
      t[2] = pow(10.,ap2 + b + c + d + e + (B - 2.)*TMath::Log10(y[0]));
    }

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B+3),x[0]) == 0.){
      t[3] = 0.;
    }else{
      Double_t ap3 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B+3),x[0])/TMath::BesselI(0,x[0]));
      t[3] = pow(10.,ap3 + b + c + d + e + (B - 3.)*TMath::Log10(y[0]));
    }

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B-1),x[0]) == 0.){
      t[4] = 0.;
    }else{
      Double_t am1 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B-1),x[0])/TMath::BesselI(0,x[0]));
      t[4] = pow(10.,am1 + b + c + d + e + (B + 1.)*TMath::Log10(y[0]));
    }

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B-2),x[0]) == 0.){
      t[5] = 0.;
    }else{
      Double_t am2 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B-2),x[0])/TMath::BesselI(0,x[0]));
      t[5] = pow(10.,am2 + b + c + d + e + (B + 2.)*TMath::Log10(y[0]));
    }

    if(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B-3),x[0]) == 0.){
      t[6] = 0.;
    }else{
      Double_t am3 = TMath::Log10(TMath::BesselI((Int_t)TMath::Abs(3*m+2*n-B-3),x[0])/TMath::BesselI(0,x[0])); 
      t[6] = pow(10.,am3 + b + c + d + e + (B + 3.)*TMath::Log10(y[0]));
    }

  }
}

//__________________________________________________________________________
Int_t TTMThermalModelSQ::PrimPartDens()
{
  // Calculates the primordial particle densities and populates the density
  // hash table. The parameters are not constrained first!. This is the 
  // function used by GenerateParticleDens() and calls Term(). 
  // Returns 1 if there are problems 0 else.
  //

  fDensTable->Delete();
  TIter next(GetParticleTable());
  TTMParticle *part;

  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = 0.;

  Int_t check = 0;

  if(fNonBaryonQStats){

    part = (TTMParticle *) next();
    
    do{ 

      if(part->GetB() == 0.){

//       TTMParameterSetBSQ pGC(fParm->GetT(),fParm->GetMuS(),0.,fParm->GetMuQ(),fParm->GetGammas()); 
        TTMParameterSetBSQ pGC(fParm->GetT(),0.,fParm->GetMuS(),fParm->GetMuQ(),fParm->GetGammas()); 

        TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part,&pGC);

        if(!ptr->ParametersAllowed()){
          check = 1;
        }
          //delete ptr; // possible memory leak B.H. 09/02/2016
      }
      
      part = (TTMParticle *) next();
      
    }while(check == 0 && part);

    next.Reset();

  }
    
  if(check){

    cout<<"WARNING: Chemical Potentials not allowed!"<<endl;

    return 1;

  }else{

    Double_t Z = 0.;

    Double_t k0 = 0., kp1 = 0., kp2 = 0., kp3 = 0., km1 = 0., km2 = 0., km3 = 0.;

    Double_t r = fParm->GetCanRadius();
    Double_t volume = 4. * TMath::Pi() / 3. * r * r * r;

    while((part = (TTMParticle *) next())){

      Double_t B = part->GetB();

      if(B == 0. && fNonBaryonQStats){

        Double_t GCPressure;

        TTMParameterSetBSQ pGC(fParm->GetT(),0.,fParm->GetMuS(),fParm->GetMuQ(),fParm->GetGammas()); 
 //       TTMParameterSetBSQ pGC(fParm->GetT(),fParm->GetMuS(),0.,fParm->GetMuQ(),fParm->GetGammas()); 

        TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part,&pGC);

        if(part->GetWidth() == 0. || !fWidth){
          GCPressure = ptr->PressureQStatNoWidth();
        }else{
          GCPressure = ptr->PressureQStatWidth();
        }	      

        delete ptr;

        Z = GCPressure * volume / fParm->GetT();

      }else{

        Double_t GCPartDens;

        TTMThermalParticleSQ *ptr = new TTMThermalParticleSQ(part, fParm, 1);

        if(part->GetWidth() == 0 || !fWidth){
          GCPartDens = ptr->DensityBoltzmannNoWidth();
        }else{
          GCPartDens = ptr->DensityBoltzmannWidth();
        }	      

        delete ptr;

        Z = GCPartDens * volume;    // only true in Boltzmann case since PV/T = nV

      }

      if (B == 0.){
        k0 += Z;     
      } else if (B == 1.) {
        kp1 += Z;
      } else if (B == 2.) {
        kp2 += Z;
      } else if (B == 3.) {
        kp3 += Z;
      } else if (B == -1.) {
        km1 += Z;
      } else if (B == -2.) {
        km2 += Z;
      } else if (B == -3.) {
        km3 += Z;
      }
    }

    Double_t x[3];

    x[0] = 2. * sqrt(kp1 * km1);
    x[1] = 2. * sqrt(kp2 * km2);
    x[2] = 2. * sqrt(kp3 * km3);

    Double_t y[3];

    y[0] = sqrt(kp1 / km1);
    y[1] = sqrt(kp2 / km2);
    y[2] = sqrt(kp3 / km3);

    Double_t xmax = 709.;		// maximum argument for Bessel fncs in ROOT

    if(x[0] > xmax || x[1] > xmax || x[2] > xmax ){

      cout<<"x's too large"<<endl;
      return 1;

    }else if(y[0] == 0 || y[1] == 0 || y[2] == 0){

      cout<<"y's zero"<<endl;
      return 1;

    }else{

      Double_t TBPartition3;
      Double_t Tfactorkm1;
      Double_t Tfactorkp1;
      Double_t Tfactorkm2;
      Double_t Tfactorkp2;
      Double_t Tfactorkm3;
      Double_t Tfactorkp3;

      Double_t term[7];

      // Include (0,0) term //

      TBPartition3 = 0.;
      Tfactorkm1 = 0.;
      Tfactorkp1 = 0.;
      Tfactorkm2 = 0.;
      Tfactorkp2 = 0.;
      Tfactorkm3 = 0.;
      Tfactorkp3 = 0.;

      Term(x,y,0,0,term);

      TBPartition3 = term[0];  
      Tfactorkm1 = term[4];
      Tfactorkp1 = term[1];
      Tfactorkm2 = term[5];
      Tfactorkp2 = term[2];
      Tfactorkm3 = term[6];
      Tfactorkp3 = term[3];

      Double_t dZ;
      Double_t dp1;
      Double_t dp2;
      Double_t dp3;
      Double_t dm1;
      Double_t dm2;
      Double_t dm3;

      Int_t m,n;
      Int_t mmax; 
      Double_t tol = 1e-40;

      // First Quandrant +m and +n axes included //

      n = 0;
      m = 1;

      do{    

        do{

          Term(x,y,m,n,term);
 
          TBPartition3 += term[0];  
          Tfactorkm1 += term[4];
          Tfactorkp1 += term[1];
          Tfactorkm2 += term[5];
          Tfactorkp2 += term[2];
          Tfactorkm3 += term[6];
          Tfactorkp3 += term[3];

          dZ = term[0]/TBPartition3;
          dp1 = term[1]/Tfactorkp1; 
          dp2 = term[2]/Tfactorkp2;
          dp3 = term[3]/Tfactorkp3;
          dm1 = term[4]/Tfactorkm1;
          dm2 = term[5]/Tfactorkm2;
          dm3 = term[6]/Tfactorkm3;
     
          m++;

        }while(dZ>tol);

        n++;
        mmax = m - 1;
        m=0;

      }while(mmax > 0);

      // Second Quandrant -m axis included

      n = 0;
      m = -1;

      do{    

        do{

          Term(x,y,m,n,term);
 
          TBPartition3 += term[0];  
          Tfactorkm1 += term[4];
          Tfactorkp1 += term[1];
          Tfactorkm2 += term[5];
          Tfactorkp2 += term[2];
          Tfactorkm3 += term[6];
          Tfactorkp3 += term[3];

          dZ = term[0]/TBPartition3;
          dp1 = term[1]/Tfactorkp1; 
          dp2 = term[2]/Tfactorkp2;
          dp3 = term[3]/Tfactorkp3;
          dm1 = term[4]/Tfactorkm1;
          dm2 = term[5]/Tfactorkm2;
          dm3 = term[6]/Tfactorkm3;
     
          m--;

        }while(dZ>tol);

        n++;
        mmax = m + 1;
        m = -1;

      }while(mmax < -1);

      // Third Quandrant -n axis included

      n = -1;
      m = 0;

      do{

        do{

          Term(x,y,m,n,term);
 
          TBPartition3 += term[0];  
          Tfactorkm1 += term[4];
          Tfactorkp1 += term[1];
          Tfactorkm2 += term[5];
          Tfactorkp2 += term[2];
          Tfactorkm3 += term[6];
          Tfactorkp3 += term[3];

          dZ = term[0]/TBPartition3;
          dp1 = term[1]/Tfactorkp1; 
          dp2 = term[2]/Tfactorkp2;
          dp3 = term[3]/Tfactorkp3;
          dm1 = term[4]/Tfactorkm1;
          dm2 = term[5]/Tfactorkm2;
          dm3 = term[6]/Tfactorkm3;
     
          m--;

        }while(dZ>tol);

        n--;
        mmax = m + 1;
        m = 0;

      }while(mmax < 0);

      // Fourth Quandrant

      n = -1;
      m = 1;

      do{    

        do{

          Term(x,y,m,n,term);
 
          TBPartition3 += term[0];  
          Tfactorkm1 += term[4];
          Tfactorkp1 += term[1];
          Tfactorkm2 += term[5];
          Tfactorkp2 += term[2];
          Tfactorkm3 += term[6];
          Tfactorkp3 += term[3];

          dZ = term[0]/TBPartition3;
          dp1 = term[1]/Tfactorkp1; 
          dp2 = term[2]/Tfactorkp2;
          dp3 = term[3]/Tfactorkp3;
          dm1 = term[4]/Tfactorkm1;
          dm2 = term[5]/Tfactorkm2;
          dm3 = term[6]/Tfactorkm3;
     
          m++;

        }while(dZ>tol);

        n--;
        mmax = m - 1;
        m = 1;

      }while(mmax > 1);

      fCorrP1 = Tfactorkp1 / TBPartition3;
      fCorrP2 = Tfactorkp2 / TBPartition3;
      fCorrP3 = Tfactorkp3 / TBPartition3;
      fCorrM1 = Tfactorkm1 / TBPartition3;
      fCorrM2 = Tfactorkm2 / TBPartition3;
      fCorrM3 = Tfactorkm3 / TBPartition3;

      flnZtot = TMath::Log(TBPartition3) + TMath::Log(TMath::BesselI(0,x[0])) + TMath::Log(TMath::BesselI(0,x[1])) + TMath::Log(TMath::BesselI(0,x[2])) + k0;

      flnZ0 = k0;

      fExactMuB = fParm->GetT() * log(fCorrP1);
   
      next.Reset();

      while ((part = (TTMParticle *) next())) {
        Double_t B = part->GetB();
        Double_t CorrFactor;

        if (B == 1) {
          CorrFactor = fCorrP1;
        } else if (B == 2) {
          CorrFactor = fCorrP2;
        } else if (B == 3) {
          CorrFactor = fCorrP3;
        } else if (B == -1) {
          CorrFactor = fCorrM1;
        } else if (B == -2) {
          CorrFactor = fCorrM2;
        } else if (B == -3) {
          CorrFactor = fCorrM3;
        } else {
          CorrFactor = 1.;
        }

        Double_t PartDens;

        if(B == 0. && fNonBaryonQStats){

          TTMParameterSetBSQ pGC(fParm->GetT(),0.,fParm->GetMuS(),fParm->GetMuQ(),fParm->GetGammas()); 

          TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part, &pGC);

          if(part->GetWidth() == 0 || !fWidth){
            PartDens = ptr->DensityQStatNoWidth();
          }else{
            PartDens = ptr->DensityQStatWidth();
          }	      

          delete ptr;

        }else{

          TTMThermalParticleSQ *ptr = new TTMThermalParticleSQ(part, fParm, CorrFactor);
      
          if(part->GetWidth() == 0 || !fWidth){
            PartDens = ptr->DensityBoltzmannNoWidth();
          }else{
            PartDens = ptr->DensityBoltzmannWidth();
          }	      	      
     
          delete ptr;

        }

        TTMDensObj *dens = new TTMDensObj(part->GetID());
        dens->SetPrimDensity(PartDens);
        fDensTable->Add(dens);

        fDensity += PartDens;

        if(part->GetS() > 0.){
          fSplus += PartDens * (part->GetS());
        }else if(part->GetB() < 0.){
          fSminus += PartDens * (part->GetS()); 
        }

        if(part->GetB() > 0.){
          fBplus += PartDens * (part->GetB()); 
        }else if(part->GetB() < 0.){
          fBminus += PartDens * (part->GetB()); 
        }

        if(part->GetQ() > 0.){
          fQplus += PartDens * (part->GetQ());
        }else if(part->GetQ() < 0.){
          fQminus += PartDens * (part->GetQ());
        }

        if(part->GetCharm() > 0.){
          fCplus += PartDens * (part->GetCharm());
        }else if(part->GetCharm() < 0.){
          fCminus += PartDens * (part->GetCharm()); 
        }

        if(part->GetBeauty() > 0.){
          fbplus += PartDens * (part->GetBeauty());
        }else if(part->GetBeauty() < 0.){
          fbminus += PartDens * (part->GetBeauty()); 
        }


      }

      fStrange = fSplus + fSminus;
      fBaryon = fBplus + fBminus;
      fCharge = fQplus + fQminus;
      fCharm = fCplus + fCminus;
      fBeauty = fbplus + fbminus;

      return 0;
    }
  }
}

//__________________________________________________________________________
Int_t TTMThermalModelSQ::GenerateParticleDens()
{
  // Calculates the Primordial particle densities and populates the density 
  // hash table after first constraining muQ if required. The 
  // Wroblewski factor and the decay contributions are also calculated 
  // provided the decays have been entered into the particle set through
  // TTMParticleSet::InputDecays(). 
  //

  Int_t check = 1;

  if (fParm->GetCorrRConstrain()) {
    fParm->SetCanRadius(fParm->GetRadius());
    fParm->GetParameter(4)->SetStatus("(Set to Fireball Radius)");
  }
  check = PrimPartDens();
  
  if(!check){
    CalcWroblewski();
    GenerateDecayPartDens();
  }else{
    cout<<"Problems"<<endl;
  }

  return check;
}

//__________________________________________________________________________
void TTMThermalModelSQ::GenerateEnergyDens()
{  
  // Iterates through the density hash table calculating the primordial
  // energy density of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to impose constraints and calculate the canonical correction
  // factors.
  //  
 	
  TIter next(fDensTable);
  TTMDensObj *dens;

  fEnergy = 0.;

  while ((dens = (TTMDensObj *) next())) {
    TTMParticle *part = fPartSet->GetParticle(dens->GetID());
    Double_t CorrFactor;

    if (part->GetB() == 1) {
      CorrFactor = fCorrP1;
    } else if (part->GetB() == 2) {
      CorrFactor = fCorrP2;
    } else if (part->GetB() == 3) {
      CorrFactor = fCorrP3;
    } else if (part->GetB() == -1) {
      CorrFactor = fCorrM1;
    } else if (part->GetB() == -2) {
      CorrFactor = fCorrM2;
    } else if (part->GetB() == -3) {
      CorrFactor = fCorrM3;
    } else {
      CorrFactor = 1.;
    }

    Double_t PartEnergy;

    if(part->GetB() == 0. && fNonBaryonQStats){

      TTMParameterSetBSQ pGC(fParm->GetT(),0.,fParm->GetMuS(),fParm->GetMuQ(),fParm->GetGammas()); 

      TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part, &pGC);

      if(part->GetWidth() == 0 || !fWidth){
        PartEnergy = ptr->EnergyQStatNoWidth();
      }else{
        PartEnergy = ptr->EnergyQStatWidth();
      }	      

      delete ptr;

    }else{

      TTMThermalParticleSQ *ptr = new TTMThermalParticleSQ(part, fParm, CorrFactor);
    
      if(part->GetWidth() == 0 || !fWidth){
        PartEnergy = ptr->EnergyBoltzmannNoWidth();
      }else{
        PartEnergy = ptr->EnergyBoltzmannWidth();
      }	    
  	
      delete ptr;

    }
      
    dens->SetPrimEnergy(PartEnergy);
    fEnergy += PartEnergy;
  }
}

//__________________________________________________________________________
void TTMThermalModelSQ::GenerateEntropyDens()
{
  // Iterates through the density hash table calculating the primordial
  // entropy density of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to impose constraints and calculate the canonical correction
  // factors. Remember that the total entropy does not split into the sum
  // of particle entropies.
  //  
   
  TIter next(fDensTable);
  TTMDensObj *dens;

  fEntropy = 0.;

  while ((dens = (TTMDensObj *) next())) {
    TTMParticle *part = fPartSet->GetParticle(dens->GetID());
    Double_t CorrFactor;

    if (part->GetB() == 1) {
      CorrFactor = fCorrP1;
    } else if (part->GetB() == 2) {
      CorrFactor = fCorrP2;
    } else if (part->GetB() == 3) {
      CorrFactor = fCorrP3;
    } else if (part->GetB() == -1) {
      CorrFactor = fCorrM1;
    } else if (part->GetB() == -2) {
      CorrFactor = fCorrM2;
    } else if (part->GetB() == -3) {
      CorrFactor = fCorrM3;
    } else {
      CorrFactor = 1.;
    }

    Double_t PartEntropy;

    if(part->GetB() == 0.){

      TTMParameterSetBSQ pGC(fParm->GetT(),0.,fParm->GetMuS(),fParm->GetMuQ(),fParm->GetGammas());
      
      TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part, &pGC);

      if(fNonBaryonQStats){

        if(part->GetWidth() == 0 || !fWidth){
          PartEntropy = ptr->EntropyQStatNoWidth();
        }else{
          PartEntropy = ptr->EntropyQStatWidth();
        }	      

      }else{

        if(part->GetWidth() == 0 || !fWidth){
          PartEntropy = ptr->EntropyBoltzmannNoWidth();
        }else{
          PartEntropy = ptr->EntropyBoltzmannWidth();
        }	      

      }

      delete ptr;

    }else{

      TTMThermalParticleSQ *ptr = new TTMThermalParticleSQ(part, fParm, CorrFactor);
    
      Double_t PartEnergy;

      if(part->GetWidth() == 0 || !fWidth){
        PartEnergy = ptr->EnergyBoltzmannNoWidth();
      }else{
        PartEnergy = ptr->EnergyBoltzmannWidth();
      }
	      	      
      delete ptr;
   
      Double_t PartDens = dens->GetPrimDensity();

      Double_t S = part->GetS();
      Double_t Q = part->GetQ();
      Double_t muS = fParm->GetMuS();
      Double_t muQ = fParm->GetMuQ();

      Double_t mu = S*muS + Q*muQ;

      PartEntropy = (PartEnergy - mu * PartDens) / fParm->GetT() ;

    }
 
    dens->SetPrimEntropy(PartEntropy);
    fEntropy += PartEntropy;

  }

  Double_t r = fParm->GetCanRadius();
  Double_t volume = 4. * TMath::Pi() / 3. * r * r * r;

  fEntropy += (flnZtot - flnZ0)/volume;

}

//__________________________________________________________________________
void TTMThermalModelSQ::GeneratePressure()
{
  // Iterates through the density hash table calculating the primordial
  // pressure of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to impose constraints and calculate the canonical correction
  // factors.
  //  
   
  TIter next(fDensTable);
  TTMDensObj *dens;

  fPressure = 0.;

  while ((dens = (TTMDensObj *) next())) {
    TTMParticle *part = fPartSet->GetParticle(dens->GetID());
    Double_t CorrFactor;

    if (part->GetB() == 1) {
      CorrFactor = fCorrP1;
    } else if (part->GetB() == 2) {
      CorrFactor = fCorrP2;
    } else if (part->GetB() == 3) {
      CorrFactor = fCorrP3;
    } else if (part->GetB() == -1) {
      CorrFactor = fCorrM1;
    } else if (part->GetB() == -2) {
      CorrFactor = fCorrM2;
    } else if (part->GetB() == -3) {
      CorrFactor = fCorrM3;
    } else {
      CorrFactor = 1.;
    }

    Double_t PartPressure;

    if(part->GetB() == 0. && fNonBaryonQStats){

      TTMParameterSetBSQ pGC(fParm->GetT(),0.,fParm->GetMuS(),fParm->GetMuQ(),fParm->GetGammas()); 

      TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part, &pGC);

      if(part->GetWidth() == 0 || !fWidth){
        PartPressure = ptr->PressureQStatNoWidth();
      }else{
        PartPressure = ptr->PressureQStatWidth();
      }	      

      delete ptr;

    }else{

      TTMThermalParticleSQ *ptr = new TTMThermalParticleSQ(part, fParm, CorrFactor);

      if(part->GetWidth() == 0 || !fWidth){
        PartPressure = ptr->PressureBoltzmannNoWidth();
      }else{
        PartPressure = ptr->PressureBoltzmannWidth();
      }	      	      

      delete ptr;
  
    }

    dens->SetPrimPressure(PartPressure);
    fPressure += PartPressure;
  }
}

//__________________________________________________________________________
void TTMThermalModelSQ::ListInfo()
{
  // List model information
  //
 	
  cout <<"  ************************************************************"
       <<"*****************"
       << endl;
  cout <<"  ***************************** Thermal Model Info ***********"
       <<"*****************"
       << endl << endl;

  cout << "\t Particle set: " << endl << "\t\t" << fPartSet->GetFilename()
       << endl << endl;

  if(fNonBaryonQStats){
    cout << "\t Quantum statistics for B-0 hadrons" << endl;
  }else{
    cout << "\t Boltzmann statistics for B-0 hadrons" << endl;
  }   
  cout << "\t Boltzmann statistics for baryons" << endl;
   
  if(fWidth){
    cout << "\t Resonance width included " << endl;
  }else{
    cout << "\t Resonance width excluded " << endl;
  }
   
  cout << endl << endl;

  fParm->List();
  cout <<"  ***************************** Thermal Quantities ***********"
       <<"****************** "
       << endl << endl;
  cout << "\t B required in canonical volume: ";
  cout.width(10);
  cout << fParm->GetB();
  cout << endl << endl;
  Double_t r = fParm->GetCanRadius();
  Double_t volume = 4. * TMath::Pi() / 3. * r * r * r;
  cout << "\t B in canonical volume (model) = ";
  cout.width(10);
  cout << fBaryon*volume;
  cout << endl << endl;
  cout << "\t B/2Q     = ";
  cout.width(10);
  cout << fBaryon / 2. / fCharge;
  cout << "\t";
  cout << endl << endl;
  cout << "\t lambda_s = ";
  cout.width(8);
  cout << fWroblewski << endl;
  cout << endl << endl;
  cout << "\t Primordial Densities: " << endl;
  if (fDensity != 0.) {
    cout << "\t\t\t\t n = " << fDensity << endl;
  }
  if (fEnergy != 0.) {
    cout << "\t\t\t\t e = " << fEnergy << endl;
  }
  if (fEntropy != 0.) {
    cout << "\t\t\t\t s = " << fEntropy << endl << endl;
  }
  ListStableParticles();
  cout << endl << endl;
  cout <<"  ***********************************************************"
       <<"*******************"
       << endl;
  cout <<"  ***********************************************************"
       <<"*******************"
       << endl << endl;
}
