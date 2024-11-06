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

// Author: Spencer Wheaton 26 April 2010 //
//__________________________________________________________________________
// Grand Canonical thermal model class.
//   

#include <TThermalModelBSQ.h>
#include <FncsConstrain.h>

ClassImp(TTMThermalModelBSQ)

//__________________________________________________________________________
TTMThermalModelBSQ::TTMThermalModelBSQ()
{
  fDescriptor = "GCanonical";
  fPartSet = (TTMParticleSet*) 0;
  fParm = (TTMParameterSetBSQ*) 0;
  fQStats = true;
  fWidth = true;
  fExclVolCorrection = false;
  fDensTable->Delete();
  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = 0.;
  fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = fPressure = 0.;
  fExclVolPressure = 0.;
  fExclVolDenominator = 1.;
}

//__________________________________________________________________________
TTMThermalModelBSQ::TTMThermalModelBSQ(TTMParticleSet *particles,
                                       TTMParameterSetBSQ *parameters, 
				       Bool_t qstats, Bool_t width)
{
  fDescriptor = "GCanonical";
  fPartSet = particles;
  fParm = parameters;
  fQStats = qstats;
  fWidth = width;
  fExclVolCorrection = false;
  fDensTable->Delete();
  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = 0.;
  fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = fPressure = 0.;
  fExclVolPressure = 0.;
  fExclVolDenominator = 1.;
}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::PrimPartDens()
{
  // Calculates the primordial particle densities and populates the density
  // hash table with these values. The parameters are not constrained first!
  // This is the function used by GenerateParticleDens().
  //

  fDensTable->Delete();
 
  THashTable *part_table = GetParticleTable();
  TIter next(part_table);
  TTMParticle *part;

  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = fPressure = 0.;

  fExclVolPressure = 0.;
  fExclVolDenominator = 1.;

  Int_t check = 0;

  if(fQStats){

    part = (TTMParticle *) next();
    
    do{ 

      TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part, fParm);   

      if(!ptr->ParametersAllowed()){
        check = 1;
      }

      part = (TTMParticle *) next();

      delete ptr;
      
    }while(check == 0 && part);

    next.Reset();

  }
    
  if(check){

    cout<<"WARNING: Chemical Potentials not allowed!"<<endl;

    return 1;

  }else{

    if(!fExclVolCorrection){

      while ((part = (TTMParticle *) next())) {

        TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(part, fParm);
      
        Double_t PartDens;
        Double_t Width = part->GetWidth();
        Int_t Stat = part->GetStat();
      
        if(Stat == 0 || !fQStats){		       // Boltzmann stats
          if(Width == 0 || !fWidth){		       // no width
            PartDens = ptr->DensityBoltzmannNoWidth();
          }else{					       // width
            PartDens = ptr->DensityBoltzmannWidth();
          }
        }else{					       // Quantum stats	      
          if(Width == 0 || !fWidth){		       // no width
            PartDens = ptr->DensityQStatNoWidth();
          }else{				    	       // width
            PartDens = ptr->DensityQStatWidth();
          }
        }
	      
        TTMDensObj *dens = new TTMDensObj(part->GetID());
        dens->SetPrimDensity(PartDens);

        fDensTable->Add(dens);
      
        fDensity += PartDens;

        if(part->GetS() > 0.){
          fSplus += PartDens * (part->GetS());
        }else if(part->GetS() < 0.){
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

        delete ptr;
      }

    }else{

      CalcExclVolPressure();

      Double_t Correction = 1.;

      TTMParticleSet CopySet = *fPartSet;

      TTMParameterSetBSQ CopyParm = *fParm;

      TIter CopyNext(CopySet.GetParticleTable());
      //TTMParticle *part;

      while ((part = (TTMParticle *) CopyNext())) {

        TTMThermalParticleBSQ *ptr = ShiftedParticle(part,&CopyParm,fExclVolPressure);
      
        Double_t PartDens;
        Double_t Width = part->GetWidth();
        Int_t Stat = part->GetStat();
      
        if(Stat == 0 || !fQStats){		       // Boltzmann stats
          if(Width == 0 || !fWidth){		       // no width
            PartDens = ptr->DensityBoltzmannNoWidth();
          }else{					       // width
            PartDens = ptr->DensityBoltzmannWidth();
          }
        }else{					       // Quantum stats	      
          if(Width == 0 || !fWidth){		       // no width
            PartDens = ptr->DensityQStatNoWidth();
          }else{				    	       // width
            PartDens = ptr->DensityQStatWidth();
          }
        }
	      
        TTMDensObj *dens = new TTMDensObj(part->GetID());
        dens->SetPrimDensity(PartDens);	// ideal gas with shifted chemical potential
        
        fDensTable->Add(dens);

        Double_t radius = part->GetRadius();
        Double_t volume = 4. * 4./3. * TMath::Pi() * pow(radius,3);

        Correction += volume*PartDens;

        delete ptr;

      }

      fExclVolDenominator = Correction;

      TIter nextDens(fDensTable);
      TTMDensObj *dens;

      while ((dens = (TTMDensObj *) nextDens())) {

        //TTMParticle *part = fPartSet->GetParticle(dens->GetID());
        part = fPartSet->GetParticle(dens->GetID());
        dens->SetPrimDensity(dens->GetPrimDensity()/fExclVolDenominator);

        Double_t PartDens = dens->GetPrimDensity();

        fDensity += PartDens;

        if(part->GetS() > 0.){
          fSplus += PartDens * (part->GetS());
        }else if(part->GetS() < 0.){
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

    }
    
    fStrange = fSplus + fSminus;
    fBaryon = fBplus + fBminus;
    fCharge = fQplus + fQminus;
    fCharm = fCplus + fCminus;
    fBeauty = fbplus + fbminus;

    return 0;
  }
}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::GenerateParticleDens()
{
  // Calculates the Primordial particle densities and populates the density 
  // hash table after first constraining muS and/or muQ and/or muC and/or 
  // mub if required. The Wroblewski factor and the decay contributions are 
  // also calculated provided the decays have been entered into the particle 
  // set through TTMParticleSet::InputDecays(). Note: the Wroblewski factor
  // is only correct in the absence of charm and beauty... 
  //

  Int_t check = 1;
	
  if (fParm->GetMuSConstrain() && fParm->GetMuQConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()) {
    check = BSQConstrainSQ(this);
    if (check) {
      cout << "Not Constrained" << endl;
    }
  } else if (fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()) {
    check = BSQConstrainS(this);     
    if (check) {
      cout << "Not Constrained" << endl;
    }
  } else if (fParm->GetMuQConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()) {
    check = BSQConstrainQ(this);
    if (check) {
      cout << "Not Constrained" << endl;
    }
  } else if (fParm->GetMuSConstrain() && fParm->GetMuQConstrain() && fParm->GetMuCConstrain() && !fParm->GetMubConstrain()) {
    check = BSQConstrainSQC(this);
    if (check) {
      cout << "Not Constrained" << endl;
    }
  } else if (fParm->GetMuSConstrain() && fParm->GetMuQConstrain() && fParm->GetMuCConstrain() && fParm->GetMubConstrain()) {
    check = BSQConstrainSQCb(this);
    if (check) {
      cout << "Not Constrained" << endl;
    }
  } else if(!fParm->GetMuSConstrain() && !fParm->GetMuQConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){
    check = PrimPartDens();
  } else {
    check = 1;
    cout << "Constraint not yet coded" << endl;
  }
   
  if(!check){
    CalcWroblewski();
    GenerateDecayPartDens();
  }else{
    cout<<"Problems"<<endl;
  }

  return check;
}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainEoverN(Double_t EoverN)
{
  // Constrains muB by the given E/N in the case where either both muS and 
  // muQ or just muS are also to be constrained while muC and mub are not to be 
  // constained. The other possibilities 
  // are yet to be coded. 1 is returned if there are problems 0 otherwise. 
  //


  Int_t check;

  if(fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){

    if(fParm->GetMuQConstrain()){

      check = BSQConstrainSQEN(this,EoverN);
      if(check){
        cout<<"Not Constrained"<<endl;
      }

    }else{

      check = BSQConstrainSEN(this,EoverN);
      if(check){
        cout<<"Not Constrained"<<endl;
      }
    }

  }else{

    check = 1;
    cout<<"Such constraint not yet coded"<<endl;
  }

  return check;

}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainBSQ(Double_t B, Double_t S, Double_t Q)
{
  // Code to find muB, muS, and muQ given B, S, Q and V (muC and mub not
  // constrained)
  //
  //

  Int_t check;

  if(!fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){
     check = BSQConstrainBSQ(this,B,S,Q);
  }else{
     check = 1;
     cout<<"Such constraint not yet coded"<<endl;
  }

  if(check){
     cout<<"Not Constrained"<<endl;
  }

  return check;

}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainQ(Double_t Q)
{
  // Code to find muQ given Q and V (muS, muC and mub not
  // constrained)
  //
  //

  Int_t check;

  if(!fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){
     check = BSQConstrainQQ(this,Q);
  }else{
     check = 1;
     cout<<"Such constraint not yet coded"<<endl;
  }

  if(check){
     cout<<"Not Constrained"<<endl;
  }

  return check;

}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainSoverT3(Double_t SoverT3)
{
  // Constrains muB by the given S/T^3 in the case where either both muS and 
  // muQ or just muS are also to be constrained while muC and mub are not. The 
  // other possibilities are yet to be coded. 1 is returned if there are 
  // problems 0 otherwise. 
  //

  Int_t check;

  if(fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){

    if(fParm->GetMuQConstrain()){

      check = BSQConstrainSQST3(this,SoverT3);
      if(check){
        cout<<"Not Constrained"<<endl;
      }

    }else{

      check = BSQConstrainSST3(this,SoverT3);
      if(check){
        cout<<"Not Constrained"<<endl;
      }
    }

  }else{

    check = 1;
    cout<<"Such constraint not yet coded"<<endl;
  }

  return check;

}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainPercolation()
{
  // Constrains muB by the percolation model. 1 is returned if there are problems 
  // 0 otherwise. 
  //

  Int_t check;

  if(fParm->GetMuQConstrain() && fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){

    check = BSQConstrainSQPercolation(this);
    if(check){
      cout<<"Not Constrained"<<endl;
    }
  }else if(fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){
   check = BSQConstrainSPercolation(this);
    if(check){
      cout<<"Not Constrained"<<endl;
    }
  }else{
      check = 1;
      cout<<"Such Constraint Not Yet Coded"<<endl;
  }

  return check;

}

//__________________________________________________________________________
void TTMThermalModelBSQ::GenerateEnergyDens()
{
  // Iterates through the density hash table calculating the primordial
  // energy density of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to impose constraints and ensure consistency between particle
  // and energy density.
  //   	
   
  TIter next(fDensTable);
  TTMDensObj *dens;
  fEnergy = 0.;
  Double_t PartEnergy;

  while ((dens = (TTMDensObj *) next())) {

    TTMThermalParticleBSQ *ptr;

    if(!fExclVolCorrection){
      
      TTMParticle *part = fPartSet->GetParticle(dens->GetID());
      ptr = new TTMThermalParticleBSQ(part,fParm);
      
      Double_t Width = part->GetWidth();
      Int_t Stat = part->GetStat();
      
      if(Stat == 0 || !fQStats){		       // Boltzmann stats
        if(Width == 0 || !fWidth){		       // no width
          PartEnergy = ptr->EnergyBoltzmannNoWidth();
        }else{					       // width
          PartEnergy = ptr->EnergyBoltzmannWidth();
        }
      }else{					       // Quantum stats	      
        if(Width == 0 || !fWidth){		       // no width
          PartEnergy = ptr->EnergyQStatNoWidth();
        }else{				    	       // width
          PartEnergy = ptr->EnergyQStatWidth();
        }
      }
   
    }else{ 

      TTMParticleSet CopySet = *fPartSet;

      TTMParameterSetBSQ CopyParm = *fParm;

      TTMParticle *CopyPart = CopySet.GetParticle(dens->GetID());

      ptr = ShiftedParticle(CopyPart,&CopyParm,fExclVolPressure);
     
      Double_t Width = CopyPart->GetWidth();
      Int_t Stat = CopyPart->GetStat();
      
      if(Stat == 0 || !fQStats){		       // Boltzmann stats
        if(Width == 0 || !fWidth){		       // no width
          PartEnergy = ptr->EnergyBoltzmannNoWidth();
        }else{					       // width
          PartEnergy = ptr->EnergyBoltzmannWidth();
        }
      }else{					       // Quantum stats	      
        if(Width == 0 || !fWidth){		       // no width
          PartEnergy = ptr->EnergyQStatNoWidth();
        }else{				    	       // width
          PartEnergy = ptr->EnergyQStatWidth();
        }
      }

      PartEnergy = PartEnergy/fExclVolDenominator;

    }

    dens->SetPrimEnergy(PartEnergy);
    fEnergy += PartEnergy;
    delete ptr;

  }
}

//__________________________________________________________________________
void TTMThermalModelBSQ::GenerateEntropyDens()
{
  // Iterates through the density hash table calculating the primordial
  // entropy density of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to impose constraints and ensure consistency between particle
  // and entropy density.
  //   	

  TIter next(fDensTable);
  TTMDensObj *dens;
  fEntropy = 0.;
  Double_t PartEntropy;

  while ((dens = (TTMDensObj *) next())) {
      
    TTMThermalParticleBSQ *ptr;

    if(!fExclVolCorrection){
      
      TTMParticle *part = fPartSet->GetParticle(dens->GetID());	   
      ptr = new TTMThermalParticleBSQ(part,fParm);
     
      Double_t Width = part->GetWidth();
      Int_t Stat = part->GetStat();
      
      if(Stat == 0 || !fQStats){		       // Boltzmann stats
        if(Width == 0 || !fWidth){		       // no width
          PartEntropy = ptr->EntropyBoltzmannNoWidth();
        }else{					       // width
          PartEntropy = ptr->EntropyBoltzmannWidth();
        }
      }else{					       // Quantum stats	      
        if(Width == 0 || !fWidth){		       // no width
          PartEntropy = ptr->EntropyQStatNoWidth();
        }else{				    	       // width
          PartEntropy = ptr->EntropyQStatWidth();
        }
      }

    }else{ 

      TTMParticleSet CopySet = *fPartSet;

      TTMParameterSetBSQ CopyParm = *fParm;

      TTMParticle *CopyPart = CopySet.GetParticle(dens->GetID());

      ptr = ShiftedParticle(CopyPart,&CopyParm,fExclVolPressure);
     
      Double_t Width = CopyPart->GetWidth();
      Int_t Stat = CopyPart->GetStat();
      
      if(Stat == 0 || !fQStats){		       // Boltzmann stats
        if(Width == 0 || !fWidth){		       // no width
          PartEntropy = ptr->EntropyBoltzmannNoWidth();
        }else{					       // width
          PartEntropy = ptr->EntropyBoltzmannWidth();
        }
      }else{					       // Quantum stats	      
        if(Width == 0 || !fWidth){		       // no width
          PartEntropy = ptr->EntropyQStatNoWidth();
        }else{				    	       // width
          PartEntropy = ptr->EntropyQStatWidth();
        }
      }

      PartEntropy = PartEntropy/fExclVolDenominator;

    }
      
    dens->SetPrimEntropy(PartEntropy);
    fEntropy += PartEntropy;
    delete ptr;

  }
}

//__________________________________________________________________________
void TTMThermalModelBSQ::GeneratePressure()
{
  // Iterates through the density hash table calculating the primordial
  // pressure of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to impose constraints and ensure consistency between particle
  // density and pressure.
  //   	

  TIter next(fDensTable);
  TTMDensObj *dens;
  fPressure = 0.;
  Double_t PartPressure;

  while ((dens = (TTMDensObj *) next())) {

    TTMThermalParticleBSQ *ptr;

    if(!fExclVolCorrection){

      TTMParticle *part = fPartSet->GetParticle(dens->GetID());	   
      ptr = new TTMThermalParticleBSQ(part,fParm);
      
      Double_t Width = part->GetWidth();
      Int_t Stat = part->GetStat();
      
      if(Stat == 0 || !fQStats){		       // Boltzmann stats
        if(Width == 0 || !fWidth){		       // no width
          PartPressure = ptr->PressureBoltzmannNoWidth();
        }else{					       // width
          PartPressure = ptr->PressureBoltzmannWidth();
        }
      }else{					       // Quantum stats	      
        if(Width == 0 || !fWidth){		       // no width
          PartPressure = ptr->PressureQStatNoWidth();
        }else{				    	       // width
          PartPressure = ptr->PressureQStatWidth();
        }
      }

    }else{ 

      TTMParticleSet CopySet = *fPartSet;

      TTMParameterSetBSQ CopyParm = *fParm;

      TTMParticle *CopyPart = CopySet.GetParticle(dens->GetID());

      ptr = ShiftedParticle(CopyPart,&CopyParm,fExclVolPressure);
     
      Double_t Width = CopyPart->GetWidth();
      Int_t Stat = CopyPart->GetStat();
      
      if(Stat == 0 || !fQStats){		       // Boltzmann stats
        if(Width == 0 || !fWidth){		       // no width
          PartPressure = ptr->PressureBoltzmannNoWidth();
        }else{					       // width
          PartPressure = ptr->PressureBoltzmannWidth();
        }
      }else{					       // Quantum stats	      
        if(Width == 0 || !fWidth){		       // no width
          PartPressure = ptr->PressureQStatNoWidth();
        }else{				    	       // width
          PartPressure = ptr->PressureQStatWidth();
        }
      }

    }
      
    dens->SetPrimPressure(PartPressure);
    fPressure += PartPressure;
    delete ptr;

  }
}

//__________________________________________________________________________
void TTMThermalModelBSQ::ListInfo()
{
  // Lists model information
  // 	
	
  cout <<"  ***************************************************************"
       <<"**************"
       << endl;
  cout <<"  ***************************** Thermal Model Info **************"
       <<"**************"
       << endl << endl;

  cout << "\t Particle set: " << endl << "\t\t" << fPartSet->GetFilename()
       << endl << endl;
   
  if(fQStats){
    cout << "\t Quantum Statistics " << endl;
  }else{
    cout << "\t Boltzmann Statistics " << endl;
  }
  if(fWidth){
    cout << "\t Resonance width included " << endl;
  }else{
    cout << "\t Resonance width excluded " << endl;
  }
  if(fExclVolCorrection){
    cout << "\t Excluded volume corrections included " << endl;
  }else{
    cout << "\t Excluded volume corrections excluded " << endl;
  }
  cout << endl << endl;
   
  fParm->List();
  cout <<"  ***************************** Thermal Quantities **************"
       <<"*************** "
       << endl << endl;

  cout << "\t S/V        = ";
  cout.width(10);
  cout << fStrange;
  cout << "\t";
  if (fParm->GetMuSConstrain()) {
    cout << "(constraint :";
    cout.width(8);
    cout << fParm->GetSDens() << ")";
  }
  cout << endl << endl;
  cout << "\t B/2Q     = ";
  cout.width(10);
  cout << fBaryon / 2. / fCharge;
  cout << "\t";
  if (fParm->GetMuQConstrain()) {
    cout << "(constraint :";
    cout.width(8);
    cout << fParm->GetB2Q() << ")";
  }
  cout << endl << endl;
  cout << "\t C/V        = ";
  cout.width(10);
  cout << fCharm;
  cout << "\t";
  if (fParm->GetMuCConstrain()) {
    cout << "(constraint :";
    cout.width(8);
    cout << fParm->GetCDens() << ")";
  }
  cout << endl << endl;
  cout << "\t b/V        = ";
  cout.width(10);
  cout << fBeauty;
  cout << "\t";
  if (fParm->GetMubConstrain()) {
    cout << "(constraint :";
    cout.width(8);
    cout << fParm->GetbDens() << ")";
  }
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
  if (fPressure != 0.) {
    cout << "\t\t\t\t P = " << fPressure << endl << endl;
  }
  ListStableParticles();
  cout << endl << endl;
  cout <<"  *************************************************************"
       <<"*****************"
       << endl;
  cout <<"  *************************************************************"
       <<"*****************"
       << endl << endl;
}

//__________________________________________________________________________
void TTMThermalModelBSQ::CalcExclVolPressure()
{
  // Calculates the total pressure in the system with excluded volume 
  // corrections 
  // 	

  TTMParticleSet CopySet = *fPartSet;

  TTMParameterSetBSQ CopyParm = *fParm;

  TTMThermalModelBSQ CopyMod(&CopySet,&CopyParm,fQStats,fWidth);

  CopyMod.SetExcludedVolume(kFALSE);
  CopyMod.PrimPartDens();
  CopyMod.GeneratePressure();

  Double_t P = CopyMod.GetPressure();

  fExclVolPressure = FindExclVolPressure(this,P);

}

//__________________________________________________________________________
Double_t TTMThermalModelBSQ::ExclVolShiftedPressure(Double_t px)
{
  // Calculates the total ideal gas pressure corresponding to chemical 
  // potentials shifted by an amount px
  //

  TTMParticleSet CopySet = *fPartSet;
  TTMParameterSetBSQ CopyParm = *fParm;

  Double_t TotalShiftedPressure = 0.;

  TIter next(CopySet.GetParticleTable());
  TTMParticle *part;

  while ((part = (TTMParticle *) next())) {

    TTMThermalParticleBSQ *ptr = ShiftedParticle(part,&CopyParm,px);

    Double_t PartPressure;
    Double_t Width = part->GetWidth();
    Int_t Stat = part->GetStat();

    if(Stat == 0 || !fQStats){		       // Boltzmann stats
      if(Width == 0 || !fWidth){		       // no width
        PartPressure = ptr->PressureBoltzmannNoWidth();
      }else{					       // width
        PartPressure = ptr->PressureBoltzmannWidth();
      }
    }else{					       // Quantum stats	      
      if(Width == 0 || !fWidth){		       // no width
        PartPressure = ptr->PressureQStatNoWidth();
      }else{				    	       // width
        PartPressure = ptr->PressureQStatWidth();
      }
    }
      
    TotalShiftedPressure += PartPressure;
    delete ptr;
  }

  return TotalShiftedPressure;

}

//__________________________________________________________________________
TTMThermalParticleBSQ* TTMThermalModelBSQ::ShiftedParticle(TTMParticle* CopyPart, 
                                                           TTMParameterSetBSQ* CopyParm, Double_t P)
{
  // Returns a pointer to a TTMThermalParticleBSQ object with a chemical 
  // potential shifted by an amount - (hadron volume) * P. Memory is allocated
  // off the heap and must be cleaned up after use.
  //

  Double_t mu = CopyPart->GetB()*fParm->GetMuB() + CopyPart->GetS()*fParm->GetMuS() + CopyPart->GetQ()*fParm->GetMuQ() + CopyPart->GetCharm()*fParm->GetMuC() + CopyPart->GetBeauty()*fParm->GetMub();

  CopyPart->SetB(1); 
  CopyPart->SetS(0); 
  CopyPart->SetQ(0); 
  CopyPart->SetCharm(0); 
  CopyPart->SetBeauty(0); 
  Double_t radius = CopyPart->GetRadius();
  Double_t volume = 4. * 4./3. * TMath::Pi() * pow(radius,3);

  Double_t shift = volume * P;

  CopyParm->SetMuB(mu-shift);
  CopyParm->SetMuS(0.);
  CopyParm->SetMuQ(0.);
  CopyParm->SetMuC(0.);
  CopyParm->SetMub(0.);

  TTMThermalParticleBSQ *ptr = new TTMThermalParticleBSQ(CopyPart,CopyParm);

  return ptr;

}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainTotalBaryonDensity(Double_t nb)
{
  // Constrains muB by the given nb in the case where either just muS or else 
  // both muS and muQ are also to be constrained while muC and mub are not 
  // constrained. The other possibilities are 
  // yet to be coded. 1 is returned if there are problems 0 otherwise. 
  //

  Int_t check;

  if(fParm->GetMuSConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){

    if(fParm->GetMuQConstrain()){
      check = BSQConstrainSQBDens(this,nb);
      if(check){
        cout<<"Not Constrained"<<endl;
      }

    }else{
      check = BSQConstrainSBDens(this,nb);
      if(check){
        cout<<"Not Constrained"<<endl;
      }
      
    }

    }else{

      check = 1;
      cout<<"Such constraint not yet coded"<<endl;
    }

      return check;

}

//__________________________________________________________________________
Int_t TTMThermalModelBSQ::ConstrainNetBaryonDensity(Double_t nb)
{
  // Constrains muB by the given nb in the case where 
  // both muS and muQ are also to be constrained while muC and mub are not 
  // to be constrained. The other possibilities are 
  // yet to be coded. 1 is returned if there are problems 0 otherwise. 
  //

  Int_t check;

  if(fParm->GetMuSConstrain() && fParm->GetMuQConstrain() && !fParm->GetMuCConstrain() && !fParm->GetMubConstrain()){

    check = BSQConstrainSQNetBDens(this,nb);
    if(check){
      cout<<"Not Constrained"<<endl;
    }

  }else{

    check = 1;
    cout<<"Such constraint not yet coded"<<endl;
  }

  return check;

}
