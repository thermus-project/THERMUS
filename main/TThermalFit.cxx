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

// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Base class for all thermal fit classes.
//   

#include <TThermalFit.h>
#include <cstdio>

void fit_function(TTMThermalFit *fit,Int_t flag = 0);
TString Int_2_String(Int_t x);

ClassImp(TTMThermalFit)

//__________________________________________________________________________
TTMThermalFit::TTMThermalFit()
{
   
  fPartSet = (TTMParticleSet *) 0;
  fChiSquare = 0;
  fQuadDev = 0;
  fYields = new TList();
  fYields->SetOwner(kTRUE);
  fMinuit = (TMinuit *) 0;
}

//__________________________________________________________________________
TTMThermalFit::~TTMThermalFit()
{
  if(fYields){
    fYields->Delete();
    delete fYields;
  }
  if(fMinuit) delete fMinuit;
}

//__________________________________________________________________________
void TTMThermalFit::AddYield(TTMYield *yield)
{

  TString name = yield->GetName();

  if((TTMYield *)fYields->FindObject(name)){
    cout << "Give Yield Unique Descriptor Before Adding" << endl; 
  }else{
    fYields->Add(yield);
  }
}

//__________________________________________________________________________
TString TTMThermalFit::GetTMName(Int_t id1,Int_t id2,TString descr)
{
  TString name;

  if(id2 == 0){
    if(id1==1){
      name = name.Append("Npart");
    }else if(id1==2){
      name = name.Append("h-");
    }else if(id1==3){
      name = name.Append("h+");
    }else if(id1==33340){
      name = name.Append("Omega + anti-Omega");
    }else{	
      name = fPartSet->GetParticle(id1)->GetPartName();
    }
  }else{
    if(id1==3130){
      name = "<";
      name = name.Append(fPartSet->GetParticle(313)->GetPartName());
      name = name.Append(">");
      name = name.Append("/");
      name = name.Append(fPartSet->GetParticle(id2)->GetPartName());
    }else if(id1==33340){
      name = name.Append("Omega + anti-Omega");
      name = name.Append("/");
      name = name.Append(fPartSet->GetParticle(id2)->GetPartName());
    }else if(id2==2){
      name = fPartSet->GetParticle(id1)->GetPartName();
      name = name.Append("/");
      name = name.Append("h-");
    }else if(id2==3){
      name = fPartSet->GetParticle(id1)->GetPartName();
      name = name.Append("/");
      name = name.Append("h+");
    }else if(id2==33340){
      name = fPartSet->GetParticle(id1)->GetPartName();
      name = name.Append("/");
      name = name.Append("Omega + anti-Omega");
    }else{
      name = fPartSet->GetParticle(id1)->GetPartName();
      name = name.Append("/");
      name = name.Append(fPartSet->GetParticle(id2)->GetPartName());
    }
  }
  if(descr!=""){
    name.Append(" ");
    name.Append(descr);
  } 
  return name;
}

//__________________________________________________________________________
void TTMThermalFit::RemoveYield(Int_t id1,Int_t id2,TString descr)
{
  TString name = GetTMName(id1,id2,descr);
  if(fYields->FindObject(name)){
    TTMYield *yield = (TTMYield *)fYields->FindObject(name);
    fYields->Remove(yield);
    delete yield;
  }else{
    cout << "Yield not in list" << endl;
  }
}

//__________________________________________________________________________
TTMYield* TTMThermalFit::GetYield(Int_t id1,Int_t id2,TString descr)
{
  TString name = GetTMName(id1,id2,descr);
 
  if(fYields->FindObject(name)){
    return (TTMYield *)fYields->FindObject(name);
  }
  else{
    cout << "Yield not in list" << endl;
    return 0; 
  }
}

//__________________________________________________________________________
void TTMThermalFit::InputExpYields(const char *file)
{
  // Inserts the experimental yields listed in the specified file 
  // in *fYields.
  //
  // The input file has the following format:
  //
  // ID \t Descriptor \t exp_yield \t exp_error\n
  // ID (of numerator) \t ID (of denominator) \t Descriptor \t exp_ratio \t exp_error\n
  // etc .... 
  //
  // In addition to all of the particle id's in fPartSet, the following are
  // also allowed: 
  // 
  // 		ID=1 	: Npart
  //		ID=2	: h-
  //		ID=3 	: h+
  //

  fYields->Delete();
  ifstream data(file);
  if (!data) {
    cout << "WARNING: Cannot open file: " << file << endl;
  } else {
    char ch;
    Int_t data_points = 0;
    Int_t tabs = 0;           // counts tabs to determine yield or ratio
    Int_t flags[100];         // records whether yield (0) or ratio (1)

    while (!data.eof()) {
      data.get(ch);
      if (ch == '#') { data.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); // BH 10/05/2014
        flags[data_points] = -1;
        data_points++;
      }
      if (ch == '\t') {
        tabs++;
      }
      if (ch == '\n') {
        if (tabs == 4) {
          flags[data_points] = 1;
        } else if (tabs == 3) {
          flags[data_points] = 0;
        } else if (!data.eof()) {
          cout << "WARNING: " << file << " has incorrect format!"
               << endl;
        }
        data_points++;
        tabs = 0;
      }
    }
    data_points = data_points - 1;
    data.close();

    ifstream measure(file);
      
    for (Int_t i = 1; i <= data_points; i++) {
      TString name, descr;
      Int_t ID_1, ID_2;
      Double_t value, error;
      if (flags[i - 1] == -1) { //Comment
        measure >> descr;
        measure.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
      } else if (flags[i - 1] == 0) {	//Yield
        measure >> ID_1;
        measure >> descr;
        measure >> value;
        measure >> error;
        name = GetTMName(ID_1,0,descr);
        TTMYield *yield = new TTMYield(name, value, error, ID_1);
        if (value < 1e-15) {
          cout << "WARNING: yield on line "<< i <<" of the input file is very probably null so it will be predicted !" << endl;
          yield->Predict();
        }
        yield->SetPartSet(fPartSet);
        AddYield(yield);
      } else if (flags[i - 1] == 1) {	//Ratio
        measure >> ID_1;
        measure >> ID_2;
        measure >> descr;
        measure >> value;
        measure >> error;
               
        name = GetTMName(ID_1,ID_2,descr);
        TTMYield *yield = new TTMYield(name, value, error, ID_1, ID_2);
        if (value < 1e-15) {
          cout << "WARNING: ratio on line "<< i <<" of the input file is very probably null so it will be predicted !" << endl;
          yield->Predict();
        }
        yield->SetPartSet(fPartSet,fPartSet);
        AddYield(yield);
      }
    }
    measure.close();
  }
}

//__________________________________________________________________________
void TTMThermalFit::FitData(Int_t flag)
{
  // flag = 0 : chi-square fit
  // flag = 1 : quad-deviation fit
  // 
  // The experimental yields must first be input using the function
  // InputExpYields(char *file).
  // 
  // The function checks which parameters are to be fitted and performs 
  // the required fit. On completion, the model predictions for each fitted 
  // yield appear in the list of TTMYield objects, while the parameter set 
  // reflects the best-fit parameters.  
  //   
  // To calculate the model predictions, GenerateYields() is run, which
  // calls GenerateParticleDens() and GenerateDecayPartDens() through an
  // intermediate model object.
  // 
    
  fit_function(this, flag);
}

//__________________________________________________________________________
void TTMThermalFit::GenerateYields()
{
  // Calculates the primordial particle densities of all particles in the 
  // set, then iterates through the list of yields, calculates 
  // their decay contributions and inserts the new model
  // predictions into the TTMYield objects, and calculates chi-squared 
  // and the quadratic deviation.  
  //             
                                                      
  fChiSquare = fQuadDev = 0.;

  TIter nextYield(fYields);
  TTMYield *next_yield;

  TTMThermalModel *model = GenerateThermalModel(fPartSet);

  model->GenerateParticleDens();
  while ((next_yield = (TTMYield *) nextYield())) {
    Int_t id_1, id_2;
    Double_t volume = GetParameterSet()->GetVolume();
    id_1 = next_yield->GetID1();
    id_2 = next_yield->GetID2();
    if (id_2 == 0)            // Yield
      {
        model->SetParticleSet(next_yield->GetPartSet1());
        if (id_1 == 1)         //Npart
          {
            next_yield->SetModelValue(model->GetBaryon() * volume);
          } else if(id_1 == 2){
            Double_t hmin = 0.;
            model->GenerateDecayPartDens();
            TIter next(model->GetDensityTable());
            TTMDensObj *dens;
	    while((dens = (TTMDensObj *) next())){
              TTMParticle *part = (next_yield->GetPartSet1())->GetParticle(dens->GetID());
              if(part->GetStable()&&part->GetQ()< 0){
                hmin += dens->GetFinalDensity();	
              }
		
            } 
            next_yield->SetModelValue(hmin * volume);
          } else if(id_1 == 3){
            Double_t hplus = 0.;
            model->GenerateDecayPartDens();
            TIter next(model->GetDensityTable());
            TTMDensObj *dens;
	    while((dens = (TTMDensObj*) next())){
              TTMParticle *part = (next_yield->GetPartSet1())->GetParticle(dens->GetID());
              if(part->GetStable()&&part->GetQ()> 0){
                hplus += dens->GetFinalDensity();	
              }
		
            } 
            next_yield->SetModelValue(hplus * volume);
          } else if(id_1 == 33340){
	    model->GenerateDecayPartDens(3334);
	    model->GenerateDecayPartDens(-3334);
            next_yield->SetModelValue((model->GetDensities(3334)->GetFinalDensity() + 
                                       model->GetDensities(-3334)->GetFinalDensity()) * volume);
          } else {
	    model->GenerateDecayPartDens(id_1);
            TTMDensObj *part_dens = model->GetDensities(id_1);
            next_yield->SetModelValue(part_dens->GetFinalDensity() *
                                      volume);
          }
      } else    // Ratio
        {
          if(id_1 == 3130){
            model->SetParticleSet(next_yield->GetPartSet1());
            model->GenerateDecayPartDens(313);
            model->GenerateDecayPartDens(-313);
	    Double_t kstar = model->GetDensities(313)->GetFinalDensity();
	    Double_t akstar = model->GetDensities(-313)->GetFinalDensity();
	    Double_t num = (kstar + akstar)/2.;
            model->SetParticleSet(next_yield->GetPartSet2());
            model->GenerateDecayPartDens(id_2);
            Double_t den = model->GetDensities(id_2)->GetFinalDensity();
            next_yield->SetModelValue(num/den);
          } else if(id_1 == 33340){
            model->SetParticleSet(next_yield->GetPartSet1());
            model->GenerateDecayPartDens(3334);
            model->GenerateDecayPartDens(-3334);
	    Double_t omega = model->GetDensities(3334)->GetFinalDensity();
	    Double_t aomega = model->GetDensities(-3334)->GetFinalDensity();
	    Double_t num = omega + aomega;
            model->SetParticleSet(next_yield->GetPartSet2());
            model->GenerateDecayPartDens(id_2);
            Double_t den = model->GetDensities(id_2)->GetFinalDensity();
            next_yield->SetModelValue(num/den);
          } else if(id_2 == 2){
            model->SetParticleSet(next_yield->GetPartSet1());
            model->GenerateDecayPartDens(id_1);
	    Double_t num = model->GetDensities(id_1)->GetFinalDensity();
            model->SetParticleSet(next_yield->GetPartSet2());
            Double_t hmin = 0.;
            model->GenerateDecayPartDens();
            TIter next(model->GetDensityTable());
            TTMDensObj *dens;
	    while((dens = (TTMDensObj*) next())){
              TTMParticle *part = (next_yield->GetPartSet2())->GetParticle(dens->GetID());
              if(part->GetStable()&&part->GetQ()< 0){
                hmin += dens->GetFinalDensity();	
              }
            } 
	    Double_t den = hmin;
            next_yield->SetModelValue(num/den);
          } else if(id_2 == 3){
            model->SetParticleSet(next_yield->GetPartSet1());
            model->GenerateDecayPartDens(id_1);
	    Double_t num = model->GetDensities(id_1)->GetFinalDensity();
            model->SetParticleSet(next_yield->GetPartSet2());
            Double_t hplus = 0.;
            model->GenerateDecayPartDens();
            TIter next(model->GetDensityTable());
            TTMDensObj *dens;
	    while((dens = (TTMDensObj*) next())){
              TTMParticle *part = (next_yield->GetPartSet2())->GetParticle(dens->GetID());
              if(part->GetStable()&&part->GetQ()> 0){
                hplus += dens->GetFinalDensity();	
              }
            } 
	    Double_t den = hplus;
            next_yield->SetModelValue(num/den);
          } else if(id_2 == 33340){
            model->SetParticleSet(next_yield->GetPartSet1());
            model->GenerateDecayPartDens(id_1);
	    Double_t num = model->GetDensities(id_1)->GetFinalDensity();
            model->SetParticleSet(next_yield->GetPartSet2());
            model->GenerateDecayPartDens(3334);
            model->GenerateDecayPartDens(-3334);
            Double_t den = model->GetDensities(3334)->GetFinalDensity() + model->GetDensities(-3334)->GetFinalDensity();
            next_yield->SetModelValue(num/den);
          } else{
            model->SetParticleSet(next_yield->GetPartSet1());
            model->GenerateDecayPartDens(id_1);
	    Double_t num = model->GetDensities(id_1)->GetFinalDensity();
            model->SetParticleSet(next_yield->GetPartSet2());
            model->GenerateDecayPartDens(id_2);
            Double_t den = model->GetDensities(id_2)->GetFinalDensity();
            next_yield->SetModelValue(num/den);
          }
        }
    if(next_yield->GetFit()){
      fChiSquare += pow(next_yield->GetStdDev(), 2);
      fQuadDev += pow(next_yield->GetQuadDev(), 2);
    }
  }
}

//__________________________________________________________________________
void TTMThermalFit::ListYields()
{
  // Lists the experimental values and model predictions 
  // 

  cout << "******************************** " << endl << endl;

  TIter next(fYields);
  TTMYield *yield;

  while ((yield = (TTMYield *) next())) {
    yield->List();
  }
  cout << "  **************************************************";
  cout << "****************************" << endl << endl;

}

//__________________________________________________________________________
void TTMThermalFit::ListMinuitInfo()
{
  // Lists the info of the Minuit fit
  // 

  if(fMinuit){
   
    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    fMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    cout << "\t\t\t" << "FCN = " << amin << endl;
    cout << "\t\t\t" << "EDM = " << edm << endl;
    cout << "\t\t\t" << "Errdef = " << errdef << endl;
   
    if(icstat == 0){cout << "\t\t" << "Covariance matrix not calculated" << endl;}
    else if(icstat == 1){cout << "\t\t" << "Covariance matrix approximated only - not accurate" << endl;}
    else if(icstat == 2){cout << "\t\t" << "Full covariance matrix calculated but forced positive definite" << endl;}
    else if(icstat == 3){cout << "\t\t" << "Full accurate covariance matrix calculated" << endl;}

    fMinuit->mnprin(2, amin);
    fMinuit->mnmatu(1);
  } else{
    cout << "Run FitData() to instantiate a TMinuit object" << endl;
  }

}
