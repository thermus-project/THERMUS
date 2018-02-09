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
// Base class for all thermal model classes.
//  

#include <TThermalModel.h>

TString Int_2_String(Int_t x);

ClassImp(TTMThermalModel)
	
//__________________________________________________________________________
TTMDensObj *TTMThermalModel::GetDensities(Int_t id)
{
  // Returns a pointer to the density object of the particle with specified 
  // ID
  //

  TString name = Int_2_String(id);

  TTMDensObj *dens = (TTMDensObj *) fDensTable->FindObject(name);

  if (!dens) {
    cout << "WARNING: Particle " << id << " not in table!" << endl;
    return 0;
  } else {
    return dens;
  }
}

//__________________________________________________________________________
void TTMThermalModel::GenerateDecayPartDens()
{
  // Iterates through the hash table of TTMDensObj objects and calculates the
  // decay contributions to the densities of the particles considered 
  // stable in the particle set.
  //

  THashTable *part_table = GetParticleTable();
  TIter next_d(part_table);
  TTMParticle *d;

  while ((d = (TTMParticle *) next_d())) {
    TTMDensObj *daughter_dens = GetDensities(d->GetID());
    daughter_dens->SetDecayDensity(0.);

    if (d->GetStable())       // stable particle
      {
        Double_t decay = 0.;
        TIter next_p(part_table);
        TTMParticle *p;

        while ((p = (TTMParticle *) next_p())) {
          if (!p->GetStable() && p->GetDecay(d->GetID())) {
            TTMDensObj *parent_dens = GetDensities(p->GetID());
            Double_t BR = p->GetDecay(d->GetID())->GetBRatio();
            decay += BR * (parent_dens->GetPrimDensity());
          }
        }
        daughter_dens->SetDecayDensity(decay);
      }
  }
}

//__________________________________________________________________________
void TTMThermalModel::GenerateDecayPartDens(Int_t id)
{
  // Calculates the decay contribution to particle with ID id
  // based on the primordial densities already in the hash table
  // provided the particle is stable!
  // (NB must first populate the hash table with primordial densities)
  //

  THashTable *part_table = GetParticleTable();
  TTMParticle *d = fPartSet->GetParticle(id);

  TTMDensObj *daughter_dens = GetDensities(id);
  daughter_dens->SetDecayDensity(0.);

  if (d->GetStable())       // stable particle
    {
      Double_t decay = 0.;
      TIter next_p(part_table);
      TTMParticle *p;

      while ((p = (TTMParticle *) next_p())) {
        if (!p->GetStable() && p->GetDecay(d->GetID())) {
          TTMDensObj *parent_dens = GetDensities(p->GetID());
          Double_t BR = p->GetDecay(d->GetID())->GetBRatio();
          decay += BR * (parent_dens->GetPrimDensity());
        }
      }
      daughter_dens->SetDecayDensity(decay);
    }
}

//__________________________________________________________________________
void TTMThermalModel::ListDecayContributions(Int_t d_id)
{
  // Lists the decay contributions (in percentage and absolute terms) to daughter
  // with ID d_id. The primordial and decay densities must already appear in the 
  // density hash table (i.e. GenerateParticleDens() should be run first).
  //

  TList* list = new TList();
  fPartSet->GetParents(list,d_id);
  TIter next(list);
  TTMDecay *decay;
  Double_t total = 0.;
  Double_t n_decay = GetDensities(d_id)->GetDecayDensity();
  cout << n_decay << endl;
  while((decay = (TTMDecay *)next())){
    Int_t p_id = decay->GetParentID();
    Double_t contr = GetDensities(p_id)->GetPrimDensity()*decay->GetBRatio();
    cout << p_id << "\t" << fPartSet->GetParticle(p_id)->GetPartName() << "\t" << contr << "\t" 
         << contr/n_decay * 100. << "%" << endl;
    total += contr;
  }
  cout << total << endl;
  delete list;
}

//__________________________________________________________________________
void TTMThermalModel::ListDecayContribution(Int_t p_id, Int_t d_id)
{
  // Lists the contribution (in percentage and absolute terms) of the decay from 
  // parent with ID p_id to daughter with ID d_id. The hash table must already 
  // contain the primordial and decay densities (i.e. first run GenerateParticleDens())
  //

  TList* list = new TList();
  fPartSet->GetParents(list,d_id);
  TIter next(list);
  TTMDecay *decay;
  Double_t n_decay = GetDensities(d_id)->GetDecayDensity();
  Bool_t found = false;
   
  while((decay = (TTMDecay *)next())){
    if(found){break;}   
    if(decay->GetParentID() == p_id){
      Double_t contr = GetDensities(p_id)->GetPrimDensity()*decay->GetBRatio();
      cout << p_id << "\t" << fPartSet->GetParticle(p_id)->GetPartName() << "\t" << contr << "\t" 
           << contr/n_decay * 100. << "%" << endl;
      found = true;
    }
  }
  if(!found){
    cout << "No such decay" << endl;
  }
  delete list;
}

//__________________________________________________________________________
void TTMThermalModel::CalcWroblewski()
{
  // Calculates the Wroblewski factor in the model from the primordial 
  // particle densities in the hash table after it has been populated 
  // by functions in the derived classes. Only correct if charm and 
  // beauty are excluded.
  //

  Double_t absstrange, str, nstrange, nstr;
  Double_t absstrangemesons, absstrangebaryons, nstrangemesons,
    nstrangebaryons;
  absstrangemesons = absstrangebaryons = nstrangemesons
    = nstrangebaryons = 0.;

  THashTable *part_table = GetParticleTable();
  TIter next(part_table);
  TTMParticle *p;

  while ((p = (TTMParticle *) next())) {
    TTMDensObj *dens = GetDensities(p->GetID());
    Double_t part_dens = dens->GetPrimDensity();

    if (p->GetB() == 0.)      //mesons
      {
        absstrangemesons += p->GetSContent()
          * part_dens;
        nstrangemesons += (2. - p->GetSContent())
          * part_dens;

      } else                    //baryons
        {
          absstrangebaryons += p->GetSContent()
            * part_dens;
          nstrangebaryons += (3. - p->GetSContent())
            * part_dens;

        }
  }

  absstrange = absstrangemesons + absstrangebaryons;
  //no. of s + sbar
  str = absstrange / 2.;
  //no. of s sbar (s = sbar since S=0)
  nstrange = nstrangemesons + nstrangebaryons - fBaryon * 3.;
  //no. of new ud
  nstr = nstrange / 2.;
  fWroblewski = 2. * str / nstr;

}

//__________________________________________________________________________
void TTMThermalModel::ListStableParticles()
{
  // Lists the stable particles in the model
  //

  fPartSet->ListStableParticles();
}

//__________________________________________________________________________
void TTMThermalModel::ListStableDensities()
{
  // Lists the densities of the stable particles in the model
  //  

  cout << "\t Densities of Stable Particles: " << endl;
  THashTable *part_table = GetParticleTable();
  TIter next(part_table);
  TTMParticle *p;

  while ((p = (TTMParticle *) next())) {
    if (p->GetStable()) {
      if (GetDensities(p->GetID())) {
        TTMDensObj *dens = GetDensities(p->GetID());
        cout << "\t ***** " << p->GetPartName() << " ***** :" << endl;
        cout << "\t\t Primordial Density:\t"
             << dens->GetPrimDensity() << endl;
        cout << "\t\t Decay Density Contribution:\t"
             << dens->GetDecayDensity() << endl << endl;
      }
    }
  }
}

//__________________________________________________________________________
TTMThermalModel::~TTMThermalModel()
{
  // return memory to heap
  //

  if(fDensTable){
    fDensTable->Delete();
    delete fDensTable;
  }
}
