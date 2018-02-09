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
// A collection of particles. A hash table is used for quick access to 
// entries. The entries each have a unique fName (based on their id's).
// When searching for entries, the id must first be converted into a 
// TString in exactly the same way as the TTMParticle fNames are set. 
// The function Int_2_String(Int_t) is used for this. 
//

#include <TTMParticleSet.h>
#include <TTMDecay.h>
#include <TTMIDObj.h>
#include <TTMDecayChannel.h>

ClassImp(TTMParticleSet)

//__________________________________________________________________________
TTMParticleSet::TTMParticleSet()
{
  fPartTable = new THashTable();
  fFilename = TString("");
  fParticleNumber = 0;
}

//__________________________________________________________________________
TTMParticleSet::TTMParticleSet(const TTMParticleSet &obj)
{
  fPartTable = new THashTable();
  fParticleNumber = 0;
  TIter next(obj.GetParticleTable());
  TTMParticle *a;
  while((a = (TTMParticle *)next())){
    TTMParticle *part = new TTMParticle();
   
    part->SetPartName(a->GetPartName());
    part->SetID(a->GetID());
    if(a->GetStable()){
      part->SetStable();
    }else{
      part->SetUnstable();
    }
    part->SetDeg(a->GetDeg());
    part->SetStat(a->GetStat());
    part->SetB(a->GetB());
    part->SetS(a->GetS());
    part->SetQ(a->GetQ());
    part->SetCharm(a->GetCharm());
    part->SetBeauty(a->GetBeauty());
    part->SetTop(a->GetTop());
    part->SetMass(a->GetMass());
    part->SetWidth(a->GetWidth());
    part->SetThreshold(a->GetThreshold());
    part->SetSContent(a->GetSContent());
    part->SetCContent(a->GetCContent());
    part->SetbContent(a->GetbContent());
    part->SetTContent(a->GetTContent());

    part->SetThresholdCalc(a->GetThresholdCalc());
    part->SetThresholdFlag(a->GetThresholdFlag());

    part->SetRadius(a->GetRadius());
   
    TIter nextch(a->GetDecayChannels());
    TTMDecayChannel *old;

    while((old=(TTMDecayChannel *)nextch())){
      TIter nextid(old->GetDaughterList());
      TList *ndlist = new TList();
      TTMIDObj *oldid;

      while((oldid = (TTMIDObj *)nextid())){
        Int_t nid = oldid->GetID();
        TTMIDObj *newid = new TTMIDObj(nid);
        ndlist->Add(newid);
      }

      TTMDecayChannel *nch = new TTMDecayChannel(old->GetBRatio(),ndlist);
      part->GetDecayChannels()->Add(nch);
    }
    part->UpdateDecaySummary();
    fPartTable->Add(part);
    fParticleNumber++;
  }   
  GenerateBRatios();
  fFilename = obj.GetFilename();
}

//__________________________________________________________________________
TTMParticleSet::TTMParticleSet(const char* file, Bool_t CB)
{
  // Populates the hash table with particles listed in the specified file.
  // This file lists only PARTICLES. If a particle listed in the file has
  // an anti-particle (i.e. if it has any non-zero quantum numbers), then
  // the anti-particle is automatically entered. In this way, particle and 
  // anti-particle are treated symmetrically.
  //
  // With CB set to true, the required file format is as follows:
  // File format: (tab separated)
  //   
  // Stable Name ID Deg. Stat. Mass S  B  Q  C  B  |S| Width Threshold (*) \n 
  // \n
  // Stable Name ID Deg. Stat. Mass S  B  Q  C  B  |S| Width Threshold (*) \n 
  // \n
  // etc ...
  //
  // With CB set to false, the file format is as above but without the 
  // C and B columns
  //
  // If the width is non-zero, the decay products in the channel which 
  // determines the threshold is listed (*).    
  //

  fPartTable = new THashTable();
  fFilename = file;
  fParticleNumber = 0;
  ifstream data(file);
  if (!data) {
    cout << "WARNING: Cannot open file " << file << endl;
  } else {
    char ch;
    Int_t lines = 0;
    data.seekg(0, ios::beg);
    while (!data.eof()) {
      data.get(ch);
    
      if (ch == '\n') {
        lines++;
      }
    }
    lines = lines - 1;
      
    data.close();
    ifstream table(file);

    for (Int_t i = 1; i <= lines / 2; i++)	// loop over entries in file
      {

        TTMParticle* part = new TTMParticle();

        Int_t temp_int;
        Double_t temp_double;
        TString temp_string;
        char temp_ch;

        table >> temp_int;
        if (temp_int) {
          part->SetStable();
        } else {
          part->SetUnstable();
        }

        table >> temp_string;
        part->SetPartName(temp_string);

        table >> temp_int;
        part->SetID(temp_int);

        table >> temp_int;
        part->SetDeg(temp_int);

        table >> temp_int;
        part->SetStat(temp_int);

        table >> temp_double;
        part->SetMass(temp_double);

        table >> temp_int;
        part->SetS(temp_int);

        table >> temp_int;
        part->SetB(temp_int);

        table >> temp_int;
        part->SetQ(temp_int);

	if(CB){
           table >> temp_int;
           part->SetCharm(temp_int);
           table >> temp_int;
           part->SetBeauty(temp_int);
        }else{
           part->SetCharm(0);
           part->SetBeauty(0);
        }

	part->SetTop(0);

        table >> temp_double;
        part->SetSContent(temp_double);

        part->SetCContent(TMath::Abs(part->GetCharm()));
        part->SetbContent(TMath::Abs(part->GetBeauty()));
        part->SetTContent(TMath::Abs(part->GetTop()));

        table >> temp_double;
        part->SetWidth(temp_double);

        table >> temp_double;
        part->SetThreshold(temp_double);

        part->SetRadius(0.);

        if (temp_double != 0) {
          table >> temp_string;
        }
        table.get(temp_ch);

        fPartTable->Add(part);
        fParticleNumber++;

        if (part->GetB() != 0 || part->GetS() != 0 || part->GetQ() != 0 || part->GetCharm() != 0 || part->GetBeauty() != 0 || part->GetTop() != 0) {
          TTMParticle* apart = new TTMParticle();
          if (part->GetStable()) {
            apart->SetStable();
          } else {
            apart->SetUnstable();
          }
          TString temp;
          temp = part->GetPartName();
          temp = temp.Prepend("anti-");
          apart->SetPartName(temp);
          apart->SetID(-part->GetID());
          apart->SetDeg(part->GetDeg());
          apart->SetStat(part->GetStat());
          apart->SetMass(part->GetMass());
          apart->SetS(-part->GetS());
          apart->SetB(-part->GetB());
          apart->SetQ(-part->GetQ());
          apart->SetCharm(-part->GetCharm());
	  apart->SetBeauty(-part->GetBeauty());
          apart->SetTop(-part->GetTop());
          apart->SetSContent(part->GetSContent());
          apart->SetCContent(part->GetCContent());
          apart->SetbContent(part->GetbContent());
          apart->SetTContent(part->GetTContent());
          apart->SetWidth(part->GetWidth());
          apart->SetThreshold(part->GetThreshold());

          apart->SetRadius(part->GetRadius());

          fPartTable->Add(apart);
          fParticleNumber++;
        }
      }
    fPartTable->Rehash(fPartTable->Capacity());
    table.close();
  }
}

//__________________________________________________________________________
TTMParticleSet::TTMParticleSet(TDatabasePDG *pdg)
{
  // Instantiates a TTMParticleSet object and populates the set with 
  // particles listed in the TDatabasePDG object pointed at by pdg.
  // Only mesons and baryons are included-- generators, sparticles,
  // quarks etc excluded. At the moment, the default ROOT file containing the 
  // particle properties does not include the charm, degeneracy (no spin!), 
  // strangeness (and hence |S|), beauty, topness or THRESHOLD of the particle! 
  // Also the lifetime, parity, isospin, and I3 are absent from the default file.
  // Anti-particles should be included in the TDatabasePDG object too as they are 
  // not automatically added in this constructor.  
  //

  fPartTable = new THashTable();
  fFilename = TString("");
  fParticleNumber = 0;
   
  const THashList *list = pdg->ParticleList();
  TIter next(list);
  TParticlePDG *pdg_part;

  while((pdg_part = (TParticlePDG *)next())){

    TString type = pdg_part->ParticleClass();

    if(type == TString("Meson") || type == TString("CharmedMeson") || type == TString("HiddenCharmMeson") || type == TString("B-Meson")){

      TTMParticle *part = new TTMParticle();
      part->SetPartName(TString(pdg_part->GetName()));
      part->SetID(pdg_part->PdgCode());

      if(pdg_part->Stable()){
        part->SetStable();
      }else{
        part->SetUnstable();
      }

      part->SetDeg(Int_t(2*pdg_part->Spin()+1));
      part->SetStat(-1);
      part->SetB(0);
      part->SetS(Int_t(pdg_part->Strangeness()));
      part->SetQ(Int_t(pdg_part->Charge()/3.));
      part->SetCharm(Int_t(pdg_part->Charm()));
      part->SetBeauty(Int_t(pdg_part->Beauty()));
      part->SetTop(Int_t(pdg_part->Top()));
      part->SetMass(pdg_part->Mass());
      part->SetWidth(pdg_part->Width());
      part->SetThreshold(0.);
      part->SetSContent(TMath::Abs(pdg_part->Strangeness()));         
      part->SetCContent(TMath::Abs(pdg_part->Charm()));         
      part->SetbContent(TMath::Abs(pdg_part->Beauty()));  
      part->SetTContent(TMath::Abs(pdg_part->Top()));  

      part->SetRadius(0.);

      fPartTable->Add(part);
      fParticleNumber++;

    }else if(type == TString("Baryon") || type == TString("CharmedBaryon") || type == TString("B-Baryon")){

      TTMParticle *part = new TTMParticle();
      part->SetPartName(TString(pdg_part->GetName()));
      part->SetID(pdg_part->PdgCode());

      if(pdg_part->Stable()){
        part->SetStable();
      }else{
        part->SetUnstable();
      }

      part->SetDeg(Int_t(2*pdg_part->Spin()+1));
      part->SetStat(+1);

      if(part->GetID()>0){
        part->SetB(1);
      }else{part->SetB(-1);}

      part->SetS(Int_t(pdg_part->Strangeness()));
      part->SetQ(Int_t(pdg_part->Charge()/3.));
      part->SetCharm(Int_t(pdg_part->Charm()));
      part->SetBeauty(Int_t(pdg_part->Beauty()));
      part->SetTop(Int_t(pdg_part->Top()));
      part->SetMass(pdg_part->Mass());
      part->SetWidth(pdg_part->Width());
      part->SetThreshold(0.);                                 
      part->SetSContent(TMath::Abs(pdg_part->Strangeness()));         
      part->SetCContent(TMath::Abs(pdg_part->Charm()));         
      part->SetbContent(TMath::Abs(pdg_part->Beauty()));         
      part->SetTContent(TMath::Abs(pdg_part->Top()));         
      part->SetRadius(0.);
      fPartTable->Add(part);
      fParticleNumber++;

    }
  }

  fPartTable->Rehash(fPartTable->Capacity());
}

//__________________________________________________________________________
TTMParticle* TTMParticleSet::GetParticle(Int_t id)
{
  // Returns a pointer to a particle object with the specified ID if
  // that particle is in the set. If not in the set, 0 is returned.
  // First id has to be converted into a TString in exactly the same 
  // way as the particle fNames are set. 
  //  	
   
  TString name = Int_2_String(id);
  
  if (TTMParticle* temp_part = (TTMParticle*) fPartTable->FindObject(name)) {
    return temp_part;
  } else {
    cout << "WARNING: Particle " << id << " not in set" << endl;
    return 0;
  }

}

//__________________________________________________________________________
void TTMParticleSet::AddParticle(TTMParticle* part)
{
  // Add new particle to set if not already in. Corresponding anti-particles 
  // are added if they exist and are not already in the set.	
  //

  if(!GetParticle(part->GetID())){
    fPartTable->Add(part);
    fParticleNumber++;
   
    if (part->GetB() != 0 || part->GetS() != 0 || part->GetQ() != 0 || part->GetCharm() != 0 || part->GetBeauty() != 0 || part->GetTop() != 0) {
      if(!GetParticle(-part->GetID())){
        TTMParticle* apart = new TTMParticle();
        if (part->GetStable()) {
          apart->SetStable();
        } else {
          apart->SetUnstable();
        }
        TString temp;
        temp = part->GetPartName();
        temp = temp.Prepend("anti-");
        apart->SetPartName(temp);
        apart->SetID(-part->GetID());
        apart->SetDeg(part->GetDeg());
        apart->SetStat(part->GetStat());
        apart->SetMass(part->GetMass());
        apart->SetS(-part->GetS());
        apart->SetB(-part->GetB());
        apart->SetQ(-part->GetQ());
        apart->SetCharm(-part->GetCharm());
        apart->SetBeauty(-part->GetBeauty());
        apart->SetTop(-part->GetTop());
        apart->SetSContent(part->GetSContent());
        apart->SetCContent(part->GetCContent());
        apart->SetbContent(part->GetbContent());
        apart->SetTContent(part->GetTContent());
        apart->SetWidth(part->GetWidth());
        apart->SetThreshold(part->GetThreshold());

        apart->SetRadius(part->GetRadius());

        fPartTable->Add(apart);
        fParticleNumber++;
      }
    }
  }
}

//__________________________________________________________________________
void TTMParticleSet::RemoveParticle(Int_t id)
{
  // Remove particle if specified particle is in the set.
  //
	
  if (TTMParticle* part = GetParticle(id)) {
    fPartTable->Remove(part);
    fParticleNumber--;
  }
}

//__________________________________________________________________________
void TTMParticleSet::ListParticle(Int_t id)
{
  // List particle with specified ID if it is in set	
  //
	
  if (TTMParticle* part = GetParticle(id)) {
    part->List();
  }
}

//__________________________________________________________________________
void TTMParticleSet::MassCut(Double_t max)
{
  // Removes all particles with mass (in GeV) greater than max. Run this 
  // function before InputDecays() to ensure that the feeding of those 
  // massive unstable particles above the mass cut are excluded from 
  // the lighter hadrons.
  //
 
  TIter next(fPartTable);
  TTMParticle *part;
   
  while((part = (TTMParticle *) next())){
    if(part->GetMass() >= max){
      fPartTable->Remove(part);
      fParticleNumber--;
    }
  }
}

//__________________________________________________________________________
void TTMParticleSet::ListStableParticles()
{ 
  // Lists stable particles (and anti-particles) in the set.
  //
	
  TIter next(fPartTable);
  TTMParticle* part;

  cout << "\t ************** STABLE PARTICLES **************" << endl;
  while ((part = (TTMParticle*) next())) {
    if (part->GetStable()) {
      cout << "\t\t" << part->GetPartName() << endl;
    }
  }
  cout << "\t **********************************************" << endl;
}

//__________________________________________________________________________
void TTMParticleSet::InputDecays(TDatabasePDG *pdg)
{
  // The function to use when inputting decays into a TTMParticleSet object
  // instantiated with the TTMParticleSet(TDatabasePDG *) constructor. The
  // function iterates over all particles in the set, inserting the decays 
  // listed in the TParticlePDG objects in the TDatabasePDG object if the
  // particle is unstable. Before inserting a decay channel it is first 
  // checked that all the daughters are included in the set (this excludes
  // gammas and neutrinos from the decay lists). Anti-particle decays are 
  // treated in the exact same way-- they are not inserted based on the
  // particle decays.
  //
  // This method populates with TTMDecay objects rather than TTMDecayChannel 
  // objects.
  //

  TIter next(fPartTable);
  TTMParticle* part;

  while ((part = (TTMParticle*) next())) {
    if (!part->GetStable() )
      {
        TObjArray *list = pdg->GetParticle(part->GetID())->DecayList();
        TDecayChannel *ch;
        TIter next_d(list);
        part->GetDecaySummary()->Clear();
        while((ch = (TDecayChannel *)next_d())){
          Bool_t check = true; 
          for(Int_t i = 0 ; i < ch->NDaughters() ; i++){
            if(!GetParticle(ch->DaughterPdgCode(i))){
              check = false;
            }  
	    if(!check){break;}
          }
	  if(check){
            Double_t BRatio = ch->BranchingRatio();
            for(Int_t i = 0 ; i < ch->NDaughters() ; i++){
              if(part->GetDecay(ch->DaughterPdgCode(i))){
                TTMDecay *d = part->GetDecay(ch->DaughterPdgCode(i));
                d->SetBRatio(d->GetBRatio() + BRatio);
              }else{
                TTMDecay *new_decay = new TTMDecay(part->GetID(),ch->DaughterPdgCode(i),BRatio);
                part->GetDecaySummary()->Add(new_decay);
              }
            }
          }
	}
      }
  }
  GenerateBRatios();
}

//__________________________________________________________________________
void TTMParticleSet::InputDecays(TString dir, Bool_t ScaleBRatios)
{
  // Iterates through the particle table and for each unstable PARTICLE it  
  // attempts to read in the decay channels listed in the particle's decay file 
  // in the specified directory ("dir/fPartName_decay.txt"). If 
  // this file is not found, the particle is set to stable. If it is found
  // the channels are input using SetDecayChannels(char *) which generates 
  // also a summary decay list. GenerateBRatios is then called to transform 
  // this summary list into one containing just the stable daughters.
  // Once the PARTICLE has been dealt with, the decays of the corresponding 
  // ANTI-PARTICLE (if it exists) are entered. 
  //  
 
  TIter next(fPartTable);
  TTMParticle* part;

  while ((part = (TTMParticle*) next())) {
    if (!part->GetStable() && part->GetID() > 0)	//unstable PARTICLES!
      {
        TString temp;
        temp = part->GetPartName();
        temp = temp.Append("_decay.txt");
        temp = temp.Prepend("/");
        temp = temp.Prepend(dir);
        part->SetDecayChannels(temp,ScaleBRatios);
      }

    if  ((part->GetB() != 0 || part->GetS() != 0 || part->GetQ() != 0 || part->GetCharm() != 0 || part->GetBeauty() != 0 || part->GetTop() != 0) && GetParticle(-part->GetID())) {
      GetParticle(-part->GetID())->GetDecayChannels()->Delete();
      TList* part_decays = part->GetDecayChannels();
      TIter nextInList(part_decays);
      TTMDecayChannel* decay;
      while ((decay = (TTMDecayChannel*) nextInList())) {

        Double_t BR = decay->GetBRatio();
        TIter nd(decay->GetDaughterList());
        TTMIDObj *d_id;

        TList *adlist = new TList();

        while((d_id=(TTMIDObj *)nd())) {
          Int_t id = d_id->GetID();
          TTMParticle* decay_part = GetParticle(id);
          if (decay_part) {
            if (decay_part->GetB() != 0 || decay_part->GetS() != 0 ||
                decay_part->GetQ() != 0 || decay_part->GetCharm() != 0 ||
                decay_part->GetBeauty() != 0 || decay_part->GetTop() != 0) {
              TTMIDObj* aid = new TTMIDObj(-id);
              adlist->Add(aid);
            } else {
              TTMIDObj* aid = new TTMIDObj(id);
              adlist->Add(aid);
            }
          }
        }
        TTMDecayChannel *apch = new TTMDecayChannel(BR,adlist);
        GetParticle(-part->GetID())->GetDecayChannels()->Add(apch);
      }
      GetParticle(-part->GetID())->UpdateDecaySummary();
    }

  }

  GenerateBRatios();

}

//__________________________________________________________________________
void TTMParticleSet::GenerateBRatios()
{
  // Iterates over the particle table updating the summary decay lists based 
  // on the channels, before calling the recursive GenerateBRatios(part) 
  // function to determine the branching ratios to the stable particles in 
  // the set. This function must be run after any changes have been made to 
  // the decay channels in the set.  
  //

  TIter next(fPartTable);
  TTMParticle *part;
 
  while ((part = (TTMParticle*) next())) {
    part->UpdateDecaySummary();
  }

  next.Reset();
  while ((part = (TTMParticle*) next())) {
    GenerateBRatios(part);
  }

}

//__________________________________________________________________________
void TTMParticleSet::GetParents(TList* parents, Int_t id)
{
  // On return parents points to a list containing the parents of particle id.
  //

  parents->Delete();
  TIter next_p(fPartTable);
  TTMParticle* p;
  while ((p = (TTMParticle*) next_p())){
    TList *parent_decays = p->GetDecaySummary();
    TIter next_decay(parent_decays);
    TTMDecay* d;
    while ((d = (TTMDecay*) next_decay())){
      if(d->GetDaughterID()==id){parents->Add(d);}
    }
  }
}

//__________________________________________________________________________
void TTMParticleSet::ListParents(Int_t id)
{
  // Lists parents of particle id.
  //

  TList* list = new TList();
  GetParents(list,id);
  TIter next(list);
  TTMDecay *next_decay;

  while((next_decay = (TTMDecay *)next())){
    cout << GetParticle(next_decay->GetParentID())->GetPartName() << " :" 
         << next_decay->GetBRatio() << endl;           
  }

  delete list;

}

//__________________________________________________________________________
void TTMParticleSet::GenerateBRatios(TTMParticle* parent)
{
  // Private recursive function used by GenerateBRatios() to convert 
  // decay summaries to lists containing just stable particles. 
  // GenerateBRatios() updates the summaries first.
  // 

  TList* parent_decays = parent->GetDecaySummary();

  if (!parent->GetStable())    //Parent unstable
    {

      TList* temp_decays = new TList();

      TIter p_next(parent_decays);
      TTMDecay* p_decay;
      Int_t flag;

      do {
        p_next.Reset();
        flag = 0;
        temp_decays->Delete();

        while ((p_decay = (TTMDecay*) p_next())) {
          TTMParticle* daughter = GetParticle(p_decay->GetDaughterID());
          if(daughter){			//if daughter is in the set
            if (daughter->GetStable())	//if daughter is stable
              {
                TIter t_next(temp_decays);
                TTMDecay* t_decay;
                Bool_t found = false;
                while ((t_decay = (TTMDecay*) t_next()) && !found) {
                  if (t_decay->GetDaughterID() == daughter->GetID())
                    //daughter already in temp decay table of parent
                    {
                      t_decay->SetBRatio(t_decay->GetBRatio()
                                         + p_decay->GetBRatio());
                      found = true;
                    }
                }
                if (!found)
                  //daughter not yet in temp decay table of parent
                  {
                    TTMDecay *decay = new TTMDecay(p_decay->GetParentID(),
                                                   p_decay->GetDaughterID(),p_decay->GetBRatio());
                    temp_decays->AddLast(decay);
                  }
              } else              // if daughter is unstable
                {
                  TList* daughter_decays = daughter->GetDecaySummary();
                  TIter d_next(daughter_decays);
                  TTMDecay* d_decay;

                  while ((d_decay = (TTMDecay*) d_next())) {
                    TTMParticle* grandaughter = GetParticle
                      (d_decay->GetDaughterID());

                    if (grandaughter->GetStable()) 
                      //if grandaughter is stable
                      {
                        TIter t_next(temp_decays);
                        TTMDecay* t_decay;
                        Bool_t found = 0;
                        while ((t_decay = (TTMDecay*) t_next()) && !found) {
                          if (t_decay->GetDaughterID() ==
                              grandaughter->GetID())
                            //grandaughter in temp decay table of parent
                            {
                              t_decay->SetBRatio(t_decay->GetBRatio() +
                                                 p_decay->GetBRatio() *
                                                 d_decay->GetBRatio());
                              found = true;
                            }
                        }
                        if (!found)
                          //grandaughter not yet in temp decay 
                          //table of parent
                          {
                            TTMDecay* decay =
                              new TTMDecay(parent->GetID(),
                                           grandaughter->GetID(),
                                           p_decay->GetBRatio() *
                                           d_decay->GetBRatio());
                            temp_decays->AddLast(decay);
                          }
                      } else
                        //if grandaughter is unstable
                        {

                          GenerateBRatios(daughter);
                          flag = 1;  //once exited we try again
                        }
                    if (flag == 1) {
                      break;
                    }
                  }
                  if (flag == 1) {
                    break;
                  }
                }
          }else{
            cout<<"WARNING: daughter "<<p_decay->GetDaughterID()<<" not in set (omitting this contribution)"<<endl;
          }
          if (flag == 1) {
            break;
          }
        }
      } while (flag == 1); 
      parent->SetDecaySummary(temp_decays);
    }
}

//__________________________________________________________________________
TTMParticleSet& TTMParticleSet::operator=(const TTMParticleSet& obj)
{
  if (this == &obj) return *this;

  fPartTable->Delete();
  fParticleNumber = 0;

  TIter next(obj.GetParticleTable());
  TTMParticle *a;
  while((a = (TTMParticle *)next())){
    TTMParticle *part = new TTMParticle();
   
    part->SetPartName(a->GetPartName());
    part->SetID(a->GetID());
    if(a->GetStable()){
      part->SetStable();
    }else{
      part->SetUnstable();
    }
    part->SetDeg(a->GetDeg());
    part->SetStat(a->GetStat());
    part->SetB(a->GetB());
    part->SetS(a->GetS());
    part->SetQ(a->GetQ());
    part->SetCharm(a->GetCharm());
    part->SetBeauty(a->GetBeauty());
    part->SetTop(a->GetTop());
    part->SetMass(a->GetMass());
    part->SetWidth(a->GetWidth());
    part->SetThreshold(a->GetThreshold());
    part->SetSContent(a->GetSContent());
    part->SetCContent(a->GetCContent());
    part->SetbContent(a->GetbContent());
    part->SetTContent(a->GetTContent());

    part->SetRadius(a->GetRadius());

    part->SetThresholdCalc(a->GetThresholdCalc());
    part->SetThresholdFlag(a->GetThresholdFlag());

    TIter nextch(a->GetDecayChannels());
    TTMDecayChannel *old;

    while((old=(TTMDecayChannel *)nextch())){
      TIter nextid(old->GetDaughterList());
      TList *ndlist = new TList();
      TTMIDObj *oldid;

      while((oldid = (TTMIDObj *)nextid())){
        Int_t nid = oldid->GetID();
        TTMIDObj *newid = new TTMIDObj(nid);
        ndlist->Add(newid);
      }

      TTMDecayChannel *nch = new TTMDecayChannel(old->GetBRatio(),ndlist);
      part->GetDecayChannels()->Add(nch);
    }
    part->UpdateDecaySummary();
    fPartTable->Add(part);
    fParticleNumber++;
  }   
  GenerateBRatios();
  fFilename = obj.GetFilename();
  return *this;
}

//__________________________________________________________________________
void TTMParticleSet::CalculateThresholds()
{
  // Calculates thresholds based on channels in lists of particles. This 
  // function iterates through the table calling 
  // CalculateThreshold(TTMParticle *).
  //

  TIter npart(fPartTable);
  TTMParticle *part;

  while((part=(TTMParticle *)npart())){
    CalculateThreshold(part);
  }

}

//__________________________________________________________________________
void TTMParticleSet::CalculateThreshold(TTMParticle* part)
{
  // Private function that calculates the threshold of the particle pointed at
  //

  Bool_t check;
  Double_t thresh;

  do{

    check = false;
    thresh = 0.;
    TIter nch(part->GetDecayChannels());
    TTMDecayChannel *ch;

    while((ch=(TTMDecayChannel *)nch())){
      TIter nid(ch->GetDaughterList());
      TTMIDObj *id;
      Double_t mch = 0.;

      while((id = (TTMIDObj *)nid())){
        TTMParticle* d = GetParticle(id->GetID());

        if(d->GetWidth()==0.){
          mch += d->GetMass();
        }else if(d->GetThresholdFlag()){
          mch += TMath::Max(d->GetMass()-2.*d->GetWidth(),d->GetThresholdCalc());
        }else{
          check = true;
          CalculateThreshold(d);
        }

      }
      if(mch>thresh){thresh=mch;}
    }

  }while(check);
  part->SetThresholdCalc(thresh);
  part->SetThresholdFlag(true);
}

//__________________________________________________________________________
void TTMParticleSet::SetDecayEfficiency(Int_t p_id, Int_t d_id, Double_t eff)
{
  // Set the decay efficiency of all decays leading from parent (p_id) to 
  // daughter (d_id) to eff %. This is reflected only in the summary decay 
  // list of the parent from which the decay contributions are calculated 
  // (i.e. this change is not reflected in the decay channel list). Therefore, 
  // running UpdateDecaySummary will erase the change brought about by calling 
  // this method.
  //

  TTMParticle *parent = GetParticle(p_id);
  TTMDecay *decay = parent->GetDecay(d_id);

  if (decay) {
    decay->SetBRatio(decay->GetBRatio() * eff / 100.);
  }

}

//__________________________________________________________________________
void TTMParticleSet::SetRadii(Double_t radius)
{
  // BH 24/04/2014
  // Iterates over the particle table and updates all particles with radius.
  
  TIter next(fPartTable);
  TTMParticle *part;
  cout<<"WARNING: updating all particle radii with radius = "<<radius<<" (to be used for excluded volume corrections)."<<endl;
  while ((part = (TTMParticle*) next())) {
      part->SetRadius(radius);
  }
  
}

//__________________________________________________________________________
TTMParticleSet::~TTMParticleSet()
{
  // Return memory to heap
  //

  if(fPartTable){	
    fPartTable->Delete();        // deletes heap-based entries in Hash Table
    delete fPartTable;           // returns memory to heap
  }
}
