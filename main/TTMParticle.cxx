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
// Stores particle properties relevant to the thermal model.
// fName is determined from the Monte Carlo particle ID using the function 
// Int_2_String(Int_t x). fName is the key used to find objects in ROOT hash 
// tables (as in the TTMParticleSet class). It is important that the 
// conversion from ID to fName is consistent with that in TTMParticleSet. 
//

#include <TTMParticle.h>
#include <TTMIDObj.h>
#include <TTMDecay.h>
#include <TTMDecayChannel.h>

ClassImp(TTMParticle)

//__________________________________________________________________________
TTMParticle::TTMParticle()
{
  // Note: the id is set through the function SetID() which sets fName too.
  // 
 	
  fPartName = "";
  SetID(0);
  fDeg = fStat = fB = fS = fQ = fCharm = fBeauty = fTop = 0;
  fStable = false;
  fMass = fWidth = fThreshold = fSContent = fCContent = fbContent = fTContent = 0.;
  fThresholdCalc = 0.;
  fThresholdFlag = false;
  fRadius = 0.;
  fDecaySummary = new TList();
  fDecayChannels = new TList();
}

//__________________________________________________________________________
TTMParticle::TTMParticle(const TTMParticle &obj)
{
  fPartName = obj.GetPartName();
  SetID(obj.GetID());
  fStable = obj.GetStable();
  fDeg = obj.GetDeg();
  fStat = obj.GetStat();
  fB = obj.GetB();
  fS = obj.GetS();
  fQ = obj.GetQ();
  fCharm = obj.GetCharm();
  fBeauty = obj.GetBeauty();
  fTop = obj.GetTop();
  fMass = obj.GetMass();
  fWidth = obj.GetWidth();
  fThreshold = obj.GetThreshold();
  fRadius = obj.GetRadius();
  fThresholdCalc = obj.GetThresholdCalc();
  fThresholdFlag = obj.GetThresholdFlag();
  fSContent = obj.GetSContent();
  fCContent = obj.GetCContent();
  fbContent = obj.GetbContent();
  fTContent = obj.GetTContent();

  fDecaySummary = new TList();
  fDecayChannels = new TList();
 
  fDecaySummary->Delete();
  fDecayChannels->Delete();

  TIter nextch(obj.GetDecayChannels());
  TTMDecayChannel *oldch;

  while((oldch=(TTMDecayChannel *)nextch())){
    Double_t BRatio = oldch->GetBRatio();
    TIter nextd(oldch->GetDaughterList());
    TTMIDObj *oldid;
    TList *dlist = new TList();

    while((oldid=(TTMIDObj *)nextd())){
      TTMIDObj *newid = new TTMIDObj(oldid->GetID());
      dlist->Add(newid);
    }

    TTMDecayChannel *dch = new TTMDecayChannel(BRatio,dlist);
    fDecayChannels->Add(dch);
  }

  UpdateDecaySummary();
}

//__________________________________________________________________________
void TTMParticle::SetID(Int_t x)
	
  // Sets fName for retrieval from Hash Table (fName based on particle ID) 
  // 
{
  fID = x;
  fName = Int_2_String(x);
}

//__________________________________________________________________________
void TTMParticle::List()
{
  cout << endl << endl;
  cout << "\t ********* LISTING FOR PARTICLE " << fPartName
       << " *********" << endl;
  cout << endl << endl << "\t";

  cout << "\t ID = " << fID << endl << endl << "\t";
  cout << "\t Deg. = " << fDeg << endl << endl << "\t";
  cout << "\t STAT = " << fStat << endl << endl << "\t";
  cout << "\t Mass \t\t= ";
  cout.width(6);
  cout << fMass;
  cout << " GeV" << endl << "\t";
  cout << "\t Width \t\t= ";
  cout.width(6);
  cout << fWidth;
  cout << " GeV" << endl << "\t";
  cout << "\t Threshold \t= ";
  cout.width(6);
  cout << fThreshold;
  cout << " GeV" << endl << endl << "\t";
  if(fThresholdFlag){
    cout << "\t Calculated Threshold \t= ";
    cout.width(6);
    cout << fThresholdCalc;
    cout << " GeV" << endl << endl << "\t";
  }

  cout << "\t Hard sphere radius = " << fRadius << endl << endl << "\t";

  cout << "\t B = " << fB << endl << "\t";
  cout << "\t S = " << fS << "\t" << "\t" << "|S| = " << fSContent << endl
       << "\t";
  cout << "\t Q = " << fQ << endl << "\t";
  cout << "\t Charm = " << fCharm << "\t" << "\t" << "|C| = " << fCContent << endl 
       << "\t";
  cout << "\t Beauty = " << fBeauty << "\t" << "\t" << "|b| = " << fbContent << endl 
       << "\t" ;
  cout << "\t Top = " << fTop << "\t" << "\t" << "|T| = " << fTContent << endl 
       << "\t" ;

  if (fStable) {
    cout << "\t\t STABLE" << endl;
  } else {
    cout << "\t\t UNSTABLE" << endl << endl;

    cout << "\t\t Decay Channels:" << endl;

    TIter nextch(fDecayChannels);
    TTMDecayChannel* dch;
    while ((dch = (TTMDecayChannel*) nextch())) {
      dch->List();
    }

    cout<<endl<<endl;

    cout << "\t\t Summary of Decays: " << endl;

    TIter next(fDecaySummary);
    TTMDecay* d;
    while ((d = (TTMDecay*) next())) {
      cout << "\t\t" << d->GetDaughterID() << "\t\t"
           << d->GetBRatio() * 100 << "%" << endl;
    }

  }
  cout << "\t ********************************************** " << endl;
}

//__________________________________________________________________________
void TTMParticle::SetDecayChannels(const char *file, Bool_t ScaleBRatios)
{
  // Inputs decay channels listed in the specified file in the decay channel 
  // list of the particle (the parent).      	
  //
  // Required file format: 
  // 				BR1(%) \t d_id1(1) \t d_id2(1) ... d_idN(1)\n
  //              		BR2(%) \t d_id1(2) \t d_id2(2) ... d_idN(2)\n
  //				etc ...              		
  //
  // If the file is not found, the particle is set to stable. Once the decay 
  // channel list has been populated, fDecaySummary is populated with each 
  // daughter appearing only once using UpdateDecaySummary.
  // 
	
  fDecaySummary->Delete();
  fDecayChannels->Delete();
  ifstream data(file);

  if (!data) {
    cout << "Cannot open file: " << file << endl;
    cout << "Setting Particle " << fID << " to stable!" << endl;
    fStable = true;                                                           
  } else {
    char ich;
    Int_t nch = 0;
    Int_t no_daughters[500];   //# of daughters in each decay channel of parent
    Double_t BRatio;
    Int_t ID;

    no_daughters[0] = 0;
    data.seekg(0, ios::beg);
    while (!data.eof()) {
      data.get(ich);
      if (ich == '\t') {
        no_daughters[nch]++;
      }
      if (ich == '\n') {
        nch++;
        no_daughters[nch] = 0;	//initialise
      }
    }
    nch = nch - 1;
    data.close();

    // ********************************************************************************** //
			    // Calculate sum of branching ratios //

    ifstream calcBRatios(file);

    Double_t TotalBRatio = 0.;
    Double_t xratio;
    Double_t t;

    for (Int_t i = 1; i <= nch; i++) {
      calcBRatios >> xratio;
      TotalBRatio += xratio;

      for (Int_t j = 1; j <= no_daughters[i - 1]; j++) // loop over daughters
        {
          calcBRatios >> t;
        }
    
    }

    calcBRatios.close();

    // ********************************************************************************** //

    ifstream decays(file);

    for (Int_t i = 1; i <= nch; i++) {




      decays >> BRatio;

      if(!ScaleBRatios){
      BRatio = BRatio / 100.;	// store fraction
      }else{
      BRatio = BRatio / 100. * 100. / TotalBRatio;	// store fraction
      }

      TList *ld = new TList();

      for (Int_t j = 1; j <= no_daughters[i - 1]; j++) // loop over daughters
        {
    
          decays >> ID;

          TTMIDObj *id = new TTMIDObj(ID);
          ld->Add(id);
        }
    
      TTMDecayChannel *ch = new TTMDecayChannel(BRatio,ld);
      fDecayChannels->Add(ch);
    }

    decays.close();

  }


  UpdateDecaySummary();

}

//__________________________________________________________________________
void TTMParticle::SetDecayChannels(TList* x)
{
  // Sets decay channel list and generates fDecaySummary using 
  // UpdateDecaySummary
  // 

  if(fDecayChannels){
    fDecayChannels->Delete();
    delete fDecayChannels;
  }
    
  fDecayChannels = x;

  UpdateDecaySummary();

}

//__________________________________________________________________________
void TTMParticle::SetDecaySummary(TList* x)
{
  // Sets fDecaySummary-- allows for incompatibility between channel and 
  // summary lists.  
  //
  
  if(fDecaySummary){
    fDecaySummary->Delete();
    delete fDecaySummary;
  }

  fDecaySummary = x;

}

//__________________________________________________________________________
void TTMParticle::UpdateDecaySummary()
{
  // Updates the fDecaySummary summary list from the fDecayChannels list by 
  // iterating over these channels. Uses GetDecay().
  //

  fDecaySummary->Delete();
  
  TIter nch(fDecayChannels);
  TTMDecayChannel *ch;

  while((ch = (TTMDecayChannel *)nch())){
    Double_t BRatio = ch->GetBRatio();
    TIter nd(ch->GetDaughterList());
    TTMIDObj *daughter;

    while((daughter = (TTMIDObj *)nd())){
      Int_t d_id = daughter->GetID();
        
      if (GetDecay(d_id))   // daughter already in summary decay list of parent
        {
          TTMDecay* decay = GetDecay(d_id);
          decay->SetBRatio(decay->GetBRatio() + BRatio);
        } else              // daughter not yet in summary decay list of parent
          {
            TTMDecay* new_decay = new TTMDecay(fID, d_id, BRatio);
            fDecaySummary->Add(new_decay);
          }
    }
  }
}

//__________________________________________________________________________
void TTMParticle::SetDecayChannelEfficiency(Int_t channel, Double_t eff)
{
  // Modifies the branching ratio of the specified channel according to the 
  // given efficiency (given as a percentage).
  //
 
  TTMDecayChannel *ch = GetDecayChannel(channel);

  if(ch){
    ch->SetBRatio(ch->GetBRatio() * eff/100.);
    UpdateDecaySummary();
  }

}

//__________________________________________________________________________
TTMDecayChannel* TTMParticle::GetDecayChannel(Int_t channel)
{
  // Provides access to the specified decay channel. Returns 0 if the 
  // specified channel is not in the list.  
  //

  TIter next(fDecayChannels);
  TTMDecayChannel *ch = 0;

  for(Int_t i = 1 ; i <= channel ; i++){

    ch = (TTMDecayChannel *)next();
  }
  if(ch){
    return ch;
  }else{
    cout<<"WARNING: Channel #"<<channel<<" does not exist!"<<endl;
    return 0;
  }

}

//__________________________________________________________________________
TTMDecay* TTMParticle::GetDecay(Int_t daughter_id)
{  
  // Access to decay object of parent involving daughter with id daughter_id.	
  // Returns 0 if daughter_id is not in the list of daughter id's.
  //

  TIter next(fDecaySummary);
  TTMDecay* decay;
  Bool_t found = false;
  while ((decay = (TTMDecay*) next())) {
    if (decay->GetDaughterID() == daughter_id) {
      found = true;
      break;
    }
  }
  if (found) {
    return decay;
  } else {
    return 0;
  }
}

//__________________________________________________________________________
TTMParticle& TTMParticle::operator=(const TTMParticle& obj)
{

  if (this == &obj) return *this;

  fPartName = obj.GetPartName();
  SetID(obj.GetID());
  fStable = obj.GetStable();
  fDeg = obj.GetDeg();
  fStat = obj.GetStat();
  fB = obj.GetB();
  fS = obj.GetS();
  fQ = obj.GetQ();
  fCharm = obj.GetCharm();
  fBeauty = obj.GetBeauty();
  fTop = obj.GetTop();
  fMass = obj.GetMass();
  fWidth = obj.GetWidth();
  fThreshold = obj.GetThreshold();
  fRadius = obj.GetRadius();
  fThresholdCalc = obj.GetThresholdCalc();
  fThresholdFlag = obj.GetThresholdFlag();
  fSContent = obj.GetSContent();
  fCContent = obj.GetCContent();
  fbContent = obj.GetbContent();
  fTContent = obj.GetTContent();

  fDecaySummary->Delete();
  fDecayChannels->Delete();

  TIter nextch(obj.GetDecayChannels());
  TTMDecayChannel *oldch;

  while((oldch=(TTMDecayChannel *)nextch())){
    Double_t BRatio = oldch->GetBRatio();
    TIter nextd(oldch->GetDaughterList());
    TTMIDObj *oldid;
    TList *dlist = new TList();

    while((oldid=(TTMIDObj *)nextd())){
      TTMIDObj *newid = new TTMIDObj(oldid->GetID());
      dlist->Add(newid);
    }

    TTMDecayChannel *dch = new TTMDecayChannel(BRatio,dlist);
    fDecayChannels->Add(dch);
  }

  UpdateDecaySummary();
 
  return *this;
}

//__________________________________________________________________________
TTMParticle::~TTMParticle()
{
  // returns the memory associated with fDecaySummary and fDecayChannels 
  // to the heap
  //

  if(fDecaySummary){	 	
    fDecaySummary->Delete();           // kills all heap-based entries
    delete fDecaySummary;              // returns memory to heap
  }

  if(fDecayChannels){	 

    TIter next(fDecayChannels);
    TTMDecayChannel *ch;
    
    while((ch = (TTMDecayChannel *)next())){
      ch->GetDaughterList()->Delete();
      delete ch->GetDaughterList();
    }

    fDecayChannels->Delete();           // kills all heap-based entries
    delete fDecayChannels;              // returns memory to heap
  }

}
