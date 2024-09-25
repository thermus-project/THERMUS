// rootlogon.C for loading THERMUS shared object libraries
// checked on 20/07/2021 with root-v6-24-02 (B.H.)

{
  TString sProgName = gProgName;
  TString gPrompt = gProgName;

  // Switch on THERMUS:
  int SwitchThermus = 1;

  TString THERMUS=gSystem->Getenv("THERMUS");
  // Start of Config aliroot
  if(SwitchThermus && ((!sProgName.CompareTo("root.exe"))||(!sProgName.CompareTo("root"))) ){
    printf(" *** Info: Root setup: \t THERMUS is included \t ***\n");
    gSystem->AddIncludePath("-I$THERMUS/particles  \
                              -I$THERMUS/main");
    gSystem->Load("$THERMUS/build/functions/libFunctions.so");
    gSystem->Load("$THERMUS/build/main/libTHERMUS.so");
    
   }
  if(!SwitchThermus && ((!sProgName.CompareTo("root.exe"))||(!sProgName.CompareTo("root"))) ){
    printf(" *** Info: Root setup: \t THERMUS is NOT included \t ***\n");
  }
   printf(" *** Info: local node: \t %s \t ***\n",gSystem->HostName());
}
