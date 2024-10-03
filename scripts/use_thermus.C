// Run this code to include Thermus in your ROOT session
{
  TString THERMUS=gSystem->Getenv("THERMUS");
  TString THERMUS_LIB=gSystem->Getenv("THERMUS_LIB");
  gSystem->AddIncludePath("-I"+THERMUS+"/includes/thermus");
  gSystem->Load(THERMUS_LIB+"/libFunctions.so");
  gSystem->Load(THERMUS_LIB+"/libTHERMUS.so");
  std::cout << " *** Info: Root setup: \t THERMUS is included ( " << THERMUS << " )\t ***\n";
}
