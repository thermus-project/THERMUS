// Run this code to include Thermus in your ROOT session
R__ADD_INCLUDE_PATH($THERMUS/include/Thermus)
R__ADD_LIBRARY_PATH($THERMUS_LIB)
R__LOAD_LIBRARY(libFunctions)
R__LOAD_LIBRARY(libTHERMUS)

TString THERMUS = gSystem->Getenv("THERMUS");
TString THERMUS_LIB = gSystem->Getenv("THERMUS_LIB");

void use_thermus()
{
   std::cout << " *** Info: Root setup: \t THERMUS is included ( " << THERMUS << " )\t ***\n";
}
