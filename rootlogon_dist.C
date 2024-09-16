// rootlogon.C for loading THERMUS shared object libraries
// checked on 20/07/2021 with root-v6-24-02 (B.H.)

{
  TString sProgName = gProgName;
  TString gPrompt = gProgName;

  gPrompt += " [%d] ";    
  ((TRint*)gROOT->GetApplication())->SetPrompt( gPrompt.Data()); 
  //   gSystem->Exec("echo Welcome, user $USER, on %OS% for a local session ");  // LINUX
  //   gSystem->Exec("echo Welcome, user %USER%, on %OS% for a local session "); // WXP
  gSystem->Exec("echo Welcome, user $USER, on $OSTYPE for a local session "); // DARWIN
  {
    TDatime start;
    int idate=start.GetDate();
    int itime=start.GetTime();
    int year=idate/10000;
    int month=(idate%10000)/100;
    int day=idate%100;
    int hh=itime/10000;
    int mm=(itime%10000)/100;
    int ss=itime%100;
    TString cmonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		        "Jul","Aug","Sep","Oct","Nov","Dec"};
    printf(" *** Start at Date : %i-%s-%i Time : %i:%i:%i \t ***\n",day, (cmonth[month-1]).Data(), year, hh, mm, ss);
  }
  gROOT->SetStyle("Plain");// Default white background for all plots
  gStyle->SetCanvasColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetPadColor(10);
   
  // Settings for statistics information
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
   
  // SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24
  gStyle->SetPaperSize(20,24); 
  int font = 42; 
  gStyle->SetDrawBorder(0);
  gStyle->SetTitleFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetTextFont(font);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTextSize(0.05);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleOffset(1.0,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetMarkerSize(1.6);
  const Int_t NRGBs = 5;
  const Int_t NCont = 20;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(1); 
  gStyle->SetOptStat(1); 
  gStyle->SetOptFit(1); 
  gStyle->SetEndErrorSize(5);   

  // Delegate proxy
  printf(" *** Info: Delegate Grid proxy \t \t \t ***\n");
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
   
  // Switch on THERMUS:
  int SwitchThermus = 1;


  // Start of Config aliroot
  if(SwitchThermus && ((!sProgName.CompareTo("root.exe"))||(!sProgName.CompareTo("root"))) ){
    printf(" *** Info: Root setup: \t THERMUS is included ( %s )\t ***\n",gEnv->GetValue("THERMUS","/usr"));
    gSystem->AddIncludePath("-I$THERMUS/includes/thermus");
    gSystem->Load("$THERMUS/lib64/libFunctions.so");
    gSystem->Load("$THERMUS/lib64/libTHERMUS.so");
    
   }
  if(!SwitchThermus && ((!sProgName.CompareTo("root.exe"))||(!sProgName.CompareTo("root"))) ){
    printf(" *** Info: Root setup: \t THERMUS is NOT included \t ***\n");
  }
   printf(" *** Info: local node: \t %s \t ***\n",gSystem->HostName());
}
