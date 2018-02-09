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

// Author: Spencer Wheaton 14 July 2004         //

#include <TROOT.h>
#include <TObject.h>
#include <TF2.h>
#include <TMath.h>

using namespace std;

Double_t Integrate2DSimpson(TF2* func, Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t h)
{
  // Integrate 2D function from ax to bx and ay to by with step-size h using Simpson's method
  //

  Double_t sum = 0.;

  Double_t x0 = ax;
  Double_t x1 = x0 + h;
  Double_t x2 = x0 + 2.*h;

  do{

    Double_t y0 = ay;
    Double_t y1 = y0 + h;
    Double_t y2 = y0 + 2.*h;

    do{

      Double_t term = h/3.*h/3.*(func->Eval(x0,y0) + 4.*func->Eval(x0,y1) + func->Eval(x0,y2) 
                                 + 4.*func->Eval(x1,y0) + 16.*func->Eval(x1,y1) + 4.*func->Eval(x1,y2) 
                                 + func->Eval(x2,y0) + 4.*func->Eval(x2,y1) + func->Eval(x2,y2));
      sum += term;

      y0 = y2;
      y1 = y0 + h;
      y2 = y0 + 2.*h;

    }while(y2 <= by);

    x0 = x2;
    x1 = x0 + h;
    x2 = x0 + 2.*h;

  }while(x2 <= bx);

  return sum;

}

Double_t Integrate2DO6(TF2* func, Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t h)
{
  // Integrate 2D-function from ax to bx and ay to by with step-size h 
  //

  Double_t sum = 0.;

  Double_t xc = ax + h; 
  Double_t dh = TMath::Sqrt(3./5.)*h;

  do{

    Double_t yc = ay + h;

    do{
      Double_t f1 = func->Eval(xc-dh,yc-dh);
      Double_t f2 = func->Eval(xc,yc-dh);
      Double_t f3 = func->Eval(xc+dh,yc-dh);
      Double_t f4 = func->Eval(xc-dh,yc);
      Double_t f5 = func->Eval(xc,yc);
      Double_t f6 = func->Eval(xc+dh,yc);
      Double_t f7 = func->Eval(xc-dh,yc+dh);
      Double_t f8 = func->Eval(xc,yc+dh);
      Double_t f9 = func->Eval(xc+dh,yc+dh);

      Double_t term = 4. * h * h *( 16./81. * f5 + 25./324. * (f1+f3+f7+f9) + 10./81. * (f2+f4+f6+f8));

      sum += term;
    
      yc = yc + 2.*h;

    }while(yc + h <= by);

    xc = xc + 2.*h;

  }while(xc + h <= bx);

  return sum;

}

Double_t Integrate2DLaguerre32Legendre32(TF2* func, Double_t ay, Double_t by)
{
  // Integrate 2D function from 0 to infinity using Gauss-Laguerre integration
  // with 32 points and from ay to by using Gauss-Legendre integration with 32 points
  //

  Double_t xleg[32];
  Double_t wleg[32];
  Double_t x[32];
  Double_t w[32];

  xleg[0] = -0.997263861849; 
  xleg[1] = -0.985611511545;
  xleg[2] = -0.964762255588; 
  xleg[3] = -0.934906075938; 
  xleg[4] = -0.896321155766; 
  xleg[5] = -0.849367613733; 
  xleg[6] = -0.794483795968; 
  xleg[7] = -0.732182118740; 
  xleg[8] = -0.663044266930; 
  xleg[9] = -0.587715757241; 
  xleg[10] = -0.506899908932; 
  xleg[11] = -0.421351276131; 
  xleg[12] = -0.331868602282; 
  xleg[13] = -0.239287362252; 
  xleg[14] = -0.144471961583; 
  xleg[15] = -0.048307665688; 
  xleg[16] = 0.048307665688; 	
  xleg[17] = 0.144471961583; 	
  xleg[18] = 0.239287362252; 	
  xleg[19] = 0.331868602282; 	
  xleg[20] = 0.421351276131; 	
  xleg[21] = 0.506899908932; 	
  xleg[22] = 0.587715757241;	
  xleg[23] = 0.663044266930; 	
  xleg[24] = 0.732182118740; 	
  xleg[25] = 0.794483795968; 	
  xleg[26] = 0.849367613733; 	
  xleg[27] = 0.896321155766; 	
  xleg[28] = 0.934906075938; 	
  xleg[29] = 0.964762255588; 	
  xleg[30] = 0.985611511545; 	
  xleg[31] = 0.997263861849; 	

  wleg[0] = 0.007018610009;
  wleg[1] = 0.016274394716;
  wleg[2] = 0.025392065309;
  wleg[3] = 0.034273862913;
  wleg[4] = 0.042835898022;
  wleg[5] = 0.050998059262;
  wleg[6] = 0.058684093479;
  wleg[7] = 0.065822222776;
  wleg[8] = 0.072345794109;
  wleg[9] = 0.078193895787;
  wleg[10] = 0.083311924227;
  wleg[11] = 0.087652093004;
  wleg[12] = 0.091173878696;
  wleg[13] = 0.093844399081;
  wleg[14] = 0.095638720079;
  wleg[15] = 0.096540088515;
  wleg[16] = 0.096540088515;
  wleg[17] = 0.095638720079;
  wleg[18] = 0.093844399081;
  wleg[19] = 0.091173878696;
  wleg[20] = 0.087652093004;
  wleg[21] = 0.083311924227;
  wleg[22] = 0.078193895787;
  wleg[23] = 0.072345794109;
  wleg[24] = 0.065822222776;
  wleg[25] = 0.058684093479;
  wleg[26] = 0.050998059262;
  wleg[27] = 0.042835898022;
  wleg[28] = 0.034273862913;
  wleg[29] = 0.025392065309;
  wleg[30] = 0.016274394716;
  wleg[31] = 0.007018610009;

  Double_t xlag[32];
  Double_t wlag[32];

  xlag[0] = 0.044489365833; 	  	
  xlag[1] = 0.234526109520;	 	
  xlag[2] = 0.576884629302; 		
  xlag[3] = 1.072448753818; 		
  xlag[4] = 1.722408776445; 		
  xlag[5] = 2.528336706426; 		
  xlag[6] = 3.492213273022; 		
  xlag[7] = 4.616456769750; 		
  xlag[8] = 5.903958504174; 		
  xlag[9] = 7.358126733186; 		
  xlag[10] = 8.982940924213; 		
  xlag[11] = 10.783018632540; 	
  xlag[12] = 12.763697986743; 	
  xlag[13] = 14.931139755523; 	
  xlag[14] = 17.292454336715; 	
  xlag[15] = 19.855860940336; 	
  xlag[16] = 22.630889013197; 	
  xlag[17] = 25.628636022459; 	
  xlag[18] = 28.862101816323; 	
  xlag[19] = 32.346629153965; 	
  xlag[20] = 36.100494805752; 	
  xlag[21] = 40.145719771539; 	
  xlag[22] = 44.509207995755; 	
  xlag[23] = 49.224394987309; 	
  xlag[24] = 54.333721333397; 	
  xlag[25] = 59.892509162134; 	
  xlag[26] = 65.975377287935; 	
  xlag[27] = 72.687628090663;	
  xlag[28] = 80.187446977914;	
  xlag[29] = 88.735340417892; 	
  xlag[30] = 98.829542868284;	
  xlag[31] = 111.751398097938; 	

  wlag[0] = 0.114187105768;
  wlag[1] = 0.266065216898;
  wlag[2] = 0.418793137325;
  wlag[3] = 0.572532846500;
  wlag[4] = 0.727648788381;
  wlag[5] = 0.884536719340;
  wlag[6] = 1.043618875892;
  wlag[7] = 1.205349274152;
  wlag[8] = 1.370221338522;
  wlag[9] = 1.538777256469;
  wlag[10] = 1.711619352686;
  wlag[11] = 1.889424063449;
  wlag[12] = 2.072959340247;
  wlag[13] = 2.263106633997;
  wlag[14] = 2.460889072488;
  wlag[15] = 2.667508126397;
  wlag[16] = 2.884392092922;
  wlag[17] = 3.113261327040;
  wlag[18] = 3.356217692596;
  wlag[19] = 3.615869856484;
  wlag[20] = 3.895513044949;
  wlag[21] = 4.199394104712;
  wlag[22] = 4.533114978534;
  wlag[23] = 4.904270287611;
  wlag[24] = 5.323500972024;
  wlag[25] = 5.806333214234;
  wlag[26] = 6.376614674160;
  wlag[27] = 7.073526580707;
  wlag[28] = 7.967693509296;
  wlag[29] = 9.205040331278;
  wlag[30] = 11.163013090768;
  wlag[31] = 15.390180415261;

  Double_t sum = 0.;

  for(Int_t i = 0 ; i < 32 ; i++){
    for(Int_t j = 0 ; j < 32 ; j++){
      x[j] = (by-ay)/2.*xleg[j] + (by+ay)/2.;
      w[j] = (by-ay)/2.*wleg[j];

      sum += wlag[i]*w[j]*func->Eval(xlag[i],x[j]);
    }
  }
  return sum;

}

Double_t Integrate2DLegendre32(TF2* func, Double_t ax, Double_t bx, Double_t ay, Double_t by)
{
  // Integrate 2D function using Gauss-Legendre integration with 32 points
  //

  Double_t xleg[32];
  Double_t wleg[32];

  Double_t x[32];
  Double_t wx[32];

  Double_t y[32];
  Double_t wy[32];

  xleg[0] = -0.997263861849; 
  xleg[1] = -0.985611511545;
  xleg[2] = -0.964762255588; 
  xleg[3] = -0.934906075938; 
  xleg[4] = -0.896321155766; 
  xleg[5] = -0.849367613733; 
  xleg[6] = -0.794483795968; 
  xleg[7] = -0.732182118740; 
  xleg[8] = -0.663044266930; 
  xleg[9] = -0.587715757241; 
  xleg[10] = -0.506899908932; 
  xleg[11] = -0.421351276131; 
  xleg[12] = -0.331868602282; 
  xleg[13] = -0.239287362252; 
  xleg[14] = -0.144471961583; 
  xleg[15] = -0.048307665688; 
  xleg[16] = 0.048307665688; 	
  xleg[17] = 0.144471961583; 	
  xleg[18] = 0.239287362252; 	
  xleg[19] = 0.331868602282; 	
  xleg[20] = 0.421351276131; 	
  xleg[21] = 0.506899908932; 	
  xleg[22] = 0.587715757241;	
  xleg[23] = 0.663044266930; 	
  xleg[24] = 0.732182118740; 	
  xleg[25] = 0.794483795968; 	
  xleg[26] = 0.849367613733; 	
  xleg[27] = 0.896321155766; 	
  xleg[28] = 0.934906075938; 	
  xleg[29] = 0.964762255588; 	
  xleg[30] = 0.985611511545; 	
  xleg[31] = 0.997263861849; 	

  wleg[0] = 0.007018610009;
  wleg[1] = 0.016274394716;
  wleg[2] = 0.025392065309;
  wleg[3] = 0.034273862913;
  wleg[4] = 0.042835898022;
  wleg[5] = 0.050998059262;
  wleg[6] = 0.058684093479;
  wleg[7] = 0.065822222776;
  wleg[8] = 0.072345794109;
  wleg[9] = 0.078193895787;
  wleg[10] = 0.083311924227;
  wleg[11] = 0.087652093004;
  wleg[12] = 0.091173878696;
  wleg[13] = 0.093844399081;
  wleg[14] = 0.095638720079;
  wleg[15] = 0.096540088515;
  wleg[16] = 0.096540088515;
  wleg[17] = 0.095638720079;
  wleg[18] = 0.093844399081;
  wleg[19] = 0.091173878696;
  wleg[20] = 0.087652093004;
  wleg[21] = 0.083311924227;
  wleg[22] = 0.078193895787;
  wleg[23] = 0.072345794109;
  wleg[24] = 0.065822222776;
  wleg[25] = 0.058684093479;
  wleg[26] = 0.050998059262;
  wleg[27] = 0.042835898022;
  wleg[28] = 0.034273862913;
  wleg[29] = 0.025392065309;
  wleg[30] = 0.016274394716;
  wleg[31] = 0.007018610009;

  Double_t sum = 0.;

  for(Int_t i = 0 ; i < 32 ; i++){
    x[i] = (bx-ax)/2.*xleg[i] + (bx+ax)/2.;
    wx[i] = (bx-ax)/2.*wleg[i];
    for(Int_t j = 0 ; j < 32 ; j++){
      y[j] = (by-ay)/2.*xleg[j] + (by+ay)/2.;
      wy[j] = (by-ay)/2.*wleg[j];

      sum += wx[i]*wy[j]*func->Eval(x[i],y[j]);
    }
  }
  return sum;

}

Double_t Integrate2DLaguerre15Legendre20(TF2* func, Double_t ay, Double_t by)
{
  // Integrate 2D function from 0 to infinity using Gauss-Laguerre integration
  // with 15 points and from ay to by using Gauss-Legendre integration with 20 points
  //

  Double_t xleg[20];
  Double_t wleg[20];
  Double_t x[20];
  Double_t w[20];

  xleg[0] = -0.993128599185;
  xleg[1] = -0.963971927278;
  xleg[2] = -0.912234428251;
  xleg[3] = -0.839116971822;
  xleg[4] = -0.746331906460;
  xleg[5] = -0.636053680727;
  xleg[6] = -0.510867001951;
  xleg[7] = -0.373706088715;
  xleg[8] = -0.227785851142;
  xleg[9] = -0.076526521133;
  xleg[10] = 0.076526521133;
  xleg[11] = 0.227785851142;
  xleg[12] = 0.373706088715;
  xleg[13] = 0.510867001951;
  xleg[14] = 0.636053680727;
  xleg[15] = 0.746331906460;
  xleg[16] = 0.839116971822;
  xleg[17] = 0.912234428251;
  xleg[18] = 0.963971927278;
  xleg[19] = 0.993128599185;

  wleg[0] = 0.017614007139;
  wleg[1] = 0.040601429768;
  wleg[2] = 0.062672048333;
  wleg[3] = 0.083276741577;
  wleg[4] = 0.101930119817;
  wleg[5] = 0.118194531962;
  wleg[6] = 0.131688638449;
  wleg[7] = 0.142096109318;
  wleg[8] = 0.149172986473;
  wleg[9] = 0.152753387131;
  wleg[10] = 0.152753387131;
  wleg[11] = 0.149172986473;
  wleg[12] = 0.142096109318;
  wleg[13] = 0.131688638449;
  wleg[14] = 0.118194531962;
  wleg[15] = 0.101930119817;
  wleg[16] = 0.083276741577;
  wleg[17] = 0.062672048333;
  wleg[18] = 0.040601429768;
  wleg[19] = 0.017614007139;

  Double_t xlag[15];
  Double_t wlag[15];

  xlag[0] = 0.093307812017;
  xlag[1] = 0.492691740302;
  xlag[2] = 1.215595412071;
  xlag[3] = 2.269949526204;
  xlag[4] = 3.667622721751;
  xlag[5] = 5.425336627414;
  xlag[6] = 7.565916226613;
  xlag[7] = 10.120228568019;
  xlag[8] = 13.130282482176;
  xlag[9] = 16.654407708330;
  xlag[10] = 20.776478899449;
  xlag[11] = 25.623894226729;
  xlag[12] = 31.407519169754;
  xlag[13] = 38.530683306486;
  xlag[14] = 48.026085572686;

  wlag[0] = 0.239578170311;
  wlag[1] = 0.560100842793;
  wlag[2] = 0.887008262919;
  wlag[3] = 1.223664402148;
  wlag[4] = 1.574448721630;
  wlag[5] = 1.944751976530;
  wlag[6] = 2.341502056636;
  wlag[7] = 2.774041926826;
  wlag[8] = 3.255643346398;
  wlag[9] = 3.806311714226;
  wlag[10] = 4.458477753837;
  wlag[11] = 5.270017784430;
  wlag[12] = 6.359563469731;
  wlag[13] = 8.031787632117;
  wlag[14] = 11.527772100941;

  Double_t sum = 0.;

  for(Int_t i = 0 ; i < 15 ; i++){
    for(Int_t j = 0 ; j < 20 ; j++){
      x[j] = (by-ay)/2.*xleg[j] + (by+ay)/2.;
      w[j] = (by-ay)/2.*wleg[j];

      sum += wlag[i]*w[j]*func->Eval(xlag[i],x[j]);
    }
  }
  return sum;

}
