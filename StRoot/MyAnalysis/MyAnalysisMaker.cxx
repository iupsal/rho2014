#include "MyAnalysisMaker.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h" 
#include <iostream>
#include <fstream>
#include <iomanip>
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TSystem.h>


//#include <StRunInfo.h>
//#include <StEventInfo.h>
//#include "StRoot/StPicoDstMaker/StPicoHelix.h"
//#include "StThreeVector.hh"
//#include "StRoot/StPicoDstMaker/StPhysicalHelixD.h"
#include "StThreeVectorF.hh"
//#include "StBTofHeader.h"


#define NumberOfTH1F      200                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTH2F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTH3F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTProfile  100                    // Number of Histograms (must be less than max in .h file)

#define MinTracks          0                    // Discard events with fewer than MinTracks (Flow with 2 tracks makes no sense :-)
#define MinEvents          0                    // Discard DST data when we have fewer than MinEvents on a file sequence

#define pi                 TMath::Pi()

#define LambdaPdgMass      1.11568
#define ProtonPdgMass      0.938272
#define PionPdgMass        0.139570
#define LambdaPdg          3122
#define ProtonPdg          2212
#define PiMinusPdg         -211
#define AntiLambdaPdg      -3122
#define AntiProtonPdg      -2212
#define PiPlusPdg          211
#define Rho0Mass           0.77549
#define Rho0Pdg            113
#define K0PdgMass          0.497648
#define K0Pdg              311
#define PhiPdgMass         1.01946
#define PhiPdg             333
#define cbins              5


#define NMAX         10000
float DataArray[NMAX][20];
int ran_map[NMAX];


ClassImp(MyAnalysisMaker)                       // Macro for CINT compatibility


MyAnalysisMaker::MyAnalysisMaker( StPicoDstMaker* maker, TString FileNameBase, int shiftnumber, int calibmode, bool epmode) : StMaker("MyAnalysisMaker")

{ // Initialize and/or zero all public/private data members here.

  mFileNameBase = FileNameBase;
  mShiftNumber = shiftnumber;
  mCalibMode1 = calibmode;
  mEpMode = epmode;


  mTree = 0;
  mrootS = -99;
  mRunYear = -99;
  // mShiftNumber = shiftnumber;
  // mCalibMode1 = calibmode;
  // mEpMode = epmode;
  for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointerOBs
    {
      histogram1D[i] = NULL    ;
    }

  for ( Int_t i = 0 ; i < NumberOfTH2F ; i++ )  // Zero the 2D histogram pointers
    {
      histogram2D[i] = NULL  ;
    }

  for ( Int_t i = 0 ; i < NumberOfTH3F ; i++ )  // Zero the 2D histogram pointers
    {
      histogram3D[i] = NULL  ;
    }

  for ( Int_t i = 0 ; i < NumberOfTProfile ; i++ )  // Zero the Profile histogram pointers
    {
      histogramTProfile[i] = NULL  ;
      histogramTProfile2D[i] = NULL  ;
    }

  mHistogramOutput  =  NULL  ;                  // Zero the pointer to histogram output file
  mFemtoDstFile = NULL;

  mPicoDstMaker       =  maker ;                  // Pass picoDst pointer to AnlysisMaker class member functions
  mPicoDst            =  mPicoDstMaker->picoDst() ;


  mEventsStarted    =  0     ;
  mEventsProcessed  =  0     ;                  // Zero the Number of Events processed by the maker
}


MyAnalysisMaker::~MyAnalysisMaker()

{
  // Destroy and/or zero out all public/private data members here.
}


Int_t MyAnalysisMaker::Init( )

{ // Do once at the start of every analysis

  //StRefMultCorr centrality def stuff
  gSystem->Load("StRefMultCorr");
  //These params are dataset dependent, it's kind of a disaster 
  // mRefmultCorrUtil = new StRefMultCorr(); 
  // mRefmultCorrUtil->init(12126101);
  mRefmultCorrUtil = new StRefMultCorr("grefmult_P16id"); 
  mRefmultCorrUtil->init(15075008);
  mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");


 

  const struct {UInt_t bins; Double_t min; Double_t max;}
  HistEp1 = {100, 0., 2.*pi},
  HistEp2 = {100, 0., pi},
  RefMultHist = {1000,-0.5,999.5},
  CentID9Hist = {10,-1.5,8.5},
  CentID16Hist = {18,-1.5,16.5},
  // HistDeltaPhi = {1, -pi, pi},
  // HistPhiMinusPsi = {1, -pi, pi},
  // HistPt = {1, 0., 1.},
  // HistMass = {1, 0., 1.},
  // K0MassHist = {1,0.4,0.6};
  HistDeltaPhi = {24, -pi, pi},
  HistPhiMinusPsi = {24, -pi, pi},
  HistDeltaPhiMinusPsi = {24, -pi, pi},
  HistPt = {400, 0., 2.},
  HistMass = {400, 0., 2.},
  K0MassHist = {400,0.4,0.6};
  
  

  mPsi1 = new Double_t[6];

  // Create Histogram output file
  mHistogramOutput = new TFile(Form("%s.histograms.root",mFileNameBase.Data()),"recreate"); 
  // Create Histograms
  //1D histograms
  histogram1D[0]    = new TH1F ( "EventCutCounter", "EventCutCounter", 20, -.5, 19.5 ) ;
  histogram1D[1]    = new TH1F ( "RefMult", "Multiplicity (RefMult)", 1000, -0.5, 999.5 ) ;
  histogram1D[2]    = new TH1F ( "Zvtx", "Zvtx", 1000,-100, 100 ) ;
  histogram1D[3]    = new TH1F ( "VpdVz", "VpdVz", 1000,-100, 100 ) ;
  histogram1D[4]    = new TH1F ( "VzDiff", "VzDiff", 1000,-20, 20 ) ;
  histogram1D[5]    = new TH1F ( "ZdcCoincidence", "ZdcCoincidence", 1000,-0.5, 9.99995e+04) ;

  histogram1D[6]    = new TH1F ( "TrackCutCounter", "TrackCutCounter", 20, -.5, 19.5 ) ;
  histogram1D[7]    = new TH1F ( "TrackPt", "TrackPt", 1000, 0., 10. ) ;
  histogram1D[8]    = new TH1F ( "TrackDca", "TrackDca", 1000, 0., 10. ) ;
  histogram1D[9]    = new TH1F ( "TrackHitsFit", "TrackHitsFit", 100, -0.5, 99.5 ) ;
  histogram1D[10]   = new TH1F ( "TrackHitsPoss", "TrackHitsPoss", 100, -0.5, 99.5 ) ;
  histogram1D[11]   = new TH1F ( "TrackFitRatio", "TrackFitRatio", 1000, 0., 1. ) ;

  histogram1D[12]    = new TH1F ( "K0sCutCounter", "K0sCutCounter", 20, -.5, 19.5 ) ;
  histogram1D[13]    = new TH1F ( "K0sPiPlusDCA", "K0sPiPlusDCA", 1000, 0., 10. ) ;
  histogram1D[14]    = new TH1F ( "K0sPiMinusDCA", "K0sPiMinusDCA", 1000, 0., 10. ) ;
  histogram1D[15]    = new TH1F ( "K0sDaughterDca", "K0sDaughterDca", 1000, 0., 10. ) ;
  histogram1D[16]    = new TH1F ( "K0sDca", "K0sDca", 2000, 0., 20. ) ;
  histogram1D[17]    = new TH1F ( "K0sDecayLength", "LambdaDecayLength", 5000, 0., 50. ) ;

  histogram1D[18]    = new TH1F ( "ProtonCutCounter", "ProtonCutCounter", 20, -.5, 19.5 ) ;
  histogram1D[19]    = new TH1F ( "ProtonDCA", "ProtonDCA", 1000, 0., 10. ) ;
  histogram1D[20]    = new TH1F ( "APiMinusDCA", "APiMinusDCA", 1000, 0., 10. ) ;
  histogram1D[21]    = new TH1F ( "ADaughterDca", "ADaughterDca", 1000, 0., 10. ) ;
  histogram1D[22]    = new TH1F ( "AntiLambdaDca", "AntiLambdaDca", 2000, 0., 20. ) ;
  histogram1D[23]    = new TH1F ( "AntiLambdaDecayLength", "AntiLambdaDecayLength", 5000, 0., 50. ) ;

  histogram1D[24]    = new TH1F ( "K0sPt", "K0sPt", 1000, 0., 10. ) ;
  histogram1D[25]    = new TH1F ( "K0sY", "K0sY", 1000, -1., 1. ) ;
  histogram1D[26]    = new TH1F ( "ProtonPt", "ProtonPt", 1000, 0., 10. ) ;
  histogram1D[27]    = new TH1F ( "ProtonY", "ProtonY", 1000, -1., 1. ) ;

  histogram1D[28]    = new TH1F ( "K0sCounter", "K0sCounter", 20, -.5, 19.5 ) ;
  histogram1D[29]    = new TH1F ( "ProtonCounter", "ProtonCounter", 20, -.5, 19.5 ) ;
  histogram1D[30]    = new TH1F ( "K0sProtonTypeCounter", "K0sProtonTypeCounter", 20, -.5, 19.5 ) ;

  histogram1D[31]    = new TH1F ( "CorrectedRefMult", "Multiplicity (RefMult)", 1000, -0.5, 999.5 ) ;
  histogram1D[32]    = new TH1F ( "RawRefMultAfterCuts", "Multiplicity (RefMult)", 1000, -0.5, 999.5 ) ;

  histogram1D[33]    = new TH1F ( "K0sHelicity", "Counts vs #hat{p}_{p}* #bullet #hat{p}_{K0s}", 1000, -1., 1. ) ;
  histogram1D[34]    = new TH1F ( "AntiLambdaHelicity", "Counts vs #hat{p}_{#hat{p}}* #bullet #hat{p}_{#bar{#Lambda}}", 1000, -1., 1. ) ;

  histogram1D[35]    = new TH1F ( "PiPlusPt", "Counts vs pi+ {p}_{T}", 1000, 0., 5. ) ;
  histogram1D[36]    = new TH1F ( "AntiProtonPt", "Counts vs anti-proton {p}_{T}", 1000, 0., 5. ) ;

  histogram1D[37]    = new TH1F ( "Psi2Raw", "Psi2Raw", HistEp2.bins, HistEp2.min, HistEp2.max ) ;
  histogram1D[38]    = new TH1F ( "Psi2Phi", "Psi2Phi", HistEp2.bins, HistEp2.min, HistEp2.max ) ;
  histogram1D[39]    = new TH1F ( "Psi2Shift", "Psi2Shift", HistEp2.bins, HistEp2.min, HistEp2.max ) ;

  histogram1D[40]    = new TH1F ( "Psi1Raw", "Psi1Raw", HistEp1.bins, HistEp1.min, HistEp1.max ) ;
  histogram1D[41]    = new TH1F ( "Psi1Phi", "Psi1Phi", HistEp1.bins, HistEp1.min, HistEp1.max ) ;
  histogram1D[42]    = new TH1F ( "Psi1Shift", "Psi1Shift", HistEp1.bins, HistEp1.min, HistEp1.max ) ;

  histogram1D[43]    = new TH1F ( "Psi1Lt10", "Psi1Lt10", HistEp1.bins, HistEp1.min, HistEp1.max ) ;
  histogram1D[44]    = new TH1F ( "Psi1Lt5", "Psi1Lt5", HistEp1.bins, HistEp1.min, HistEp1.max ) ;






  //2D histograms
  histogram2D[0]    = new TH2F ( "Vy_vs_Vx","Vy vs Vx", 300,-3.,3.,300,-3.,3.) ;
  histogram2D[1]    = new TH2F ( "BbcAdcSumEvsW","BbcAdcSumEvsW", 500,0.,5.e4,500,0.,5.e4) ;
  histogram2D[2]    = new TH2F ( "TpcSumEvsW","TpcSumEvsW", 500,0.,1.e2,500,0.,1.e2) ;

  histogram2D[3]    = new TH2F ( "PhiVsEta","PhiVsEta", 500,0.,2.*pi,500,-1.5,1.5) ;
  histogram2D[4]    = new TH2F ( "PVsdEdx","q*p vs dEdx", 500,-2.,2.,500,1.5,10.) ;
  histogram2D[5]    = new TH2F ( "PVsTofBeta","p vs 1/beta", 500,0.,3.5,500,0.5,3.) ;
  histogram2D[6]    = new TH2F ( "PVsTofMass","p vs mass sqr", 500,0.,3.5,500,-0.2,1.2) ;

  histogram2D[7]    = new TH2F ( "PhiVsEta_AC","PhiVsEta (AC)", 500,0.,2.*pi,500,-1.5,1.5) ;
  histogram2D[8]    = new TH2F ( "PVsdEdx_AC","q*p vs dEdx (AC)", 500,-2.,2.,500,1.5,10.) ;
  histogram2D[9]    = new TH2F ( "PVsTofBeta_AC","p vs 1/beta (AC)", 500,0.,3.5,500,0.5,3.) ;
  histogram2D[10]   = new TH2F ( "PVsTofMass_AC","p vs mass sqr (AC)", 500,0.,3.5,500,-0.2,1.2) ;

  // histogram2D[7]    = new TH2F ( "K0sMass","K0sMass", 10, -0.5, 9.5, K0MassHist.bins, K0MassHist.min, K0MassHist.max) ;
  // histogram2D[8]    = new TH2F ( "AntiK0sMass","AntiK0sMass", 10, -0.5, 9.5, K0MassHist.bins, K0MassHist.min, K0MassHist.max) ;

  // histogram2D[9]    = new TH2F ( "K0sPhiVsEta","K0sPhiVsEta", 500,0.,2.*pi,500,-1.5,1.5) ;
  // histogram2D[10]   = new TH2F ( "ProtonPhiVsEta","ProtonPhiVsEta", 500,0.,2.*pi,1000,-1.5,1.5) ;

  histogram2D[11]   = new TH2F ( "MultVsCentId9","MultVsCentId9", 10, -1.5, 8.5, 1000, -0.5, 999.5) ;
  histogram2D[12]   = new TH2F ( "MultVsCentId16","MultVsCentId16", 16, -1.5, 15.5, 1000, -0.5, 999.5) ;

  histogram2D[13]   = new TH2F ( "LambdaVertexXY","LambdaVertexXY", 500, -20., 20., 500, -20., 20.) ;
  histogram2D[14]   = new TH2F ( "AntiLambdaVertexXY","AntiLambdaVertexXY", 500, -20., 20., 500, -20., 20.) ;

  histogram2D[15]   = new TH2F ( "PrimaryEtaVsPhi","PrimaryEtaVsPhi", 64,0.,2.*pi,64,-1.,1.) ;

  histogram2D[16]   = new TH2F ( "Psi2EvsWRaw","Psi2EvsWRaw", HistEp2.bins, HistEp2.min, HistEp2.max, HistEp2.bins, HistEp2.min, HistEp2.max) ;
  histogram2D[17]   = new TH2F ( "Psi2EvsWPhi","Psi2EvsWPhi", HistEp2.bins, HistEp2.min, HistEp2.max, HistEp2.bins, HistEp2.min, HistEp2.max) ;
  histogram2D[18]   = new TH2F ( "Psi2EvsWShift","Psi2EvsWShift", HistEp2.bins, HistEp2.min, HistEp2.max, HistEp2.bins, HistEp2.min, HistEp2.max) ;

  histogram2D[19]   = new TH2F ( "Psi1EvsWRaw","Psi1EvsWRaw", HistEp1.bins, HistEp1.min, HistEp1.max, HistEp1.bins, HistEp1.min, HistEp1.max) ;
  histogram2D[20]   = new TH2F ( "Psi1EvsWPhi","Psi1EvsWPhi", HistEp1.bins, HistEp1.min, HistEp1.max, HistEp1.bins, HistEp1.min, HistEp1.max) ;
  histogram2D[21]   = new TH2F ( "Psi1EvsWShift","Psi1EvsWShift", HistEp1.bins, HistEp1.min, HistEp1.max, HistEp1.bins, HistEp1.min, HistEp1.max) ;

  histogram2D[22]   = new TH2F ( "PrimaryTrackNSigmas","PrimaryTrackNSigmas", 8, 0.5, 8.5, 200, -5., 5.) ;

  histogram2D[23]   = new TH2F ( "RefMultVsgRefMult", "RefMult Vs gRefMult", 1000, -0.5, 999.5, 1000, -0.5, 999.5 ) ;
  histogram2D[24]   = new TH2F ( "RefMultVsCentId9","RefMultVsCentId9", 10, -1.5, 8.5, 1000, -0.5, 999.5) ;
  histogram2D[25]   = new TH2F ( "RefMultVsCentId16","RefMultVsCentId16", 10, -1.5, 8.5, 1000, -0.5, 999.5) ;

  histogram2D[26]   = new TH2F ( "CentIDvsPsi1","CentIDvsPsi1", HistEp1.bins, HistEp1.min, HistEp1.max, 10, -1.5, 8.5) ;

  //3D histograms
  histogram3D[0+0*cbins]   = new TH3F ( "PtVsMassVsDeltaPhi_5060","PtVsMassVsDeltaPhi_5060", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[1+0*cbins]   = new TH3F ( "PtVsMassVsDeltaPhi_6070","PtVsMassVsDeltaPhi_6070", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[2+0*cbins]   = new TH3F ( "PtVsMassVsDeltaPhi_7080","PtVsMassVsDeltaPhi_7080", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[3+0*cbins]   = new TH3F ( "PtVsMassVsDeltaPhi_80100","PtVsMassVsDeltaPhi_80100", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[4+0*cbins]   = new TH3F ( "PtVsMassVsDeltaPhi_LT5","PtVsMassVsDeltaPhi_LT5", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;

  histogram3D[0+1*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiPP_5060","PtVsMassVsDeltaPhiPP_5060", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[1+1*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiPP_6070","PtVsMassVsDeltaPhiPP_6070", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[2+1*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiPP_7080","PtVsMassVsDeltaPhiPP_7080", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[3+1*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiPP_80100","PtVsMassVsDeltaPhiPP_80100", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[4+1*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiPP_LT5","PtVsMassVsDeltaPhiPP_LT5", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;

  histogram3D[0+2*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiNN_5060","PtVsMassVsDeltaPhiNN_5060", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[1+2*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiNN_6070","PtVsMassVsDeltaPhiNN_6070", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[2+2*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiNN_7080","PtVsMassVsDeltaPhiNN_7080", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[3+2*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiNN_80100","PtVsMassVsDeltaPhiNN_80100", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;
  histogram3D[4+2*cbins]   = new TH3F ( "PtVsMassVsDeltaPhiNN_LT5","PtVsMassVsDeltaPhiNN_LT5", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistDeltaPhi.bins, HistDeltaPhi.min, HistDeltaPhi.max) ;





  histogram3D[0+3*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsi_5060","PtVsMassVsPhiMinusPsi_5060", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[1+3*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsi_6070","PtVsMassVsPhiMinusPsi_6070", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[2+3*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsi_7080","PtVsMassVsPhiMinusPsi_7080", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[3+3*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsi_80100","PtVsMassVsPhiMinusPsi_80100", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[4+3*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsi_LT5","PtVsMassVsPhiMinusPsi_LT5", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;

  histogram3D[0+4*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiPP_5060","PtVsMassVsPhiMinusPsiPP_5060", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[1+4*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiPP_6070","PtVsMassVsPhiMinusPsiPP_6070", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[2+4*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiPP_7080","PtVsMassVsPhiMinusPsiPP_7080", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[3+4*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiPP_80100","PtVsMassVsPhiMinusPsiPP_80100", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[4+4*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiPP_LT5","PtVsMassVsPhiMinusPsiPP_LT5", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;

  histogram3D[0+5*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiNN_5060","PtVsMassVsPhiMinusPsiNN_5060", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[1+5*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiNN_6070","PtVsMassVsPhiMinusPsiNN_6070", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[2+5*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiNN_7080","PtVsMassVsPhiMinusPsiNN_7080", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[3+5*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiNN_80100","PtVsMassVsPhiMinusPsiNN_80100", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;
  histogram3D[4+5*cbins]   = new TH3F ( "PtVsMassVsPhiMinusPsiNN_LT5","PtVsMassVsPhiMinusPsiNN_LT5", HistMass.bins, HistMass.min, HistMass.max, HistPt.bins, HistPt.min, HistPt.max, HistPhiMinusPsi.bins, HistPhiMinusPsi.min, HistPhiMinusPsi.max) ;










  //TProfiles
  histogramTProfile[0]  = new TProfile ( "ResCorr2", "ResCorr2", 10, -1.5, 8.5) ;
  histogramTProfile[1]  = new TProfile ( "ResCorr2Raw", "ResCorr2Raw", 10, -1.5, 8.5) ;
  histogramTProfile[2]  = new TProfile ( "ResCorr2Phi", "ResCorr2Phi", 10, -1.5, 8.5) ;

  histogramTProfile[3]  = new TProfile ( "ResCorr1", "ResCorr1", 10, -1.5, 8.5) ;
  histogramTProfile[4]  = new TProfile ( "ResCorr1Raw", "ResCorr1Raw", 10, -1.5, 8.5) ;


  //2D Profiles
  //TPC EP2 shift correction
  histogramTProfile2D[0] = new TProfile2D("ShiftTpc2SinEast","ShiftTpc2SinEast",1,0,1,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[1] = new TProfile2D("ShiftTpc2CosEast","ShiftTpc2CosEast",1,0,1,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[2] = new TProfile2D("ShiftTpc2SinWest","ShiftTpc2SinWest",1,0,1,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[3] = new TProfile2D("ShiftTpc2CosWest","ShiftTpc2CosWest",1,0,1,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[4] = new TProfile2D("ShiftTpc2Sin","ShiftTpc2Sin",1,0,1,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[5] = new TProfile2D("ShiftTpc2Cos","ShiftTpc2Cos",1,0,1,10,0,10.0,-1.0,1.0);

  //ZDC beam center correction
  histogramTProfile2D[6] = new TProfile2D("ZDCSMDBeamCenter","ZDCSMDBeamCenter", 4, 0.5, 4.5, 100, -0.5, 99.5, -20, 20, "");

  //ZDC EP1 shift correction
  histogramTProfile2D[7] = new TProfile2D("ShiftZdc1SinEast","ShiftZdc1SinEast",50,0,50,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[8] = new TProfile2D("ShiftZdc1CosEast","ShiftZdc1CosEast",50,0,50,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[9] = new TProfile2D("ShiftZdc1SinWest","ShiftZdc1SinWest",50,0,50,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[10] = new TProfile2D("ShiftZdc1CosWest","ShiftZdc1CosWest",50,0,50,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[11] = new TProfile2D("ShiftZdc1Sin","ShiftZdc1Sin",50,0,50,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[12] = new TProfile2D("ShiftZdc1Cos","ShiftZdc1Cos",50,0,50,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[13] = new TProfile2D("Q1Centrality","Q1Centrality", 4, 0.5, 4.5, 10, -1.5, 8.5, -20, 20, "");

  //TPC EP2 shift after correction
  histogramTProfile2D[14] = new TProfile2D("AfterShiftTpc2SinEast","AfterShiftTpc2SinEast",1,0,1,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[15] = new TProfile2D("AfterShiftTpc2CosEast","AfterShiftTpc2CosEast",1,0,1,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[16] = new TProfile2D("AfterShiftTpc2SinWest","AfterShiftTpc2SinWest",1,0,1,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[17] = new TProfile2D("AfterShiftTpc2CosWest","AfterShiftTpc2CosWest",1,0,1,10,0,10.0,-1.0,1.0);

  histogramTProfile2D[18] = new TProfile2D("AfterShiftTpc2Sin","AfterShiftTpc2Sin",1,0,1,10,0,10.0,-1.0,1.0);
  histogramTProfile2D[19] = new TProfile2D("AfterShiftTpc2Cos","AfterShiftTpc2Cos",1,0,1,10,0,10.0,-1.0,1.0);



  //Book the (Isaac's) femto dst code
  // mTreeEvent = new IEvent;
  // BookTree();

  SetEpCorrections();
  //Random = new TRandom() ;

  return kStOK;

}

//________________________________________________________________________________
void MyAnalysisMaker::SetEpCorrections(){
  //mCalibMode1 = 0;  //EP correction level already performed: 0, 1, 2 is Raw (no correction), centering correction, Shift correction 
  //mCalibMode2 = 0;  //EP correction level already performed: 0, 1, 2 is Raw (no correction), Phi weight correction, Shift correction 


  TString FName = "./Corrections/Shift";
  FName += (mShiftNumber);
  FName += ".root";
  cout << FName << endl;

  // ZDC-based EP1 corrections:
  if(mCalibMode1 > 0){  

    mCorrectionZdcBeamCenter = new TFile(FName,"READ"); //ZDC beam center correction
    //mCorrectionZdcBeamCenter = new TFile("./Corrections/Zdc/Subhash/BeamCenter.root","READ");
    //mCorrectionZdcBeamCenter = new TFile("./centering.root","READ");

    if(mCalibMode1 > 1){
      mCorrectionFileEp1 = new TFile(FName,"READ"); //ZDC shift center correction
      //mCorrectionFileEp1 = new TFile("./Corrections/Zdc/Subhash/shift107.root","READ");
    }
  }

  //set up ZDC EP finder code, adapted from Subhash
  //if you set mCalibMode1 too high the code works fine
  mZdcEP = new ZdcEventPlane(mCalibMode1, mCorrectionFileEp1, mCorrectionZdcBeamCenter, true); 


  //TPC-based EP2 corrections:
  // float entries2D[NumberOfCorrectionTH2F];
  // if(mCalibMode2 > 0){
  //   mCorrectionFileEp2 = new TFile(FName,"READ");
  //   histogram2DCorr[0] = (TH2F*)mCorrectionFileEp2->Get("PrimaryEtaVsPhi");
  //   entries2D[0] = histogram2DCorr[0]->GetEntries();
  //   histogram2DCorr[0]->Scale(1./entries2D[0]);

  //   if(mCalibMode2 > 1){
  //     TProfile2D* mTpc2_shiftWest_c = (TProfile2D*)mCorrectionFileEp2->Get("ShiftTpc2CosWest");
  //     TProfile2D* mTpc2_shiftWest_s = (TProfile2D*)mCorrectionFileEp2->Get("ShiftTpc2SinWest");
  //     TProfile2D* mTpc2_shiftEast_c = (TProfile2D*)mCorrectionFileEp2->Get("ShiftTpc2CosEast");
  //     TProfile2D* mTpc2_shiftEast_s = (TProfile2D*)mCorrectionFileEp2->Get("ShiftTpc2SinEast");
  //     TProfile2D* mTpc2_shiftFull_c = (TProfile2D*)mCorrectionFileEp2->Get("ShiftTpc2Cos");
  //     TProfile2D* mTpc2_shiftFull_s = (TProfile2D*)mCorrectionFileEp2->Get("ShiftTpc2Sin");

  //     for(int cent=0; cent<9; cent++){
  //       for(int i=1;i<TpcSBin+1;i++) {
  //         mTpc2ShiftEast_cos[cent][i-1] = mTpc2_shiftEast_c->GetBinContent(i, cent+2);
  //         mTpc2ShiftEast_sin[cent][i-1] = mTpc2_shiftEast_s->GetBinContent(i, cent+2);
  //         mTpc2ShiftWest_cos[cent][i-1] = mTpc2_shiftWest_c->GetBinContent(i, cent+2);
  //         mTpc2ShiftWest_sin[cent][i-1] = mTpc2_shiftWest_s->GetBinContent(i, cent+2);
  //         mTpc2ShiftFull_cos[cent][i-1] = mTpc2_shiftFull_c->GetBinContent(i, cent+2);
  //         mTpc2ShiftFull_sin[cent][i-1] = mTpc2_shiftFull_s->GetBinContent(i, cent+2);
  //       }
  //     }
  //   }
  // }


}
//_____________________________________________________________________________
Int_t MyAnalysisMaker::Make( )


{

  if(!mPicoDst) return kStOK;
  StPicoEvent* picoEvent =  mPicoDst->event()        ;  if ( !picoEvent ) return kStOK ;


  Int_t RunId = picoEvent->runId();
  mRunYear = floor( RunId/pow(10,6) );  //actual year data was taken + 1
  int RunDay = floor( (RunId - mRunYear*pow(10,6))/pow(10,3) );  
  int RunDayBin = -10;
  if(15 == mRunYear) RunDayBin = RunDay - 74;
  else if(12 == mRunYear) RunDayBin = RunDay - 122;
  else if(17 == mRunYear) RunDayBin = RunDay - 38;
  int EventId = picoEvent->eventId();
  float MagneticField = picoEvent->bField();
  mrootS = 200; // it doesn't look like I can get this from the picoDsts  
  Int_t RefMult = picoEvent->refMult() ;
  Int_t RefMultPos = picoEvent->refMultPos() ;
  Int_t RefMultNeg = picoEvent->refMultNeg() ;

  TVector3 PrimaryVertex = picoEvent->primaryVertex();
  Double_t Vx = PrimaryVertex.X(); 
  Double_t Vy = PrimaryVertex.Y(); 
  Double_t Vz = PrimaryVertex.Z(); 
  mVpdVz = -999;
  if(mrootS > 25){
    mVpdVz = picoEvent->vzVpd();
  }     
  Double_t VzDiff = Vz-mVpdVz;
  Double_t ZdcCoincidenceRate = 0.;
  if(mrootS > 60) ZdcCoincidenceRate = picoEvent->ZDCx();

  histogram1D[1]->Fill(picoEvent->grefMult());
  histogram1D[2]->Fill(Vz);
  histogram1D[3]->Fill(mVpdVz);
  histogram1D[4]->Fill(VzDiff);
  histogram1D[5]->Fill(ZdcCoincidenceRate); 
  histogram2D[0]->Fill( Vx, Vy ) ;

  mBbcSumEast = 0.; mBbcSumWest = 0.; //sums don't seem automatic in the PicoDsts
  for(int ibbc = 0; ibbc < 24; ibbc++){
    mBbcSumEast += picoEvent->bbcAdcEast(ibbc);
    mBbcSumWest += picoEvent->bbcAdcWest(ibbc);
  }

  histogram2D[1]->Fill(mBbcSumEast, mBbcSumWest);

  mEventsStarted++ ;


  histogram1D[0]->Fill(0);
  //if(!RejectRunNumbers(picoEvent)) return kStOk;
  histogram1D[0]->Fill(1);
  if ( !AcceptTrigger(picoEvent) ) return kStOK ;                  // Skip this event if it doens't pass the trigger cuts
  histogram1D[0]->Fill(2);
  if ( !AcceptEvent(picoEvent) )   return kStOK ;                  // Skip this event if it doens't pass the 'event' cuts
  histogram1D[0]->Fill(10);
  if (!EtaSymmetryCut())  return kStOK ;   // Skip this event 'etaSymetry' cuts
  histogram1D[0]->Fill(11);


  Int_t CentralityID9 = -99;
  Int_t CentralityID16 = -99;
  Float_t ReweightGrefmult = -99.;
  Float_t gRefMultCorr = -99.;
  Int_t gRefMult=picoEvent->grefMult(); 

  mRefmultCorrUtil->init((Int_t)picoEvent->runId());


  if(200 == mrootS && 15 == mRunYear){
    
    //mRefmultCorrUtil->init(RunId);  
    mRefmultCorrUtil->initEvent(gRefMult,Vz,ZdcCoincidenceRate);
    //float refMultCorr = mRefmultCorrUtil->getRefMultCorr() ;  
    //
    //Doing-8-bin analysis
    CentralityID9 = mRefmultCorrUtil->getCentralityBin9(); 
    CentralityID16 = mRefmultCorrUtil->getCentralityBin16();
    //Reject for bad centrality and Run number 
    //if(mRefMultCorr->isBadRun(picoEvent->runId())) continue;

    ReweightGrefmult = mRefmultCorrUtil->getWeight() ;
    gRefMultCorr = mRefmultCorrUtil->getRefMultCorr() ;
  }
  if(200 == mrootS && 17 == mRunYear){

    mRefmultCorrUtil->initEvent(gRefMult,Vz,ZdcCoincidenceRate);
    ReweightGrefmult = mRefmultCorrUtil->getWeight() ;
    gRefMultCorr = mRefmultCorrUtil->getRefMultCorr() ;
    CentralityID9 = mRefmultCorrUtil->getCentralityBin9(); 
    CentralityID16 = mRefmultCorrUtil->getCentralityBin16();
  }
  else{
    mRefmultCorrUtil->initEvent(RefMult, Vz, ZdcCoincidenceRate);
    CentralityID9 = mRefmultCorrUtil->getCentralityBin9();
    CentralityID16 = mRefmultCorrUtil->getCentralityBin16();

    ReweightGrefmult = mRefmultCorrUtil->getWeight() ;
    gRefMultCorr = mRefmultCorrUtil->getRefMultCorr() ;
  }

  
  //if (CentralityID9 < 0) return kStOK;
  histogram1D[0]->Fill(12);

  bool badRefMultCorr = mRefmultCorrUtil->isBadRun((Int_t)picoEvent->runId());
  if (badRefMultCorr) return kStOK;
  histogram1D[0]->Fill(13);
  

  histogram1D[31]->Fill(gRefMultCorr);
  histogram1D[32]->Fill(RefMult);
  histogram2D[11]->Fill(CentralityID9, gRefMultCorr);
  histogram2D[12]->Fill(CentralityID16, gRefMultCorr);

  histogram2D[23]->Fill(RefMult, gRefMultCorr);
  histogram2D[24]->Fill(CentralityID9, RefMult);
  histogram2D[25]->Fill(CentralityID16, RefMult);

  if(!CentralityCut(CentralityID9, gRefMultCorr)) return kStOK;
  histogram1D[0]->Fill(12);

  //cout << "RefMult " << RefMult << endl;
  //cout << "gRefMultCorr " << gRefMultCorr << endl;


  //cout << "CentralityID9 " << CentralityID9 << endl;
  //cout << "gRefMultCorr " << gRefMultCorr << endl;
  float Psi1ERaw(0.), Psi1WRaw(0.), Psi1Raw(0.), Psi1E(0.), Psi1W(0.), Psi1(0.);
  for(int i=0;i<6;i++) mPsi1[i] = kFALSE;  
  mZdcEP->AnalyzeZdcEvent(picoEvent, CentralityID9);
  mPsi1 = mZdcEP->GetPsi(); 
  Double_t* GetQ = mZdcEP->GetQCenter(); 
  Double_t* GetQCorr = mZdcEP->GetQCorrected(); 
  bool ShiftIsOkay = true;
  for(int i =0; i<3; i++){if (mPsi1[i] == kFALSE)  return kStOK;}
  for(int i =0; i<6; i++){if (mPsi1[i] == kFALSE)  ShiftIsOkay = false;}
  histogram1D[0]->Fill(12);
  if(ShiftIsOkay == true){
    Psi1E = mPsi1[3];
    Psi1W = mPsi1[4];   
    Psi1 = mPsi1[5];
  }
  if(!ShiftIsOkay && mCalibMode1==2) return kStOK;
  histogram1D[0]->Fill(13);
  Psi1ERaw = mPsi1[0];
  Psi1WRaw = mPsi1[1];
  Psi1Raw = mPsi1[2];
  // cout << "Psi1ERaw " << Psi1ERaw << ", Psi1WRaw " << Psi1WRaw << endl;

  // cout << "GetQ[0] " << GetQ[0] << endl;
  // zdc beam center correction
  for(int i=0; i<4; i++){
    histogramTProfile2D[6]->Fill(i+1, RunDayBin, GetQ[i]);
    histogramTProfile2D[13]->Fill(i+1, CentralityID9, GetQCorr[i]);
  }

  histogram1D[40]->Fill(Psi1Raw);
  histogram1D[42]->Fill(Psi1);

  if(gRefMultCorr <= 5) histogram1D[44]->Fill(Psi1);

  histogram2D[19]->Fill(Psi1WRaw, Psi1ERaw);
  //histogram2D[20]->Fill(Psi1WRaw, Psi1ERaw);
  histogram2D[21]->Fill(Psi1W, Psi1E);

  float EpCorr1 = -cos(Psi1E-Psi1W);
  float EpCorr1Raw = -cos(Psi1ERaw-Psi1WRaw);

  histogramTProfile[3]->Fill(CentralityID9, EpCorr1) ;
  histogramTProfile[4]->Fill(CentralityID9, EpCorr1Raw) ;

  //fill shift correction histograms for EP1
  if(mCalibMode1 > 0){
    for (int nfill=1;nfill<=50;nfill++){      
      histogramTProfile2D[7]->Fill(nfill-1, CentralityID9+1, sin(nfill*Psi1ERaw));
      histogramTProfile2D[8]->Fill(nfill-1, CentralityID9+1, cos(nfill*Psi1ERaw));
      
      histogramTProfile2D[9]->Fill(nfill-1, CentralityID9+1, sin(nfill*Psi1WRaw));
      histogramTProfile2D[10]->Fill(nfill-1, CentralityID9+1, cos(nfill*Psi1WRaw));
      
      histogramTProfile2D[11]->Fill(nfill-1, CentralityID9+1, sin(nfill*Psi1Raw));
      histogramTProfile2D[12]->Fill(nfill-1, CentralityID9+1, cos(nfill*Psi1Raw));
    }
  }
  //okay, reset useful EP1 to maximally corrected EP
  if(mCalibMode1 < 2) Psi1 = Psi1Raw;
  histogram2D[26]->Fill(Psi1, CentralityID9);
  if(1 == mEpMode) return kStOK; //use this for only psi1 corrections 




  TObjArray* PiPlusTrackArray = new TObjArray();
  TObjArray* PiMinusTrackArray = new TObjArray();

  //find max global track index
  int maxGBTrackIndex = -1;
  for(unsigned int iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++) 
  {
    StPicoTrack *gTrack = mPicoDst->track(iTrack);
    if (! gTrack) continue;
    int index = gTrack->id();
    if(index > maxGBTrackIndex)
      maxGBTrackIndex = index;

    histogram1D[6]->Fill(14);
    if(!(gTrack->isPrimary())) continue;
    //track quantities
    TVector3 trackMomentum = gTrack->pMom();
    float trackNSigmaPion = gTrack->nSigmaPion();
    float trackCharge= gTrack->charge();
    float trackP = trackMomentum.Mag();
    float trackBeta = -999.;
    int trackTofIndex = gTrack->bTofPidTraitsIndex();
    if(trackTofIndex >= 0) trackBeta=mPicoDst->btofPidTraits(trackTofIndex)->btofBeta();
    else trackBeta = -999; 
    float trackdEdx=gTrack->dEdx();
    float trackDCA = gTrack->gDCA(Vx,Vy,Vz);
    float MassSquare = 0;
    if(trackBeta > 0.){
      MassSquare = trackP*trackP*(1.0/(trackBeta*trackBeta)-1.0);
    }
    float trackPhi = trackMomentum.Phi();
    if(trackPhi < 0.) trackPhi += 2.*pi;

    histogram2D[3]->Fill(trackPhi, trackMomentum.Eta());
    histogram2D[4]->Fill(trackCharge*trackP, trackdEdx);
    histogram2D[5]->Fill(trackP, 1./trackBeta);
    histogram2D[6]->Fill(trackP, MassSquare);
    histogram1D[7]->Fill(trackMomentum.Pt());
    histogram1D[8]->Fill(trackDCA);
    
    histogram1D[6]->Fill(15);
    if(trackDCA > 1.) continue;

    //track cuts
    histogram1D[6]->Fill(16);
    if(!AcceptTrack(gTrack)) continue;
    histogram1D[6]->Fill(17);
    if(!PionCheck(trackNSigmaPion)) continue;
    histogram1D[6]->Fill(18);
    if(!PionMCheck(trackBeta,trackP)) continue;
    histogram1D[6]->Fill(19);
    histogram2D[7]->Fill(trackPhi, trackMomentum.Eta());
    histogram2D[8]->Fill(trackCharge*trackP, trackdEdx);
    histogram2D[9]->Fill(trackP, 1./trackBeta);
    histogram2D[10]->Fill(trackP, MassSquare);
    //if(fabs(-999. - trackBeta) < 0.1) cout << "track beta is bad?!" << endl;

    if(trackCharge > 0) PiPlusTrackArray->Add(gTrack);
    else if(trackCharge < 0) PiMinusTrackArray->Add(gTrack);
  }  //end global track loop

  histogram1D[0]->Fill(14);

  int centralityBin = -1;
  float centralityWeight = 0.;





  //_________________________________________________________________
  //Set up an IEvent that contains information relavent to making a tree
  //write in event quantities - the lambda data are added later.
  // TVector3 PrimVTXTV3(Vx,Vy,Vz);

  // mTreeEvent->SetRunNumber(RunId);
  // mTreeEvent->SetRefMult(gRefMultCorr);
  // mTreeEvent->SetCentralityID9(CentralityID9);
  // mTreeEvent->SetCentralityID16(CentralityID16);
  // mTreeEvent->SetBfield(MagneticField);
  // mTreeEvent->SetVpdVz(mVpdVz);
  // mTreeEvent->SetPrimaryVertex(PrimVTXTV3);
  // mTreeEvent->SetEp1(Psi1);
  //look at the end where Psi2 is added
  //_____________________________________________________________  



  

  int eventId = -1;
  int runId = -1;
  

  centralityWeight = 1;
  Int_t nevent = 100000;




  //make a map from trackID to track index in global track array
  //Do not do this step unless you want some information which is not in the KFParticle
  vector<int> trackMap;
  trackMap.resize(maxGBTrackIndex+1, -1);
  TLorentzVector Rho04Momentum, PiPlus4Momentum, PiMinus4Momentum;
  TLorentzVector PvPiPlus4Momentum, PvPiMinus4Momentum;
  TVector3 Rho03Momentum;
  StPicoPhysicalHelix PiPlusHelix, PiMinusHelix;
  float PiPlusDCA(-999.), PiMinusDCA(-999);


  TVector3 PiPlus3Momentum, PiMinus3Momentum;
  TVector3 PiPlus3Momentum1, PiMinus3Momentum1;
  TVector3 PiPlus3Momentum2, PiMinus3Momentum2;
  TLorentzVector PiPlus4Momentum1, PiMinus4Momentum1;
  TLorentzVector PiPlus4Momentum2, PiMinus4Momentum2;


  Float_t RhoDeltaPhi(0.), RhoPt(0.), RhoInvMass(0.), RhoDeltaPhiMinusPsi(0.), RhoPhiMinusPsi(0.);

  int K0sCounter = 0;


  //code for manually looping over global tracks.
  //pi+ and pi-
  for(int pimiter=0 ; pimiter < PiMinusTrackArray->GetEntriesFast(); pimiter++){

    StPicoTrack* PiMinusTrack = ((StPicoTrack*)PiMinusTrackArray->At(pimiter));  
    //PiMinus3Momentum = PiMinusTrack->gMom();
    PiMinus3Momentum = PiMinusTrack->pMom();
    PiMinus4Momentum.SetXYZM(PiMinus3Momentum.X(), PiMinus3Momentum.Y(), PiMinus3Momentum.Z(), PionPdgMass);

    for(int pipiter=0 ; pipiter < PiPlusTrackArray->GetEntriesFast(); pipiter++){

      StPicoTrack* PiPlusTrack  = ((StPicoTrack*)PiPlusTrackArray->At(pipiter));    //second tracks are pion tracks
      //PiPlus3Momentum = PiPlusTrack->gMom();
      PiPlus3Momentum = PiPlusTrack->pMom();
      PiPlus4Momentum.SetXYZM(PiPlus3Momentum.X(), PiPlus3Momentum.Y(), PiPlus3Momentum.Z(), PionPdgMass);

      Rho04Momentum = PiMinus4Momentum + PiPlus4Momentum;

      RhoInvMass = Rho04Momentum.M();
      RhoPt = Rho04Momentum.Pt();
      RhoDeltaPhi = GetDeltaPhi(PiMinus4Momentum.Vect(), PiPlus4Momentum.Vect());

      RhoDeltaPhiMinusPsi = RhoDeltaPhi - Psi1;
      if(RhoDeltaPhiMinusPsi < -pi) RhoDeltaPhiMinusPsi += 2.*pi;
      else if(RhoDeltaPhiMinusPsi > pi) RhoDeltaPhiMinusPsi -= 2.*pi;

      RhoPhiMinusPsi = Rho04Momentum.Phi() - Psi1;
      if(RhoPhiMinusPsi < -pi) RhoPhiMinusPsi += 2.*pi;
      else if(RhoPhiMinusPsi > pi) RhoPhiMinusPsi -= 2.*pi;


      if(2 == CentralityID9) histogram3D[0+0*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(1 == CentralityID9) histogram3D[1+0*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(0 == CentralityID9) histogram3D[2+0*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(-1 == CentralityID9) histogram3D[3+0*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi);
      if(gRefMultCorr <= 5) histogram3D[4+0*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 

      if(2 == CentralityID9) histogram3D[0+3*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(1 == CentralityID9) histogram3D[1+3*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(0 == CentralityID9) histogram3D[2+3*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(-1 == CentralityID9) histogram3D[3+3*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(gRefMultCorr <= 5) histogram3D[4+3*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 

    }
  }


  //pi+ and pi+
  for(int pipiter1=0 ; pipiter1 < PiPlusTrackArray->GetEntriesFast()-1; pipiter1++){

    StPicoTrack* PiPlusTrack1 = ((StPicoTrack*)PiPlusTrackArray->At(pipiter1));  
    //PiPlus3Momentum1 = PiPlusTrack1->gMom();
    PiPlus3Momentum1 = PiPlusTrack1->pMom();
    PiPlus4Momentum1.SetXYZM(PiPlus3Momentum1.X(), PiPlus3Momentum1.Y(), PiPlus3Momentum1.Z(), PionPdgMass);

    for(int pipiter2=pipiter1+1 ; pipiter2 < PiPlusTrackArray->GetEntriesFast(); pipiter2++){

      StPicoTrack* PiPlusTrack2  = ((StPicoTrack*)PiPlusTrackArray->At(pipiter2));    //second tracks are pion tracks
      //PiPlus3Momentum2 = PiPlusTrack2->gMom();
      PiPlus3Momentum2 = PiPlusTrack2->pMom();
      PiPlus4Momentum2.SetXYZM(PiPlus3Momentum2.X(), PiPlus3Momentum2.Y(), PiPlus3Momentum2.Z(), PionPdgMass);

      Rho04Momentum = PiPlus4Momentum1 + PiPlus4Momentum2;

      RhoInvMass = Rho04Momentum.M();
      RhoPt = Rho04Momentum.Pt();
      RhoDeltaPhi = GetDeltaPhi(PiPlus4Momentum1.Vect(), PiPlus4Momentum2.Vect());

      RhoDeltaPhiMinusPsi = RhoDeltaPhi - Psi1;
      if(RhoDeltaPhiMinusPsi < -pi) RhoDeltaPhiMinusPsi += 2.*pi;
      else if(RhoDeltaPhiMinusPsi > pi) RhoDeltaPhiMinusPsi -= 2.*pi;

      RhoPhiMinusPsi = Rho04Momentum.Phi() - Psi1;
      if(RhoPhiMinusPsi < -pi) RhoPhiMinusPsi += 2.*pi;
      else if(RhoPhiMinusPsi > pi) RhoPhiMinusPsi -= 2.*pi;


      if(2 == CentralityID9) histogram3D[0+1*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(1 == CentralityID9) histogram3D[1+1*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(0 == CentralityID9) histogram3D[2+1*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(-1 == CentralityID9) histogram3D[3+1*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi);
      if(gRefMultCorr <= 5) histogram3D[4+1*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 

      if(2 == CentralityID9) histogram3D[0+4*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(1 == CentralityID9) histogram3D[1+4*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(0 == CentralityID9) histogram3D[2+4*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi);
      if(-1 == CentralityID9) histogram3D[3+4*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(gRefMultCorr <= 5) histogram3D[4+4*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 

    }
  }


  //pi- and pi-
  for(int pimiter1=0 ; pimiter1 < PiMinusTrackArray->GetEntriesFast()-1; pimiter1++){

    StPicoTrack* PiMinusTrack1 = ((StPicoTrack*)PiMinusTrackArray->At(pimiter1));  
    //PiMinus3Momentum1 = PiMinusTrack1->gMom();
    PiMinus3Momentum1 = PiMinusTrack1->pMom();
    PiMinus4Momentum1.SetXYZM(PiMinus3Momentum1.X(), PiMinus3Momentum1.Y(), PiMinus3Momentum1.Z(), PionPdgMass);

    for(int pimiter2=pimiter1+1 ; pimiter2 < PiMinusTrackArray->GetEntriesFast(); pimiter2++){

      StPicoTrack* PiMinusTrack2  = ((StPicoTrack*)PiMinusTrackArray->At(pimiter2));    //second tracks are pion tracks
      //PiMinus3Momentum2 = PiMinusTrack2->gMom();
      PiMinus3Momentum2 = PiMinusTrack2->pMom();
      PiMinus4Momentum2.SetXYZM(PiMinus3Momentum2.X(), PiMinus3Momentum2.Y(), PiMinus3Momentum2.Z(), PionPdgMass);

      Rho04Momentum = PiMinus4Momentum1 + PiMinus4Momentum2;

      RhoInvMass = Rho04Momentum.M();
      RhoPt = Rho04Momentum.Pt();
      RhoDeltaPhi = GetDeltaPhi(PiMinus4Momentum1.Vect(), PiMinus4Momentum2.Vect());


      RhoDeltaPhiMinusPsi = RhoDeltaPhi - Psi1;
      if(RhoDeltaPhiMinusPsi < -pi) RhoDeltaPhiMinusPsi += 2.*pi;
      else if(RhoDeltaPhiMinusPsi > pi) RhoDeltaPhiMinusPsi -= 2.*pi;

      RhoPhiMinusPsi = Rho04Momentum.Phi() - Psi1;
      if(RhoPhiMinusPsi < -pi) RhoPhiMinusPsi += 2.*pi;
      else if(RhoPhiMinusPsi > pi) RhoPhiMinusPsi -= 2.*pi;


      if(2 == CentralityID9) histogram3D[0+2*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(1 == CentralityID9) histogram3D[1+2*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(0 == CentralityID9) histogram3D[2+2*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(-1 == CentralityID9) histogram3D[3+2*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 
      if(gRefMultCorr <= 5) histogram3D[4+2*cbins]->Fill(RhoInvMass, RhoPt, RhoDeltaPhi); 

      if(2 == CentralityID9) histogram3D[0+5*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(1 == CentralityID9) histogram3D[1+5*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(0 == CentralityID9) histogram3D[2+5*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(-1 == CentralityID9) histogram3D[3+5*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 
      if(gRefMultCorr <= 5) histogram3D[4+5*cbins]->Fill(RhoInvMass, RhoPt, RhoPhiMinusPsi); 

    }
  }






  mEventsProcessed++;
  return kStOK;
}
//------------------------------------------------------------
Int_t MyAnalysisMaker::Finish(){
  // Do final update of histograms
  // Write histograms to disk, output miscellaneous other information
  cout << "Finishing" << endl;
  mHistogramOutput->Write();
  mHistogramOutput->Close();

  cout << endl << endl << "Events Started " << mEventsStarted << "  Events Fully Processed " << mEventsProcessed << endl << endl << endl ;
  return kStOK ;
}
//------------------------------------------------------------
bool MyAnalysisMaker::AcceptTrigger(StPicoEvent* picoEvent){
  
  //2011
  if (12 == mRunYear && 200 == mrootS && !(picoEvent->isTrigger(350003) || picoEvent->isTrigger(350023) || picoEvent->isTrigger(350033) || picoEvent->isTrigger(350043) ))  return false ; 
  
  //2014
  if (15 == mRunYear && 200 == mrootS && !(picoEvent->isTrigger(450050) || picoEvent->isTrigger(450060) || picoEvent->isTrigger(450005) || picoEvent->isTrigger(450015) || picoEvent->isTrigger(450025) ))  return false ; 
  //if ( ! ( picoEvent->triggerIdCollection().nominal().isTrigger(280001) ))  return false ;// ... this is from http://www.star.bnl.gov/protected/bulkcorr/hmasui/2010/centrality_39GeV_Aug09/hiroshi_centrality_AuAu39GeV_Run10_Aug10_2010.pdf

  //2016
  if (17 == mRunYear && 200 == mrootS && !(picoEvent->isTrigger(520001) || picoEvent->isTrigger(520011) || picoEvent->isTrigger(520021) || picoEvent->isTrigger(520031) || picoEvent->isTrigger(520041) || picoEvent->isTrigger(520051) ))  return false ; 
  
  return true ;
}
//----------------------------------------------------------------
bool MyAnalysisMaker::RejectRunNumbers(StPicoEvent* picoEvent) //cut to reject runs with "bad run lists" from Alex/Hiroshi
{
//Badrunlist from the hist
//Bad runs for 2011, Au+Au 19 GeV Data from from Hiroshi
//
// MUST END THE ARRAYS WITH ZERO!! - this is how they're run over
//

const Int_t null_list[0+1] = {0};
const Int_t bad_run_list_7GeV[328+1]  = {11114084,11114085,11114086,11114088,11114089,11114094,11114095,11114100,11114109,11115005,11115007,11115013,11115019,11115025,11115027,11115028,11115030,11115032,11115051,11115062,11115064,11115069,11115072,11115078,11115079,11115080,11115086,11115088,11115094,11116001,11116002,11116005,11116006,11116010,11116014,11116020,11116023,11116028,11116060,11116061,11116062,11116064,11116068,11116070,11116072,11116073,11116075,11117002,11117006,11117031,11117033,11117034,11117036,11117039,11117044,11117045,11117046,11117052,11117055,11117063,11117064,11117071,11117075,11117085,11117088,11117089,11117090,11117093,11117094,11117095,11117098,11117100,11117103,11117104,11117107,11118007,11118008,11118016,11118024,11118025,11118026,11118039,11118044,11119001,11119003,11119006,11119007,11119009,11119012,11119013,11119015,11119016,11119017,11119022,11119024,11119026,11119029,11119030,11119056,11119057,11119060,11119062,11119067,11119069,11119070,11119071,11119074,11119075,11119077,11119079,11119081,11119090,11119091,11119100,11119101,11120003,11120006,11120008,11120011,11120014,11120019,11120023,11120030,11120034,11120037,11120039,11120040,11120045,11120052,11120057,11120062,11120063,11120069,11120070,11120071,11120074,11120077,11120078,11120084,11120092,11121006,11121014,11121015,11121019,11121029,11121030,11121034,11121035,11121043,11121044,11121054,11121058,11121066,11121067,11121070,11121075,11121082,11122001,11122007,11122008,11122010,11122017,11122024,11122037,11122038,11122047,11122048,11122049,11122050,11122053,11122058,11122062,11122069,11122073,11122078,11122085,11122097,11123003,11123004,11123015,11123026,11123028,11123040,11123044,11123055,11123057,11123058,11123059,11123067,11123075,11123076,11123077,11123079,11123081,11123084,11123086,11123088,11123089,11123093,11123094,11123095,11123100,11123101,11123102,11123104,11124001,11124005,11124007,11124008,11124015,11124016,11124018,11124041,11124046,11124050,11124051,11124052,11124053,11124058,11124060,11124061,11124062,11124063,11124064,11124065,11124066,11124069,11125002,11125003,11125004,11125005,11125006,11125008,11125012,11125013,11125014,11125015,11125016,11125017,11125020,11125021,11125022,11125023,11125073,11125081,11125089,11125090,11125096,11125097,11126005,11126006,11126007,11126016,11126018,11126022,11126023,11127001,11127002,11127043,11128005,11128012,11128018,11128050,11128056,11128072,11129018,11129022,11129028,11129051,11130027,11130034,11130057,11131038,11131062,11132013,11132070,11133006,11133019,11134053,11134060,11134067,11134076,11135068,11136003,11136005,11136006,11136007,11136008,11136012,11136013,11136014,11136061,11136076,11136101,11136130,11136160,11136163,11137019,11138027,11138049,11138086,11138124,11139014,11140076,11140086,11141063,11142117,11143026,11143028,11144001,11144009,11144031,11144033,11144040,11144043,11144052,11145008,11145028,11145035,11146061,11146076,11146079,11147004,11147006,11147014,11147017,11147021,11147023, 0};
const Int_t bad_run_list_11GeV[27+1]  = {11148039,11148045,11149001,11149008,11149010,11149011,11149015,11149047,11150016,11150025,11150028,11151036,11151040,11151050,11152016,11152036,11152078,11153032,11153042,11155001,11155009,11156003,11156009,11157012,11158006,11158022,11158024, 0};
const Int_t bad_run_list_14GeV[253] = {15046073,15046089,15046094,15046096,15046102,15046103,15046104,15046105,15046106,15046107,15046108,15046109,15046110,15046111,15047004,15047015,15047016,15047019,15047021,15047023,15047024,15047026,15047027,15047028,15047029,15047030,15047039,15047040,15047041,15047044,15047047,15047050,15047052,15047053,15047056,15047057,15047061,15047062,15047063,15047064,15047065,15047068,15047069,15047070,15047071,15047072,15047074,15047075,15047082,15047085,15047086,15047087,15047093,15047096,15047097,15047098,15047100,15047102,15047104,15047106,15048003,15048004,15048012,15048013,15048014,15048016,15048017,15048018,15048019,15048020,15048021,15048023,15048024,15048025,15048026,15048028,15048029,15048030,15048031,15048033,15048034,15048074,15048075,15048076,15048077,15048078,15048079,15048080,15048081,15048082,15048083,15048084,15048085,15048086,15048087,15048088,15048089,15048091,15048092,15048093,15048094,15048095,15048096,15048097,15048098,15049002,15049003,15049009,15049013,15049014,15049015,15049016,15049017,15049018,15049019,15049020,15049021,15049022,15049023,15049025,15049026,15049027,15049028,15049030,15049031,15049032,15049033,15049037,15049038,15049039,15049040,15049041,15049074,15049077,15049083,15049084,15049085,15049086,15049087,15049088,15049089,15049090,15049091,15049092,15049093,15049094,15049096,15049097,15049098,15049099,15050001,15050002,15050003,15050004,15050005,15050006,15050010,15050011,15050012,15050013,15050014,15050015,15050016,15051131,15051132,15051133,15051134,15051137,15051141,15051144,15051146,15051147,15051148,15051149,15051156,15051157,15051159,15051160,15052001,15052004,15052005,15052006,15052007,15052008,15052009,15052010,15052011,15052014,15052015,15052016,15052017,15052018,15052019,15052020,15052021,15052022,15052023,15052024,15052025,15052026,15052040,15052041,15052042,15052043,15052060,15052061,15052062,15052063,15052064,15052065,15052066,15052067,15052068,15052069,15052070,15052073,15052074,15052075, 15053027,15053028,15053029,15053034,15053035,15053052,15053054,15053055,15054053,15054054,15055018,15055137,15056117,15057055,15057059,15058006,15058011,15058021,15059057,15059058,15061001,15061009,15062006,15062069,15065012,15065014,15066070,15068013,15068014,15068016,15068018,15069036,15070008,15070009,15070010
};
const Int_t bad_run_list_19GeV[35+1]  = {12113091,12114007,12114035,12114078,12114092,12114116,12115009,12115014,12115015,12115016,12115018,12115019,12115020,12115022,12115023,12115062,12115073,12115093,12115094,12116012,12116054,12117010,12117016,12117020,12117065,12119040,12119042,12120017,12120026,12121017,12121022,12121034,12121050,12121067,12122019, 0};
const Int_t bad_run_list_27GeV[34+1]  = {12172050,12172051,12172055,12173030,12173031,12173032,12173033,12173034,12174067,12174085,12175062,12175087,12175113,12175114,12175115,12176001,12176044,12176054,12176071,12177015,12177061,12177092,12177099,12177101,12177106,12177107,12177108,12178003,12178004,12178005,12178006,12178013,12178099,12178120, 0};
const Int_t bad_run_list_39GeV[38+1]  = {11199124,11100002,11100045,11101046,11102012,11102051,11102052,11102053,11102054,11102055,11102058,11103035,11103056,11103058,11103092,11103093,11105052,11105053,11105054,11105055,11107007,11107042,11107057,11107061,11107065,11107074,11108101,11109013,11109077,11109088,11109090,11109127,11110013,11110034,11110073,11110076,11111084,11111085, 0};
const Int_t bad_run_list_62GeV[105+1] = {11080072,11081023,11081025,11082012,11082013,11082046,11082056,11082057,11084009,11084011,11084012,11084013,11084020,11084021,11084035,11084044,11084064,11085015,11085025,11085030,11085046,11085055,11085056,11085057,11086005,11086007,11087001,11087002,11087003,11087004,11088013,11089026,11089028,11089029,11089055,11089068,11089072,11091007,11091015,11091021,11091078,11092010,11092011,11092012,11092032,11092033,11092034,11092067,11092096,11093001,11094016,11094017,11094018,11094019,11094020,11094021,11094022,11094023,11094024,11094027,11094028,11094042,11094044,11094045,11094046,11094047,11094048,11094050,11094051,11094052,11094053,11094054,11094055,11094074,11094075,11094077,11095001,11095002,11095003,11095004,11095005,11095006,11095009,11095010,11095011,11095012,11095013,11095014,11095015,11095022,11095040,11095048,11095050,11095051,11095061,11095062,11095063,11095064,11095082,11095087,11096024,11096039,11096043,11096044,11097093, 0};
const Int_t bad_run_list_200GeV_2010[219+1]  = { 11002120, 11002121, 11002126, 11002127, 11002129, 11003010, 11003011, 11003101, 11003102, 11004007, 11004008, 11004009, 11004010, 11004011, 11004012, 11004013, 11004014, 11004015, 11004016, 11004018, 11004020, 11004021, 11004023, 11004024, 11004025, 11004026, 11004028, 11004029, 11004030, 11004032, 11004033, 11004034, 11004035, 11004037, 11004038, 11005042, 11006004, 11006005, 11006008, 11007015, 11010031, 11011019, 11011053, 11015069, 11015071, 11016024, 11017006, 11018003, 11018007, 11018008, 11018036, 11019001, 11019080, 11019081, 11021027, 11021028, 11021031, 11023048, 11025034, 11025038, 11025054, 11025067, 11025069, 11026005, 11026008, 11026021, 11026022, 11026023, 11026025, 11026067, 11026068, 11028004, 11028005, 11028006, 11028007, 11028008, 11028009, 11028010, 11028011, 11028012, 11028013, 11028018, 11028019, 11028020, 11028021, 11028022, 11028023, 11028024, 11028025, 11028026, 11028027, 11030041, 11030080, 11031061, 11031064, 11035008, 11035009, 11035072, 11036026, 11037035, 11037037, 11037060, 11037066, 11037067, 11038048, 11038049, 11038050, 11039047, 11039067, 11040078, 11040083, 11041022, 11041023, 11041040, 11041041, 11042001, 11042002, 11042003, 11042004, 11042005, 11042006, 11042007, 11042008, 11042011, 11042012, 11042018, 11042019, 11042020, 11042021, 11042022, 11042023, 11042024, 11042025, 11042026, 11042027, 11042042, 11042043, 11042044, 11042045, 11042046, 11042047, 11042048, 11042049, 11044029, 11047059, 11047065, 11047066, 11047067, 11048037, 11049001, 11049002, 11049005, 11049023, 11051038, 11051049, 11051051, 11051055, 11051063, 11051064, 11051068, 11052011, 11053057, 11054021, 11054022, 11054024, 11054059, 11054062, 11054066, 11057012, 11057035, 11057036, 11058005, 11058050, 11058083, 11059043, 11059055, 11059060, 11059075, 11059076, 11059077, 11060008, 11060049, 11060059, 11060069, 11060076, 11061008, 11061009, 11061021, 11061034, 11061037, 11061038, 11061095, 11063006, 11063007, 11063008, 11063011, 11063013, 11063014, 11063015, 11063016, 11063017, 11063036, 11063083, 11064003, 11064023, 11065038, 11066024, 11066045, 11071056, 11072032, 11072044, 11072045, 11073001, 11073002, 11073003, 11073049, 11075039, 11075045, 11075048, 0 };

const Int_t bad_run_list_200GeV_2011[134+1]  = {
12113091, 12114007, 12114035, 12114078, 12114092, 12114116, 12115009, 12115014, 12115015, 12115016, 12115018, 12115019, 12115020, 12115022, 12115023, 12115062,
12115073, 12115093, 12115094, 12116012, 12116054, 12117010, 12117016, 12117020, 12117065, 12119040, 12119042, 12120017, 12120026, 12121017, 12121022, 12121034,
12121050, 12121067, 12122019, 12127003, 12127010, 12127011, 12127017, 12127018, 12127032, 12128025, 12132043, 12132061, 12133018, 12134023, 12136005, 12136006,
12136014, 12136017, 12136022, 12136023, 12136024, 12136025, 12136027, 12136028, 12136029, 12136030, 12136031, 12136034, 12136054, 12138005, 12138017, 12138021,
12146004, 12146006, 12146007, 12146008, 12151035, 12153002, 12153004, 12153007, 12153013, 12157038, 12157051, 12158040, 12158041, 12158054, 12158056, 12158057,
12162055, 12162056, 12162057, 12162058, 12164037, 12164078, 12164079, 12166002, 12166003, 12167015, 12167024, 12167052, 12168002, 12168009, 12168022, 12168077,
12170044, 12170045, 12170054, 12170056, 12172050, 12172051, 12172055, 12173030, 12173031, 12173032, 12173033, 12173034, 12174067, 12174085, 12175062, 12175087,
12175113, 12175114, 12175115, 12176001, 12176044, 12176054, 12176071, 12177015, 12177061, 12177092, 12177099, 12177101, 12177106, 12177107, 12177108, 12178003,
12178004, 12178005, 12178006, 12178013, 12178099, 12178120, 0 };

const Int_t bad_run_list_200GeV_2014[95+1]  = {
15075007, 15075065, 15075073, 15075079, 15076101, 15076102, 15076105, 15076109, 15077080, 15078075, 15078110, 15078111, 15079042, 15079046, 15079047, 15079048, 15079050, 15079052, 15079057, 15079058, 15079060, 15079061, 15079063, 15080003, 15080004, 15080005, 15080006, 15080007, 15080008, 15080013, 15080014, 15080015, 15080016, 15080024, 15080029, 15080035, 15080036, 15080037, 15080038, 15080039, 15080042, 15080044, 15080045, 15080053, 15080054, 15080055, 15080056, 15080057, 15080058, 15080059, 15080061, 15080063, 15081001, 15081003, 15081015, 15081017, 15081024, 15081025, 15082016, 15084030, 15086076, 15097032, 15097034, 15099001, 15100100, 15100101, 15100102, 15100103, 15102024, 15104016, 15104018, 15108018, 15108019, 15110032, 15121062, 15122045, 15131052, 15131053, 15151041, 15151042, 15156008, 15161051, 15166014, 15166015, 15166016, 15166017, 15097006, 15109039, 15121062, 15122045, 15124004, 15124031, 15124033, 15131052, 15131053, 0};

  const Int_t *use_run;
  if(8 == mrootS) use_run = bad_run_list_7GeV; 
  else if(11 == mrootS) use_run = bad_run_list_11GeV;
  else if(15 == mrootS) use_run = bad_run_list_14GeV;
  else if(20 == mrootS) use_run = bad_run_list_19GeV;
  else if(27 == mrootS) use_run = bad_run_list_27GeV;
  else if(39 == mrootS) use_run = bad_run_list_39GeV;
  else if(62 == mrootS) use_run = bad_run_list_62GeV;
  else if(11 == mRunYear && 200 == mrootS) use_run = bad_run_list_200GeV_2010;
  else if(12 == mRunYear && 200 == mrootS) use_run = bad_run_list_200GeV_2011;
  else if(15 == mRunYear && 200 == mrootS) use_run = bad_run_list_200GeV_2014;
  else {use_run = null_list; cout << "There is no bad run list!!" << endl;}

  for(int i=0; use_run[i] != 0; i++){
    if(picoEvent->runId()== use_run[i]) return false;
  }
  return true;
}
//----------------------------------------------------------------------------
bool MyAnalysisMaker::AcceptEvent(StPicoEvent* picoEvent) {
  // Cut parameters for each event

  Float_t VertexZMin  = 0.0 ;  //distance event occurs along beam direction from center (cm)
  if(mrootS == 11) VertexZMin  = -50.0 ;
  else if(mrootS == 39) VertexZMin  = -40.0 ;
  else if(mrootS != 39 && mrootS != 11 && mrootS < 50) VertexZMin  = -70.0 ;
  if(200 == mrootS && 12 == mRunYear) VertexZMin = -30.;
  if(200 == mrootS && 15 == mRunYear) VertexZMin = -6.;
  if(200 == mrootS && 17 == mRunYear) VertexZMin = -6.; 
  if(VertexZMin == 0.0) return false;
  histogram1D[0]->Fill(3);

  const Float_t VertexZMax  =  - VertexZMin ;  // cm
  Float_t VertexDiff = 1.e6;
  if(200 == mrootS) VertexDiff = 3.;

  const Int_t   MultMax     = 1000 ;  //max multiplicity

  Int_t Multiplicity =  picoEvent->refMult() ;  // Reference Multiplicity for Centrality determination
  if ( /*Multiplicity < MultMin     ||*/ Multiplicity >= MultMax ){     //BOTH are used by Yadav
             return false;}
 // if ( Multiplicity < MultMin     || Multiplicity >= MultMax ) return false;
  histogram1D[0]->Fill(4);

  Int_t tofMult = picoEvent->btofTrayMultiplicity ();
  if ( tofMult < 2 ) return false;
  histogram1D[0]->Fill(5);

  Float_t vertex[3] = {0} ;
  vertex[0] = picoEvent->primaryVertex().X() ;
  vertex[1] = picoEvent->primaryVertex().Y() ;
  vertex[2] = picoEvent->primaryVertex().Z() ;

  float epsilon = 1.e-6;
  if ( vertex[0] < epsilon && vertex[1] < epsilon && vertex[2] < epsilon ) return false;  // Skip events without a primary vertex
  histogram1D[0]->Fill(6);
  // Cut on Vertex location
  if ( vertex[2] < VertexZMin || vertex[2] > VertexZMax ) return false;

  //cut for Beam Pipe events
  if(mrootS != 15) {
    if ( sqrt(vertex[0]*vertex[0]+ vertex[1]*vertex[1]) > 2.0){
      return false;}
  }
  else if(mrootS == 15) {
    if ( sqrt(vertex[0]*vertex[0]+ (vertex[1]+0.89)*(vertex[1]+0.89)) > 1.0){
      return false;}
  }
  histogram1D[0]->Fill(7);

  if(fabs(vertex[2] - mVpdVz) > VertexDiff) return false;
  histogram1D[0]->Fill(8);

  //Cut for BBC Adc Saturation
  if(mBbcSumEast<75 || mBbcSumWest<75) return false;
  histogram1D[0]->Fill(9);

  return true ;}
//----------------------------------------------------
bool MyAnalysisMaker::EtaSymmetryCut()
//enforce eta symmetry - is this necessary?
{
  int mEtaSymPosTpcN = 0;
  int mEtaSymNegTpcN = 0;
  for(unsigned int itrack1 = 0; itrack1 < mPicoDst->numberOfTracks(); itrack1++) {
    StPicoTrack* track = mPicoDst->track(itrack1);
    float eta = track->gMom().Eta();
    if (eta > 0.) mEtaSymPosTpcN++;
    else mEtaSymNegTpcN++;    
  }
  if ((mEtaSymPosTpcN + mEtaSymNegTpcN)== 0.0 ) return false;
  histogram2D[2]->Fill(mEtaSymPosTpcN, mEtaSymNegTpcN);
  float etaSymTpc =(mEtaSymPosTpcN - mEtaSymNegTpcN)/sqrt(mEtaSymPosTpcN + mEtaSymNegTpcN);
  if (etaSymTpc > 5.0 || etaSymTpc < -5.0 ) return false;

  return true ;
}
//------------------------------------------------------------------

bool MyAnalysisMaker::AcceptTrack(StPicoTrack* track)

{

  // Cut Parameters for individual tracks

  //const Float_t dcaCut =   0.1  ; // cm
  const Float_t PtMin  =   0.15 ; // GeV
  const Float_t PtMax  =   10.0  ; // GeV
  const Float_t EtaMin =  -1.0  ;
  const Float_t EtaMax =   1.0  ;
  const Float_t FitRatio =  0.52  ;      // Number of hits over number of hits possible
  const Int_t   nHitMin  =    15  ;      // 15 is typical but sometimes goes as high as 25
  const Int_t   nHitMax  =   100  ;      // 45 pad rows in the TPC and so anything bigger than 45+Silicon is infinite
  const Int_t   nHitPossMin =  5  ;      // Don't bother to fit tracks if # possible hits is too low, also protects / 0.

  TVector3 trackMomentum = track->gMom();
  
  histogram1D[6]->Fill(0);

  //for some reason StPicoTrack doesn't have "flag" so I just got rid of this cut
  //if ( track->flag() <  0  )     return false ;          // Track quality
  histogram1D[6]->Fill(1);
  //if ( track->dcaGlobal().mag() < dcaCut  )  return false ;   // 3D DCA for global tracks

  if ( track->nHitsMax() <  nHitPossMin )  return false ;    // Minimum number of Possible hits, see above.
  histogram1D[6]->Fill(2);

  if ( trackMomentum.Eta()       <  EtaMin  || trackMomentum.Eta()      >  EtaMax  ) return false ;
  histogram1D[6]->Fill(3);

  if ( trackMomentum.Pt()        <  PtMin   || trackMomentum.Pt()       >  PtMax   ) return false ;
  histogram1D[6]->Fill(4);

  if ( track->nHitsFit()  <  nHitMin || track->nHitsFit() >  nHitMax ) return false;
  histogram1D[6]->Fill(5);

  //if ( ( (float)track->nHitsFit() / (float)track->nHitsMax() ) < FitRatio )  return false ;
  histogram1D[6]->Fill(6);

  return true ;

}
//-------------------------------------------------------------------------
bool MyAnalysisMaker::ProtonCheck(Float_t proton){
  if (proton <-3.0 || proton > 3.0) return false;   //checking nsigma for proton
  histogram1D[6]->Fill(9);
  return true ;
}
//------------------------------------------------------
bool MyAnalysisMaker::ProtonMCheck(Float_t beta, Float_t p){
  //checking the proton mass
  if(beta == -999) return true;
  float masssqr = p*p*(1.0/(beta*beta)-1.0);//mass squre

  if (masssqr < 0.5 || masssqr > 1.5) return false ;//Alex www.star.bnl.gov/protected/heavy/aschmah/BES_PID_v2_paper/technical_note/Technical_note_BES_PID_v2_V0_03.pdf
  histogram1D[6]->Fill(10);
  return true ;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::ProtonPtCut(bool isProton, float ProtonPt){
  //minimum pt cut for protons and antiprotons to remove knock-outs
  //isProton is true for protons, false for antiprotons
  if(isProton && ProtonPt < 0.5) return false;
  if(!isProton && ProtonPt < 0.3) return false;

  return true;  
}
//------------------------------------------------------
bool MyAnalysisMaker::PionCheck(Float_t pion){
  //nsigma
  if (pion <-3.0 || pion > 3.0) return false;
  histogram1D[6]->Fill(7);
  return true;
}
//------------------------------------------------------
bool MyAnalysisMaker::PionMCheck(Float_t beta, Float_t p){

  if(beta == -999) return true;
  float masssqr = p*p*(1.0/(beta*beta)-1.0);//mass squre
  if (masssqr < 0.017-0.013*p || masssqr > 0.04) return false ;
  histogram1D[6]->Fill(8);
  return true;

}
//--------------------------------------------------------------
bool MyAnalysisMaker::V0DirectionCut(Float_t pLxLxPV){
  //pLxLxPV is the product (of vectors)        p_Lambda . (x_Lambda - x_PrimaryVertex)     a cut on the direction of the lambda
  if(pLxLxPV < 0) return false;
  return true;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::DaughterDcaCut(Float_t DaughterDCA, float LambdaPt){
  //cut on the DCA of proton to pion (daughters)
  float cutfactor = 1.; //loosen relative to noninal cuts
  if(DaughterDCA > 0.75*cutfactor) return false;

  return true;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::LambdaDcaCut(Float_t LambdaDCA, float LambdaPt){

  float cutfactor = 1.;
  if(LambdaPt < 0.8 && LambdaDCA > 0.7*cutfactor) return false;
  else if(LambdaPt > 0.8 && LambdaPt < 3.6 && LambdaDCA > 0.75*cutfactor) return false;
  else if(LambdaPt > 3.6 && LambdaDCA > 0.4*cutfactor) return false;

  return true;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::LambdaDecayLengthCut(Float_t DecayLength, float LambdaPt){
  //LPVdiff is the magnitude of the difference between the lambda's decay point and the PV
  //cut on the DCA of proton to pion (daughters), the lambda, and direction of momentum
  float cutfactor = 1.;

  if(LambdaPt < 0.8 && DecayLength < 4./cutfactor) return false;
  else if(LambdaPt > 0.8 && LambdaPt < 3.6 && DecayLength < 4./cutfactor) return false;
  else if(LambdaPt > 3.6) if(DecayLength < 10./cutfactor || DecayLength > 125.) return false;

  return true;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::ProtonDcaCut(Float_t ProtonDCA, float LambdaPt){
  //DCAs are to primary vertex, ProTOF is the pathlength in the tof (just a quantity to test is there is TOF data)
  float cutfactor = 1.;

  if(LambdaPt < 0.8 && ProtonDCA < 1.0/cutfactor) return false;
  else if(LambdaPt > 0.8 && LambdaPt < 3.6 && ProtonDCA < 0.75/cutfactor) return false;
  else if(LambdaPt > 3.6 && ProtonDCA < 0./cutfactor) return false;

  return true;  
}
//--------------------------------------------------------------
bool MyAnalysisMaker::PionDcaCut(Float_t PionDCA, float LambdaPt){
  //pion DCAs are to primary vertex, ProTOF is the pathlength in the tof (just a quantity to test is there is TOF data)
  //float LooseCutPimDCA = 1./3.;
  float cutfactor = 1.; 
  if(LambdaPt < 0.8 && PionDCA < 2.5/cutfactor) return false;
  else if(LambdaPt > 0.8 && LambdaPt < 3.6 && PionDCA < 2.0/cutfactor) return false;
  else if(LambdaPt > 3.6 && PionDCA < 1.0/cutfactor) return false;  

  return true;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::LambdaMassCut(Float_t LambdaMass){
  //this is the distance between the decay vertex and the primary vertex
  //it is NOT the distance from the decay vertex to the DCA to the PV
  float cutfactor = 1.; //loosen relative to noninal cuts

  if(LambdaMass < 1.115683 - 0.005*cutfactor) return false;
  else if(LambdaMass > 1.115683 + 0.005*cutfactor) return false;

  return true;
}
//--------------------------------------------------------------
double MyAnalysisMaker::CalculateHelicity(TLorentzVector PLam4, TLorentzVector PPro4){
  //function to calculate helicity, p_{proton}^{*} . p_{Lambda}

  PPro4.Boost(-PLam4.BoostVector());
  TVector3 ProtonInLambdaFrame = PPro4.Vect();
  ProtonInLambdaFrame = ProtonInLambdaFrame*(1./(ProtonInLambdaFrame.Mag()));
  TVector3 LambdaInLabFrame = PLam4.Vect();
  LambdaInLabFrame = LambdaInLabFrame*(1./(LambdaInLabFrame.Mag()));

  return ProtonInLambdaFrame.Dot(LambdaInLabFrame);
}
//--------------------------------------------------------------
float MyAnalysisMaker::GetDeltaPhi(TVector3 p1, TVector3 p2){
  float phi1 = (p1+p2).Phi();   
  
  TVector3 pminus;
  TRandom3 rand;
  rand.SetSeed(0);
  float rng = rand.Rndm();
  if(rng > 0.5) pminus = p1-p2;
  else pminus = p2-p1;
  float phi2 = pminus.Phi();

  float deltaphi = phi1 - phi2;
  if(deltaphi < -pi) deltaphi += 2.*pi;

  else if(deltaphi > pi) deltaphi -= 2.*pi;
  // cout << "rng " << rng << endl;
  // cout << "phi2 " << phi2 << endl;
  // cout << "deltaphi " << deltaphi << endl;
  return deltaphi;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::CentralityCut(int CentralityID, float RefMult){
  if (CentralityID <= 2 || RefMult <= 5) return true;
  return false;
}
//--------------------------------------------------------------
bool MyAnalysisMaker::RhoPtCut(float pt){
  if (pt > 0.2) return false;
  return true;
}
//--------------------------------------------------------------

