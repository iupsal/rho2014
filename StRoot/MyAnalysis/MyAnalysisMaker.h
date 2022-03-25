// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#ifndef MyAnalysisMaker_def
#define MyAnalysisMaker_def

#include "StMaker.h"
#include "TString.h"
#include <TTree.h>
#include "StThreeVectorF.hh"
#include "ZdcEventPlane.h"

#include "StRefMultCorr/CentralityMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h" 

class StPicoDstMaker ;
class StPicoDst ;
class StPicoEvent    ;
class StPicoTrack    ;

class TFile        ;
class TH1F         ;
class TH2F         ;
class TH3F         ;
class TProfile     ;
class TProfile2D   ;
class TRandom      ;
class TRandom3     ;
class TVector3     ;
class TLorentzVector     ;

#define MaxNumberOfTH1F      200
#define MaxNumberOfTH2F      100
#define MaxNumberOfTH3F      100
#define MaxNumberOfTProfile  200

class MyAnalysisMaker : public StMaker
{
  
 private:
  TFile* mFemtoDstFile;
  TFile* mHistogramOutput             ;     //  Histograms outputfile pointer
  

  StPicoDstMaker* mPicoDstMaker ;     
  StPicoDst*      mPicoDst;                          //!
                   //  Make PicoDst pointer available to member functions
  TTree* mTree; //!   This is when I finally get around to making a TTree.  It is written out in Finish()
  StRefMultCorr*  mRefmultCorrUtil;

  Double_t* mPsi1 ;

  TH1F* histogram1D[MaxNumberOfTH1F]   ;     //  1D Histograms
  TH2F* histogram2D[MaxNumberOfTH2F] ;     //  2D Histograms
  TH3F* histogram3D[MaxNumberOfTH3F] ;     //  3D Histograms
  TProfile* histogramTProfile[MaxNumberOfTProfile] ;     //  Profile Histograms
  TProfile2D* histogramTProfile2D[MaxNumberOfTProfile] ;     //  Profile Histograms

  TFile* histogram_output             ;     //  Histograms outputfile pointer
  TFile* histogram_output_broken             ;     //  Histograms w/ files that don't pass minEvents

  ZdcEventPlane* mZdcEP;


  int mrootS;
  int mRunYear;
  float mVpdVz;
  float mBbcSumEast;
  float mBbcSumWest;
  int mCalibMode1;
  int mCalibMode2;
  int mLambdaTotalCounter;
  int mShiftNumber;
  bool mEpMode;


  
  ULong_t mEventsStarted               ;     //  Number of Events read
  ULong_t mEventsProcessed             ;     //  Number of Events processed and analyzed

  void SetEpCorrections();
  bool AcceptTrack(StPicoTrack*)   ;     //  Function to make cuts on track quality (ie. nhits or eta)

  bool ProtonCheck(Float_t nsigma);
  bool ProtonMCheck(Float_t beta,Float_t p);
  bool ProtonPtCut(bool isProton, float ProtonPt);
  bool PionCheck(Float_t nsigma);
  bool PionMCheck(Float_t beta,Float_t p);

  bool V0DirectionCut(Float_t pLxLxPV); 
  bool ProtonDcaCut(Float_t dca, Float_t pt); 
  bool PionDcaCut(Float_t dca, Float_t pt); 
  bool DaughterDcaCut(Float_t dca, Float_t pt); 
  bool LambdaDcaCut(Float_t dca, Float_t pt); 
  bool LambdaDecayLengthCut(Float_t decaylength, Float_t pt); 
  bool LambdaMassCut(Float_t LambdaMass); 


  bool AcceptEvent   (StPicoEvent*)   ;     //  Function to make cuts on event quanitites (ie. vertex position)
  bool AcceptTrigger (StPicoEvent*)   ;     //  Function to make cuts on the trigger words (Cusomize each year)

  bool EtaSymmetryCut()   ;       // new cut introduced flow low energy
  bool RejectRunNumbers(StPicoEvent*);      // Function to reject bad runnumbers

  double CalculateHelicity(TLorentzVector PLam4, TLorentzVector PPro4);
  float GetDeltaPhi(TVector3 p1, TVector3 p2);
  bool CentralityCut(int CentralityID, float RefMult);
  bool RhoPtCut(float RhoPt);

  TString mFileNameBase;  // malisa - all filenames will be based on this string.
  TFile*  mCorrectionZdcBeamCenter;
  TString mHistogramOutputFileName     ;     //  Name of the histogram output file 
  TFile*  mCorrectionFileEp1;

 protected:


 public:
  //MyAnalysisMaker(StPicoDstMaker* maker, TString FileNameBase="LocalOut");    //  Constructor
  MyAnalysisMaker(StPicoDstMaker* maker, TString FileNameBase="LocalOut", int shiftnumber = 0, int calibmode = 2, bool epmode = 0);    //  Constructor
  virtual          ~MyAnalysisMaker( )   ;          //  Destructor

  Int_t Init    ( ) ;                               //  Initiliaze the analysis tools ... done once
  Int_t Make    ( ) ;                               //  The main analysis that is done on each event
  Int_t Finish  ( ) ;                               //  Finish the analysis, close files, and clean up.

  ClassDef(MyAnalysisMaker,1)                     //  Macro for CINT compatability
    
};

#endif




