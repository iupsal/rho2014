

#include "TMath.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "ZdcEventPlane.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

// ZDCSMD gain factor
// 200GeV AuAu  run11
Double_t zdcsmdGainFac[2][2][8] = { 
  //2011 weights
  // { 
  //   { 1,1.2567,1.12124,1.10505,0.884667,1.00386,1.15781,1},
  //   { 1,1.2691,1.40142,1.37204,1.11616,1.36999,1.41652,1.51231}
  // },
  // {
  //   { 1,1.24905,1.46427,1.22947,1.09283,1.06398,1.09407,1},
  //   { 1,1.05491,1.2956,1.07312,1.15465,1.50818,1.15796,1.03961}
  // }

  //2014 weights from Subhash
  {
    {1,1.02411,0.902382,0.882074,0.767194,0.891737,1.13113,1},
    {1,1.20689,1.18671,1.15595,1.0111,1.17598,1.29584,1.61315}
  },
  {
    {1,1.07303,1.14271,1.03752,0.98926,0.969333,1.12912,1.1703},
    {1,0.841123,0.907441,0.733425,0.820226,1.00744,0.93808,0.966097}
  }  
};

//float pi = TMath::Pi();




ClassImp(ZdcEventPlane)


ZdcEventPlane::~ZdcEventPlane(){}
//-----------------------------------------------------------------------------------------------//
ZdcEventPlane::ZdcEventPlane(int rCalibMode, TFile* rShiftFile, TFile* rZDCSMDConstFile, bool CanWander = false){
  mCalibMode = rCalibMode;
  mShiftFile = rShiftFile;
  mZDCSMDConstFile = rZDCSMDConstFile;
  mCorrectionWander = CanWander;  
  //let correction histogram wander away from correct time-bin to find some empty bin
  mWanderStep = 0; 
  GetHists();
}
//-----------------------------------------------------------------------------------------------//
void ZdcEventPlane::GetHists(){
  if(mCalibMode > 0){
    mZDCSMDBeamCenter2D = (TProfile2D *)mZDCSMDConstFile->Get("ZDCSMDBeamCenter");
    mProjectedY = mZDCSMDBeamCenter2D->ProfileY();
  }
  if(2 == mCalibMode){
    //ZDC psi shift   
    mZDC_shiftWest_c = (TProfile2D*)mShiftFile->Get("ShiftZdc1CosWest");
    mZDC_shiftWest_s = (TProfile2D*)mShiftFile->Get("ShiftZdc1SinWest");
    mZDC_shiftEast_c = (TProfile2D*)mShiftFile->Get("ShiftZdc1CosEast");
    mZDC_shiftEast_s = (TProfile2D*)mShiftFile->Get("ShiftZdc1SinEast");
    mZDC_shiftFull_c = (TProfile2D*)mShiftFile->Get("ShiftZdc1Cos");
    mZDC_shiftFull_s = (TProfile2D*)mShiftFile->Get("ShiftZdc1Sin");
    // mZDC_shiftWest_c = (TProfile2D*)mShiftFile->Get("mZDC_shiftwest_cos");
    // mZDC_shiftWest_s = (TProfile2D*)mShiftFile->Get("mZDC_shiftwest_sin");
    // mZDC_shiftEast_c = (TProfile2D*)mShiftFile->Get("mZDC_shifteast_cos");
    // mZDC_shiftEast_s = (TProfile2D*)mShiftFile->Get("mZDC_shifteast_sin");
    // mZDC_shiftFull_c = (TProfile2D*)mShiftFile->Get("mZDC_shiftfull_cos");
    // mZDC_shiftFull_s = (TProfile2D*)mShiftFile->Get("mZDC_shiftfull_sin");

  }
}
//-----------------------------------------------------------------------------------------------//
void ZdcEventPlane::AnalyzeZdcEvent(StPicoEvent* picoEvent, int rCentralityID){

  mCentralityID = rCentralityID;

  mRunID = picoEvent->runId(); 
  mRunYear = floor( mRunID/pow(10,6) );
  mRunDay = floor( (mRunID - mRunYear*pow(10,6))/pow(10,3) );  
  mRunDayBin = -10;
  if(15 == mRunYear) mRunDayBin = mRunDay - 74;

  //cout << "mRunID " << mRunID << endl;
  Double_t fZdcSmdEastHorizontal[8];
  Double_t fZdcSmdWestHorizontal[8];
  Double_t fZdcSmdEastVertical[8];
  Double_t fZdcSmdWestVertical[8];
   //added subhash
    //zdcsmd loop
  for(int strip=0; strip<8; strip++){
      fZdcSmdEastHorizontal[strip] = picoEvent->ZdcSmdEastHorizontal(strip);
      fZdcSmdWestHorizontal[strip] = picoEvent->ZdcSmdWestHorizontal(strip);
      fZdcSmdEastVertical[strip] = picoEvent->ZdcSmdEastVertical(strip);
      fZdcSmdWestVertical[strip] = picoEvent->ZdcSmdWestVertical(strip);
      //cout << strip << "  " << fZdcSmdEastHorizontal[strip]<<"  "<< fZdcSmdWestHorizontal[strip]  <<"  "<< fZdcSmdEastVertical[strip]<< "  "<< fZdcSmdWestVertical[strip]<<endl;
  }
      
  //[ east=0, west=1] [ 0=vertical, 1=Horizontal] [ strips]
  //Double_t mZDCSMD[0][1][8];// = 0.0;
  //Double_t mZDCSMD[1][1][8];// = 0.0;
  //Double_t mZDCSMD[0][0][8];// = 0.0;
  //Double_t mZDCSMD[1][0][8];// = 0.0;

  for(int strip=0; strip<8; strip++)
    {

    if(fZdcSmdEastHorizontal[strip]) mZDCSMD[0][1][strip] = fZdcSmdEastHorizontal[strip]*zdcsmdGainFac[0][1][strip];//east,hor
    if(fZdcSmdWestHorizontal[strip]) mZDCSMD[1][1][strip] = fZdcSmdWestHorizontal[strip]*zdcsmdGainFac[1][1][strip];//west,hor
    if(fZdcSmdEastVertical[strip]) mZDCSMD[0][0][strip] = fZdcSmdEastVertical[strip]*zdcsmdGainFac[0][0][strip];//east,ver
    if(fZdcSmdWestVertical[strip]) mZDCSMD[1][0][strip] = fZdcSmdWestVertical[strip]*zdcsmdGainFac[1][0][strip];//west, ver

    //cout << "strip " << strip << fZdcSmdEastHorizontal[strip] << ", mZDCSMD[0][1][strip] " << mZDCSMD[0][1][strip] << ", fZdcSmdEastHorizontal[strip] " << fZdcSmdEastHorizontal[strip] << ", zdcsmdGainFac[0][1][strip] " << zdcsmdGainFac[0][1][strip] << endl;
    //cout << mZDCSMD[0][1][strip] << " "<<mZDCSMD[1][1][strip]<<" "<< mZDCSMD[0][0][strip] <<" "<<mZDCSMD[1][1][strip]<<endl;
    }

  
  //ZDC East, West, Full 
  TVector2 mQE;
  TVector2 mQW;
  TVector2 mQ;
  
  Double_t mQEx=0., mQEy=0.;
  Double_t mQWx=0., mQWy=0.;
  Double_t mQx=0.,  mQy=0.;
  
  Float_t eXsum=0.,  eYsum=0.,  eXWgt=0., eYWgt=0., psi_e=0., psi_e_s=0.;
  Float_t wXsum=0.,  wYsum=0.,  wXWgt=0., wYWgt=0., psi_w=0., psi_w_s=0.;
  Float_t eXFsum=0., wXFsum=0., eYFsum=0., wYFsum=0., eWFgt=0., wWFgt=0., psi_f=0., psi_f_s=0.;
  
  for(int strip=0;strip<8;strip++){
    eYsum += ZDC_GetPosition(0,1,strip) * mZDCSMD[0][1][strip];  //cout << "ZDC_GetPosition(0,1,strip) " << ZDC_GetPosition(0,1,strip) << ", mZDCSMD[0][1][strip] " << mZDCSMD[0][1][strip] << endl;
    wYsum += ZDC_GetPosition(1,1,strip) * mZDCSMD[1][1][strip];  //cout << "ZDC_GetPosition(1,1,strip, mRunID) " << ZDC_GetPosition(1,1,strip, mRunID) << ", mZDCSMD[1][1][strip] " << mZDCSMD[1][1][strip] << endl;
    eYWgt += mZDCSMD[0][1][strip];
    wYWgt += mZDCSMD[1][1][strip];
    //if(strip>7) continue;
    if(strip>6) continue;
    eXsum += ZDC_GetPosition(0,0,strip) * mZDCSMD[0][0][strip];  //cout << "ZDC_GetPosition(0,0,strip, mRunID) " << ZDC_GetPosition(0,0,strip, mRunID) << ", mZDCSMD[0][0][strip] " << mZDCSMD[0][0][strip] << endl;
    wXsum += ZDC_GetPosition(1,0,strip) * mZDCSMD[1][0][strip];  //cout << "ZDC_GetPosition(1,0,strip, mRunID) " << ZDC_GetPosition(1,0,strip, mRunID) << ", mZDCSMD[1][0][strip] " << mZDCSMD[1][0][strip] << endl;
    eXWgt += mZDCSMD[0][0][strip];
    wXWgt += mZDCSMD[1][0][strip];
    
  }
  //cout << "eYWgt " << eYWgt << ", wYWgt " << wYWgt << ", eXWgt " << eXWgt << ", wXWgt " << wXWgt << endl;
  // East    
  mQEx= (eXWgt>0.0) ? eXsum/eXWgt:0.0;
  mQEy= (eYWgt>0.0) ? eYsum/eYWgt:0.0;
  
  // West   
  mQWx= (wXWgt>0.0) ? wXsum/wXWgt:0.0;
  mQWy= (wYWgt>0.0) ? wYsum/wYWgt:0.0;
  
  // Full  
  mQx=(eXWgt>0. && wXWgt>0.) ? wXsum/wXWgt - eXsum/eXWgt:0.;
  mQy=(eYWgt>0. && wYWgt>0.) ? wYsum/wYWgt - eYsum/eYWgt:0.;
  
  mQE.Set(mQEx,mQEy);
  mQW.Set(mQWx,mQWy);
  mQ.Set(mQx,mQy);
  
  //cout << "mQEx " << mQEx << ", mQEy " << mQEy << ", mQWx " << mQWx << ", mQWy " << mQWy << endl;

  if(mQE.Mod() && mQW.Mod() && mQ.Mod() ){
    psi_e=mQE.Phi();
    psi_w=mQW.Phi();
    psi_f= mQ.Phi();
    if(psi_e<0.0) psi_e +=2.*pi;
    if(psi_w<0.0) psi_w +=2.*pi;
    if(psi_f<0.0) psi_f +=2.*pi;
  }
  else
    {
      psi_e=-9999;
      psi_w=-9999;
      psi_f=-9999;
    }
  
  //------------------------------------------------------------------------------------------------------------------------//
  //------------------------------------------------------------------------------------------------------------------------//
  //I want to get good calibrated Q vectors even if calibration mode is > 0, so that my correction histograms are always good
  //ZDC East, West, Full 
  TVector2 mCalQE;
  TVector2 mCalQW;
  
  Double_t mCalQEx=0., mCalQEy=0.;
  Double_t mCalQWx=0., mCalQWy=0.;
  
  Float_t CaleXsum=0.,  CaleYsum=0.,  CaleXWgt=0., CaleYWgt=0.;
  Float_t CalwXsum=0.,  CalwYsum=0.,  CalwXWgt=0., CalwYWgt=0.;
  
  for(int strip=0;strip<8;strip++){
    CaleYsum += ZDC_GetPositionCal(0,1,strip) * mZDCSMD[0][1][strip];  //cout << "ZDC_GetPosition(0,1,strip) " << ZDC_GetPosition(0,1,strip) << ", mZDCSMD[0][1][strip] " << mZDCSMD[0][1][strip] << endl;
    CalwYsum += ZDC_GetPositionCal(1,1,strip) * mZDCSMD[1][1][strip];  //cout << "ZDC_GetPosition(1,1,strip, mRunID) " << ZDC_GetPosition(1,1,strip, mRunID) << ", mZDCSMD[1][1][strip] " << mZDCSMD[1][1][strip] << endl;
    CaleYWgt += mZDCSMD[0][1][strip];
    CalwYWgt += mZDCSMD[1][1][strip];
    //if(strip>7) continue;
    if(strip>6) continue;
    CaleXsum += ZDC_GetPositionCal(0,0,strip) * mZDCSMD[0][0][strip];  //cout << "ZDC_GetPosition(0,0,strip, mRunID) " << ZDC_GetPosition(0,0,strip, mRunID) << ", mZDCSMD[0][0][strip] " << mZDCSMD[0][0][strip] << endl;
    CalwXsum += ZDC_GetPositionCal(1,0,strip) * mZDCSMD[1][0][strip];  //cout << "ZDC_GetPosition(1,0,strip, mRunID) " << ZDC_GetPosition(1,0,strip, mRunID) << ", mZDCSMD[1][0][strip] " << mZDCSMD[1][0][strip] << endl;
    CaleXWgt += mZDCSMD[0][0][strip];
    CalwXWgt += mZDCSMD[1][0][strip];
  }
  //cout << "CaleYWgt " << CaleYWgt << ", CalwYWgt " << CalwYWgt << ", CaleXWgt " << CaleXWgt << ", CalwXWgt " << CalwXWgt << endl;
  // East    
  mCalQEx= (CaleXWgt>0.0) ? CaleXsum/CaleXWgt:0.0;
  mCalQEy= (CaleYWgt>0.0) ? CaleYsum/CaleYWgt:0.0;
  
  // West   
  mCalQWx= (CalwXWgt>0.0) ? CalwXsum/CalwXWgt:0.0;
  mCalQWy= (CalwYWgt>0.0) ? CalwYsum/CalwYWgt:0.0;
   
  mCalQE.Set(mCalQEx,mCalQEy);
  mCalQW.Set(mCalQWx,mCalQWy);
  //------------------------------------------------------------------------------------------------------------------------//
  //------------------------------------------------------------------------------------------------------------------------//

  if(0 == mCalibMode){
    if(fabs(mQE.X() - mCalQE.X()) > 1e-5 || fabs(mQE.Y() - mCalQE.Y()) > 1e-5){cout << "Isaac your ZDC code is wrong (East)!!!" << endl; mCalQE = mQE;}   
    if(fabs(mQW.X() - mCalQW.X()) > 1e-5 || fabs(mQW.Y() - mCalQW.Y()) > 1e-5){cout << "Isaac your ZDC code is wrong (West)!!!" << endl; mCalQW = mQW;}   
  }

  double mReactionPlaneQEx = mCalQE.X();
  double mReactionPlaneQEy = mCalQE.Y();
  double mReactionPlaneQWx = mCalQW.X();
  double mReactionPlaneQWy = mCalQW.Y();

  

  if(fabs(mReactionPlaneQEx)<10  && fabs(mReactionPlaneQEy)<10 && fabs(mReactionPlaneQWx)<10  && fabs(mReactionPlaneQWy)<10 ) {
    mQCenter[0] = mReactionPlaneQEx;
    mQCenter[1] = mReactionPlaneQEy;
    mQCenter[2] = mReactionPlaneQWx;
    mQCenter[3] = mReactionPlaneQWy;
    if(mCalibMode > 0){
      mQCorrected[0] = mQE.X();
      mQCorrected[1] = mQE.Y();
      mQCorrected[2] = mQW.X();
      mQCorrected[3] = mQW.Y();     
    }
  }
  //cout << "mQCenter[0] " << mQCenter[0] << ", mQCenter[1] " << mQCenter[1] << ", mQCenter[2] " << mQCenter[2] << ", mQCenter[3] " << mQCenter[3] << endl;

  
  if(mCalibMode==2)
    {
      //Read shift file 
      //TFile *mShiftFile;
      //mShiftFile = new TFile("./Correction/shift.root","READ");
      //mShiftFile = new TFile("/global/homes/s/subhash/nsm/for_shbhs/test/read_v0_tree/make_d0/shift/shift.root","READ");
      
      float  mZDCshiftEast_cos[50];
      float  mZDCshiftEast_sin[50];
      float  mZDCshiftWest_cos[50];
      float  mZDCshiftWest_sin[50];
      float  mZDCshiftFull_cos[50];
      float  mZDCshiftFull_sin[50];

      
      for(int i=1;i<51;i++) {
        mZDCshiftEast_cos[i-1] = mZDC_shiftEast_c->GetBinContent(i, mCentralityID+2);
        mZDCshiftEast_sin[i-1] = mZDC_shiftEast_s->GetBinContent(i, mCentralityID+2);
        mZDCshiftWest_cos[i-1] = mZDC_shiftWest_c->GetBinContent(i, mCentralityID+2);
        mZDCshiftWest_sin[i-1] = mZDC_shiftWest_s->GetBinContent(i, mCentralityID+2);
        mZDCshiftFull_cos[i-1] = mZDC_shiftFull_c->GetBinContent(i, mCentralityID+2);
        mZDCshiftFull_sin[i-1] = mZDC_shiftFull_s->GetBinContent(i, mCentralityID+2);
      }
      
      
      //Initilize the value
      psi_e_s=psi_e;
      psi_w_s=psi_w;
      psi_f_s=psi_f;
      
      for(int j=1;j<51;j++){
        psi_e_s += 2.*(-mZDCshiftEast_sin[j-1]*cos(j*psi_e)+mZDCshiftEast_cos[j-1]*sin(j*psi_e))/(float)j ;
        psi_w_s += 2.*(-mZDCshiftWest_sin[j-1]*cos(j*psi_w)+mZDCshiftWest_cos[j-1]*sin(j*psi_w))/(float)j ;
        psi_f_s += 2.*(-mZDCshiftFull_sin[j-1]*cos(j*psi_f)+mZDCshiftFull_cos[j-1]*sin(j*psi_f))/(float)j ;
      }
      
      while(psi_e_s<0.0)   psi_e_s +=2.*pi;
      while(psi_e_s>2.*pi) psi_e_s -=2.*pi;
      
      while(psi_w_s<0.0)   psi_w_s +=2.*pi;
      while(psi_w_s>2.*pi) psi_w_s -=2.*pi;
      
      while(psi_f_s<0.0)   psi_f_s +=2.*pi;
      while(psi_f_s>2.*pi) psi_f_s -=2.*pi;


    }//shift calib mode=2
  
  
  mPsi[0] = psi_e; 
  mPsi[1] = psi_w;
  mPsi[2] = psi_f;
  mPsi[3] = psi_e_s;
  mPsi[4] = psi_w_s;
  mPsi[5] = psi_f_s;



  //for test
  float tqx = (mQ.Mod())*cos(psi_f_s);
  float tqy = (mQ.Mod())*sin(psi_f_s);
  TVector2 tQ;
  tQ.Set(tqx, tqy);
  mReactionPlaneFS1 = tQ.Phi();
}
//-----------------------------------------------------------------------------------------------//
Double_t* ZdcEventPlane::GetPsi() {
  return mPsi;
}
//-----------------------------------------------------------------------------------------------//
Double_t* ZdcEventPlane::GetQCenter() {
  return mQCenter;
}
//-----------------------------------------------------------------------------------------------//
Double_t* ZdcEventPlane::GetQCorrected() {
  return mQCorrected;
}
//-----------------------------------------------------------------------------------------------//
Double_t ZdcEventPlane::ZDC_GetPosition(int ew, int vh, int strip) {


  Double_t mZDCCenterEX = 0.0;
  Double_t mZDCCenterEY = 0.0;
  Double_t mZDCCenterWX = 0.0;
  Double_t mZDCCenterWY = 0.0;
  
  if(mCalibMode==1 || mCalibMode==2){
    
    //Read and init ZDC-SMD contant file  and Get the BeamCenter
    //Apply Beam center Correction
    
    int EffectiveDayBin = mRunDayBin + mWanderStep;
    //EffectiveDayBin = 4700-1; //4700 is special for testing Subhash's correction file
    int HistEntries = mProjectedY->GetBinEntries(EffectiveDayBin+1);
    //cout << "before entries " << HistEntries << ", mWanderStep " << mWanderStep << endl;
    if(0 == HistEntries && true == mCorrectionWander){
      //this code is massively stupid. My brain is too foggy to fix right now
      //I'm trying to make it so that the code uses the nearest-to-date nonzero correction

      //skip if done before
      if(0 == HistEntries && abs(mWanderStep) > 0){
        EffectiveDayBin += mWanderStep;
        HistEntries = mProjectedY->GetBinEntries(EffectiveDayBin+1);
      } 
      int fbin = 1;
      while(0 == HistEntries){ 
        int BinLo = mProjectedY->GetBinEntries(EffectiveDayBin-fbin+1);
        int BinHi = mProjectedY->GetBinEntries(EffectiveDayBin+fbin+1);
        HistEntries = max(BinLo, BinHi);
        if (HistEntries && BinLo > BinHi){
          fbin = -1.*fbin;
          break;
        }
        if (fbin < 0) fbin = fbin -1;
        else fbin ++;
        if(abs(fbin) > 1000){
          cout << "no entries within 1000 bins!!! BREAK OUT" << endl;
          fbin = 0;
          break;
        }        
      }
      if(fbin > 0) fbin -= 1;
      EffectiveDayBin = EffectiveDayBin + fbin;
      mWanderStep = fbin;
      //cout << "after EffectiveDayBin " << EffectiveDayBin << ", fbin " << fbin << ", mWanderStep " << mWanderStep << endl;
    }
    //cout << "mWanderStep " << mWanderStep << ", mRunDayBin " << mRunDayBin << endl;

    mZDCCenterEX = mZDCSMDBeamCenter2D->GetBinContent(1,EffectiveDayBin+1);
    mZDCCenterEY = mZDCSMDBeamCenter2D->GetBinContent(2,EffectiveDayBin+1);
    mZDCCenterWX = -1.*mZDCSMDBeamCenter2D->GetBinContent(3,EffectiveDayBin+1);
    mZDCCenterWY = mZDCSMDBeamCenter2D->GetBinContent(4,EffectiveDayBin+1);
    // mZDCCenterEX = mZDCSMDBeamCenter2D->GetBinContent(1,4700);  
    // mZDCCenterEY = mZDCSMDBeamCenter2D->GetBinContent(2,4700);
    // mZDCCenterWX = -1.*mZDCSMDBeamCenter2D->GetBinContent(3,4700);
    // mZDCCenterWY = mZDCSMDBeamCenter2D->GetBinContent(4,4700);



    //cout << "beam center values: ex, ey, wx, wy"<<endl;
    //cout << mRunDayBin+1 << " " << mZDCCenterEX << " "<< mZDCCenterEY <<" "<<mZDCCenterWX <<" "<< mZDCCenterWY<<endl;
  }

  Float_t zdcsmd_x[7] = { 0.5,  2,    3.5,  5,    6.5,   8,     9.5         };
  Float_t zdcsmd_y[8] = { 1.25, 3.25, 5.25, 7.25, 9.25, 11.25, 13.25, 15.25 };
  
  
  
  if(mCalibMode==0){

    //with gain correction
    if(ew==0 && vh==0) return zdcsmd_x[strip];
    if(ew==1 && vh==0) return -zdcsmd_x[strip];
    if(ew==0 && vh==1) return zdcsmd_y[strip]/sqrt(2.0);
    if(ew==1 && vh==1) return zdcsmd_y[strip]/sqrt(2.0);
  }
  else {
    //default with beam center correction  
    if(ew==0 && vh==0) return zdcsmd_x[strip] - mZDCCenterEX;     
    if(ew==1 && vh==0) return mZDCCenterWX - zdcsmd_x[strip]; 
    if(ew==0 && vh==1) return zdcsmd_y[strip]/sqrt(2.0) - mZDCCenterEY; 
    if(ew==1 && vh==1) return zdcsmd_y[strip]/sqrt(2.0) - mZDCCenterWY;
  }
  return 0;
}

//-----------------------------------------------------------------------------------------------//
Double_t ZdcEventPlane::ZDC_GetPositionCal(int ew, int vh, int strip) {
  //I want meaningful centering corrections to be written no matter what the calibration mode is

  Double_t mZDCCenterEX = 0.0;
  Double_t mZDCCenterEY = 0.0;
  Double_t mZDCCenterWX = 0.0;
  Double_t mZDCCenterWY = 0.0;


  Float_t zdcsmd_x[7] = { 0.5,  2,    3.5,  5,    6.5,   8,     9.5         };
  Float_t zdcsmd_y[8] = { 1.25, 3.25, 5.25, 7.25, 9.25, 11.25, 13.25, 15.25 };
  
  
  if(ew==0 && vh==0) return zdcsmd_x[strip];
  if(ew==1 && vh==0) return -zdcsmd_x[strip];
  if(ew==0 && vh==1) return zdcsmd_y[strip]/sqrt(2.0);
  if(ew==1 && vh==1) return zdcsmd_y[strip]/sqrt(2.0);

  return 0;
}






