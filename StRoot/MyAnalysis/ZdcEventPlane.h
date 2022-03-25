#ifndef ZdcEventPlane_hh
#define ZdcEventPlane_hh

#include "TObject.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StThreeVectorF.hh"
#include "TFile.h"
#include "TProfile2D.h"


class StPicoTrack;
class StPicoEvent;


class ZdcEventPlane : public TObject {


  public:
    ZdcEventPlane(int CalibMode, TFile* ShiftCorrection, TFile* ZdcBeamCenter, bool CanWander);
    ~ZdcEventPlane();
    void GetHists();
    void AnalyzeZdcEvent(StPicoEvent* picoEvent, int CentralityID);
    Double_t* GetPsi();
    Double_t* GetQCenter();
    Double_t* GetQCorrected();
    int GetWanderSteps(){return mWanderStep;}

  private:
    int mCentralityID;
    Float_t mZDCSMD[2][2][8];
    Double_t mZdc_phi;
    Double_t mPsi_full;
    Double_t mPsi_e;
    Double_t mPsi_w;	
    Double_t ZDC_GetPosition(int , int , int);// 
    Double_t ZDC_GetPositionCal(int , int , int);// 
    Double_t mPsi[6];  // return value of get phi
    Double_t mQCenter[4];
    Double_t mQCorrected[4]; //This is the corrected Q vector
    Double_t mReactionPlaneFS1; //I don't know what this is
    int mRunYear;
    int mRunID;
    int mRunDay;
    int mRunDayBin;
    int mCalibMode;
    bool mCorrectionWander;
    int mWanderStep;

    float  mZDCshiftEast_cos[50];
    float  mZDCshiftEast_sin[50];
    float  mZDCshiftWest_cos[50];
    float  mZDCshiftWest_sin[50];
    float  mZDCshiftFull_cos[50];
    float  mZDCshiftFull_sin[50];

    //run11 200 GeV
    Double_t  mZDCCenterEX;// 0.0 ;//4.225;   //zdcsmd_wx0
    Double_t  mZDCCenterEY;// 5.89;//4.727;   //zdcsmd_ex0
    Double_t  mZDCCenterWX;// 4.78;//4.967;   //zdcsmd_wy0
    Double_t  mZDCCenterWY;// 5.611;//5.621;   //zdcsmd_ey0


    TFile* mShiftFile;
    TFile* mZDCSMDConstFile;
    TProfile2D *mZDCSMDBeamCenter; //I don't know why two of these exist, just copying Subhash
    TProfile2D* mZDCSMDBeamCenter2D;
    TProfile* mProjectedY;

    TProfile2D* mZDC_shiftWest_c;
    TProfile2D* mZDC_shiftWest_s;
    TProfile2D* mZDC_shiftEast_c;
    TProfile2D* mZDC_shiftEast_s;
    TProfile2D* mZDC_shiftFull_c;
    TProfile2D* mZDC_shiftFull_s;



    ClassDef(ZdcEventPlane,1)
};
#endif
