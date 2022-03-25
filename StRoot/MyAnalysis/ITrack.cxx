#include "ITrack.h"
#include "TVector3.h"
#include "TLorentzVector.h"

ClassImp(ITrack)

//___________________
ITrack::ITrack(){
  mTrackID = -99;
  mNhits = -99;
  mMomentum.SetXYZM(-99,-99,-99,-99);
  mTopologyMap0 = 999;
  mTopologyMap1 = 999;
}
//_________________________
ITrack::~ITrack(){
  mTrackID = -99;
  mNhits = -99;
  mMomentum.SetXYZM(-99,-99,-99,-99);
  mTopologyMap0 = 999;
  mTopologyMap1 = 999;
}
//_________________________
void ITrack::SetTopologyMap(int a, unsigned long long topMap){
  if(0 == a) mTopologyMap0 = topMap;
  else if (1 == a) mTopologyMap1 = topMap;
  else{
    printf("WRONG ASSIGNMENT OF TopologyMap integer, nothing saved!!! \n");
  }
}
//_________________________
unsigned long long ITrack::TopologyMap(int a){
  if(0 == a) return mTopologyMap0;
  else if (1 == a) return mTopologyMap1;
  else{
    printf("WRONG ASSIGNMENT OF TopologyMap integer!!! \n");
    return 999;
  }
}
