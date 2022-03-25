// Tiny DST class for Lambda mixing

#ifndef ITRACK_H
#define ITRACK_H

#include "TLorentzVector.h"
class IEvent;

class ITrack : public TObject {

 protected:

  Int_t mTrackID;
  Int_t mNhits;
  TLorentzVector mMomentum; //use momentum at the decay point of the Lambda, not the track momentum
  unsigned long long mTopologyMap0;
  unsigned long long mTopologyMap1;


 public:

  ITrack();
  ~ITrack();

  //  void SetParentEvent(DEvent* p){mParentEvent=p;}
  void SetNhits(Int_t n){mNhits=n;}
  void SetTrackID(Int_t id){mTrackID=id;}
  void SetMomentum(TLorentzVector a){mMomentum = a;}
  void SetTopologyMap(int a, unsigned long long topMap);

  Int_t TrackID(){return mTrackID;}
  Int_t Nhits(){return mNhits;}
  TLorentzVector Momentum(){return mMomentum;}
  unsigned long long TopologyMap(int a);

  ClassDef(ITrack,1)  // my track
};

#endif
