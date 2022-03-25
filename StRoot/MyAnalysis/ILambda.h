// Tiny DST class for Lambda mixing

#ifndef ILAMBDA_H
#define ILAMBDA_H

#include "ITrack.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class IEvent;    // forward declaration
class ILambda : public TObject{

 protected:

  TVector3 mThreeMomentum;
  TVector3 mDecayPoint;
  Float_t mCalcMass;
  Int_t mIndex;	//necessary for adding list of Lambdas to Event


  ITrack mPosDaughter;
  ITrack mNegDaughter;

 public:
  ILambda();
  ~ILambda();

  void SetThreeMomentum(TVector3 m){mThreeMomentum=m;}
  void SetDecayPoint(TVector3 m){mDecayPoint=m;}
  void SetCalcMass(Float_t b){mCalcMass=b;}
  void SetPosDaughter(ITrack p){mPosDaughter=p;}
  void SetNegDaughter(ITrack p){mNegDaughter=p;}
  void SetIndex(Int_t i){mIndex=i;}

  TVector3 ThreeMomentum(){return mThreeMomentum;}
  TVector3 DecayPoint(){return mDecayPoint;}
  Float_t CalcMass(){return mCalcMass;}
  ITrack PosDaughter(){return mPosDaughter;}
  ITrack NegDaughter(){return mNegDaughter;}
  Float_t GetIndex(){return mIndex;}

  TLorentzVector RealMomentum(); //using PDG mass
  TLorentzVector CalcMomentum(); //using mass calculated from daughters

  friend class IEvent;  // so that IEvent can set the daughter keys without going through interface functions.

  ClassDef(ILambda,1)  // my Lambda

};
    
#endif
