#include "ILambda.h"
#include "IEvent.h"
#include "TVector3.h"
#include "TLorentzVector.h"

ClassImp(ILambda)

//______________
ILambda::ILambda(){
  mThreeMomentum.SetXYZ(-99.0,-99.0,-99.0);
  mDecayPoint.SetXYZ(-99.0,-99.0,-99.0);
  mCalcMass = 99.;
  mIndex = -99;
}
//__________________________
ILambda::~ILambda(){
  mThreeMomentum.SetXYZ(-99.0,-99.0,-99.0);
  mDecayPoint.SetXYZ(-99.0,-99.0,-99.0);
  mCalcMass = 99.;
  mIndex = -99;
}
//__________________________
TLorentzVector ILambda::RealMomentum(){
  float mass = 1.115683;
  TLorentzVector RealMomentum;
  RealMomentum.SetXYZM(mThreeMomentum.X(),mThreeMomentum.Y(),mThreeMomentum.Z(),mass);  
  return RealMomentum;
}
//__________________________
TLorentzVector ILambda::CalcMomentum(){
  float mass = mCalcMass;
  TLorentzVector CalcMomentum;
  CalcMomentum.SetXYZM(mThreeMomentum.X(),mThreeMomentum.Y(),mThreeMomentum.Z(),mass);  
  return CalcMomentum;
}
