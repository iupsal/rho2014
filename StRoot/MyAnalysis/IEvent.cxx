#include "IEvent.h"
#include "stdio.h"

#include "ITrack.h"
#include "ILambda.h"


ClassImp(IEvent)

//_________________
IEvent::IEvent(){

  fLambdaCollection = new TClonesArray("ILambda");
  fAntiLambdaCollection = new TClonesArray("ILambda");

  mRunNumber = -99;
  mRefMult = -99;
  mCentralityID9 = -99;
  mCentralityID16 = -99;
  mPrimaryVertex.SetXYZ(-99.0,-99.0,-99.0);
  mBfield = -99.;
  mVpdVz = -99.;
}
//__________________
IEvent::~IEvent(){
  ClearEvent();
  if (fLambdaCollection)     {delete fLambdaCollection;     fLambdaCollection=0;}
  if (fAntiLambdaCollection)    {delete fAntiLambdaCollection;    fAntiLambdaCollection=0;}
}
//_________________
void IEvent::Init(){
  //SetUpPointers();  // allows the Lambdas to point to their daughters, which are sitting in the TClonesArray of this DEvent
}
//_________________
void IEvent::ClearEvent(){

  mRunNumber = -99;
  mRefMult = -99;
  mCentralityID9 = -99;
  mCentralityID16 = -99;
  mPrimaryVertex.SetXYZ(-99.0,-99.0,-99.0);
  mBfield = -99.;
  mVpdVz = -99.;

  fLambdaCollection->Delete();
  fAntiLambdaCollection->Delete();
}
//_________________
void IEvent::AddLambdaVertex(TClonesArray* ca, ILambda* dt){
  Int_t i = (Int_t)(ca->GetLast());
  i++;
  TClonesArray &ar = *ca;
//  new(ar[i]) ILambda(*dt);
  ar[i] = new ILambda(*dt);
  ((ILambda*)(ca->At(i)))->SetIndex(i);
}
//_________________
ILambda* IEvent::GetLambdaVertex(TClonesArray* ca, Int_t index) {
  ILambda* t = (ILambda*)(ca->At(index));
  if (!t) {
    Error("GetLambda","no Lambda found");
    return NULL;
  }
  else {
    return t;
  }
}
//__________________________
ILambda* IEvent::GetLambda(Int_t index){
  return GetLambdaVertex(fLambdaCollection,index);
}
//__________________________
void IEvent::AddLambda(ILambda* lam){
  AddLambdaVertex(fLambdaCollection,lam);
}

//__________________________
ILambda* IEvent::GetAntiLambda(Int_t index){
  return GetLambdaVertex(fAntiLambdaCollection,index);
}
//__________________________
void IEvent::AddAntiLambda(ILambda* lam){
  AddLambdaVertex(fAntiLambdaCollection,lam);
}
//__________________________


