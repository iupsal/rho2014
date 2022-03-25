// Tiny DST class for Lambda mixing

#ifndef IEVENT_H
#define IEVENT_H

#include "TVector3.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <vector>


class ITrack;
class ILambda;

/*
  Store simple event summary information
*/

class IEvent : public TObject {

 protected:

  Int_t mRunNumber;
  Int_t mRefMult;
  Int_t mCentralityID9;
  Int_t mCentralityID16;
  TVector3 mPrimaryVertex;
  Float_t mBfield;
  Float_t mVpdVz;

  TClonesArray* fLambdaCollection;
  TClonesArray* fAntiLambdaCollection;

  void AddLambdaVertex(TClonesArray*, ILambda*);
  ILambda* GetLambdaVertex(TClonesArray*, Int_t index);


 public:

  IEvent();
  ~IEvent();

  void ClearEvent();
  void Init();


  Int_t   NLambda(){return fLambdaCollection->GetEntriesFast();}
  ILambda*    GetLambda(Int_t index);
  void    AddLambda(ILambda* lam);

  Int_t   NAntiLambda(){return fAntiLambdaCollection->GetEntriesFast();}
  ILambda*    GetAntiLambda(Int_t index);
  void    AddAntiLambda(ILambda* lam);

  Int_t GetRunNumber(){return mRunNumber;}
  Int_t GetRefMult(){return mRefMult;}
  Int_t GetCentralityID9(){return mCentralityID9;}
  Int_t GetCentralityID16(){return mCentralityID16;}
  TVector3 PrimaryVertex() {return mPrimaryVertex;}
  Float_t Bfield(){return mBfield;}
  Float_t GetVpdVz(){return mVpdVz;}

  void SetRunNumber(Int_t rn){mRunNumber=rn;}
  void SetRefMult(Int_t rm){mRefMult = rm;}
  void SetCentralityID9(Int_t cent){mCentralityID9=cent;}
  void SetCentralityID16(Int_t cent){mCentralityID16=cent;}
  void SetPrimaryVertex(TVector3 pv){mPrimaryVertex = pv;}
  void SetBfield(Double_t f){mBfield=f;}
  void SetVpdVz(Double_t f){mVpdVz=f;}

  ClassDef(IEvent,1)  // my event

};

#endif


