#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;


StChain *chain;


void Analysis(Int_t nEvents, TString InputFileList, TString OutputDir, TString JobIdName, int ShiftWeek = 0, int CalibMode = 2, bool EpMode = 0)
{
  // ShiftWeek is the week the data was taken. First week = 1, etc.
  // Calibmode is 0, 1, 2 for Gain, Recenter, Shift EP corrections
  // EpMode is for running the EP corrections, true means run only EP corrections

  if ( nEvents == 0 )  nEvents = 100000000 ;       // Take all events in nFiles if nEvents = 0



  //root -l Analysis.C\(10,\"tinyfile200.list\",\"./\",\"LocalOut\"\)
  //root -l Analysis.C\(50,\"root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_low_2014/ReversedFullField/P16id.SL18f/2014/095/15095020/st_physics_15095020_raw_1000013.picoDst.root\",\"./\",\"LocalOut\",4,2,0\)
  

  // SL18f
  //Load all the System libraries
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("MyAnalysis.so");

  chain                            =  new StChain();
  StPicoDstMaker*    picoMaker     =  new StPicoDstMaker(2,InputFileList,"picoDst");
  MyAnalysisMaker*   AnalysisCode  =  new MyAnalysisMaker(picoMaker,JobIdName,ShiftWeek,CalibMode,EpMode) ;


  if( chain->Init()==kStErr ){ 
    cout<<"chain->Init();"<<endl;
    return;
  }

  int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++){

    if(i%1000==0)
    cout << "Working on eventNumber " << i << endl;
    
    chain->Clear();
    int iret = chain->Make(i);
    
    if (iret) { cout << "Bad return code!" << iret << endl; break;}

    total++;
  }
  
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
  
  delete AnalysisCode;
  delete picoMaker;
  delete chain;
  
}
