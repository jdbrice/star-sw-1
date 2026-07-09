void pico3(int n=0, int opt=11, int run=23030048, int set=0) {
  //void pico3(int n=0, int opt=11, int run=1, int set=0) {
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");
  gSystem->Load("StEvent");
  gSystem->Load("StDbLib");
  gSystem->Load("StFcsDbMaker");
  gSystem->Load("StPicoEvent");

  gROOT->ProcessLine(".L picoMatch.C++");
  gROOT->ProcessLine(".L picoDilepton.C++");
  gROOT->ProcessLine(Form(".x readpico3.C++(%d,%d,%d,%d)",n,opt,run,set));
}
