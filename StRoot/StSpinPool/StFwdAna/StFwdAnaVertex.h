/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpoe of this class is to check and set the vertex information and fill histograms related to the vertex QA

  DESCRIPTION
  Grabs the primary vertex information first from the Mudst, then VPD, then EPD vertex from #StMuFcsAnaEpdQaAndVert, and finally BBC vertex. It will then save all the vertex information to #FcsEventInfo and set #StMuFcsAnaData::mFoundVertex and #StMuFcsAnaData::mUseVertex accordingly based on the vertex information it found. The #StMuFcsAnaData::mUseVertex is the z vertex that can be used for analysis.

  LOG
  @[January 15, 2026] > First instance where relevant functionality was copied from #StMuFcsTreeMaker
  @[July 1, 2026] > Changed name from StMuFcsAnaVertex to StFwdAnaVertex
  @[July 24, 2026] > Implemented #LoadGraphs() and #FillGraphs() which is for doing run by run qa on the vertex histograms  
  @[August 11, 2026] > Adjusted header loading to use minimum required
*/


#ifndef STFWDANA_STFWDANAVERTEX_HH
#define STFWDANA_STFWDANAVERTEX_HH

#include "StEvent/StTriggerData.h"
#include "StEvent/StTriggerId.h"

#include "StFwdAnaVirtual.h"

class StFwdAnaVertex : public StFwdAnaVirtual
{
public:
  StFwdAnaVertex();
  ~StFwdAnaVertex();

  virtual UInt_t LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata);
  virtual Int_t DoMake(StFwdAnaData* anadata);
  virtual UInt_t LoadGraphs(TFile* file, TObjArray* graphs, StFwdAnaData* anadata );
  virtual void FillGraphs(Int_t irun);

  void DrawVertex(TCanvas* canv, const char* savename) const;
  void DrawVertexCorrelation(TCanvas* canv, const char* savename) const;
  void DrawVertexCorrelationNoZdc(TCanvas* canv, const char* savename) const;

  void DrawVertexGraphs(TCanvas* canv, const char* savename="testGraphVertex.png") const;
  
protected:
  TH1* mH2F_VertexPrim_yVx = 0;       ///< Vertex histogram from primary vertex y vs. x
  TH1* mH1F_VertexPrimZ = 0;          ///< Vertex histograms from Primary Vertex
  TH1* mH1F_VertexVpd = 0;            ///< Vertex histograms from VPD
  TH1* mH1F_VertexBbc = 0;            ///< Vertex histograms from BBC
  TH1* mH1F_BbcTimeDiff = 0;          ///< BBC Time difference used to compute the vertex
  TH1* mH1F_VertexZdc = 0;            ///< Vertex from ZDC
  TH1* mH1F_VertexEpd = 0;            ///< Vertex from EPD

  TH1* mH2F_VertexZ_vpdVprim = 0;     ///< Correlation histogram between computed VPD z vertex vs. primary vertex z
  TH1* mH2F_VertexZ_epdVprim = 0;     ///< Correlation histogram between computed EPD z vertex vs. primary vertex z
  TH1* mH2F_VertexZ_bbcVprim = 0;     ///< Correlation histogram between computed BBC z vertex vs. primary vertex z
  TH1* mH2F_VertexZ_vpdVepd = 0;      ///< Correlation histogram between VPD z vertex vs. computed EPD z vertex
  TH1* mH2F_VertexZ_zdcVepd = 0;      ///< Correlation histogram between ZDC z vertex vs. computed EPD z vertex
  TH1* mH2F_VertexZ_bbcVepd = 0;      ///< Correlation histogram between BBC z vertex vs. computed EPD z vertex
  TH1* mH2F_VertexZ_vpdVbbc = 0;      ///< Correlation histogram between VPD z vertex vs. the BBC z vertex
  TH1* mH2F_VertexZ_vpdVzdc = 0;      ///< Correlation histogram between VPD z vertex vs. the ZDC z vertex
  TH1* mH2F_VertexZ_zdcVbbc = 0;      ///< Correlation histogram between ZDC z vertex vs. the BBC z vertex
  
  TH1* mH2F_foundVvertex = 0;           ///< found vertex bit vs. Vertex

  TGraphErrors* mGE_VertexPrimZ = 0;     ///< Graph for Mean primary Z vertex and Err as RMS vs. Run Index
  TGraphErrors* mGE_VertexVpd = 0;       ///< Graph for Mean VPD vertex and Err as RMS vs. Run Index
  TGraphErrors* mGE_VertexBbc = 0;       ///< Graph for Mean BBC vertex and Err as RMS vs. Run Index
  TGraphErrors* mGE_VertexEpd = 0;       ///< Graph for Mean EPD vertex and Err as RMS vs. Run Index
  
  ClassDef(StFwdAnaVertex,1)
};

#endif

