#include "StEvent/StEventTypes.h"
#include "StEvent/StFstEvtCollection.h"
#include "StEvent/StFstRawHit.h"

#include "StFwdAnaFstRun22Qa.h"

StFwdAnaFstRun22Qa::StFwdAnaFstRun22Qa()
{
  memset(mH2S_FstRawHitStrip_rVphi,  0, sizeof(mH2S_FstRawHitStrip_rVphi));
  memset(mH2S_FstHitStripMean_rVphi, 0, sizeof(mH2S_FstHitStripMean_rVphi));
  memset(mH2S_FstHit_rVphi,          0, sizeof(mH2S_FstHit_rVphi));
  memset(mH2S_FstHit_apvVgeoid,      0, sizeof(mH2S_FstHit_apvVgeoid));
  memset(mH2S_FstHitGlobal_yVx,      0, sizeof(mH2S_FstHitGlobal_yVx));
  memset(mH2S_FstHitGlobal_rVphi,    0, sizeof(mH2S_FstHitGlobal_rVphi));

  memset(mH2S_FstRawHit_adcVgeoid,   0, sizeof(mH2S_FstRawHit_adcVgeoid));
}

StFwdAnaFstRun22Qa::~StFwdAnaFstRun22Qa()
{}

UInt_t StFwdAnaFstRun22Qa::LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata)
{
  UInt_t nloaded = 0;
  std::stringstream ss_histname;
  std::stringstream ss_histtitle;
  for( int iwedge=0; iwedge<kFstNumWedges; ++iwedge ){
    for( int isensor=0; isensor<kFstNumSensorsPerWedge; ++isensor ){
      ss_histname.str("");
      ss_histtitle.str("");
      if( mRawHitOn ){
	ss_histname << "H2S_FstRawHitStrip_rVphi_"<< iwedge*3+isensor+1;
	ss_histtitle << "Raw Hit rstrip vs. #phi strip: Wedge "<<iwedge+1<< " Sensor " <<  isensor << ";#phi Strip;r Strip";
	nloaded += histman->AddH2S(file,mH2S_FstRawHitStrip_rVphi[iwedge*3+isensor],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumPhiSegPerWedge,0,kFstNumPhiSegPerWedge, kFstNumRStripsPerWedge,0,kFstNumRStripsPerWedge);
      }
      //ss_histname.str("");
      //ss_histtitle.str("");
      //sprintf(buffer,"numOfRawHitsVsEventId_Sensor%d", iwedge*3+isensor+1);
      //sprintf(histname, "Number of raw hits vs. EventID: Wedge %d Sensor %d", iwedge+1, isensor);
      //numOfRawHits_EventId[iwedge*3+isensor] = new TProfile(buffer, histname, 10000, 0, 10000);
      //numOfRawHits_EventId[iwedge*3+isensor]->GetXaxis()->SetTitle("EventID (Time)");
      //numOfRawHits_EventId[iwedge*3+isensor]->GetYaxis()->SetTitle("<rawhits>");

      if( mHitOn ){
	ss_histname.str("");
	ss_histtitle.str("");
	ss_histname << "H2S_FstHitStripMean_rVphi_"<<iwedge*3+isensor+1;
	ss_histtitle << "Hit mean rstrip vs. #phi strip: Wedge "<<iwedge+1<< " Sensor " << isensor <<";Mean #phi strip;Mean r strip";
	nloaded += histman->AddH2S(file,mH2S_FstHitStripMean_rVphi[iwedge*3+isensor],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumPhiSegPerWedge,0,kFstNumPhiSegPerWedge, kFstNumRStripsPerWedge,0,kFstNumRStripsPerWedge);
      }
    }
  }

  if( mHitOn ){
    for( int idisk = 0; idisk<kFstNumDisk; ++idisk ){
      ss_histname.str("");
      ss_histtitle.str("");
      ss_histname << "H2S_FstHit_rVphi_Disk"<<idisk+1;
      ss_histtitle << "FST Hit Map for Disk "<<idisk+1 << " r vs #phi;#phi strip;r strip";
      nloaded += histman->AddH2S(file,mH2S_FstHit_rVphi[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumPhiSegPerWedge*kFstNumWedgePerDisk,0,kFstNumPhiSegPerWedge*kFstNumWedgePerDisk, kFstNumRStripsPerWedge,0,kFstNumRStripsPerWedge);

      ss_histname.str("");
      ss_histtitle.str("");
      ss_histname << "H2S_FstHit_apvVgeoid_Disk"<<idisk;
      ss_histtitle << "FST hit map in APV vs Wedge Idx for Disk "<<idisk <<";WedgeIdx [wedge-(disk-1)*WedgePerDisk+1];APV Geometry ID";
      nloaded += histman->AddH2S(file,mH2S_FstHit_apvVgeoid[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumWedgePerDisk,1,kFstNumWedgePerDisk+1, kFstApvsPerWedge,0,kFstApvsPerWedge);

      ss_histname.str("");
      ss_histtitle.str("");
      ss_histname << "H2S_FstHitGlobal_yVx_Disk"<<idisk;
      ss_histtitle << "FST Hit Map in STAR coordinates y vs. x for Disk "<<idisk+1 <<";STAR X (cm);STAR Y(cm)";
      nloaded += histman->AddH2S(file,mH2S_FstHitGlobal_yVx[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(), 140,-35,35, 140,-35,35);

      ss_histname.str("");
      ss_histtitle.str("");
      ss_histname << "H2S_FstHitGlobal_rVphi_Disk"<<idisk;
      ss_histtitle << "FST Hit Map in STAR coordinates r vs. #phi for Disk "<<idisk+1 <<";STAR #phi (rad.);STAR r (cm)";
      nloaded += histman->AddH2S(file,mH2S_FstHitGlobal_rVphi[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(),kFstNumPhiSegPerWedge*kFstNumWedgePerDisk/8, -TMath::Pi(), TMath::Pi(), kFstNumRStripsPerWedge, 5, 28);
    }
  }

  if( mRawHitOn ){
    for( int itb=0; itb<kFstNumTimeBins; ++itb ){
      ss_histname.str("");
      ss_histtitle.str("");
      ss_histname << "H2S_FstRawHit_adcVgeoid_tb"<<itb;
      ss_histtitle << "ADC of raw hits at timebin "<<itb << "vs. channel geometry ID;Channel Geometry ID;ADC of Raw Hits";
      nloaded += histman->AddH2S(file,mH2S_FstRawHit_adcVgeoid[itb],ss_histname.str().c_str(),ss_histtitle.str().c_str(),288, 0, 36864, 512, 0, kFstMaxAdc);    
    }
  }
    
  if( mRawHitOn ){
    nloaded += histman->AddH2S(file,mH2S_FstRawHit_adcerrVgeoid,"H2S_FstRawHit_adcerrVgeoid","RMS noise of raw hits vs. channel geometry ID;Channel Geometry ID;RMS noise of Raw Hits",288, 0, 36864, 128, 0, 64);
    nloaded += histman->AddH2S(file,mH2S_Fst_nrawhitsVid,"H2S_Fst_nrawhitsVid","Number of Raw Hits vs. Sensor ID;Sensor ID;Number of Raw Hits", kFstNumSensors,0,kFstNumSensors, 128,0,128);
  }
  if( mHitOn ){
    nloaded += histman->AddH2S(file,mH2S_FstRawHit_maxtbVapv,"H2S_FstHit_maxtbVapv","Max time bin of hit vs. APV ID;APV ID [48*(RDO-1)+16*ARM+APV];Max Time Bin Index",kFstNumApvs, 0, kFstNumApvs, kFstNumTimeBins, 0, kFstNumTimeBins);
    nloaded += histman->AddH2S(file,mH2S_FstHit_adcVid,"H2S_FstHit_adcVid","ADC of hits vs. Sensor ID;Sensor ID;ADC of Hits",kFstNumSensors,0,kFstNumSensors, 512,0,kFstMaxAdc);
    nloaded += histman->AddH2S(file,mH2S_FstHit_adcerrVid,"H2S_FstHit_adcerrVid","RMS noise of hits vs Sensor ID;Sensor ID;RMS noise of Hits", kFstNumSensors,0,kFstNumSensors, 128,0,64);
    nloaded += histman->AddH2S(file,mH2S_FstHit_maxtbVid,"H2S_FstHit_maxtbVid","Max time bin of hits vs. Sensor ID;Sensor ID;Max Time Bin Index",kFstNumSensors,0,kFstNumSensors, kFstNumTimeBins,0,kFstNumTimeBins);

    nloaded += histman->AddH2S(file,mH2S_Fst_nhitsVid,"H2S_Fst_nhitsVid","The number of hits vs. Sensor ID;Sensor ID,Number of Hits", kFstNumSensors,0,kFstNumSensors, 128, 0, 128);
    nloaded += histman->AddH2S(file,mH2S_FstHit_nrawhitsVid,"H2S_FstHit_nrawhitsVid","Number of Raw Hits in a Cluster(Hit) vs. Sensor ID;Sensor ID;Number of Raw Hits in Cluster",kFstNumSensors,0,kFstNumSensors, 20,0,20);
    nloaded += histman->AddH2S(file,mH2S_FstHit_nrawhitsrVid,"H2S_FstHit_nrawhitsrVid","Hit size in R of raw hits vs Sensor ID;Sensor ID;Hit Size in R of Raw Hits",kFstNumSensors,0,kFstNumSensors, 20,0,20);
    nloaded += histman->AddH2S(file,mH2S_FstHit_nrawhitsphiVid,"H2S_FstHit_nrawhitsphiVid","Hit size in #phi of raw hits vs. Sensor ID;Sensor ID;Hit Size in #phi of Raw Hits", kFstNumSensors,0,kFstNumSensors, 20,0,20);
  }

  return nloaded;
}


Int_t StFwdAnaFstRun22Qa::DoMake(StFwdAnaData* anadata)
{
  //StMuFstCollection* mufstcoll = anadata->muFstColl();
  TClonesArray* tc_mufstraw = anadata->muFstRawHitColl();  //Collection exists but no raw hits stored
  TClonesArray* tc_mufsthit = anadata->muFstHitColl();
  if( tc_mufstraw==0 ){
    LOG_WARN << "StFwdAnaFstRun22Qa::No FST Raw Hit Collection" << endm;
    return kStWarn;
  }
  if( tc_mufsthit==0 ){
    LOG_WARN << "StFwdAnaFstRun22Qa::No FST Hit Collection" << endm;
    return kStWarn;
  }

  //loop over raw hits
  unsigned int rawcounter[kFstNumSensors]; //raw hit multiplicity per sensor per event
  memset(rawcounter,0,sizeof(rawcounter));
  if( mRawHitOn ){
    //StFstCollection* fstcoll = anadata->fstColl();
    //StFstEvtCollection* fstevtcoll = anadata->event()->fstEvtCollection();
    //std::cout << "|fstevtcoll:"<<fstevtcoll << std::endl;
    //StSPtrVecFstRawHit& rawhits = fstevtcoll->rawHits();

    //std::cout << "|tc_mufstraw:"<<tc_mufstraw->GetEntriesFast() <<"|raw:"<<rawhits.size()<< "|tc_mufsthit:"<<tc_mufsthit->GetEntriesFast() << std::endl;
    for( int iraw=0; iraw<tc_mufstraw->GetEntriesFast(); ++iraw ){
      //for( unsigned int iraw = 0; iraw<rawhits.size(); ++iraw ){
      StMuFstRawHit* rawhit = (StMuFstRawHit*)tc_mufstraw->At(iraw);
      //for(int wedgeIdx=0; wedgeIdx<kFstNumWedges; ++wedgeIdx ){
      //StFstRawHitCollection *rawhitcoll = fstcoll->getRawHitCollection( wedgeIdx );
      //if( rawhitcoll ){
      //std::vector<StFstRawHit*>& rawhitvec = rawhitcoll->getRawHitVec();
      //std::vector< StFstRawHit* >::iterator rawhititer;
      //for( rawhititer = rawhitvec.begin(); rawhititer != rawhitvec.end(); ++rawhititer ){
      //StFstRawHit* rawhit = *rawhititer;
    
      //StFstRawHit* rawhit = (StFstRawHit*)rawhits.at(iraw);
      if( rawhit==0 ){ continue; }
      //Cast to unsigned int to avoid bad char to int conversions
      unsigned int wedge = rawhit->getWedge();
      unsigned int sensor = rawhit->getSensor();
      unsigned int maxTimeBin = rawhit->getMaxTimeBin();
      unsigned int sensorId = (wedge-1)*static_cast<unsigned int>(kFstNumSensorsPerWedge) + sensor;
      ++rawcounter[sensorId];

      for( unsigned char timeBin = 0; timeBin < kFstNumTimeBins; ++timeBin ) {
	mH2S_FstRawHit_adcVgeoid[timeBin]->Fill(rawhit->getGeoId(), (int)rawhit->getCharge( timeBin ));
      }
      mH2S_FstRawHit_adcerrVgeoid->Fill(rawhit->getGeoId(), (int)(rawhit->getChargeErr( maxTimeBin )+0.5));
      mH2S_FstRawHitStrip_rVphi[sensorId]->Fill((int)rawhit->getPhiStrip(), (int)rawhit->getRStrip());
      mH2S_FstRawHit_maxtbVapv->Fill(((int)rawhit->getRdo()-1)*48+(int)rawhit->getArm()*16+(int)rawhit->getApv(), (int)maxTimeBin);
    }
  }

  unsigned int hitcounter[kFstNumSensors];  //hit multiplicity per sensor
  if( mHitOn ){
    memset(hitcounter,0,sizeof(hitcounter));
    for( int ihit = 0; ihit<tc_mufsthit->GetEntriesFast(); ++ihit ){
      StMuFstHit* hit = (StMuFstHit*)tc_mufsthit->At(ihit);
      if( hit!=0 ){
	const TVector3& hitpos = hit->xyz();
	int sensorIdxTemp = ((int)hit->getWedge()-1)*kFstNumSensorsPerWedge + (int)hit->getSensor(); // 0-107
	++hitcounter[sensorIdxTemp];
	int diskIdxTemp = ((int)hit->getWedge()-1)/kFstNumWedgePerDisk + 1; // 1-3
	int wedgeIdxTemp = (int)hit->getWedge() - (diskIdxTemp-1)*kFstNumWedgePerDisk; // 1-12
	int phiIdxTemp = (wedgeIdxTemp-1)*kFstNumPhiSegPerWedge+(int)(hit->getMeanPhiStrip()+0.5);
	int rIdxTemp = (int)(hit->getMeanRStrip()+0.5);

	mH2S_FstHitStripMean_rVphi[sensorIdxTemp]->Fill((int)(hit->getMeanPhiStrip()+0.5), (int)(hit->getMeanRStrip()+0.5));
	mH2S_FstHit_rVphi[diskIdxTemp-1]->Fill(phiIdxTemp, rIdxTemp);
	mH2S_FstHit_apvVgeoid[diskIdxTemp-1]->Fill(wedgeIdxTemp, (int)hit->getApv()+(wedgeIdxTemp%2-1)*kFstApvsPerWedge);
	mH2S_FstHitGlobal_yVx[diskIdxTemp-1]->Fill((float)hitpos.X(), (float)hitpos.Y());
	mH2S_FstHitGlobal_rVphi[diskIdxTemp-1]->Fill((float)hitpos.Phi(), (float)hitpos.Perp());

	mH2S_FstHit_adcVid->Fill(sensorIdxTemp, (int)hit->getCharge());
	mH2S_FstHit_adcerrVid->Fill(sensorIdxTemp, (int)(hit->getChargeErr()+0.5));
	mH2S_FstHit_maxtbVid->Fill(sensorIdxTemp, (int)hit->getMaxTimeBin());
	mH2S_FstHit_nrawhitsVid->Fill(sensorIdxTemp, (int)hit->getNRawHits());
	mH2S_FstHit_nrawhitsrVid->Fill(sensorIdxTemp, (int)hit->getNRawHitsR());
	mH2S_FstHit_nrawhitsphiVid->Fill(sensorIdxTemp, (int)hit->getNRawHitsPhi());
      }
    }
  }

  for( int iS=0; iS<kFstNumSensors; ++iS ){
    if( mRawHitOn){ mH2S_Fst_nrawhitsVid->Fill(iS, rawcounter[iS]); }
    if( mHitOn ){ mH2S_Fst_nhitsVid->Fill(iS, hitcounter[iS]); }
    //numOfRawHits_EventId[iS]->Fill((int)eventPtr->id()/100+1, counter[iS]);
  }
  
  return kStOk;
}

void StFwdAnaFstRun22Qa::Paint_FstRawHitStrip_rVphi(TCanvas* canvas, int disk, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(6,6); //kFstNumSensor=3*36
  for( int i=36*disk; i<36*disk+36; ++i ){
    TVirtualPad* pad = canvas->cd(i-36*disk+1);
    pad->SetLogz();
    if( mH2S_FstRawHitStrip_rVphi[i] ){ mH2S_FstRawHitStrip_rVphi[i]->Draw("colz"); }
  }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHitStripMean_rVphi(TCanvas* canvas, int disk, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(6,6); //kFstNumSensor=3*36
  for( int i=36*disk; i<36*disk+36; i++ ){
    TVirtualPad* pad = canvas->cd(i-36*disk+1);
    pad->SetLogz();
    if( mH2S_FstHitStripMean_rVphi[i] ){ mH2S_FstHitStripMean_rVphi[i]->Draw("colz"); }
  }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_rVphi(TCanvas* canvas, int disk, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_rVphi[disk] ){ mH2S_FstHit_rVphi[disk]->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_apvVgeoid(TCanvas* canvas, int disk, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_apvVgeoid[disk] ){ mH2S_FstHit_apvVgeoid[disk]->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHitGlobal_yVx(TCanvas* canvas, int disk, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHitGlobal_yVx[disk] ){ mH2S_FstHitGlobal_yVx[disk]->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHitGlobal_rVphi(TCanvas* canvas, int disk, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHitGlobal_rVphi[disk] ){ mH2S_FstHitGlobal_rVphi[disk]->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstRawHit_adcVgeoid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(3,3);
  for( int i=0; i<kFstNumTimeBins; i++ ){
    TVirtualPad* pad = canvas->cd(i+1);
    pad->SetLogz();
    if( mH2S_FstRawHit_adcVgeoid[i] ){ mH2S_FstRawHit_adcVgeoid[i]->Draw("colz"); }
  }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstRawHit_adcerrVgeoid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstRawHit_adcerrVgeoid ){ mH2S_FstRawHit_adcerrVgeoid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstRawHit_maxtbVapv(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstRawHit_maxtbVapv ){ mH2S_FstRawHit_maxtbVapv->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_adcVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_adcVid ){ mH2S_FstHit_adcVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_adcerrVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_adcerrVid ){ mH2S_FstHit_adcerrVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_maxtbVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_maxtbVid ){ mH2S_FstHit_maxtbVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_Fst_nrawhitsVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_Fst_nrawhitsVid ){ mH2S_Fst_nrawhitsVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_Fst_nhitsVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_Fst_nhitsVid ){ mH2S_Fst_nhitsVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_nrawhitsVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_nrawhitsVid ){ mH2S_FstHit_nrawhitsVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_nrawhitsrVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_nrawhitsrVid ){ mH2S_FstHit_nrawhitsrVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::Paint_FstHit_nrawhitsphiVid(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  TVirtualPad* pad = canvas->cd();
  pad->SetLogz();
  if( mH2S_FstHit_nrawhitsphiVid ){ mH2S_FstHit_nrawhitsphiVid->Draw("colz"); }
  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::PaintFstRawHitQa(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(2,2);

  canvas->cd(1);
  if( mH2S_FstRawHit_adcerrVgeoid ){ mH2S_FstRawHit_adcerrVgeoid->Draw("colz"); }
  canvas->cd(2);
  if( mH2S_FstRawHit_maxtbVapv ){ mH2S_FstRawHit_maxtbVapv->Draw("colz"); }
  canvas->cd(3);
  if( mH2S_Fst_nrawhitsVid ){ mH2S_Fst_nrawhitsVid->Draw("colz"); }

  canvas->Print(savename);
}

void StFwdAnaFstRun22Qa::PaintFstHitQa(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(3,3);

  canvas->cd(1);
  if( mH2S_FstHit_adcVid ){ mH2S_FstHit_adcVid->Draw("colz"); }
  canvas->cd(2);
  if( mH2S_FstHit_adcerrVid ){ mH2S_FstHit_adcerrVid->Draw("colz"); }
  canvas->cd(3);
  if( mH2S_FstHit_maxtbVid ){ mH2S_FstHit_maxtbVid->Draw("colz"); }

  canvas->cd(4);
  if( mH2S_Fst_nhitsVid ){ mH2S_Fst_nhitsVid->Draw("colz"); }

  canvas->cd(5);
  if( mH2S_FstHit_nrawhitsVid ){ mH2S_FstHit_nrawhitsVid->Draw("colz"); }
  canvas->cd(6);
  if( mH2S_FstHit_nrawhitsrVid ){ mH2S_FstHit_nrawhitsrVid->Draw("colz"); }
  canvas->cd(7);
  if( mH2S_FstHit_nrawhitsphiVid ){ mH2S_FstHit_nrawhitsphiVid->Draw("colz"); }

  canvas->Print(savename);
}

