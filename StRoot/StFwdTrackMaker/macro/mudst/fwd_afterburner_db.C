
//usr/bin/env root4star -l root -l -q  $0; exit $?
//usr/bin/env root4star -l -b -q $0'("'${1:-st_physics_23055058_raw_1500001.MuDst.root}'",'${2:-100}')'; exit $?
// that is a valid shebang to run script as executable, but with only one arg

#include <typeinfo.h>
#include <fstream>
#include <string>

// Fast fwd tracking without DB
// bool runDb = false;
// bool runFttChain = true;
// bool runFcsChain = false; 
// bool runFwdChain = true;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = false;
// bool runPico = true;

// For EPD QA only
// bool runDb = false;
// bool runFttChain = false;
// bool runFcsChain = true;
// bool runFwdChain = false;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = true;

// Tracking with FCS (with DB)
bool runDb = true;
bool runFttChain = true;
bool runFcsChain = true;
bool runFwdChain = true;
bool refillMuDst = false;
bool runFwdQa = false;
bool runFitQa = false;
bool runPico = true;
bool runFwdResidual = true;  // FST/FTT alignment residuals (StFwdResidualMaker); MuDst-only, no file written
// Unbiased (hit-removed) FST/FTT residuals (StFwdAlignmentMaker) -- see
// proposal_alignment_path.txt. Off by default: costs ~(num FST+FTT planes)
// extra refits per Global track, not meant for routine running. Completely
// independent of runFwdResidual/StFwdResidualMaker above -- separate maker,
// separate output file, doesn't touch StEvent/MuDst/PicoDst. Set at runtime
// via fwd_afterburner_db()'s enableAlignment parameter (see signature below),
// not by editing this default -- this is what condor batch jobs (runab) use
// to turn it on for a specific production campaign without a source edit.
bool runFwdAlignment = false;

// Memory Baseline
// bool runDb = false;
// bool runFttChain = true;
// bool runFcsChain = false;
// bool runFwdChain = false;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = false;
// bool runPico = true;

#include "StMemStat.h"

void loadLibs();
//void fwd_afterburner(const Char_t * fileList = "st_physics_23037002_raw_1000064.MuDst.root",size_t nEvents = 100, int debug=0){
// residualTrackType is no longer used -- StFwdResidualMaker now runs Global+BLC+BLCVtx
// together in one pass (see fwdResiduals[] below) -- kept as a parameter only so any
// existing 4-arg call sites don't break.
// extraFileList: optional path to a plain text file, one MuDst.root path per
// line, to Add() onto the same TChain as fileList -- a guaranteed-to-work way
// to combine multiple MuDst files for more statistics (StMuDstMaker's own
// constructor does NOT auto-detect a .list/.lis fileList argument the way
// some other STAR IO makers do -- tried that directly, it silently processed
// zero events instead of chaining). Empty (default) = old single-file
// behavior, unchanged.
// enableAlignment: turns on StFwdAlignmentMaker for THIS invocation only, no
// source edit/rebuild needed -- this is how a condor batch (runab) requests
// the expensive unbiased-residual ntuple for a specific production campaign
// while the compile-time default (runFwdAlignment above) stays off for
// routine/interactive use.
// applyFstWedgeAlignment: turns on FstWedgeAligner's per-wedge FST phi
// correction (StFwdHitLoader.h) for THIS invocation only -- see
// bugreport_StFstHitMaker.txt. Off by default so existing production
// behavior doesn't silently change; turn on to validate the correction or
// once it's ready to use routinely.
// magField: where the track fitter's magnetic field comes from. There is
// nothing to set per dataset -- mode 1 reads it out of the file.
//
//    1 = FROM DATA (default). Construct StarMagField with the scale implied by
//        the MuDst's own magneticField() stamp. Correct for field-on AND
//        field-off runs: measured -4.98702 kG for the 2022 ReversedFullField
//        production, +0.02994 kG for the 2022 zeroFieldAlignment runs.
//
//    0 = FORCE ZERO. genfit::ConstField(0,0,0) via TrackFitter:zeroB.
//        Diagnostic only.
//
//   -1 = LEGACY. Construct nothing, reproducing every afterburner production
//        before 2026-09-09. StBFChain::SetDbOptions normally creates
//        StarMagField from the bfc FieldOn/FieldOff option, but this macro
//        builds its chain by hand and never did; StarFieldAdaptor guards its
//        whole body with "if (StarMagField::Instance())" and so silently
//        returned B={0,0,0} at EVERY point -- FST, FTT and beyond, not just
//        inside the uniform-field shortcut. Keep for A/B against existing
//        output; do not use for new production.
//
// Note mode -1 also makes the STARField.h uniform-field shortcut unreachable,
// which is why restoring old behaviour needs no second switch in that header.
//
// At B=0 momentum is unconstrained (a straight track has no curvature), so pT
// is meaningless and the pT-based guards in FwdTracker.h (addFttHits /
// addFstHits / addEpdHits, "blown-up state") reject nearly every track. That
// is why modes 0 and -1 yield few or no FTT residuals.
void fwd_afterburner_db(const Char_t * fileList = "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_pp500_2022/ReversedFullField/P24ia/2022/108/23108014/st_fwd_23108014_raw_2000026.MuDst.root",size_t nEvents = 100, int debug=0, int residualTrackType=0, const char* extraFileList="", bool enableAlignment=false, bool applyFstWedgeAlignment=false, bool applyFstGapFix=false, bool applyFstMirror=false, int magField=1, bool applyFttXYMirror=false, int fttTimeCutMode=2, bool fieldFstConstBz=true, double fieldZMax=450., int fttTimeCutLow=-40, int fttTimeCutHigh=100, int fttDiagType=-1, bool fttNoAdd=false, bool fttDiagMix=false, bool useDbGeometry=true){
	cout << "FileList: " << fileList << endl;
	cout << "nEvents: " << nEvents << endl;
	cout << "enableAlignment: " << enableAlignment << endl;
	cout << "applyFstWedgeAlignment: " << applyFstWedgeAlignment << endl;
	cout << "applyFstGapFix: " << applyFstGapFix << endl;
	cout << "applyFstMirror: " << applyFstMirror << endl;
	cout << "magField: " << magField
	     << (magField==1 ? "  (from MuDst)" : (magField==0 ? "  (forced zero)" : "  (LEGACY: no field at all)")) << endl;
	cout << "applyFttXYMirror: " << applyFttXYMirror << endl;
	cout << "fieldFstConstBz: " << fieldFstConstBz << "  fieldZMax: " << fieldZMax << endl;
	cout << "fttDiagType: " << fttDiagType << (fttDiagType<0 ? "  (all track types)" : "  (one type only)")
	     << "  fttNoAdd: " << fttNoAdd << (fttNoAdd ? "  (sTGC hits NEVER added to tracks)" : "")
	     << "  fttDiagMix: " << fttDiagMix << (fttDiagMix ? "  (diagnostic vs PREVIOUS event)" : "") << endl;
	cout << "fttTimeCutMode: " << fttTimeCutMode
	     << (fttTimeCutMode==2 ? "  (calibrated time)" : (fttTimeCutMode==1 ? "  (AcceptAll)" : "  (other)"))
	     << "  window [" << fttTimeCutLow << "," << fttTimeCutHigh << "] dbcid ticks" << endl;
	runFwdAlignment = enableAlignment;

	// First load some shared libraries we need
	loadLibs();

	// create the chain
	StChain *chain  = new StChain("StChain");

	const char* inMuDstFile = fileList;
	// create the StMuDstMaker
	StMuDstMaker *muDstMaker = new StMuDstMaker(  	0,
							0,
							"",
							inMuDstFile,
							"MuDst.root",
							1
							);
	TChain& muDstChain = *muDstMaker.chain();
	if (extraFileList && strlen(extraFileList) > 0) {
		std::ifstream extraIn(extraFileList);
		std::string extraLine;
		int nAdded = 0;
		while (std::getline(extraIn, extraLine)) {
			if (extraLine.empty()) continue;
			muDstChain.Add(extraLine.c_str());
			nAdded++;
		}
		cout << "Added " << nAdded << " extra MuDst files from " << extraFileList
		     << " -- chain now has " << muDstChain.GetEntries() << " events total" << endl;
	}
	printf( "MuDst file has %d events available in tree\n", muDstChain.GetEntries());
	
	/*******************************************************************************************/
	// Initialize the database
		if (runDb){
			cout << endl << "============  Data Base =========" << endl;
			St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
			// dbMk->SetDateTime(20220225, 0); see below on reading run# and time from 1st evennt and InitRun
			// things will run fine without a timestamp set, but FCS DB will give bad values ...
		}
	/*******************************************************************************************/
	

	/*******************************************************************************************/
	// Create the StMuDst2StEventMaker
		StMuDst2StEventMaker * mu2ev = new StMuDst2StEventMaker();
	/*******************************************************************************************/

	/*******************************************************************************************/
	// Setup Fcs Database if needed
        //if ( (runFcsChain && runDb) || runFitQa){
	if ( runDb || runFitQa){
		StFcsDbMaker * fcsDb = new StFcsDbMaker();
		chain->AddMaker(fcsDb);
		// fcsDb->SetDebug();
	}
	/*******************************************************************************************/
	

	/*******************************************************************************************/
	// FTT chain
	if (runFttChain){
		gSystem->Load("libStFttDbMaker.so");
		StFttDbMaker * fttDbMk = new StFttDbMaker();
		// useDbGeometry=true (default): per-quadrant sTGC offsets = AGML nominal + the
		// Geometry/stgc survey tables.  false restores StFttDb's original hardcoded
		// numbers exactly, which is the A/B for whether the tables help.
		fttDbMk->setUseDbGeometry( useDbGeometry );
		printf("fwd_afterburner_db: StFttDb geometry source = %s\n",
		       useDbGeometry ? "DB survey tables" : "hardcoded");   // LOG_INFO is not in CINT scope here
		chain->AddMaker(fttDbMk);
		StFttHitCalibMaker * ftthcm = new StFttHitCalibMaker();
		StFttClusterMaker * fttclu = new StFttClusterMaker();
		// fttTimeCutMode: 2 = kTimeCutModeCalibratedTime (production default since
		// 2026-07, window from Run22-Run24 online QA), 1 = kTimeCutModeAcceptAll.
		// Use 1 for SPARSE data such as the zeroFieldAlignment physics-stream runs.
		// Mode 2 cuts on the per-VMM calibrated time from StFttHitCalibMaker, and a
		// VMM only gets a calibrated time after HitCalibHelper::ready(), which needs
		// >200 DISTINCT dbcid values -- until then every hit is time=-4097 and is
		// rejected. Measured 2026-09-14: field-on fwd stream 23081049 reaches it in
		// <100 events (366/368 VMMs ready), but ZF run 23063028 has 31x fewer FTT
		// hits/event and after 300 events only 4/373 VMMs are ready (H 0/128). The
		// few that do pass are the NOISIEST, since distinct-dbcid count rewards noise.
		// Real fix is per-VMM anchors from a dense run instead of per-file warm-up.
		fttclu->SetTimeCut(fttTimeCutMode, fttTimeCutLow, fttTimeCutHigh);
		StFttClusterPointMaker *fttCP = new StFttClusterPointMaker();
		fttCP->setApplyXYMirror( applyFttXYMirror ); // TEMP diagnostic, see signature note
		// StFttPointMaker * fttpoint = new StFttPointMaker();
	}
	/*******************************************************************************************/

	/*******************************************************************************************/
	// FCS Chain
	if (runFcsChain){
		gSystem->Load("libStFcsWaveformFitMaker.so");
		gSystem->Load("libStFcsClusterMaker.so");
		
		StFcsWaveformFitMaker *fcsWFF = new StFcsWaveformFitMaker();
		// This should only be used for simulated data, for real data this done in the database
		//fcsWFF->setEnergySelect(0);
		// This skips waveform analysis, and only apply new gain from DB
		fcsWFF->setAnaWaveform(false);
		StFcsClusterMaker *fcsclu = new StFcsClusterMaker();
	}
	/*******************************************************************************************/

	/*******************************************************************************************/
	// FwdTrackMaker Chain
	StFwdTrackMaker *fwdTrack = NULL;
	const int kNResidualTypes = 3;
	int residualTypesToRun[kNResidualTypes] = {0, 1, 4}; // Global, BLC, BLCVtx
	StFwdResidualMaker *fwdResiduals[kNResidualTypes] = {NULL, NULL, NULL};
	if (runFwdChain){
		// FwdTrackMaker
		fwdTrack = new StFwdTrackMaker();
		fwdTrack->SetDebug(debug);
		fwdTrack->setGeoCache( "fGeom.root" );
		fwdTrack->setSeedFindingWithFst();
		fwdTrack->setTrackRefit( true );

		// Fitter Options
		fwdTrack->setFitDebugLvl( 0 );
		fwdTrack->setFitMinIterations( 40 );
		fwdTrack->setFitMaxIterations( 100 );
		
		// fwdTrack->setDeltaPval( 1e-9 );
		// fwdTrack->setRelChi2Change( 1e-9 );

		// fwdTrack->setSeedFindingOff();
		// fwdTrack->setTrackFittingOff();
		fwdTrack->setFstHitSource( 2 /* = MUDST */);
		fwdTrack->setFttHitSource( 1 /* = STEVENT */);
		fwdTrack->setApplyFstWedgeAlignment( applyFstWedgeAlignment ); // see bugreport_StFstHitMaker.txt
		fwdTrack->setApplyFstGapFix( applyFstGapFix ); // outer-sensor kFstStripGapPhi sign, see StFwdHitLoader.h
		fwdTrack->setApplyFstMirror( applyFstMirror ); // DIAGNOSTIC only, see StFwdHitLoader.h

		// magField==0 only; see the note above the signature. Modes 1 and -1 are
		// handled where the field is built (search "MAGNETIC FIELD" below).
		if (magField == 0) fwdTrack->setZeroB( true );
		// StarFieldAdaptor options (STARField.h): uniform sign-correct Bz over the FST box,
		// and |z| beyond which B = 0 (<= 0: no cut)
		fwdTrack->setConfigKeyValue( "TrackFitter:fieldFstConstBz", (bool)fieldFstConstBz );
		fwdTrack->setConfigKeyValue( "TrackFitter:fieldZMax", (double)fieldZMax );
		// FST-blind diagnostic controls -- see the comment above addFttHits in FwdTracker.h
		fwdTrack->setConfigKeyValue( "TrackFitter:fttDiagType", (int)fttDiagType );
		fwdTrack->setConfigKeyValue( "TrackFitter:fttNoAdd",    (bool)fttNoAdd );
		fwdTrack->setConfigKeyValue( "TrackFitter:fttDiagMix",  (bool)fttDiagMix );

		if (runDb) fwdTrack->setUseBeamlineFromDB( true ); // use measured beamline for BLC; off for MC

		// fwdTrack->setConfigKeyValue("TrackFitter:doBeamlineTrackFitting", false);
        // fwdTrack->setConfigKeyValue("TrackFitter:doPrimaryTrackFitting", false);
        // fwdTrack->setConfigKeyValue("TrackFitter:doSecondaryTrackFitting", false);
        // skip finding fwd vertices
	}



		if (runFcsChain){
			// FwdTrack and FcsCluster assciation
			gSystem->Load("StFcsTrackMatchMaker");
			StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
			match->setMaxDistance(6,10);
			match->setFileName("fcstrk.root");
		}

		
		
		if (runFwdQa){
			StFwdQAMaker *fwdQA = new StFwdQAMaker();
			fwdQA->SetDebug(debug);
			TString fwdqaname( gSystem->BaseName(inMuDstFile) );
			fwdqaname.ReplaceAll(".MuDst.root", ".FwdTree.root");
			cout << fwdqaname.Data() << endl;
			fwdQA->setTreeFilename(fwdqaname);

			gSystem->Load("StFwdUtils.so");
			StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
			fwdAna->setMuDstInput();
		}

		if (runFwdResidual && runFwdChain){
			// Per-plane FST/FTT alignment residuals -- reads StEvent::fwdTrackCollection()
			// directly (same source StPicoDstMaker::fillFwdTracks() uses), NOT
			// muDst->muFwdTrackCollection(): the input StMuDstMaker here only reads,
			// it never calls the (protected) fillFwdTrack(StEvent*), so the MuDst-level
			// collection would otherwise be stale production-time data, not this
			// event's freshly refit tracks.
			//
			// One StFwdResidualMaker instance per track type, all added to the same
			// chain, all fed by the same single event loop -- the tracker already
			// fits every track type each event, so there's no need to re-run the
			// whole afterburner once per type (that used to cost 3x the wall time
			// for no reason: StMaker's own name defaulted the same for every
			// instance, which StChain doesn't tolerate, hence one type at a time).
			const char* residualTypeName[6] = {"Global","BLC","Primary","FwdVtx","BLCVtx","FCSConstrained"};
			for (int i = 0; i < kNResidualTypes; i++) {
				int rt = residualTypesToRun[i];
				TString residualName( gSystem->BaseName(inMuDstFile) );
				residualName.ReplaceAll(".MuDst.root", Form(".FwdDetResidual_%s.root", residualTypeName[rt]));
				fwdResiduals[i] = new StFwdResidualMaker(residualName, Form("fwdResidual_%s", residualTypeName[rt]));
				fwdResiduals[i]->setTrackType((UChar_t)rt);
			}
		}

		const int kNAlignTypes = 1;
		int alignTypesToRun[kNAlignTypes] = {4}; // BLCVtx only -- this is what the current
		                                          // production campaign is testing (the FST1/FST2
		                                          // translation lead found in data/BLCVtx). Add 0
		                                          // (Global) back to this array for a clean-baseline
		                                          // comparison sample too, at ~2x the alignment
		                                          // wall-time cost per job.
		StFwdAlignmentMaker *fwdAlignments[kNAlignTypes] = {NULL};
		if (runFwdAlignment && runFwdChain){
			// Unbiased (hit-removed) residuals -- see proposal_alignment_path.txt.
			// Completely separate output/maker from runFwdResidual above; does not
			// read or write anything it touches. One instance per track type, same
			// pattern as the fwdResiduals[] block above -- BLCVtx added to directly
			// test the FST1/FST2 translation signal StFwdResidualMaker's sine-fit
			// found there (see debug/residual summary discussion); Global stays in
			// as the clean baseline (no vertex-constraint leakage caveat -- see
			// StFwdAlignmentMaker.h).
			const char* alignTypeName[6] = {"Global","BLC","Primary","FwdVtx","BLCVtx","FCSConstrained"};
			for (int i = 0; i < kNAlignTypes; i++) {
				int rt = alignTypesToRun[i];
				TString alignName( gSystem->BaseName(inMuDstFile) );
				alignName.ReplaceAll(".MuDst.root", Form(".FwdAlignment_%s.root", alignTypeName[rt]));
				fwdAlignments[i] = new StFwdAlignmentMaker(alignName, Form("fwdAlignment_%s", alignTypeName[rt]));
				fwdAlignments[i]->setTrackMaker(fwdTrack);
				fwdAlignments[i]->setTrackType((UChar_t)rt);
			}
		}


	// The PicoDst
	if (runPico){
		gSystem->Load("libStPicoEvent");
		gSystem->Load("libStPicoDstMaker");
		StPicoDstMaker *picoMk = (StMaker*) (new StPicoDstMaker(StPicoDstMaker::IoWrite, inMuDstFile, "picoDst"));
		cout << "picoMk = " << picoMk << endl;
		picoMk->setVtxMode(StPicoDstMaker::Vtxless);
	}

	if ( runFitQa && runFwdChain){
		StFwdFitQAMaker *fwdFitQA = new StFwdFitQAMaker();
		fwdFitQA->SetDebug(debug);
		TString fitqaoutname(gSystem->BaseName(inMuDstFile));
		fitqaoutname.ReplaceAll(".MuDst.root", ".FwdFitQA.root");
		fwdFitQA->setOutputFilename( fitqaoutname );
	}
	/*******************************************************************************************/

	// gMessMgr->MemoryOff();

	/*******************************************************************************************/
	// Initialize chain
	chain->SetDebug(kError+debug);
	Int_t iInit = chain->Init();
	chain->SetDebug(kError+debug);
	cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
	// ensure that the chain initializes

	if ( iInit )
		chain->Fatal(iInit,"on init");
	
	// print the chain status
	chain->PrintInfo();

	// Read 1st event from MuDst to get run number and event time 
	if (runDb) {
	  muDstChain.GetEntry(0);
	  int run = muDstMaker->muDst()->event()->runNumber();
	  time_t tt = (time_t)muDstMaker->muDst()->event()->eventInfo()->time();
	  printf("Run number from 1st event: %d time: %d\n", run, tt);
	  struct tm* gmt = gmtime(&tt);
	  int date = (gmt->tm_year + 1900) * 10000 + (gmt->tm_mon + 1) * 100 + gmt->tm_mday;
	  int itime = gmt->tm_hour * 10000 + gmt->tm_min * 100 + gmt->tm_sec;
	  printf("GMT date: %d time: %d\n", date, itime);
	  dbMk->SetDateTime(date, itime);
	  chain->InitRun(run);
	  fcsDb->InitRun(run); //not sure why I need to call this separately...
	  StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");
	  //make sure we get good values
	  for(int d=0; d<4; d++){
	    StThreeVectorD off=fcsdb->getDetectorOffset(d);
	    printf("FCS Offset d=%1d %8.2f  %8.2f  %8.2f\n",d,off.x(),off.y(),off.z());
	  }
	  printf("FCS Gain 352 R17c0 = %8.4f\n",fcsdb->getGainCorrection(0,352));
	  printf("FCS Gain 374 R18c0 = %8.4f\n",fcsdb->getGainCorrection(0,374));
	}

	/*******************************************************************************************/
	// MAGNETIC FIELD.
	//
	// StarMagField is normally constructed by StBFChain::SetDbOptions from the
	// bfc FieldOn/FieldOff/HalfField/ReverseField option (StBFChain.cxx:1913).
	// This macro builds its chain BY HAND and so never created one. STARField.h's
	// StarFieldAdaptor guards every lookup with "if (StarMagField::Instance())"
	// and silently leaves B={0,0,0} when it is absent -- so every afterburner
	// production to date has tracked at ZERO FIELD regardless of the data.
	// Verified 2026-09-09: Instance() is NULL after loadLibs(), after chain->Init()
	// and after chain->InitRun(); and running the same file with zeroB on and off
	// gave bit-identical fitted momenta.
	//
	// Take the field from the MuDst rather than from an option: the production
	// stamps it per event, so there is nothing to set and no way to get it wrong
	// for a given file. Measured: -4.98700 kG for the 2022 ReversedFullField
	// production, +0.02994 kG for the 2022 zeroFieldAlignment runs.
	//
	// This must come AFTER chain->Init() (which is where TrackFitter picks its
	// AbsBField) but that is fine -- StarFieldAdaptor resolves Instance() at each
	// lookup, not once at init.
	{
	  // Only read an entry if the runDb block above did not already do it.
	  // Calling GetEntry(0) a second time is NOT neutral: it perturbs the
	  // MuDst buffers enough to change the fit (57 -> 55 tracks over 5 events
	  // of run 23063028), which would have made magField=-1 fail to reproduce
	  // the legacy output it exists to reproduce.
	  if (!runDb) muDstChain.GetEntry(0);
	  double bz = muDstMaker->muDst()->event()->magneticField(); // kGauss, signed
	  double scale = bz / 4.98;                                  // StarMagField nominal full field
	  printf("MuDst magneticField = %.5f kGauss (implies StarMagField scale %.5f)\n", bz, scale);
	  if (magField == 1) {
	    if (!StarMagField::Instance()) {
	      new StarMagField( StarMagField::kMapped, scale, kTRUE );
	      printf("magField=1: created StarMagField(kMapped, %.5f, locked)\n", scale);
	    } else {
	      printf("magField=1: StarMagField instance already exists -- left alone\n");
	    }
	  } else {
	    printf("magField=%d: NOT creating StarMagField."
	           " Instance()=%p -> StarFieldAdaptor returns B={0,0,0} EVERYWHERE."
	           " This is the pre-2026-09-09 behaviour; momenta are meaningless.\n",
	           magField, (void*)StarMagField::Instance());
	  }
	}
	/*******************************************************************************************/

	StMemStat stmem;
	stmem.PrintMem("BEFORE Event Loop");
	/*******************************************************************************************/
	// MAIN EVENT LOOP
	/*******************************************************************************************/
	size_t nEntries = muDstChain.GetEntries();
	if (nEntries > nEvents && nEvents > 0) {
		nEntries = nEvents;
		cout << "Limiting to " << nEntries << " events." << endl;
	}
	size_t numProcessed = 0;
	for (int i = 0; i < nEntries; i++) {
		printf("Processing event %d of %d\n", i, nEntries);
		if (i > 0) // skip first event to make it consistent
			stmem.Start();
		chain->Clear();
		if (fwdTrack)
			fwdTrack->SetDebug(debug);
		
        if (kStOK != chain->Make())
            break;

		if (refillMuDst){
			StEvent *mStEvent = static_cast<StEvent *>(muDstMaker->GetInputDS("StEvent"));
			// muDstMaker->fillFwdTrack( mStEvent);
			fwdQA->Make();
		}
		stmem.PrintMem(TString::Format("After Event %d:", i).Data());	
		if (i > 0)
			stmem.Stop();
		// MipMaker->Make();
		// picoMk->Make();
        cout << "EVENT #" << i << " COMPLETED" << endl; 
	}
	stmem.PrintMem("After Event Loop");
	stmem.Summary();
	/*******************************************************************************************/

	// Finish every maker explicitly, while the chain is still intact, instead of relying
	// on teardown at process exit: closes the picoDst and fcstrk.root, flushes the
	// residual/alignment ntuples, and prints StFttHitCalibMaker's per-job count of hits
	// timed from the DB anchor vs on the fly. StFwdResidualMaker/StFwdAlignmentMaker
	// Finish() are idempotent, so the exit-time pass that follows is harmless.
	// (Re-enabled 2026-09-14; had been commented out since 37b3744cb8.)
	chain->Finish();

	// delete chain;
}

void loadLibs(){	
	gSystem->Load("libStarClassLibrary.so");
	gSystem->Load("libStarRoot.so");
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	
	gSystem->Load("StarMagField");
	gSystem->Load("StMagF");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StTpcDb");
	gSystem->Load("StDaqLib");
	gSystem->Load("StDbBroker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("St_db_Maker");

	gSystem->Load("StEvent");
	gSystem->Load("StEventMaker");
	gSystem->Load("StarMagField");
 
	gSystem->Load("libGeom");
	gSystem->Load("St_g2t");
	
	// Added for Run16 And beyond
	gSystem->Load("libGeom.so");
	
	gSystem->Load("St_base.so");
	gSystem->Load("StUtilities.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("StarAgmlUtil.so");
	gSystem->Load("StarAgmlLib.so");
	gSystem->Load("libStarGeometry.so");
	gSystem->Load("libGeometry.so");
	
	gSystem->Load("xgeometry");
 
	gSystem->Load("St_geant_Maker");


	// needed since I use the StMuTrack
	gSystem->Load("StarClassLibrary");
	gSystem->Load("StStrangeMuDstMaker");
	gSystem->Load("StMuDSTMaker");
	gSystem->Load("StBTofCalibMaker");
	gSystem->Load("StVpdCalibMaker");
	gSystem->Load("StBTofMatchMaker");
	gSystem->Load("StFcsDbMaker");	

	/*******************************************************************************************/
	// loading libraries
	gSystem->Load("StFcsDbMaker");
	gSystem->Load( "StFttDbMaker" );
	gSystem->Load( "StFttHitCalibMaker" );
	gSystem->Load( "StFttClusterMaker" );
	gSystem->Load( "StFttClusterPointMaker" );
	gSystem->Load( "StFttPointMaker" );
	gSystem->Load("libStarGeneratorUtil.so");
	gSystem->Load("libgenfit2");
	gSystem->Load("libKiTrack");
	gSystem->Load("libXMLIO.so");
	gSystem->Load( "StFwdTrackMaker.so" );
	gSystem->Load( "StFwdResidualMaker.so" );
	gSystem->Load( "StFwdAlignmentMaker.so" );
	gSystem->Load( "StFwdUtils.so" );
	gSystem->Load("libStEpdUtil.so");
	gSystem->Load("StStarLogger.so");
	
	/*******************************************************************************************/


}
