//usr/bin/env root4star -l -b -q  $0; exit $?

#include "TTree.h"
#include "TClonesArray.h"

void qa(){

    // setup and make sure libraries are loaded
    gSystem->Load( "libStarRoot.so" );
	gSystem->Load("libStarClassLibrary.so");
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("libgenfit2.so");
	gSystem->Load("libKiTrack.so");
	gSystem->Load( "libStFwdTrackMaker.so" );

    // now open our data file
    TFile *f = new TFile("fwdtree.root");
    TTree *t = (TTree*)f->Get("fwd");

    // create the readers for each branch
	TClonesArray *mcTracks = new TClonesArray("StMuMcTrack");
   	t->GetBranch("mcTracks")->SetAutoDelete(kFALSE);
   	t->SetBranchAddress("mcTracks",&mcTracks);

	TClonesArray *fttPoints = new TClonesArray("StMuFttPoint");
	t->GetBranch("fttPoints")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("fttPoints",&fttPoints);

	TClonesArray *fttClusters = new TClonesArray("StMuFttCluster");
	t->GetBranch("fttClusters")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("fttClusters",&fttClusters);

	TClonesArray *fstPoints = new TClonesArray("StMuFstHit");
	t->GetBranch("fstHits")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("fstHits",&fstPoints);

	TClonesArray *wcal = new TClonesArray("FcsClusterWithStarXYZ");
	t->GetBranch("wcalClusters")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("wcalClusters",&wcal);

	TClonesArray *wcalHits = new TClonesArray("FcsHitWithStarXYZ");
	t->GetBranch("wcalHits")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("wcalHits",&wcalHits);

	TClonesArray *hcal = new TClonesArray("FcsClusterWithStarXYZ");
	t->GetBranch("hcalClusters")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("hcalClusters",&hcal);

	TClonesArray *hcalHits = new TClonesArray("FcsHitWithStarXYZ");
	t->GetBranch("hcalHits")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("hcalHits",&hcalHits);

	TClonesArray *epdHits = new TClonesArray("FcsHitWithStarXYZ");
	t->GetBranch("epdHits")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("epdHits",&epdHits);

	TClonesArray *fwdTracks = new TClonesArray("StMuFwdTrack");
	t->GetBranch("reco")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("reco",&fwdTracks);

	TClonesArray *seeds = new TClonesArray("StMuFwdTrackSeedPoint");
	t->GetBranch("seeds")->SetAutoDelete(kFALSE);
	t->SetBranchAddress("seeds",&seeds);

	
    //loop over the events
    for ( int i = 0; i < t->GetEntries(); i++ ){
        t->GetEntry(i);
        for ( int j = 0; j < fwdTracks->GetEntries(); j++ ){
            StMuFwdTrack *track = fwdTracks->At(j);
            printf("Track %d: pt=%f, eta=%f, phi=%f\n", j, track->momentum().Pt(), track->momentum().Eta(), track->momentum().Phi());

			StMuFwdTrackProjection projWCAL;
			track->getProjectionFor(kFcsWcalId, projWCAL);
			printf("Projection @ WCAL: det=%d, x=%f, y=%f, z=%f\n", projWCAL.mDetId, projWCAL.mXYZ.X(), projWCAL.mXYZ.Y(), projWCAL.mXYZ.Z());


			StMuFwdTrackProjection projHCAL;
			track->getProjectionFor(kFcsHcalId, projHCAL);
			printf("Projection @ HCAL: det=%d, x=%f, y=%f, z=%f\n", projHCAL.mDetId, projHCAL.mXYZ.X(), projHCAL.mXYZ.Y(), projHCAL.mXYZ.Z());
        
			// loop over WCAL clusters
			for ( int k = 0; k < wcal->GetEntries(); k++ ){
				FcsClusterWithStarXYZ *cluster = (FcsClusterWithStarXYZ*)wcal->At(k);
				
				printf("WCAL Cluster %d: x=%f, y=%f, z=%f\n", k, cluster->mXYZ.X(), cluster->mXYZ.Y(), cluster->mXYZ.Z());

			}

        
			// loop over WCAL clusters
			for ( int k = 0; k < wcal->GetEntries(); k++ ){
				FcsClusterWithStarXYZ *cluster = (FcsClusterWithStarXYZ*)wcal->At(k);
				
				printf("WCAL Cluster %d: x=%f, y=%f, z=%f\n", k, cluster->mXYZ.X(), cluster->mXYZ.Y(), cluster->mXYZ.Z());

			}
		
		}
    }
	cout << "Processed: " << t->GetEntries() << " entries" << endl;
}