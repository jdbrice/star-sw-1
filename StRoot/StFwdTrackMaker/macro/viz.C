//usr/bin/env root -l -b -q  $0'('$1')'; exit $?
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"


TFile * fData;
TTree * fwd;
TH2 * hFrame;
TCanvas *gCan;
TPad *padRZ, *padXY, *padStat;

float LegendX, LegendY;
float lineScale = 1.5;
int mTotalSeeds = 0;

enum ProjectionType { kXY, kRZ, kRZSigned, kXZ, kYZ };

float xx( float x, float y, float z, ProjectionType proj = kRZ ){

    if ( proj == kRZ || proj == kRZSigned){
        return z;//(TMath::ATan2( y, x ) + 2*3.1415926 )/ (2*3.14159) * 360;
    } else if ( proj == kXY ){
        return x;
    } else if ( proj == kXZ ){
        return z;
    } else if ( proj == kYZ ){
        return z;
    }

    return x;
}

float yy( float x, float y, float z, ProjectionType proj = kRZ ){

    if ( proj == kRZ ){
        float r = sqrt( pow(x, 2) + pow(y, 2) );
        return r;
    } else if ( proj == kXY ){
        return y;
    } else if ( proj == kRZSigned ){
        float r = sqrt( pow(x, 2) + pow(y, 2) );
        if ( y == 0 ) return r;
        r *= y / fabs(y);
        return r;
    } else if ( proj == kXZ ){
        return x;
    } else if ( proj == kYZ ){
        return y;
    }

    return y;
}

void viz_ftt_clusters( int eventIndex, ProjectionType projType ){
    TLine ll;
    ll.SetLineWidth( 3 );
    auto cmd = "fttClusters.pos.fX:fttClusters.pos.fY:-fttClusters.pos.fZ:fttClusters.mOrientation:fttClusters.mQuadrant";
    auto cond = "fttClusters.mNStrips > 1 && fttClusters.mSumAdc > 0";
    fwd->Draw( cmd, cond, "goff", 1, eventIndex );
    int N = fwd->GetSelectedRows();
    printf( "%s : has %d results \n", cmd, N );
    printf( "Projection Mode : %d \n", projType );

    auto x = fwd->GetV1();
    auto y = fwd->GetV2();
    auto z = fwd->GetV3();
    auto d = fwd->GetV4();
    auto q = fwd->GetVal(4);


    float x0, y0, x1, y1;
    for ( int i = 0; i < N; i++ ){

        x0 = xx( x[i], y[i], z[i], projType );
        y0 = yy( x[i], y[i], z[i], projType );
        float alpha = 0.3;
        // if (q[i]!=3) continue;

        float my = 1;
        float mx = 1;
        if ( q[i] == 1 ) {my = -1; }
        if ( q[i] == 2 ) {my = -1; mx = -1;}
        if ( q[i] == 3 ) {mx = -1;}

        if ( d[i] == 0 ){
            ll.SetLineColorAlpha( kRed, 0.3 );
            x1 = xx( x[i] + 60*mx, y[i], z[i], projType );
            y1 = yy( x[i], y[i], z[i], projType );
        } else {
            ll.SetLineColorAlpha( kBlue, 0.3 );
            x1 = xx( x[i], y[i], z[i], projType );
            y1 = yy( x[i], y[i] + 60*my, z[i], projType );
        }

        ll.DrawLine( x0, y0, x1, y1 );

    }


}

void viz_points( const char* name, const char* cmd,
                 int color, int eventIndex,
                 ProjectionType projType,
                 const char *cond = ""
){

    TLine ll;
    fwd->Draw( cmd, cond, "goff", 1, eventIndex );
    int N = fwd->GetSelectedRows();
    printf( "%s : has %d results \n", cmd, N );
    printf( "Projection Mode : %d \n", projType );

    auto cmdX = fwd->GetV1();
    auto cmdY = fwd->GetV2();
    auto cmdZ = fwd->GetV3();
    auto cmdE = fwd->GetV4();
    if ( cmdE != nullptr ){
        printf( "TOWERS\n" );
    }
    float vizX[5000];
    float vizY[5000];

    TText *t = new TText(.5,.5,"Hello World !");
    // t->SetTextAlign(22);
    t->SetTextColor(kBlack);
    t->SetTextFont(43);
    t->SetTextSize(20);

    int zColorStep = 90;
    int slc = color;
    int zColors[50];
    float zSizes[] = {2.5, 2.5, 2.0, 1.5, 2.5, 2.5, 2.5};
    for ( int i = 0; i < 50; i++ )
        zColors[i] = TColor::GetColorPalette(i*zColorStep % 255 );


    bool lgZ = false;
    float alpha = 0.6;
    for ( int i = 0; i < N; i++ ){

        vizX[i] = xx( cmdX[i], cmdY[i], cmdZ[i], projType );
        vizY[i] = yy( cmdX[i], cmdY[i], cmdZ[i], projType );
        // printf( "\tpoint at (%f, %f, %f) -> (%f, %f)\n", cmdX[i], cmdY[i], cmdZ[i], vizX[i], vizY[i] );

        int zIndex = 0;
        if ( fabs( cmdZ[i] - 151.75) < 2.5 ) zIndex = 1;
        if ( fabs( cmdZ[i] - 165.25) < 2.5 ) zIndex = 2;
        if ( fabs( cmdZ[i] - 178.75) < 2.5 ) zIndex = 3;


        TMarker *mk = new TMarker( vizX[i], vizY[i], 20 );
        mk->SetMarkerSize( 2.5 );
        if (zIndex >= 1 && zIndex < 50){
            slc = zColors[zIndex];
        }mk->SetMarkerSize( zSizes[zIndex] );



        // mk->SetMarkerSize( (float)(zIndex) *  0.5 + 0.5 );

        alpha = 0.6;
        if ( cmdE != nullptr ){
            mk->SetMarkerStyle( 21 );
            mk->SetMarkerSize( 0.5 + 0.5 * cmdE[i] );
            alpha = (cmdE[i] / 10.0);
            if (alpha>=1) alpha = 1;
        }

        // printf( "\tzIndex = %d -> color = %d \n", zIndex, slc );
        mk->SetMarkerColorAlpha( slc, alpha );
        if ( zIndex >= 1 ){
            mk->SetMarkerColorAlpha( slc, alpha );
            lgZ = true;
        }

        mk->Draw("same");

        if ( "FTC" == string(name) ){
            ll.SetLineColor( slc );
            ll.SetLineWidth( 2 );
            ll.DrawLine( vizX[i], vizY[i], 0, 0 );
        }
    }

    if ( lgZ ){
        for ( int i = 1; i < 4; i++){
            TMarker *mk1 = new TMarker( LegendX, LegendY, 20 );
            mk1->SetMarkerSize( 2.5 );
            mk1->SetMarkerColorAlpha( zColors[i], 0.5 );
            mk1->Draw("same");
            t->DrawText( LegendX + 2, LegendY - 0.5, TString::Format( "%s: %d", name, i ) );

            LegendY -= 5;
        }
    } else {
        TMarker *mk1 = new TMarker( LegendX, LegendY, 20 );
        mk1->SetMarkerSize( 2.5 );
        mk1->SetMarkerColor( color );
        mk1->Draw("same");
        t->DrawText( LegendX + 2, LegendY - 0.5, TString::Format( "%s:", name ) );

        LegendY -= 5;
    }
}

void viz_seeds(int nTrk, int eventIndex, ProjectionType projType){
    TLine ll;
    ll.SetLineWidth(lineScale);
    ll.SetLineColor(kGreen);
    // ll.DrawLine( 150, 10, 250, 20 );
    // Tracks
    mTotalSeeds = 0;
    fwd->Draw( "nSeedTracks", "", "goff", 1, eventIndex );
    mTotalSeeds = fwd->GetV1()[0];

    // if (nTrk == 0) nTrk = 500;
    for ( int i = 0; i < mTotalSeeds; i++ ){
        printf( "Seeds for trkId=%d\n", i );
        fwd->Draw( "seeds.pos.fX:seeds.pos.fY:seeds.pos.fZ", TString::Format("seeds.trackId == %d", i), "goff", 1, eventIndex );
        auto seedX = fwd->GetV1();
        auto seedY = fwd->GetV2();
        auto seedZ = fwd->GetV3();

        int slc = TColor::GetColorPalette(i*100 % 255);
        ll.SetLineColor(slc);
        printf( "---Found %d seed points", fwd->GetSelectedRows()  );
        for ( int j = 0; j < fwd->GetSelectedRows()-1; j++  ){

            float seedX1 = xx( seedX[j], seedY[j], seedZ[j], projType );
            float seedY1 = yy( seedX[j], seedY[j], seedZ[j], projType );

            printf( "seed(x=%f, y=%f, z=%)->(xx=%f, yy=%f)\n", seedX[j], seedY[j], seedZ[j], seedX1, seedY1 );

            float seedX2 = xx( seedX[j+1], seedY[j+1], seedZ[j+1], projType );
            float seedY2 = yy( seedX[j+1], seedY[j+1], seedZ[j+1], projType );

            ll.DrawLine( seedX1, seedY1, seedX2, seedY2 );
        } // end loop j
    } // end loop i

} // end viz_seeds


void viz_tracks(int nTrk, int eventIndex, ProjectionType projType, bool seeds = false, int iTrack = -1, bool filter = false){
    TLine ll;
    ll.SetLineWidth(lineScale);

    // ll.DrawLine( 150, 10, 250, 20 );
    // Tracks
    int NumTracksFound = 0;
    for ( int i = 0; i < nTrk; i++ ){
        if ( iTrack >= 0 && i != iTrack ) continue;

        // fwd->Draw( TString::Format("reco[%d].projs.mXYZ.fX:reco[%d].projs.mXYZ.fY:reco[%d].projs.mXYZ.fZ", i, i, i), TString::Format("reco[%d].status>=1 && fabs(reco[%d].mChi2) > 0.5", i, i), "goff", 1, eventIndex );
        fwd->Draw( TString::Format("reco[%d].projs.mXYZ.fX:reco[%d].projs.mXYZ.fY:reco[%d].projs.mXYZ.fZ:reco[%d].mChi2:reco[%d].nFailedPoints:reco[%d].q", i, i, i, i, i, i), "", "goff", 1, eventIndex );
        // fwd->Draw( TString::Format("0:5:reco[%d].projs.mXYZ.fZ", i, i, i), "", "goff", 1, eventIndex );
        auto trkX = fwd->GetV1();
        auto trkY = fwd->GetV2();
        auto trkZ = fwd->GetV3();
        auto trkChi2 = fwd->GetV4();
        auto trkConv = fwd->GetVal(4);
        auto trkQ = fwd->GetVal(5);

        TText text;
        text.SetTextFont(43);
        text.SetTextSize(36);
        if (iTrack >= 0){
            text.DrawTextNDC( 0.05, 0.7, TString::Format( "chi2=%f", trkChi2[0] ) );
            text.DrawTextNDC( 0.05, 0.65, TString::Format( "converge=%d", trkConv[0] ) );
        }else {
            // if ( trkChi2[0] > 100 ) continue;
        }

        // if ( trkChi2[0] > 0.05 ) ll.SetLineColor(kBlue);

        

        if ( fwd->GetSelectedRows() > 0 ){
            NumTracksFound++;
        }
        printf( "Track has %d projections -> projType = %d\n", fwd->GetSelectedRows(), projType );
        for ( int j = 0; j < fwd->GetSelectedRows()-1; j++  ){
            // if (j==0) continue;
            float trkX1 = xx( trkX[j], trkY[j], trkZ[j], projType );
            float trkY1 = yy( trkX[j], trkY[j], trkZ[j], projType );

            float trkX2 = xx( trkX[j+1], trkY[j+1], trkZ[j+1], projType );
            float trkY2 = yy( trkX[j+1], trkY[j+1], trkZ[j+1], projType );

            printf( "(%f, %f, %f) -> (%f, %f, %f)\n", trkX[j], trkY[j], trkZ[j], trkX[j+1], trkY[j+1], trkZ[j+1] );
            printf( "(%f, %f) -> (%f, %f)\n", trkX1, trkY1, trkX2, trkY2 );

            if ( true){
                if ( trkQ[0] > 0 )
                    ll.SetLineColor(kRed);
                else 
                    ll.SetLineColor(kBlue);
                ll.DrawLine( trkX1, trkY1, trkX2, trkY2 );
                
            }
        }
    }

    if (seeds == false) return;

    for ( int i = 0; i < nTrk; i++ ){
        if ( iTrack >= 0 && i != iTrack ) continue;

        fwd->Draw( TString::Format("reco[%d].seeds.pos.fX:reco[%d].seeds.pos.fY:reco[%d].seeds.pos.fZ", i, i, i), TString::Format("reco[%d].nFailedPoints==0", i), "goff", 1, eventIndex );
        auto seedX = fwd->GetV1();
        auto seedY = fwd->GetV2();
        auto seedZ = fwd->GetV3();

        // printf( "Found %d seeds for track %d\n", fwd->GetSelectedRows(), i );
        // int slc = TColor::GetColorPalette(i*100 % 255);
        ll.SetLineColor(kGreen);

        for ( int j = 0; j < fwd->GetSelectedRows()-1; j++  ){

            float seedX1 = xx( seedX[j], seedY[j], seedZ[j], projType );
            float seedY1 = yy( seedX[j], seedY[j], seedZ[j], projType );

            // printf( "seed(x=%f, y=%f, z=%)->(xx=%f, yy=%f)\n", seedX[j], seedY[j], seedZ[j], seedX1, seedY1 );

            float seedX2 = xx( seedX[j+1], seedY[j+1], seedZ[j+1], projType );
            float seedY2 = yy( seedX[j+1], seedY[j+1], seedZ[j+1], projType );

            // printf( "(%f, %f) -> (%f, %f)\n", seedX1, seedY1, seedX2, seedY2 );
            ll.DrawLine( seedX1, seedY1, seedX2, seedY2 );
            
            TMarker *mk1 = new TMarker( seedX1, seedY1, 20 );
            mk1->SetMarkerSize( 2.5 );
            mk1->SetMarkerColor(kBlue);
            mk1->Draw("same");

            TMarker *mk2 = new TMarker( seedX2, seedY2, 20 );
            mk2->SetMarkerSize( 2.5 );
            mk2->SetMarkerColor(kBlue);
            mk2->Draw("same");
        } // end loop j
    } // end loop i


} // viz Tracks

float statTextY = 0.97;
void n() { statTextY -= 0.05; }
void viz_stats( int eventIndex ){
    statTextY = 0.97;
    TText text;
    text.SetTextFont(43);
    text.SetTextSize(36);
    text.DrawTextNDC( 0.05, statTextY, TString::Format("Event : %d", eventIndex) ); n();

    fwd->Draw( "reco.mom.fX", "", "goff", 1, eventIndex );
    text.DrawTextNDC( 0.05, statTextY, TString::Format("#Tracks : %d", fwd->GetSelectedRows()) ); n();

    fwd->Draw( "reco.mom.fX", "reco.mChi2<100", "goff", 1, eventIndex );
    text.DrawTextNDC( 0.05, statTextY, TString::Format("#Tracks (good) : %d", fwd->GetSelectedRows()) ); n();

    fwd->Draw( "reco.mom.fX", "reco.mChi2<100 && reco.q==1", "goff", 1, eventIndex );
    text.DrawTextNDC( 0.05, statTextY, TString::Format("#Pos Tracks (good) : %d", fwd->GetSelectedRows()) ); n();
    fwd->Draw( "reco.mom.fX", "reco.mChi2<100 && reco.q==-1", "goff", 1, eventIndex );
    text.DrawTextNDC( 0.05, statTextY, TString::Format("#Neg Tracks (good) : %d", fwd->GetSelectedRows()) ); n();

    // fwd->Draw( "seeds.trackId", "", "goff", 1, eventIndex );
    // fwd->Draw( "nSeedTracks", "", "goff", 1, eventIndex );
    // mTotalSeeds = fwd->GetV1()[0];
    text.DrawTextNDC( 0.05, statTextY, TString::Format("#Seeds : %d", mTotalSeeds ) ); n();
}


int viz_event( int eventIndex, ProjectionType projType = kRZSigned, bool frame = false ){

    if ( projType == kRZSigned || projType == kXZ || projType == kYZ ){
        hFrame = new TH2F( "hFrame", ";z;R", 500, 0, 900, 120, -60, 60 );
        hFrame->SetTitle( "Event Visualization (RZ Signed)" );
        LegendX = 10;
        LegendY = 60;
    } else if ( projType == kRZ ){
        hFrame = new TH2F( "hFrame", ";z;R", 500, 0, 900, 60, 0, 60 );
        hFrame->SetTitle( "Event Visualization (RZ Signed)" );
        LegendX = 10;
        LegendY = 60;
    } else if ( projType == kXY ){
        hFrame = new TH2F( "hFrame", ";x;y", 5, -65, 65, 5, -65, 65 );
        hFrame->SetTitle( "Event Visualization (XY)" );
        LegendX = -40;
        LegendY = 40;
    }

    if ( frame ){
        hFrame->Draw("colz");
        return 1;
    }

    printf( "Visualizing Event %d \n", eventIndex );

    fwd->Draw( "reco.mom.fX", "", "goff", 1, eventIndex );
    int nTrk = fwd->GetSelectedRows();
    printf( "Event has %lld Tracks \n", nTrk );


    hFrame->Draw("colz");

    viz_points( "FTT", "ftt.pos.fX:ftt.pos.fY:ftt.pos.fZ", kRed, eventIndex, projType );
    // viz_points( "FTC", "fttClusters.pos.fX:fttClusters.pos.fY:-fttClusters.pos.fZ", kGreen, eventIndex, projType, "fttClusters.mNStrips>2" );
    viz_points( "FST", "fst.pos.fX:fst.pos.fY:fst.pos.fZ", kRed, eventIndex, projType );
    viz_points( "FCS", "wcalClusters.pos.fX:wcalClusters.pos.Y():wcalClusters.pos.Z()+715.0:wcalClusters.mEnergy", kGray, eventIndex, projType );

    // viz_ftt_clusters(eventIndex, projType);

    
    viz_tracks(nTrk, eventIndex, projType, false);
    viz_seeds(nTrk, eventIndex, projType);

    printf( "DONE with event!\n" );
    return 1;
}



void viz( int mode = 1, int maxEvents = 10, TString fn = "fwdtree.root") {

    fData = new TFile( fn );
    fwd = (TTree*)fData->Get( "fwd" );

    gStyle->SetOptStat(0);

    float canWidth = 19 * 100;
    float canHeight = 16 * 100;
    gCan = new TCanvas( "g", "", canWidth, canHeight );
    gCan->SetMargin( 0, 0, 0, 0);
    gCan->cd();
    gCan->Draw();

    padRZ = new TPad( "padRZ", "", 0.0, 0.5, 0.95, 0.99 );
    padRZ->SetMargin( .05,.01,.05,.01 );
    padRZ->Draw("same");
    padRZ->cd();

    gCan->cd();
    padXY = new TPad( "padXY", "", 0.0, 0.0, 0.5, 0.5 );
    padXY->SetMargin( .1,.02,.05,.01 );
    padXY->Draw("same");
    padXY->cd();

    gCan->cd();
    padStat = new TPad( "padStat", "", 0.5, 0.0, 1.0, 0.5 );
    padStat->SetMargin( .1,.02,.05,.01 );
    padStat->Draw("same");
    padStat->cd();

    // gPad->SetMargin(0.1, 0.05, 0.15, 0.05);
    int nEvents = fwd->GetEntries();
    if (nEvents > maxEvents) nEvents = maxEvents;
    // Viz event by event - all tracks hits etc. together for this event
    if ( mode == 0 ){
        
        // nEvents = 5;
        for ( int iEvent = 0; iEvent < nEvents; iEvent ++ ){

            printf( "Event: %d\n", iEvent );
            padRZ->cd();
            int nTrk = viz_event( iEvent, kRZSigned );
            padXY->cd();
            viz_event( iEvent, kXY );
            if (nTrk > -1){
                padRZ->Update();
                padXY->Update();

                padStat->cd();
                padStat->Clear();
                viz_stats( iEvent );
                padStat->Update();
                gCan->Update();
                gCan->Print( TString::Format( "out_event%d.pdf", iEvent ) );
            }
            hFrame->Reset();
        }
    }

    if ( mode != 1 ) return;
    // visualize the event one track at a time
    for ( int inEvent = 0; inEvent < nEvents; inEvent++ ){     
        fwd->Draw( "reco.mom.fX", "", "goff", 1, inEvent );
        int nTrk = fwd->GetSelectedRows();
        printf( "Event %d has %lld Tracks \n", inEvent, nTrk );

        for ( int iTrack = 0; iTrack < nTrk; iTrack ++ ){

            printf( "Track: %d\n", iTrack );

            padRZ->cd();
            // int nTrk = viz_event( iEvent, kRZSigned );
            viz_event( inEvent, kRZSigned, true );
            viz_tracks(nTrk, inEvent, kRZSigned, true, iTrack);
            padXY->cd();
            viz_event( inEvent, kXY, true );
            viz_tracks(nTrk, inEvent, kXY, true, iTrack);
            if (nTrk > -1){
                padRZ->Update();
                padXY->Update();

                padStat->cd();
                padStat->Clear();
                // viz_stats( iEvent );
                padStat->Update();
                gCan->Update();
                gCan->Print( TString::Format( "out_event%d_track%d.pdf", inEvent, iTrack ) );
            }
            hFrame->Reset();
        }
    }

}
