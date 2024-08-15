int n;
vector<double> x;
vector<double> y;
vector<double> sigmaX;
vector<double> sigmaY;
vector<double> sigmaXY;
vector<double> orientation;

int orr0_n=0;
vector<double> orr0_x;
vector<double> orr0_y;
vector<double> orr0_sigmaX;
vector<double> orr0_sigmaY;
vector<double> orr0_sigmaXY;
vector<double> orr0_orientation;

int orr1_n=0;
vector<double> orr1_x;
vector<double> orr1_y;
vector<double> orr1_sigmaX;
vector<double> orr1_sigmaY;
vector<double> orr1_sigmaXY;
vector<double> orr1_orientation;

int orr2_n=0;
vector<double> orr2_x;
vector<double> orr2_y;
vector<double> orr2_sigmaX;
vector<double> orr2_sigmaY;
vector<double> orr2_sigmaXY;
vector<double> orr2_orientation;

int orr3_n=0;
vector<double> orr3_x;
vector<double> orr3_y;
vector<double> orr3_sigmaX;
vector<double> orr3_sigmaY;
vector<double> orr3_sigmaXY;
vector<double> orr3_orientation;

void makeCanvas() {
    static int ican2 = 0; // Static variable to keep track of canvas count
    TCanvas *can = new TCanvas(TString::Format("can%d", ican2++), "", 900, 900);
    can->SetBottomMargin(0.1);
    can->SetLeftMargin(0.15);
    can->SetTopMargin(0.1);
    can->SetRightMargin(0.1);

    can->SetTickx();
    can->SetTicky();
}

void parse_points(){
    for (int i=0; i<n; i++){
        if (orientation[i]==0 && !TMath::IsNaN(sigmaXY[i]) ){
            orr0_n++;
            orr0_x.push_back(x[i]);
            orr0_y.push_back(y[i]);
            orr0_sigmaX.push_back(sigmaX[i]);
            orr0_sigmaY.push_back(sigmaY[i]);
            orr0_sigmaXY.push_back(sigmaXY[i]);
            orr0_orientation.push_back(orientation[i]);
        } else if (orientation[i]==1 && !TMath::IsNaN(sigmaXY[i]) ){
            orr1_n++;
            orr1_x.push_back(x[i]);
            orr1_y.push_back(y[i]);
            orr1_sigmaX.push_back(sigmaX[i]);
            orr1_sigmaY.push_back(sigmaY[i]);
            orr1_sigmaXY.push_back(sigmaXY[i]);
            orr1_orientation.push_back(orientation[i]);
        } else if (orientation[i]==2 && !TMath::IsNaN(sigmaXY[i]) ){
            orr2_n++;
            orr2_x.push_back(x[i]);
            orr2_y.push_back(y[i]);
            orr2_sigmaX.push_back(sigmaX[i]);
            orr2_sigmaY.push_back(sigmaY[i]);
            orr2_sigmaXY.push_back(sigmaXY[i]);
            orr2_orientation.push_back(orientation[i]);
        } else if (orientation[i]==3 && !TMath::IsNaN(sigmaXY[i]) ){
            orr3_n++;
            orr3_x.push_back(x[i]);
            orr3_y.push_back(y[i]);
            orr3_sigmaX.push_back(sigmaX[i]);
            orr3_sigmaY.push_back(sigmaY[i]);
            orr3_sigmaXY.push_back(sigmaXY[i]);
            orr3_orientation.push_back(orientation[i]);
        }
    }
}

void plot_with_ellipse(vector<double> xs, vector<double> ys, vector<double> sigmaxs, vector<double> sigmays, vector<double> sigmaxys, int color){
    for (size_t i = 0; i < xs.size(); ++i) {
        // Extract covariance matrix
        TMatrixD cov(2, 2);
        cov(0,0) = sigmaxs[i] * sigmaxs[i];
        cov(1,1) = sigmays[i] * sigmays[i];
        cov(0,1) = sigmaxys[i];
        cov(1,0) = sigmaxys[i];

        // Calculate eigenvalues and eigenvectors of the covariance matrix
        TMatrixDEigen eig(cov);
        TMatrixD eigValues = eig.GetEigenValues();
        TMatrixD eigVectors = eig.GetEigenVectors();

        // Semi-major and semi-minor axes (square root of eigenvalues)
        double a = TMath::Sqrt(eigValues(0,0));
        double b = TMath::Sqrt(eigValues(1,1));

        // Rotation angle (angle of the eigenvector corresponding to the largest eigenvalue)
        double angle = TMath::ATan2(eigVectors(1,0), eigVectors(0,0));

        if (TMath::IsNaN(angle) && !TMath::IsNaN(a) && !TMath::IsNaN(b) ) { continue; }

        // Create the ellipse
        TEllipse* ellipse = new TEllipse(xs[i], ys[i], a, b, 0, 360, angle * 180.0 / TMath::Pi());
        ellipse->SetLineColor(color);
        ellipse->SetLineWidth(2);
        ellipse->SetFillStyle(0);  // No fill
        ellipse->Draw("same");
    }

    gPad->Update();
}

void vis() {
    int iEvent = 0;
    int nEvents = 10;

    TFile *file = TFile::Open("st_physics_23072003_raw_3000004.MuDst.root");

    TTree *t = NULL;
    file->GetObject("MuDst", t);

    //FttPoint.mPlane == 0 to select a particular plane
    t->Draw("FttPoint.mXYZ.y():FttPoint.mXYZ.x():FttPoint.msigma_Y:FttPoint.msigma_X", "FttPoint.mPlane == 0", "goff", nEvents, iEvent);

    n = t->GetSelectedRows();

    // Resize vectors to fit the number of selected rows
    x.resize(n);
    y.resize(n);
    sigmaX.resize(n);
    sigmaY.resize(n);

    for (int i=0; i<n; i++) {
        x[i] = t->GetV2()[i];
        y[i] = t->GetV1()[i];
        sigmaX[i] = t->GetV4()[i];
        sigmaY[i] = t->GetV3()[i];
    }

    t->Draw("FttPoint.msigma_XY:FttPoint.mOrientation", "FttPoint.mPlane == 0", "goff", nEvents, iEvent);

    sigmaXY.resize(n);
    orientation.resize(n);

    for (int i=0; i<n; i++) {
        sigmaXY[i] = t->GetV1()[i];
        orientation[i] = t->GetV2()[i];
    }

    parse_points();

    TGraph * orr0_graph = new TGraph(orr0_n, &orr0_x[0], &orr0_y[0]);
    TGraph * orr1_graph = new TGraph(orr1_n, &orr1_x[0], &orr1_y[0]);
    TGraph * orr2_graph = new TGraph(orr2_n, &orr2_x[0], &orr2_y[0]);
    TGraph * orr3_graph = new TGraph(orr3_n, &orr3_x[0], &orr3_y[0]);

    cout << "sigma x: " << orr0_sigmaX[2] << " sigma_y" << orr0_sigmaY[2] << endl;

//drawing
TCanvas* c1 = new TCanvas("c1", "Ellipse Plot", 800, 600);
orr0_graph->GetXaxis()->SetLimits(-700, 700);
orr0_graph->GetYaxis()->SetLimits(-700, 700);
orr0_graph->SetMarkerColor(1);
orr1_graph->SetMarkerColor(2);
orr2_graph->SetMarkerColor(3);
orr3_graph->SetMarkerColor(4);
orr0_graph->Draw("AP");
orr1_graph->Draw("SAME P");
orr2_graph->Draw("SAME P");
orr3_graph->Draw("SAME P");
c1->Draw();
plot_with_ellipse(orr0_x, orr0_y, orr0_sigmaX, orr0_sigmaY, orr0_sigmaXY, 1);
plot_with_ellipse(orr1_x, orr1_y, orr1_sigmaX, orr1_sigmaY, orr1_sigmaXY, 2);
plot_with_ellipse(orr2_x, orr2_y, orr2_sigmaX, orr2_sigmaY, orr2_sigmaXY, 3);
plot_with_ellipse(orr3_x, orr3_y, orr3_sigmaX, orr3_sigmaY, orr3_sigmaXY, 4);
gPad->Print("StFttClusterPointMaker_visualization.png");
}
