// Run 23 sTGC = .000002: the Run 22 iteration-3 constants.
//
// Run 23 has never been analysed for sTGC alignment (AuAu, worth doing at some point),
// so this lane carries the Run 22 numbers as the best guess, exactly as .000001 did.
// Nothing here is a Run 23 measurement.
//
// WHY. Reco built the cluster position from the strip's row-band coordinate instead
// of the strip centre, dropping 17.55 mm from every FTT hit. Commit 825d1fe21e fixed
// reco and the slow sim together, which left every table that had been FITTED against
// the old hits absorbing that constant. The correction is
//     t0 -= 1.755*sx,  t1 -= 1.755*sy   [cm]
// with StFttDb's per-quadrant signs A (+,+) B (+,-) C (-,-) D (-,+), applied with
// script/stgcTableAddDelta.py.
//
// MEASURED, not assumed: the 90-job Run-22 ZF campaign
// condor/submit_pico_iter3_20261024.txt -> blinddiag_iter3.root, read with
// script/stgcQuadOffsets.C(file,0,2,5.9), finds the per-quadrant residual moved by
// exactly 1.755*(sx,sy) over 25 independently fitted coordinates: mean +0.040 mm,
// rms 0.44 mm, max deviation 1.26 mm, at 4.8-175 sigma. Because the delta is a
// universal constant at that level, it is applied analytically to all 16 cells
// rather than cell by cell -- a per-cell fit would have left the 5 dx cells under
// stgcQuadOffsets' 5-sigma floor HELD at a value known to be 17.55 mm wrong.
//
// stgcOnTpc and stationOnStgc are NOT touched and need no version bump: the global
// average of the per-quadrant constant is exactly zero (x signs +,+,-,- and y signs
// +,-,-,+ both sum to zero), so neither ever absorbed any of it. The iteration-3
// measurement agrees -- its overall mean is +0.040 mm.
//
// t2 (per-quadrant z) is untouched: the strip-frame constant lies in the plane.
//
// MATCHED PAIR. pentOnStation is read by StFttDb AND by AgPosition when the geometry
// is built with "misalign", so this entry and its fGeom cache go together. A job that
// mixes one with the other is silently 17.55 mm wrong.
TDataSet *CreateTable() {
    if (!TClass::GetClass("St_Survey")) return 0;
Survey_st row;
St_Survey *tableSet = new St_Survey("pentOnStation",16);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 1;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 0.34345;
        row.t1   = 2.85035;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=1(Quad=A,x>0,y>0) = .000003 + (-1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 2;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 0.62397;
        row.t1   = 2.79728;
        row.t2   = -0.24400;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=2(Quad=D,x<0,y>0) = .000003 + (+1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 3;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = -3.16708;
        row.t1   = 1.90912;
        row.t2   = -0.02500;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=3(Quad=C,x<0,y<0) = .000003 + (+1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 4;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 4.16125;
        row.t1   = 2.28101;
        row.t2   = -0.17700;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=4(Quad=B,x>0,y<0) = .000003 + (-1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 5;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = -0.92100;
        row.t1   = 1.77800;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=1(Quad=A,x>0,y>0) = .000003 + (-1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 6;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 0.57787;
        row.t1   = 2.23337;
        row.t2   = -0.22200;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=2(Quad=D,x<0,y>0) = .000003 + (+1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 7;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = -3.23763;
        row.t1   = 1.94777;
        row.t2   = -0.23000;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=3(Quad=C,x<0,y<0) = .000003 + (+1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 8;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 4.29533;
        row.t1   = 2.22879;
        row.t2   = 0.14100;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=4(Quad=B,x>0,y<0) = .000003 + (-1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 9;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 0.00087;
        row.t1   = 2.84784;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=1(Quad=A,x>0,y>0) = .000003 + (-1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 10;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 0.64276;
        row.t1   = 2.38015;
        row.t2   = -0.19300;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=2(Quad=D,x<0,y>0) = .000003 + (+1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 11;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = -4.10728;
        row.t1   = 2.14889;
        row.t2   = -0.29100;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=3(Quad=C,x<0,y<0) = .000003 + (+1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 12;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 4.28568;
        row.t1   = 2.37934;
        row.t2   = -0.09400;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=4(Quad=B,x>0,y<0) = .000003 + (-1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 13;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = -0.01015;
        row.t1   = 2.78199;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=1(Quad=A,x>0,y>0) = .000003 + (-1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 14;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 0.76490;
        row.t1   = 2.00331;
        row.t2   = 0.03300;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=2(Quad=D,x<0,y>0) = .000003 + (+1.755,-1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 15;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = -3.60010;
        row.t1   = 2.03931;
        row.t2   = -0.11100;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=3(Quad=C,x<0,y<0) = .000003 + (+1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

    memset(&row,0,tableSet->GetRowSize());
        row.Id   = 16;
        row.r00  = 1.0;
        row.r01  = 0.0;
        row.r02  = 0.0;
        row.r10  = 0.0;
        row.r11  = 1.0;
        row.r12  = 0.0;
        row.r20  = 0.0;
        row.r21  = 0.0;
        row.r22  = 1.0;
        row.t0   = 3.82266;
        row.t1   = 2.43605;
        row.t2   = -0.01300;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=4(Quad=B,x>0,y<0) = .000003 + (-1.755,+1.755) cm, the strip-frame constant\x00",1);
    tableSet->AddAt(&row);

return (TDataSet *)tableSet;
}
