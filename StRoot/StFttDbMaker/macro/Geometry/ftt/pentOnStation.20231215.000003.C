// Run 24 sTGC = .000003: iteration 1 with the strip-frame constant removed.
//
// .000002 was the Run 24 ZF alignment, fitted against the pre-825d1fe21e hits, so it
// absorbed the 17.55 mm once -- the Run 24 measurement is preserved, only the bug is
// taken out. Byte-equal to the local entry pentOnStation.20231210.000002.
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
        row.t0   = 0.19294;
        row.t1   = 4.00062;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=1(Quad=A,x>0,y>0) geo(+0.809,+3.634) + align(+1.331,+1.794)\x00",1);
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
        row.t0   = -0.23634;
        row.t1   = 2.79728;
        row.t2   = -0.24400;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=2(Quad=D,x<0,y>0) geo(-0.369,+3.670) + align(-1.266,+1.911)\x00",1);
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
        row.t0   = -4.54018;
        row.t1   = 2.22204;
        row.t2   = -0.02500;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=3(Quad=C,x<0,y<0) geo(-4.251,+2.460) + align(-1.138,-1.960)\x00",1);
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
        row.t0   = 4.10278;
        row.t1   = 2.73039;
        row.t2   = -0.17700;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=4(Quad=B,x>0,y<0) geo(+4.774,+2.524) + align(+1.397,-1.685)\x00",1);
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
        row.t0   = 0.13637;
        row.t1   = 3.85474;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=1(Quad=A,x>0,y>0) geo(+0.834,+3.533) ; NO V HITS - geometry term only\x00",1);
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
        row.t0   = -0.21967;
        row.t1   = 3.80460;
        row.t2   = -0.22200;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=2(Quad=D,x<0,y>0) geo(-0.496,+3.540) + align(-1.205,+1.765)\x00",1);
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
        row.t0   = -4.15966;
        row.t1   = 2.25091;
        row.t2   = -0.23000;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=3(Quad=C,x<0,y<0) geo(-4.322,+2.442) + align(-1.046,-1.903)\x00",1);
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
    memcpy(&row.comment,"Station=2(Disk=1),Pent=4(Quad=B,x>0,y<0) geo(+4.714,+2.437) + align(+1.539,-1.666)\x00",1);
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
        row.t0   = -0.02393;
        row.t1   = 4.40812;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=1(Quad=A,x>0,y>0) geo(+0.662,+3.703) + align(+1.139,+1.795)\x00",1);
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
        row.t0   = -0.18366;
        row.t1   = 4.00752;
        row.t2   = -0.19300;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=2(Quad=D,x<0,y>0) geo(-0.436,+3.655) + align(-1.193,+1.312)\x00",1);
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
        row.t1   = 2.51969;
        row.t2   = -0.29100;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=3(Quad=C,x<0,y<0) geo(-4.487,+2.537) + align(-1.375,-1.785)\x00",1);
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
        row.t0   = 4.38535;
        row.t1   = 2.90358;
        row.t2   = -0.09400;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=4(Quad=B,x>0,y<0) geo(+4.830,+2.461) + align(+1.482,-1.481)\x00",1);
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
        row.t0   = -0.17104;
        row.t1   = 4.09087;
        row.t2   = 0.00000;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=1(Quad=A,x>0,y>0) geo(+0.754,+3.601) + align(+1.079,+1.846)\x00",1);
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
        row.t0   = -0.15113;
        row.t1   = 3.95077;
        row.t2   = 0.03300;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=2(Quad=D,x<0,y>0) geo(-0.375,+3.516) + align(-1.139,+1.410)\x00",1);
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
        row.t1   = 2.42812;
        row.t2   = -0.11100;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=3(Quad=C,x<0,y<0) geo(-4.391,+2.381) + align(-0.964,-1.739)\x00",1);
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
        row.t1   = 3.05787;
        row.t2   = -0.01300;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=4(Quad=B,x>0,y<0) geo(+4.849,+2.481) + align(+1.361,-1.476)\x00",1);
    tableSet->AddAt(&row);

return (TDataSet *)tableSet;
}
