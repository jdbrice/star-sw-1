// Run 22 sTGC survey, as it was hardcoded in StFttDb (X_shift_Quad*/Y_shift_Quad*/
// idealPlaneZLocations_Quad*). Those numbers are survey of the built detector, not
// design values -- see the header of StFttDb.h. Written here as the misalignment on
// top of the AGML placeholder, which is what the DB path composes:
//     t = hardcoded origin - AGML placeholder,  placeholder = (phx[pent], 5.9, zstation)
//     phx = {0, 0, -6.5, +6.5} cm by pent, zstation = {312.342, 329.953, 347.637, 365.422}
// so placeholder + this table reproduces the hardcoded origin exactly, to the digit.
// NO alignment is folded in -- this is the survey before any alignment iteration.
// row = 4*station + pent, pent order A,D,C,B (script/fttQuadConvention.C).
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
        row.t0   = 0.809;
        row.t1   = 3.634;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=1(Quad=A) survey, no alignment: dx=+0.809 dy=+3.634 dz=+0.000 cm\x00",1);
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
        row.t0   = -0.369;
        row.t1   = 3.67;
        row.t2   = -0.244;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=2(Quad=D) survey, no alignment: dx=-0.369 dy=+3.670 dz=-0.244 cm\x00",1);
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
        row.t0   = -4.251;
        row.t1   = 2.46;
        row.t2   = -0.025;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=3(Quad=C) survey, no alignment: dx=-4.251 dy=+2.460 dz=-0.025 cm\x00",1);
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
        row.t0   = 4.774;
        row.t1   = 2.524;
        row.t2   = -0.177;
    memcpy(&row.comment,"Station=1(Disk=0),Pent=4(Quad=B) survey, no alignment: dx=+4.774 dy=+2.524 dz=-0.177 cm\x00",1);
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
        row.t0   = 0.834;
        row.t1   = 3.533;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=1(Quad=A) survey, no alignment: dx=+0.834 dy=+3.533 dz=+0.000 cm\x00",1);
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
        row.t0   = -0.496;
        row.t1   = 3.54;
        row.t2   = -0.222;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=2(Quad=D) survey, no alignment: dx=-0.496 dy=+3.540 dz=-0.222 cm\x00",1);
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
        row.t0   = -4.322;
        row.t1   = 2.442;
        row.t2   = -0.23;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=3(Quad=C) survey, no alignment: dx=-4.322 dy=+2.442 dz=-0.230 cm\x00",1);
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
        row.t0   = 4.714;
        row.t1   = 2.437;
        row.t2   = 0.141;
    memcpy(&row.comment,"Station=2(Disk=1),Pent=4(Quad=B) survey, no alignment: dx=+4.714 dy=+2.437 dz=+0.141 cm\x00",1);
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
        row.t0   = 0.662;
        row.t1   = 3.703;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=1(Quad=A) survey, no alignment: dx=+0.662 dy=+3.703 dz=+0.000 cm\x00",1);
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
        row.t0   = -0.436;
        row.t1   = 3.655;
        row.t2   = -0.193;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=2(Quad=D) survey, no alignment: dx=-0.436 dy=+3.655 dz=-0.193 cm\x00",1);
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
        row.t0   = -4.487;
        row.t1   = 2.537;
        row.t2   = -0.291;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=3(Quad=C) survey, no alignment: dx=-4.487 dy=+2.537 dz=-0.291 cm\x00",1);
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
        row.t0   = 4.83;
        row.t1   = 2.461;
        row.t2   = -0.094;
    memcpy(&row.comment,"Station=3(Disk=2),Pent=4(Quad=B) survey, no alignment: dx=+4.830 dy=+2.461 dz=-0.094 cm\x00",1);
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
        row.t0   = 0.754;
        row.t1   = 3.601;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=1(Quad=A) survey, no alignment: dx=+0.754 dy=+3.601 dz=+0.000 cm\x00",1);
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
        row.t0   = -0.375;
        row.t1   = 3.516;
        row.t2   = 0.033;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=2(Quad=D) survey, no alignment: dx=-0.375 dy=+3.516 dz=+0.033 cm\x00",1);
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
        row.t0   = -4.391;
        row.t1   = 2.381;
        row.t2   = -0.111;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=3(Quad=C) survey, no alignment: dx=-4.391 dy=+2.381 dz=-0.111 cm\x00",1);
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
        row.t0   = 4.849;
        row.t1   = 2.481;
        row.t2   = -0.013;
    memcpy(&row.comment,"Station=4(Disk=3),Pent=4(Quad=B) survey, no alignment: dx=+4.849 dy=+2.481 dz=-0.013 cm\x00",1);
    tableSet->AddAt(&row);

    return (TDataSet *)tableSet;
}
