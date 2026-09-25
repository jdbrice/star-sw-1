// per-station: IDENTITY. Measured per-disk common terms are <=0.08 cm in x and <=0.04 cm in y, an order of magnitude below the per-quadrant terms, so nothing is put here.
TDataSet *CreateTable() {
    if (!TClass::GetClass("St_Survey")) return 0;
Survey_st row;
St_Survey *tableSet = new St_Survey("stationOnStgc",4);

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
        row.t0   = 0.00000;
        row.t1   = 0.00000;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=1(Disk=0) identity\x00",1);
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
        row.t0   = 0.00000;
        row.t1   = 0.00000;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=2(Disk=1) identity\x00",1);
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
        row.t0   = 0.00000;
        row.t1   = 0.00000;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=3(Disk=2) identity\x00",1);
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
        row.t0   = 0.00000;
        row.t1   = 0.00000;
        row.t2   = 0.0;
    memcpy(&row.comment,"Station=4(Disk=3) identity\x00",1);
    tableSet->AddAt(&row);

return (TDataSet *)tableSet;
}
