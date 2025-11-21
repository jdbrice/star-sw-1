class StFttDb;


void readFttDataWindows(){

    //-- load dBase and Table definition libraries
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("StUtilities");
  gSystem->Load("St_Tables.so");

  gSystem->Load("StDbLib.so");
  gSystem->Load("libStDb_Tables.so");

  //-- get the singleton manager
  StDbManager* dbManager = StDbManager::Instance();

  //-- connect to the db & get an empty container
  StDbConfigNode* configNode = dbManager->initConfig("Calibrations_ftt");
  string ZReadTime = "2022-12-31 23:59:59";
  dbManager->setRequestTime(ZReadTime.c_str());

  StDbTable* dataset = configNode->addDbTable("fttDataWindowsB");
  dbManager->fetchDbTable(dataset);

  Int_t rows = dataset->GetNRows();
  printf( "rows = %d\n", rows);

  fttDataWindowsB_st *st = static_cast<fttDataWindowsB_st*>(dataset);
  printf("st = %p\n", st);

  for (int i = 0; i < 384; i++){
    // printf( "uuid[%d] = %d\n", i, (int)st[0].anchor[i] );
    cout << setw(6) << (int)(st[0].uuid[i]) << setw(6) << (int)st[0].mode[i] << endl;
  }



}