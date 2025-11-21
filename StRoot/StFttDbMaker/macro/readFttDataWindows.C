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

  // Get pointer to the raw table data
  auto table = dataset->GetTable();
  printf("table = %p\n", table);

  // Access as fttDataWindowsB_st pointer
  fttDataWindowsB_st *st = reinterpret_cast<fttDataWindowsB_st*>(table);

  cout << setw(6) << "UUID"
       << setw(6) << "Mode"
       << setw(8) << "Anchor"
       << setw(8) << "Min"
       << setw(8) << "Max"
       << endl;
  cout << "==============================================" << endl;
  for (int i = 0; i < 384; i++){
    cout << setw(6) << (int)(st[0].uuid[i])
         << setw(6) << (int)st[0].mode[i]
         << setw(8) << (int)st[0].anchor[i]
         << setw(8) << (int)st[0].min[i]
         << setw(8) << (int)st[0].max[i]
         << endl;
  }



}