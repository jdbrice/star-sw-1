#include <string.h>
#include <math.h>
struct rotm_t {
  char name[5];
  //  float Thet1   Phi1  Thet2   Phi2  Thet3   Phi3;
  float r[6];
};
// Name      Thet1    Phi1   Thet2    Phi2   Thet3    Phi3
static rotm_t rotm[] = {
  {"000D",{   90.0,    0.0,   90.0,   90.0,    0.0,    0.0}}, //    (x,y,z) => ( x, y, z)
  {"XYZD",{  -90.0,    0.0,   90.0,   90.0,    0.0,    0.0}}, //    (x,y,z) => (-x, y, z) 
  {"XZXD",{   90.0,    0.0,  -90.0,   90.0,    0.0,    0.0}}, //    (x,y,z) => ( x,-y, z) 
  {"YZXD",{  -90.0,    0.0,  -90.0,   90.0,    0.0,    0.0}}, //    (x,y,z) => (-x,-y, z) 
  {"180X",{   90.0,    0.0,  -90.0,   90.0,  180.0,    0.0}}, //    (x,y,z) => ( x,-y,-z) 
  {"180R",{   90.0,    0.0,   90.0,   90.0,  180.0,    0.0}}, //    (x,y,z) => ( x, y,-z) 
  {"180D",{  -90.0,    0.0,   90.0,   90.0,  180.0,    0.0}}, //    (x,y,z) => (-x, y,-z) 
  {"180Y",{   90.0,   90.0,  -90.0,    0.0,  180.0,    0.0}}, //    (x,y,z) => ( y,-x,-z) 
  {"180Z",{   90.0,   90.0,  -90.0,    0.0,    0.0,    0.0}}, //    (x,y,z) => ( y,-x, z) 
  {"Z180",{  -90.0,   90.0,   90.0,    0.0,    0.0,    0.0}}, //    (x,y,z) => (-y, x, z) 
  {"ZSXD",{   90.0,   90.0,  180.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => ( y,-z, x) 
  {"90ZD",{  -90.0,    0.0,    0.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => (-x, z, y) 
  {"90XD",{   90.0,   90.0,    0.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => ( y, z, x) 
  {"90DX",{  -90.0,   90.0,  180.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => (-y,-z, x) 
  {"20XD",{   90.0,   70.0,    0.0,    0.0,   90.0,  -20.0}}, //    (x,y,z) => ( ?, z, ?) 
  {"60XD",{   90.0,  110.0,    0.0,    0.0,   90.0,   20.0}}, //    (x,y,z) => ( ?, z, ?) 
  {"90YZ",{  -90.0,   90.0,  180.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => (-y,-z, x) 
  {"90YD",{    0.0,    0.0,   90.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => ( z, x, y) 
  {"90XY",{    0.0,    0.0,  -90.0,   90.0,   90.0,    0.0}}, //    (x,y,z) => ( z,-y, x) 
  {"90YX",{  180.0,    0.0,   90.0,   90.0,   90.0,    0.0}}, //    (x,y,z) => (-z, y, x) 
  {"YXZ2",{   90.0,   90.0,    0.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => ( y, z, x) 
  {"YXZ3",{  -90.0,   90.0,   90.0,    0.0,  180.0,    0.0}}, //    (x,y,z) => (-y, x,-z) 
  {"YXZ4",{  -90.0,   90.0,   90.0,    0.0,    0.0,    0.0}}, //    (x,y,z) => (-y, x, z) 
  {"YXZ5",{  -90.0,   90.0,  -90.0,    0.0,  180.0,    0.0}}, //    (x,y,z) => (-y,-x,-z) 
  {"YXZ6",{   90.0,   90.0,   90.0,    0.0,    0.0,    0.0}}, //    (x,y,z) => ( y, x, z) 
  {"R000",{   90.0,    0.0,   90.0,   90.0,    0.0,    0.0}}, //    (x,y,z) => ( x, y, z) 
  {"R005",{   90.0,    5.0,   90.0,   95.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R010",{   90.0,   10.0,   90.0,  100.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R015",{   90.0,   15.0,   90.0,  105.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R020",{   90.0,   20.0,   90.0,  110.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R025",{   90.0,   25.0,   90.0,  115.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R030",{   90.0,   30.0,   90.0,  120.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R035",{   90.0,   35.0,   90.0,  125.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R040",{   90.0,   40.0,   90.0,  130.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R045",{   90.0,   45.0,   90.0,  135.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R050",{   90.0,   50.0,   90.0,  140.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R055",{   90.0,   55.0,   90.0,  145.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R060",{   90.0,   60.0,   90.0,  150.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R065",{   90.0,   65.0,   90.0,  155.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R070",{   90.0,   70.0,   90.0,  160.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R075",{   90.0,   75.0,   90.0,  165.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R080",{   90.0,   80.0,   90.0,  170.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R085",{   90.0,   85.0,   90.0,  175.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R090",{   90.0,   90.0,   90.0,  180.0,    0.0,    0.0}}, //    (x,y,z) => ( y,-x, z) 
  {"R095",{   90.0,   95.0,   90.0,  185.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R100",{   90.0,  100.0,   90.0,  190.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R105",{   90.0,  105.0,   90.0,  195.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R110",{   90.0,  110.0,   90.0,  200.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R115",{   90.0,  115.0,   90.0,  205.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R120",{   90.0,  120.0,   90.0,  210.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R125",{   90.0,  125.0,   90.0,  215.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R130",{   90.0,  130.0,   90.0,  220.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R135",{   90.0,  135.0,   90.0,  225.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R140",{   90.0,  140.0,   90.0,  230.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R145",{   90.0,  145.0,   90.0,  235.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R150",{   90.0,  150.0,   90.0,  240.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R155",{   90.0,  155.0,   90.0,  245.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R160",{   90.0,  160.0,   90.0,  250.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R165",{   90.0,  165.0,   90.0,  255.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R170",{   90.0,  170.0,   90.0,  260.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R175",{   90.0,  175.0,   90.0,  265.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R180",{   90.0,  180.0,   90.0,  270.0,    0.0,    0.0}}, //    (x,y,z) => (-x,-y, z) 
  {"R185",{   90.0,  185.0,   90.0,  275.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R190",{   90.0,  190.0,   90.0,  280.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R195",{   90.0,  195.0,   90.0,  285.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R200",{   90.0,  200.0,   90.0,  290.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R205",{   90.0,  205.0,   90.0,  295.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R210",{   90.0,  210.0,   90.0,  300.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R215",{   90.0,  215.0,   90.0,  305.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R220",{   90.0,  220.0,   90.0,  310.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R225",{   90.0,  225.0,   90.0,  315.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R230",{   90.0,  230.0,   90.0,  320.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R235",{   90.0,  235.0,   90.0,  325.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R240",{   90.0,  240.0,   90.0,  330.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R245",{   90.0,  245.0,   90.0,  335.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R250",{   90.0,  250.0,   90.0,  340.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R255",{   90.0,  255.0,   90.0,  345.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R260",{   90.0,  260.0,   90.0,  350.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R265",{   90.0,  265.0,   90.0,  355.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R270",{   90.0,  270.0,   90.0,    0.0,    0.0,    0.0}}, //    (x,y,z) => (-y, x, z) 
  {"R275",{   90.0,  275.0,   90.0,    5.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R280",{   90.0,  280.0,   90.0,   10.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R285",{   90.0,  285.0,   90.0,   15.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R290",{   90.0,  290.0,   90.0,   20.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R295",{   90.0,  295.0,   90.0,   25.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R300",{   90.0,  300.0,   90.0,   30.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R305",{   90.0,  305.0,   90.0,   35.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R310",{   90.0,  310.0,   90.0,   40.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R315",{   90.0,  315.0,   90.0,   45.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R320",{   90.0,  320.0,   90.0,   50.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R325",{   90.0,  325.0,   90.0,   55.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R330",{   90.0,  330.0,   90.0,   60.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R335",{   90.0,  335.0,   90.0,   65.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R340",{   90.0,  340.0,   90.0,   70.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R345",{   90.0,  345.0,   90.0,   75.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R350",{   90.0,  350.0,   90.0,   80.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R355",{   90.0,  355.0,   90.0,   85.0,    0.0,    0.0}}, //    (x,y,z) => ( ?, ?, z) 
  {"R360",{   90.0,  360.0,   90.0,   90.0,    0.0,    0.0}}, //    (x,y,z) => ( x, y, z) 
  {"000T",{   90.0,    0.0,   90.0,   90.0,  180.0,    0.0}}, //    (x,y,z) => ( x, y,-z) 
  {"005T",{   90.0,    5.0,   90.0,   95.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"010T",{   90.0,   10.0,   90.0,  100.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"015T",{   90.0,   15.0,   90.0,  105.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"020T",{   90.0,   20.0,   90.0,  110.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"025T",{   90.0,   25.0,   90.0,  115.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"030T",{   90.0,   30.0,   90.0,  120.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"035T",{   90.0,   35.0,   90.0,  125.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"040T",{   90.0,   40.0,   90.0,  130.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"045T",{   90.0,   45.0,   90.0,  135.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"050T",{   90.0,   50.0,   90.0,  140.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"055T",{   90.0,   55.0,   90.0,  145.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"060T",{   90.0,   60.0,   90.0,  150.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"065T",{   90.0,   65.0,   90.0,  155.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"070T",{   90.0,   70.0,   90.0,  160.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"075T",{   90.0,   75.0,   90.0,  165.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"080T",{   90.0,   80.0,   90.0,  170.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"085T",{   90.0,   85.0,   90.0,  175.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"090T",{   90.0,   90.0,   90.0,  180.0,  180.0,    0.0}}, //    (x,y,z) => ( y,-x,-z) 
  {"095T",{   90.0,   95.0,   90.0,  185.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"100T",{   90.0,  100.0,   90.0,  190.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"105T",{   90.0,  105.0,   90.0,  195.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"110T",{   90.0,  110.0,   90.0,  200.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"115T",{   90.0,  115.0,   90.0,  205.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"120T",{   90.0,  120.0,   90.0,  210.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"125T",{   90.0,  125.0,   90.0,  215.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"130T",{   90.0,  130.0,   90.0,  220.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"135T",{   90.0,  135.0,   90.0,  225.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"140T",{   90.0,  140.0,   90.0,  230.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"145T",{   90.0,  145.0,   90.0,  235.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"150T",{   90.0,  150.0,   90.0,  240.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"155T",{   90.0,  155.0,   90.0,  245.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"160T",{   90.0,  160.0,   90.0,  250.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"165T",{   90.0,  165.0,   90.0,  255.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"170T",{   90.0,  170.0,   90.0,  260.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"175T",{   90.0,  175.0,   90.0,  265.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"180T",{   90.0,  180.0,   90.0,  270.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"185T",{   90.0,  185.0,   90.0,  275.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"190T",{   90.0,  190.0,   90.0,  280.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"195T",{   90.0,  195.0,   90.0,  285.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"200T",{   90.0,  200.0,   90.0,  290.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"205T",{   90.0,  205.0,   90.0,  295.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"210T",{   90.0,  210.0,   90.0,  300.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"215T",{   90.0,  215.0,   90.0,  305.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"220T",{   90.0,  220.0,   90.0,  310.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"225T",{   90.0,  225.0,   90.0,  315.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"230T",{   90.0,  230.0,   90.0,  320.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"235T",{   90.0,  235.0,   90.0,  325.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"240T",{   90.0,  240.0,   90.0,  330.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"245T",{   90.0,  245.0,   90.0,  335.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"250T",{   90.0,  250.0,   90.0,  340.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"255T",{   90.0,  255.0,   90.0,  345.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"260T",{   90.0,  260.0,   90.0,  350.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"265T",{   90.0,  265.0,   90.0,  355.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"270T",{   90.0,  270.0,   90.0,    0.0,  180.0,    0.0}}, //    (x,y,z) => (-y, x,-z) 
  {"275T",{   90.0,  275.0,   90.0,    5.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"280T",{   90.0,  280.0,   90.0,   10.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"285T",{   90.0,  285.0,   90.0,   15.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"290T",{   90.0,  290.0,   90.0,   20.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"295T",{   90.0,  295.0,   90.0,   25.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"300T",{   90.0,  300.0,   90.0,   30.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"305T",{   90.0,  305.0,   90.0,   35.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"310T",{   90.0,  310.0,   90.0,   40.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"315T",{   90.0,  315.0,   90.0,   45.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"320T",{   90.0,  320.0,   90.0,   50.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"325T",{   90.0,  325.0,   90.0,   55.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"330T",{   90.0,  330.0,   90.0,   60.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"335T",{   90.0,  335.0,   90.0,   65.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"340T",{   90.0,  340.0,   90.0,   70.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"345T",{   90.0,  345.0,   90.0,   75.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"350T",{   90.0,  350.0,   90.0,   80.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"355T",{   90.0,  355.0,   90.0,   85.0,  180.0,    0.0}}, //    (x,y,z) => ( ?, ?,-z) 
  {"360T",{   90.0,  360.0,   90.0,   90.0,  180.0,    0.0}}, //    (x,y,z) => ( x, y,-z) 
  {"C001",{   90.0,   11.25,  0.0,    0.0,   90.0,  281.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C002",{   90.0,   22.50,  0.0,    0.0,   90.0,  292.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C003",{   90.0,   33.75,  0.0,    0.0,   90.0,  303.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C004",{   90.0,   45.00,  0.0,    0.0,   90.0,  315.0}}, //    (x,y,z) => ( ?, z, ?) 
  {"C005",{   90.0,   56.25,  0.0,    0.0,   90.0,  326.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C006",{   90.0,   67.50,  0.0,    0.0,   90.0,  337.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C007",{   90.0,   78.75,  0.0,    0.0,   90.0,  348.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C008",{   90.0,  101.25,  0.0,    0.0,   90.0,   11.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C009",{   90.0,  112.50,  0.0,    0.0,   90.0,   22.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C010",{   90.0,  123.75,  0.0,    0.0,   90.0,   33.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C011",{   90.0,  135.00,  0.0,    0.0,   90.0,   45.0}}, //    (x,y,z) => ( ?, z, ?) 
  {"C012",{   90.0,  146.25,  0.0,    0.0,   90.0,   56.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C013",{   90.0,  157.50,  0.0,    0.0,   90.0,   67.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C014",{   90.0,  168.75,  0.0,    0.0,   90.0,   78.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C015",{   90.0,  180.00,  0.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => (-x, z, y) 
  {"C016",{   90.0,  191.25,  0.0,    0.0,   90.0,  101.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C017",{   90.0,  202.50,  0.0,    0.0,   90.0,  112.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C018",{   90.0,  213.75,  0.0,    0.0,   90.0,  123.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C019",{   90.0,  225.00,  0.0,    0.0,   90.0,  135.0}}, //    (x,y,z) => ( ?, z, ?) 
  {"C020",{   90.0,  236.25,  0.0,    0.0,   90.0,  146.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C021",{   90.0,  247.50,  0.0,    0.0,   90.0,  157.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C022",{   90.0,  258.75,  0.0,    0.0,   90.0,  168.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C023",{   90.0,  270.00,  0.0,    0.0,   90.0,  180.0}}, //    (x,y,z) => (-y, z,-x) 
  {"C024",{   90.0,  281.25,  0.0,    0.0,   90.0,  191.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C025",{   90.0,  292.50,  0.0,    0.0,   90.0,  202.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C026",{   90.0,  303.75,  0.0,    0.0,   90.0,  213.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C027",{   90.0,  315.00,  0.0,    0.0,   90.0,  225.0}}, //    (x,y,z) => ( ?, z, ?) 
  {"C028",{   90.0,  337.50,  0.0,    0.0,   90.0,  247.5}}, //    (x,y,z) => ( ?, z, ?) 
  {"C029",{   90.0,  326.25,  0.0,    0.0,   90.0,  236.25}}, //    (x,y,z) => ( ?, z, ?) 
  {"C030",{   90.0,  348.75,  0.0,    0.0,   90.0,  258.75}}, //    (x,y,z) => ( ?, z, ?) 
  {"C031",{   90.0,  360.00,180.0,    0.0,   90.0,  270.0}}, //    (x,y,z) => ( x,-z,-y) 
  {"C032",{   90.0,   11.25,180.0,    0.0,   90.0,  101.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C033",{   90.0,   22.50,180.0,    0.0,   90.0,  112.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C034",{   90.0,   33.75,180.0,    0.0,   90.0,  123.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C035",{   90.0,   45.00,180.0,    0.0,   90.0,  135.0}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C036",{   90.0,   56.25,180.0,    0.0,   90.0,  146.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C037",{   90.0,   67.50,180.0,    0.0,   90.0,  157.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C038",{   90.0,   78.75,180.0,    0.0,   90.0,  168.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C039",{   90.0,   90.0,  180.0,    0.0,   90.0,  180.0}}, //    (x,y,z) => ( y,-z,-x) 
  {"C040",{   90.0,  101.25,180.0,    0.0,   90.0,  191.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C041",{   90.0,  112.50,180.0,    0.0,   90.0,  202.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C042",{   90.0,  123.75,180.0,    0.0,   90.0,  213.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C043",{   90.0,  135.00,180.0,    0.0,   90.0,  225.0}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C044",{   90.0,  146.25,180.0,    0.0,   90.0,  236.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C045",{   90.0,  157.50,180.0,    0.0,   90.0,  247.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C046",{   90.0,  168.75,180.0,    0.0,   90.0,  258.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C047",{   90.0,  180.00,180.0,    0.0,   90.0,  270.0}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C048",{   90.0,  191.25,180.0,    0.0,   90.0,  281.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C049",{   90.0,  202.50,180.0,    0.0,   90.0,  292.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C050",{   90.0,  213.75,180.0,    0.0,   90.0,  303.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C051",{   90.0,  225.00,180.0,    0.0,   90.0,  315.0}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C052",{   90.0,  236.25,180.0,    0.0,   90.0,  326.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C053",{   90.0,  247.50,180.0,    0.0,   90.0,  337.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C054",{   90.0,  258.75,180.0,    0.0,   90.0,  348.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C055",{   90.0,  270.00,180.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => (-y,-z, x) 
  {"C056",{   90.0,  281.25,180.0,    0.0,   90.0,   11.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C057",{   90.0,  292.50,180.0,    0.0,   90.0,   22.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C058",{   90.0,  303.75,180.0,    0.0,   90.0,   33.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C059",{   90.0,  315.00,180.0,    0.0,   90.0,   45.0}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C060",{   90.0,  326.25,180.0,    0.0,   90.0,   56.25}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C061",{   90.0,  337.50,180.0,    0.0,   90.0,   67.5}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C062",{   90.0,  348.75,180.0,    0.0,   90.0,   78.75}}, //    (x,y,z) => ( ?,-z, ?) 
  {"C063",{   90.0,    0.0,   0.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => ( x, z, y) 
  {"Z1YX",{   90.0,   90.0,   0.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => ( y, z, x) 
  {"Z2XY",{   90.0,    0.0, 180.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => ( x,-z, y) 
  {"Z3XY",{  -90.0,    0.0,   0.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => (-x, z, y) 
  {"Z4YX",{   90.0,   90.0, 180.0,    0.0,  -90.0,    0.0}}, //    (x,y,z) => ( y,-z,-x) 
  {"Z5YX",{  -90.0,   90.0,   0.0,    0.0,  -90.0,    0.0}}, //    (x,y,z) => (-y, z,-x) 
  {"Z6XY",{  -90.0,    0.0, 180.0,    0.0,  -90.0,   90.0}}, //    (x,y,z) => (-x,-z,-y) 
  {"Z7XY",{   90.0,    0.0,   0.0,    0.0,  -90.0,   90.0}}, //    (x,y,z) => ( x, z,-y) 
  {"Z8YX",{  -90.0,   90.0, 180.0,    0.0,   90.0,    0.0}}, //    (x,y,z) => (-y,-z, x) 
  {"XXZY",{   90.0,    0.0, 180.0,    0.0,   90.0,   90.0}}, //    (x,y,z) => ( x,-z, y) 
  {"XZXY",{   90.0,  180.0, 180.0,    0.0,   90.0,  -90.0}}, //    (x,y,z) => (-x,-z,-y) 
  {"YZXZ",{   90.0,   90.0,  90.0,    0.0,  180.0,    0.0}}, //    (x,y,z) => ( y, x,-z) 
  {"XXYZ",{   90.0,  180.0,  90.0,  -90.0,    0.0,    0.0}}  //    (x,y,z) => (-x,-y, z) 
};
static int nrotm = sizeof(rotm)/sizeof(rotm_t);
static int n = 5;
extern "C" int itorotm_(float r[6], char s[5], int nc) {
  static const char name[5] = "rotm";
  nc = n;
  for (int i = 0; i < nrotm; i++) {
    for (int j = 0; j < 6; j++) {
      if (fabs(r[j] - rotm[i].r[j])>1.e-7) goto NEXT;
    }
    memcpy(s,rotm[i].name,n); return i+1;
  NEXT:
    continue;
  }
  memcpy(s,name,n);
  return 0;
}