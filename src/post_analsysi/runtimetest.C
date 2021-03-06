// ROOT includes
//#include </usr/include/mysql/mysql.h>
#include </state/partition2/Data/2010analysis/includes/TMySQLResult.h>
#include </state/partition2/Data/2010analysis/includes/TMySQLServer.h>
#include </state/partition2/Data/2010analysis/includes/TMySQLRow.h>
#include <TApplication.h>
// C includes
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "iostream"

int main()
{
  using namespace std;
  TApplication theApp("App",0,0);

   // read in the analysis database location.
  Float_t rtime_e,rtime_w;
  Int_t runnum = 16242;
  Int_t flipperOn = 0;
  char db[100];
  fstream f;
  TMySQLResult *res;
  TMySQLRow *row;
  char buffer[200];

  // open file with current database loction
  cout << "opening" << endl;	
  f.open("mysql_db.txt",fstream::in);
  f >> db;
  f.close();
  
  cout << "Attempting to connect " << endl;
  // Connect to Caltech database
  TMySQLServer *sql = new TMySQLServer(db,"ucn","UCNBetAs");  
    
  // get run time
 
  cout << "Connected" << endl;
//  sprintf(buffer,"select live_time_e,live_time_w,live_time from analysis where run_number between  %d and %d"
//	  ,runnum,runnum+20);
  sprintf(buffer,"select unix_timestamp(end_time)-unix_timestamp(start_time) from run where run_number between %d and %d",runnum,runnum+20);

  cout << buffer << endl;

  res = (TMySQLResult*)sql->Query(buffer);
  cout << "here" << endl;
  if(res->GetRowCount() != 0){
    while((row = (TMySQLRow*)res->Next())){
      rtime_e = atof(row->GetField(0));
      rtime_w = atof(row->GetField(0));
      cout << rtime_e << "\t" << rtime_w << endl;
      delete row;
    }
  }
  cout << "here 2" << endl; 
  // get flipper status

  sprintf(buffer,"select flipper from run where run_number = %d ",runnum);

  res = (TMySQLResult*)sql->Query(buffer);
  flipperOn = 0;

  if(res->GetRowCount() != 0){
    while((row = (TMySQLRow*)res->Next())){
      if(!strncmp(row->GetField(0),"On",5)){
        flipperOn = 1;
      }
      delete row;
    }
  }
  
  // 

  cout << "at end " << endl;
  delete res;
  delete sql;

  theApp.Run();

  return runnum;
}
