#ifndef __EXTRACT_RUNLOG_CC__
#define __EXTRACT_RUNLOG_CC__

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

using namespace std;

struct linedata_t {
        char Date[12];
	int Run_Number;
        time_t  StartTime;
        time_t  StopTime;
	char DataSize[10];
        char Comment[80];
	int runtype;  // 0 - nothing , 1 - octet , 2-calibration
};

int getruntype(char *com);

int main(int argc,char *argv[])
{
  // Get the runlog file from the data 2011 data location
  FILE *fRunLog; 
  fstream fListout;
  fListout.open("parsed_output_list.txt",fstream::out);
  char *path = getenv("UCNARAWDATA2011");
  char runlog[500];
  sprintf(runlog,"%s/runlog.txt",path);
  // Open the file
  fRunLog = fopen(runlog,"r");
  // Create a structure to hold the parsed line data.
  linedata_t dataparse;
  // loop through file
  do{	
	int nseg = 0;
	char line[600];
	fgets(line,600,fRunLog);
	char *seg;
	seg = strtok(line,"\t");
	while ( seg != NULL){
		//cout << "Seg " << seg << endl;
		switch (nseg) {
			case 0 :
				sprintf(dataparse.Date,"%s",seg);
				break;
			case 1 :
				dataparse.Run_Number = atoi(seg);
				break;
			case 2 :
				cout << "Time is : " << seg << endl;
				break;
			case 3 : 
				cout << "End time is : " << seg << endl;
				break;
			case 4 : 
				sprintf(dataparse.DataSize,"%s",seg);
				break;
			case 5 :
				sprintf(dataparse.Comment,"%s",seg);
				dataparse.runtype = getruntype(dataparse.Comment);
				break;
			default :
				break;
		};
		nseg++;
		seg = strtok(NULL,"\t");
	}
	cout << "run type is " << dataparse.runtype << "\t" << dataparse.Comment << endl;
	if(dataparse.runtype == 1) fListout << dataparse.Run_Number << endl;
  }while(!(feof(fRunLog)));

  fListout.close();
  
  return -1;
}

int getruntype(char *com)
{
	string comments (com);
	// list the octet types.
	string octets[] = {"A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12",
		           "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12"};

	string calib[] = {"Xeon","Bi","Sn","Tin","Sr","Ce","Calibration"};

	for(int i = 0 ;i < 24 ; i++){	
		int nsize = octets[i].size();
		unsigned found = comments.find(octets[i]);
		if(found != std::string::npos && found < 100) return 1;
	}
	for(int i = 0 ; i < 7 ; i++){
		unsigned found = comments.find(calib[i]);
                if(found!=std::string::npos && found < 100 ) return 2;
        }

        return 0;
}

#endif
