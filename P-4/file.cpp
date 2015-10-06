#include<iostream>
#include <dirent.h>
#include <fstream>
#include<string.h>

using namespace std;

int main(){

DIR *pDIR;
        struct dirent *entry;
        if( pDIR=opendir("/u/anramesh/anramesh-gpoornac-p4/attraction_images") ){
                while(entry = readdir(pDIR)){
                        if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ){
                        	char* name = entry->d_name;
				string path = "./attraction_images/";
				string full = path + name;
				cout <<full.c_str()<<endl;
			}
                }
                closedir(pDIR);
        }

return 0;
}
