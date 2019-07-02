#include "CLCT.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>

using namespace cw;

int main(int argc, char * argv[]) {

	if(argc <= 1) {
        	std::cout << "Error! No filename specified!" << std::endl;
        	return -1;
	}

	std::string file_name(argv[1]);
	std::vector<CLCT> clcts;

	std::cout << "Reading from File prefix:  " << file_name << std::endl;
	ReadTxt(file_name, clcts);
	std::cout << "Successfully read " << clcts.size() << " Muons" << std::endl;

	for(int i=0;i<clcts.size();i++){
		PrintCLCT(clcts.at(i));
	}
	
	std::cout << std::endl;
	
	std::cout << "Generating " << file_name << "(.pat) files" << std::endl;
	
	if(WritePat(file_name, clcts)){ 
		return 0; 
	}
	else{
		return -1;
	}
/*

	std::string pref = "SamplePattern";
	std::string title = "TestRPattern";
	std::vector<CLCT> pat;
	std::vector<Hit> raw;

	int N = ReadTxt(pref, pat);
	ExtractHits(pat, raw);

	std::cout << pat.size() << std::endl;
	std::cout << raw.size() << std::endl;
	
	std::cout << std::endl;

	for(int i=0; i < pat.size(); i++){
		PrintCLCT(pat.at(i));
	}

	WriteTxt(title, pat);

	bool ek = WritePat(pat, title);

	char ch = 'a';
	while (ch == 'a') 
	{
		std::cout << "Press Any Key to Quit:	";
		std::cin >> ch;
	}

*/
	

	return 0;
}
