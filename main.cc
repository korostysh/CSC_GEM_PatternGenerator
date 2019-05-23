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
	
	std::string pref = "SamplePattern";
	std::string title = "Assbutt";
	std::vector<CLCT> pat;
	std::vector<Hit> raw;

	int N = ReadTxt(pref, pat);
	ExtractHits(pat, raw);

	std::cout << pat.size() << std::endl;
	std::cout << raw.size() << std::endl;

	WriteTxt(title, pat);

	bool ek = WritePat(pat, title);

	char ch = 'a';
	while (ch == 'a') 
	{
		std::cout << "Press Any Key to Quit:	";
		std::cin >> ch;
	}


	/*
	hits.push_back(hij);
	// fill a Hit vector
	for(int j=0; j < 4; j++){
		Hit hj(1+j,1+(j*10),2+j);
		hits.push_back(hj);
	}
	// output
	for(int j=0; j < hits.size(); j++){
		
		if(hit == hits.at(j) || true){ 
			hits.at(j).print_cmd(); 
		}
	}	
		


	std::cout << "Number of hits : " << hits.size() << '\n' << '\n';
	std::sort(hits.begin(), hits.end(), CompareHit);

	for(int i=0; i < hits.size(); i++){
		hits.at(i).set_layer(0);
		hits.at(i).print_cmd();
	}	

	std::vector<unsigned char> bytes;
	for(int i=0; i < hits.size(); i++){
		hits.at(i).get_triad(bytes);
		if(hits.at(i).set_layer(0)){
		std::bitset<8> A(bytes.at(0));
		std::bitset<8> B(bytes.at(1));
		std::bitset<8> C(bytes.at(2));
		std::cout << A << ' ' << B << ' ' << C << '\n';
	}}

	SetFlags(hits);
	for(int i=0; i < hits.size(); i++){
		hits.at(i).print_cmd();
	}
	std::cout << '\n';	
	Clean(hits);
	for(int i=0; i < hits.size(); i++){
		hits.at(i).print_cmd();
	}

	std::cout << '\n' << '\n';
	int xx;
	Hit AA(0,30,0);
	Hit BB(2,31,0);
	xx = IsConflict( AA, BB );
	std::cout << xx << '\n';
	AA = Hit(0,30,0);
	BB = Hit(0,3,0);
	xx = IsConflict( AA, BB );
	std::cout << xx << '\n';
	AA = Hit(3,0,0);
	BB = Hit(3,3,0);
	xx = IsConflict( AA, BB );
	std::cout << xx << '\n';

	CLCT CD;
	CD.print_cmd();
	CLCT CF(69, 69, 4, 6);
	CF.print_cmd();
	CLCT CG = CF;
	CG.print_cmd();


	std::vector<CLCT> clcts, clcfs;
	std::string str = "cam";
	ReadFile(str, clcts);
	
	clcts.at(0).print_cmd();

	std::cout << " \n";

	clcfs.push_back(CF);
	clcfs.push_back(CLCT(71, 70, 4, 4));
	str.at(0) = 'C';
	WriteFile(str,clcfs);

	std::vector<Hit> jit;
	clcfs.at(1).get_hits(jit);
	for(int j=0; j < jit.size(); j++){
		jit.at(j).print_cmd();
	}


	for(int i=0; i <clcts.size(); i++){
		clcts.at(i).print_cmd();
	}
	*/

	return 0;
}
