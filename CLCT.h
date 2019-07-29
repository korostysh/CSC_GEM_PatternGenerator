#include <cstdlib>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>
//#include "TaoFunc.h"
//#include "CSCConstants.h"
/// TODO
//		--.pat output function
//		--Add Overlapping Hit Checking

#define COMPILE_TYPE 0xc
#define LAYERS 6
#define NUM_CFEB 7
#define STRIPS_CFEB 32			// number halfstrips in cfeb
#define RAMPAGE_SIZE 4*1024
#define LAYER_SET 1



namespace cw {


	//	Structures & Classes
	struct Hit
	{
	public:
		// Data Members
		int bx;
		int hs;
		int lay;

		// Constructors
		Hit(void);
		Hit(int Bx, int Hs, int Layer);
                friend std::ostream& operator << (std::ostream& os, const Hit& hit);
	};

	class CLCT
	{
	public:
		// Data Members
		int bx;
		int hs;
		int pat;
		int nhits;
		std::vector<Hit> hits;

		// Constructors
		CLCT(void);
		CLCT(int, int, int, int);
		CLCT(const CLCT& clct);

		// Operators
		//CLCT& operator = (const CLCT&);
		friend std::ostream& operator << (std::ostream&, const CLCT&);
		friend std::istream& operator >> (std::istream&, CLCT&);
	};

	class Group
	{
	public:
		std::vector<unsigned char> hexdata;

		Group(void);
		Group(std::vector<Hit>&, int);

		void addHit(Hit&, int);
		friend std::ostream& operator << (std::ostream&, const Group&);
	};

	class Cluster
	{
	public:
		// Data Members	  Range
		int bx;		// 0-500
		int roll;	// 1-8
		int pad;	// 0-383
		int size;	// 1-8
		int layer;	// 1-2
		
		Cluster(void);
		Cluster(int, int, int, int, int);
		Cluster(const Cluster&);

		friend std::ostream& operator<<(std::ostream&, const Cluster&);
		friend std::istream& operator>>(std::istream&, Cluster&);
	};
	


	
	//	Functions
	int GetCFEB(int hs);		//	Out:	Cfeb given half strip
	int GetLocal(int hs);		//	Out:	Halfstrip relative to CFEB

	// Tao Style (.txt) Pattern files
	int ReadTxt(std::string&, std::vector<CLCT>&);					// input : file prefix ONLY
	std::string ReadTxt(std::string&, std::vector<CLCT>&, std::vector<Cluster>&);	// input : (file prefix) + ".txt"	GEM Capable!!
	void WriteTxt(std::string&, std::vector<CLCT>&);
	// Writes Patterns (.pat) to be loaded to EmuBoard
	bool WritePat(std::string&, std::vector<CLCT>&);

	void ExtractHits(std::vector<CLCT>& clcts, std::vector<Hit>& hits, int feb = -1);

	// Print object data to console
	void DumpGroup(Group grp, int Bx);
	void DumpHit(Hit&);

	void PrintCLCT(CLCT&);
	
	void TMB_Check(std::vector<CLCT>&, std::string&);

//	int IsConflict(Hit&, Hit&);		// EVEN => NO CONFLICT
//	int IsConflict(CLCT&, CLCT&);


}
