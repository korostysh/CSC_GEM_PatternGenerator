#include "CLCT.h"
#include "TaoFunc.h"

namespace cw {
	int layerorder_all[6][6] = { {3,1,5,6,4,2},//layerorder1
						  {4,2,5,6,3,1},
						  {1,3,6,5,2,4},
						  {2,4,6,5,1,3},
						  {4,6,2,1,3,5},
						  {1,2,3,4,5,6}
	};

//~~~~~~~~~ Support Functions ~~~~~~~~~~~~
bool EdgeCheck(int key, int hs)						// used to check if generated hit is within CSC bounds
	{
		if (hs < 0 || hs > CSCConstants::NUM_HS) // falls outside hs range
			return false;
		else if (key <= 127 && hs > 127)	// crosses Edge of ME1/1b
			return false;
		else if (key >= 128 && hs < 128)	// crosses Edge of ME1/1a
			return false;
		else
			return true;
	}

bool Hits_Generator(int Bx, int Key, int Pat, int Nhits, std::vector<Hit>& hits)	// used to Fill Hits vector in CLCT constructor
	{
		std::vector<bool> usedlayers(CSCConstants::NUM_LAYERS, false);				//  used layer accounting

		// Pattern validity check
		if (Pat < 2 || Pat >= CSCConstants::NUM_PATTERNS)
			return true;  				// invalid pattern types
		if (Nhits < 1 || Nhits > CSCConstants::NUM_LAYERS)
			return true;				// invalid number of hits

		std::srand(std::time(0));

		int n = 0;
		while (n != Nhits)
		{
			int layer = std::abs(std::rand()) % CSCConstants::NUM_LAYERS;				// randomly select a layer
			if (usedlayers.at(layer)) continue;
			else usedlayers.at(layer) = true;			

			int hs = Key + CSCConstants::CLCTPatterns[Pat][layer];      				// selects halfstrip hit from pattern

			if (EdgeCheck(Key, hs))	hits.push_back(Hit(Bx, hs, layer));

			n++;
		}

		return false;
	}

int  GetCFEB(int hs)
	{ 
		return GetInputCFEBByX<STRIPS_CFEB>(hs, COMPILE_TYPE); 
	}

int  GetLocal(int hs)
	{ 
		return GetInputXStrip<STRIPS_CFEB>(hs, COMPILE_TYPE); 
	}

/// Check
unsigned char GetTriadSeg(Hit hit, int i)
	{	
		int localhs = GetLocal(hit.hs);
		unsigned char n = 1 << (localhs / 4);
		bool leftstrip = (localhs % 4 < 2);
		bool lefths = (localhs % 2 == 0);
		
		switch (i) 
		{
		case 0:
			return n;
			break;
		case 1:
			return (leftstrip ? 0x00 : n);
			break;
		case 2:
			return (lefths ? 0x00 : n);
			break;
		default:
			// Throw some kind of error?
			return (0x00);
			break;
		}
	}

//~~~~~~~~~~ Class Structures ~~~~~~~~~~~
Hit::Hit() : bx(0), hs(0), lay(0) {}

Hit::Hit(int Bx, int Hs, int Layer) : bx(Bx), hs(Hs), lay(Layer) {}


CLCT::CLCT(void) : bx(0), hs(0), pat(0), nhits(0) {}

CLCT::CLCT(int Bx, int Key, int Pat, int Nhit) :
	bx(Bx),
	hs(Key),
	pat(Pat),
	nhits(Nhit)
	{
		bool bad = Hits_Generator(bx, hs, pat, nhits, hits);
	}

CLCT::CLCT(const CLCT& cl) :
	bx(cl.bx),
	hs(cl.hs),
	pat(cl.pat),
	nhits(cl.nhits)
	{
		hits.insert( hits.end(), cl.hits.begin(), cl.hits.end() );
	}

std::ostream& operator<<(std::ostream& os, const CLCT& cl) 
	{
		os << cl.bx << '\t' << cl.hs << '\t'
			<< cl.pat << '\t' << cl.nhits;

		return os;
	}

std::istream& operator>>(std::istream& is, CLCT& cl) 
	{
		int B, K, P, N;
		is >> B >> K >> P >> N;

		cl = CLCT(B, K, P, N);
		return is;
	}

/// Check
Group::Group(void)// : 
	//hexdata( std::vector<unsigned char>(CSCConstants::NUM_LAYERS, unsigned char(0)) )
	{
		for(int i=0; i < CSCConstants::NUM_LAYERS; i++)
			{
				hexdata.push_back(char(0));
			}
		//hexdata = std::vector<unsigned char> (CSCConstants::NUM_LAYERS, unsigned char(0));
	}

/// Check
Group::Group(std::vector<Hit>& hits, int Bx) //:
	//hexdata(std::vector<unsigned char>(CSCConstants::NUM_LAYERS, unsigned char(0)))
{
	for(int i=0;i < CSCConstants::NUM_LAYERS; i++)
		{
			hexdata.push_back(char(0));
		}
	
	for (int i=0; i < hits.size(); i++) 
	{
		int delta_t = Bx - hits[i].bx;
		if(delta_t >= 0 && delta_t < 3)
		{
			hexdata[hits.at(i).lay] = hexdata[hits.at(i).lay] | GetTriadSeg(hits[i], delta_t);
		}
	}
}

/// Check
std::ostream& operator<<(std::ostream& os, const Group& grp)
	{
		for (int i = 0; i < grp.hexdata.size(); i++)
			{
				os << grp.hexdata[ layerorder_all[i][LAYER_SET] - 1 ];	// This replaces Tao's Shuffle Layers Function hope
			}

		os << char(0) << char(0);
		return os;
	}

void DumpGroup(Group grp, int Bx)
	{
		for (int i = 0; i < grp.hexdata.size(); i++)
			{
				std::cout << std::bitset<8>(grp.hexdata[layerorder_all[i][LAYER_SET] - 1]) << std::endl;
			}
		std::cout << std::bitset<8>(char(0)) << std::endl;
		std::cout << std::bitset<8>(char(0)) << "   Bx: " << Bx << std::endl << std::endl;

		return;
	}

void DumpHit(Hit& hit, int N)
	{
		if(N == -1)	std::cout << "Hit : \n";
		else		std::cout << "Hit : " << N;

		std::cout << '\t' << "bx: " << hit.bx
			  << '\t' << "hs: " << hit.hs
			  << '\t' << "layer: " << hit.lay << '\n';

		for(int i=0; i < 3; i++)	std::cout << "  " << std::bitset<8>(GetTriadSeg(hit,i)) << '\n';
		
		return;
	}


//~~~~~~~~ Functions ~~~~~~~~~~~~~~~~~~~

void PrintCLCT(CLCT& clct)
	{
		std::string empty_layer = "-----------";
		std::vector<std::string> layers(6,empty_layer);
		
		for(int i=0; i < clct.hits.size(); i++){
			layers.at(clct.hits.at(i).lay).at(5 + (clct.hits.at(i).hs - clct.hs)) = 'X';
			std::cout << "Hit (" << i << ')' << "   Bx: " << clct.hits.at(i).bx << "   Hs: " << clct.hits.at(i).hs << "   Layer: " << clct.hits.at(i).lay << "   CFEB: " << GetCFEB(clct.hits.at(i).hs) << std::endl;
		}

		std::cout << "bx keystrip pattern nhits" << std::endl;
                std::cout << clct << std::endl;


		for(int ly=0; ly < 6; ly++)	// for each layer
		{
			std::cout << "ly" << ly;
			if(ly == 2) std::cout << "  key  ";
			else 	    std::cout << "       ";
			
			std::cout << layers.at(ly) << std::endl;	
			
		}
	}

int ReadTxt(std::string& prefix, std::vector<CLCT>& clcts) 
	{
		CLCT cl;
		std::string str; 										// String to be used for parsing purposes
		std::fstream ifs((prefix + ".txt"), std::ios_base::in); // Open File to be read;

		// First two lines are header	
		std::getline(ifs, prefix);
		std::getline(ifs, str);
		

		///prefix = prefix.substr(7, prefix.length() - 7);

		int n = 0;
		while (ifs >> cl) {
			clcts.push_back(cl);
			n++;
			std::cout << clcts.at(n - 1) << " size :" << n << '\n';
		}

		return n;
	}

void WriteTxt(std::string& str, std::vector<CLCT>& clcts) 
	{
		std::fstream text_file((str + ".txt"), std::ios_base::out); // Create output file => (input string).txt
		//std::sort(clcts.begin(), clcts.end(), CompareCLCT);	// Sort the CLCT vector in case it's asorted

		// Add header stuff to pattern(.txt) file
		text_file << "prefix:" << str << std::endl;
		text_file << "bx keystrip pattern nhits" << std::endl;

		// print each CLCT to file + "\n" unless last CLCT then No "\n"
		for (int i = 0; i < clcts.size(); i++) {
			text_file << clcts.at(i);
			if (i < (clcts.size() - 1)) {
				text_file << '\n';
			}
		}

		text_file.close();
		return;
	}
/// Check
void ExtractHits(std::vector<CLCT>& clcts, std::vector<Hit>& hits, int feb)
	{
		if (feb == -1)
		{
			for (int i = 0; i < clcts.size(); i++)
			{
				hits.insert(hits.end(), clcts.at(i).hits.begin(), clcts.at(i).hits.end());
			}
		}
		/// DOUBLE CHECK THE REST OF THIS FUNCTION	**BAD PRACTICE**
		else if(feb >= 0 && feb < CSCConstants::NUM_DCFEBS)
		{
			for(int i=0; i < clcts.size(); i++)
			{
				for(int j=0; j < clcts.at(i).hits.size(); j++)
				{
					if (GetCFEB(clcts.at(i).hits.at(j).hs) == feb) { hits.push_back(clcts.at(i).hits.at(j)); }
				}
			}
		}
		return;
	}


void fillnbytes(std::vector<std::fstream*> oss, unsigned int n, unsigned int thisfeb) {

	for (unsigned int j = 0; j < n; j++)
		(*(oss[thisfeb])) << char(0);

}

void fillnbytes_allfebs(std::vector<std::fstream*> oss, unsigned int n) {

	for (unsigned int i = 0; i < oss.size(); i++)
		for (unsigned int j = 0; j < n; j++)
			(*(oss[i])) << char(0);

}

/// DOUBLE Check
bool WritePat(std::vector<CLCT>& clcts, std::string& prefix) 
	{
		// Prepare output file streams
		std::vector<std::fstream*> oss;
		char tmbtype = COMPILE_TYPE - 0xa + 'a';
		for (int i = 0; i < CSCConstants::NUM_DCFEBS; i++)
		{
			std::stringstream ss;
			ss << prefix << "_ClctPattern_CFEB" << i << "_tmb" << tmbtype << ".pat";
			//std::cout << ss.str().c_str() << std::endl;
			oss.push_back(new std::fstream(ss.str().c_str(), std::ios_base::out) );
		}

		for (int i = 0; i < CSCConstants::NUM_DCFEBS; i++)
		{
			std::cout << "Writing (.pat) file number: " << i << std::endl;	// debug

			std::vector<Hit> hits;
			std::vector<int> times;

			ExtractHits(clcts, hits, i);
			
			// DEBUG PURPOSES~~~~~
			for(int j=0; j < hits.size(); j++){
				DumpHit(hits.at(j), j);
			}
			// ~~~~~~~~~

			if(hits.size() == 0)
			{
				fillnbytes(oss, RAMPAGE_SIZE, i);
				(*(oss.at(i))) << std::flush;		// Fill & Close file
				delete (oss.at(i));
			}
			else
			{
				// Get times
				for (int j = 0; j < hits.size(); j++)
				{
					bool tf = false;
					for (int k = 0; k < times.size(); k++)
					{
						if (times.at(k) == hits.at(j).bx) tf = true;
					}
					if (!tf) times.push_back(hits.at(j).bx);
				}
				std::sort(times.begin(), times.end());

				int Bx = 0;
				int q = 0;

				if (times.at(0) != 0)
				{
					Bx = times.at(0);
					fillnbytes(oss, times.at(0)*8, i);
				}


				while (Bx < (times[times.size() -1] + 3) && (RAMPAGE_SIZE - Bx*8) > 0)
				{
					if (Bx == (times.at(q) + 3) && (i + 1) == times.size())
					{
						fillnbytes(oss, RAMPAGE_SIZE - (Bx * 8), i);		/// Double check range here
						break;
					}
					else if (Bx == (times.at(q) + 3))
					{
						q++;
						if(times.at(q) > Bx)
						{
							fillnbytes(oss, (times.at(q) - Bx)*8, i);
							Bx = times.at(q);
							// write Group @ this->Bx
							(*(oss[i])) << Group(hits,Bx);
							DumpGroup(Group(hits, Bx), Bx);	// debug
							Bx++;
						}
						else
						{
							//Write Group @ this->Bx
							(*(oss[i])) << Group(hits, Bx);
							DumpGroup(Group(hits, Bx), Bx);	// debug
							Bx++;
						}
					}
					else
					{
						//Write Group @ this->Bx
						(*(oss[i])) << Group(hits,Bx);
						DumpGroup(Group(hits, Bx), Bx);	// debug
						Bx++;
					}
				}
				if((RAMPAGE_SIZE - Bx*8) > 0)
				{
					fillnbytes(oss, RAMPAGE_SIZE - (Bx * 8), i);		/// Double check range here
				}
				
				*(oss.at(i)) << std::flush;
				delete (oss.at(i));

			}
			
			
		}
		return true;
	}


}
