#include"Metaplot.h"
#include"Chromosome.h"
#include"MetaplotRegion.h"
#include"UserOpts.h"
using namespace std;

int main(int argc, char * argv[])
{
	try{
		UserOpts input(argc, argv);
	
		cout << "Begin splitting bed files ..." << endl;
		vector<Bed *> * bedByChrs = input.splitBedFiles();
		cout << "\tdone." << endl;
		cout << "Begin splitting wig file ..." << endl;
		vector<Wig *> wigByChrs = input.splitWigFile();
		cout << "\tdone." << endl;

		vector<string> commonChrs = input.commonChroms();

		MetaplotRegion * region = new MetaplotRegion[input.getBedNumber()];

		for (int i = 0; i < input.getBedNumber(); i++)
			region[i] = MetaplotRegion(input.getMaxWindow());
		
		if (!input.monteCarlo())
		{
//			cerr << "DEBUG not MonteCarlo" << endl;
//			cerr << "DEBUG before calculateMetaplot" << endl;
			calculateMetaplot(input.getBedNumber(), bedByChrs, wigByChrs, commonChrs, region);
//			cerr << "DEBUG after calculateMetaplot" << endl;
//			cerr << "Input bed num is " << input.getBedNumber() << endl;
//			cerr << "Input name str is " << input.getNameString() << endl;
//			cerr << "Input name str R is " << input.getNameStringR() << endl;
//			cerr << "Input max window is " << input.getMaxWindow() << endl;

//			cerr << "region bed 0 01 is " << region[0].basePairs[0][1] << endl;
			string namestr = input.getNameString();
			string namestrR = input.getNameStringR();
// IMPORTANT: DO NOT DELETE
// for some reason, using getNameString() and getNameStringR() in a function call results in getNameString() not being called an dgetNameStringR() being empty. Assign to strings before hand and pass those strings in -> no segfault ...
//			debug(region, input.getBedNumber(), input.getNameString(), input.getNameStringR()); 
			printResults(region, input.getBedNumber(), namestr, namestrR, input.getMaxWindow()); 
//			cerr << "DEBUG after printResults" << endl;
		}
		else
		{
			// for some reason control ends up in HERE ...
//			cerr << "DEBUG in MonteCarlo" << endl;
			cout << endl << "Monte Carlo Metaplot" << endl;
			calculateMetaplot(input.getBedNumber(), bedByChrs, wigByChrs, commonChrs, region);
		
			string namestr = input.getNameString();
			string namestrR = input.getNameStringR();
			
			string outfile = printResults(region, input.getBedNumber(), namestr, namestrR, input.getMaxWindow()); 
			monteCarloMetaplot(outfile, input.getBedNumber());
		}

	} // try block end paren	
	catch(int e)
	{
		if (e == 1)
		{
			cerr << "Exiting!" << endl;
			return 1;
		}
	}
	return 0;
}

