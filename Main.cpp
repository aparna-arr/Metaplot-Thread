#include"MetaplotRegionThread.h"
#include"UserOptsThread.h"
#include"ChromosomeThread.h"
#include"MetaplotThread.h"
using namespace std;

int main(int argc, char * argv[])
{
	try 
	{ 
		UserOpts * input = new UserOpts(argc, argv);
		
//		print("main: before preprocess()");

		print("- Preprocessing files");
		input->preprocess();


//		print("main: after preprocess");

		vector<Bed*> * bedsByChr = new vector<Bed *>[input->getBedNum()];
		vector<Wig*> wigsByChr;
//		print("main: after making bedsByChr wigsByChr");

		print("- Splitting files");
		input->split(bedsByChr, wigsByChr);
//		print("main: after split");

		map < string, vector<int> * > * chrs = input->getCommonChrs();

//		print("main: after getCommonChrs");

		ThreadInfo * threads = input->getThreadInfo();

//		cerr << threads->chromThreads << endl;
		
//		print("main: after getThreadInfo");

//		MetaplotRegion * result = calculate(input->getBedNum(), input->getMaxWindow(), bedsByChr, wigsByChr, threads, chrs);
//		string outfile = printResults(result, input->getBedNum(), input->getNameStr(), input->getNameStrR(), input->getMaxWindow()); 

//		cerr << "shift is " << input->getShift() << endl;

		print("- Initialize Process Object");
		Process metaplot(input->getBedNum(), input->getMaxWindow(), input->getShift(), input->getWindow(), bedsByChr, wigsByChr, threads, chrs);

//		print("main: after init Process");

		print("- Beginning metaplot calculation");
		metaplot.calculate();

//		print("main: after calculate()");

		print("- Printing results");
		string outfile = metaplot.printResults(input->getNameStr(), input->getNameStrR());

		if (input->isMonteCarlo())
			monteCarlo(outfile, input->getBedNum());
	} // try
	catch (int e)
	{		
		if (e == 0)
			return 0;
		if (e == 1)
			return 1;
		if (e == 2)
		{
			string err = "INVALID ARG: Opt requires INT argument!";
			print(err);
			return 1;
		}
	} // catch
}
