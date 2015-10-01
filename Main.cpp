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

		input->preprocess();

		vector<Bed*> * bedsByChr = new vector<Bed *>[input->getBedNum()];
		vector<Wig*> wigsByChr;

		input->split(bedsByChr, wigsByChr);

		map < string, vector<int> * > * chrs = input->getCommonChrs();

		ThreadInfo * threads = input->getThreadInfo();

		MetaplotRegion * result = calculate(input->getBedNum(), input->getMaxWindow(), bedsByChr, wigsByChr, threads, chrs);
		string outfile = printResults(result, input->getBedNum(), input->getNameStr(), input->getNameStrR(), input->getMaxWindow()); 

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
