#include"MetaplotThread.h"
using namespace std;

void print(string msg)
{
	lock_guard<mutex> lock(printMutex);
	cout << msg << endl; 
} 

/* RunningAvg functions */

RunningAvg::RunningAvg()
{
	
}

RunningAvg::RunningAvg(int window, int shift)
{
	// FIXME expect window % shift == 0
	numShifts = window / shift;
	queue = new double[numShifts];

	for (int i = 0; i < numShifts; i++)
		queue[i] = 0;
	nextEmpty = 0;
	first = 0;
}

// 1st window 
void RunningAvg::populate(double value)
{
	for (int i = 0; i < numShifts; i++)	
		queue[i] = value;

	nextEmpty = 0;
	first = 0;
}

void RunningAvg::queue(double val)
{
	currentAvg = currentAvg - queue[nextEmpty] + val;

	queue[nextEmpty] = val;
	nextEmpty = numShifts % (nextEmpty + 1);
}

// if circular, no need
void RunningAvg::dequeue(void)
{

}

double RunningAvg::getAvg(void)
{
	return currentAvg;
}

/* end RunningAvg functions */

// Maybe this should be another class: Processing?

MetaplotRegion * calculate(int bedNum, int maxWindow, std::vector<Bed *> * bedsByChr, std::vector<Wig *> wigsByChr, ThreadInfo * threads, std::map< std::string, std::vector<int> * > * chrs)
{
	// FIXME create chromVec at same time as populating the map, so we don't have this extra O(numChroms) step
	vector<string> chromVec;
	
	for (map<string, vector<int> *>::iterator iter = chrs->begin(); iter != chrs->end(); iter++)
		chromVec.push_back(iter->first);

	if (threads->chromThreads == 0) // do not thread chroms
	{
		chromThread(bedNum, maxWindow, bedsByChr, wigsByChr, threads, chrs, chromVec, 0, chrs->size());
	}
	else
	{
		// thread chroms
		thread threadAr[threads->chromThreads];
		int chromsPerThread = chrs->size() / threads->chromThreads;
		
		for (int i = 0; i < threads->chromThreads; i++)
			threadAr[i] = thread(chromThread, bedNum, maxWindow, bedsByChr, wigsByChr, threads, chrs, chromVec, chromsPerThread * i, chromsPerThread);

		for (int j = 0; j < threads->chromThreads; j++)
			threadAr[j].join();

		// join metaplot regions
	}			
	print("Done with calculate()");
}

void chromThread(int bedNum, int maxWindow, std::vector<Bed *> * bedsByChr, std::vector<Wig *> wigsByChr, ThreadInfo * threads, map< string, vector<int> *> *chrs, vector<string> chromVec, int startChrom, int chromsPerThread)
{
	if (threads->bedThreads == 0) // do not thread beds
	{
		for (vector<string>::iterator iter = chromVec.begin() + startChrom; iter != chromVec.begin() + startChrom + chromsPerThread; iter++)
		{
			wigsByChr[*((chrs->find(*iter)->second)->begin())]->readPeaks();
			Bed ** bedAr = new Bed*[bedNum];

			for (int i = 0; i < bedNum; i++)
				for (vector<int>::iterator it = chrs->find(*iter)->second->begin() + 1; it != chrs->find(*iter)->second->end(); it++)
					bedAr[it - chrs->find(*iter)->second->begin()] = *(bedsByChr[i].begin() + *it);

			bedThread(maxWindow, bedAr, wigsByChr[*((chrs->find(*iter)->second)->begin())], threads->divThreads, 0, bedNum);	
		}
	}
	else 
	{
		// thread beds
		for (vector<string>::iterator iter = chromVec.begin() + startChrom; iter != chromVec.begin() + startChrom + chromsPerThread; iter++)
		{
			wigsByChr[*((chrs->find(*iter)->second)->begin())]->readPeaks();

			Bed ** bedAr = new Bed*[bedNum];
			
			for (int i = 0; i < bedNum; i++)
				for (vector<int>::iterator it = chrs->find(*iter)->second->begin() + 1; it != chrs->find(*iter)->second->end(); it++)
					bedAr[it - chrs->find(*iter)->second->begin()] = *(bedsByChr[i].begin() + *it);
				
			thread threadAr[threads->bedThreads];
			int bedsPerThread = bedNum / threads->bedThreads;
	
			for (int i = 0; i < threads->bedThreads; i++)
				threadAr[i] = thread(bedThread, maxWindow, bedAr, wigsByChr[*((chrs->find(*iter)->second)->begin())], threads->divThreads, bedsPerThread * i, bedsPerThread);
	
			for (int j = 0; j < threads->bedThreads; j++)
				threadAr[j].join();
	
			// join metaplot regions
		}
	}
	print("Done with chromThread");
}

void bedThread(int maxWindow, Bed ** bedAr, Wig * wig, int divThreads, int startBed, int bedsPerThread)
{
	if (divThreads == 0) // do not thread divs
	{
		for (int i = startBed; i < startBed + bedsPerThread; i++)
		{
			bedAr[i]->readPeaks();	
			vector<Peak>::iterator startBed, endBed, startWig, endWig;
			bedAr[i]->getPeakDiv(1, 0, startBed, endBed); // FIXME will this work or cut off the last peak?
			wig->getPeakDiv(startBed->start, endBed->start, startWig, endWig);
			divThread(maxWindow, startBed, endBed, startWig, endWig);
		}
	}
	else
	{
		for (int i = startBed; i < startBed + bedsPerThread; i++)
		{
			bedAr[i]->readPeaks();

			thread threadAr[divThreads];	
			// need Chromosome to have a vector of peaks, not a stack, so we can have random acces -> done
			for (int i = 0; i < divThreads; i++)
			{
				vector<Peak>::iterator startBed, endBed, startWig, endWig;
				bedAr[i]->getPeakDiv(divThreads, i, startBed, endBed);
				wig->getPeakDiv(startBed->start, endBed->start, startWig, endWig);
				threadAr[i] = thread(divThread, maxWindow, startBed, endBed, startWig, endWig);
			}
			
			for (int j = 0; j < divThreads; j++)
				threadAr[j].join();
		}	
	}
	print("Done with bedThread");
} 

void divThread(int maxWindow, vector<Peak>::iterator bedStart, vector<Peak>::iterator bedEnd, vector<Peak>::iterator wigStart, vector<Peak>::iterator wigEnd)
{
	for (vector<Peak>::iterator iter = bedStart; iter != bedEnd; iter++)
	{
		// map wigs & smooth
		// make list of wig peaks that fall on this bed
		// running avg across peaks in-place?x
		// addSignal
	} // for each bed peak	
} 

MetaplotRegion * map(int maxWindow, vector<Peak>::iterator wigStart, vector<Peak>::iterator wigEnd, int bedStart, char bedStrand, int window, int shift)
{
	MetaplotRegion * region;
	
	for(vector<Peak>::iterator iter = wigStart; iter != wigEnd; iter++)
	{
		RunningAvg avg(window, shift);
		int pos = iter->start - (window - 1);
		
		// this is wrong
		avg.populate( (double)iter->value );
		
		vector<Peak> backIter = iter;

		if (iter != wigStart)
		{
			while( backIter != wigStart && backIter->end > pos )
				backIter--;

			if (backIter->end < pos)
				backIter++;
		}

		while (backIter != iter)
		{
			if (backIter->start < pos)
			{
				for (int i = 0; i < (backIter->end - pos) / shift; i++)
					avg.queue(backIter->value * shift);	
			}
		}

		while (pos < iter->start)
		{
		} // while
			
					
	} // for
} 

string printResults(MetaplotRegion * result, int bedNum, std::string nameStr, std::string nameStrR, int maxWindow)
{

}

void monteCarlo(std::string outfile, int bedNum)
{

}
