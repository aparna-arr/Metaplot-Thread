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
	queueAr = new double[numShifts];

	for (int i = 0; i < numShifts; i++)
		queueAr[i] = 0;
	nextEmpty = 0;
}

// 1st window 
void RunningAvg::populate(double value)
{
	for (int i = 0; i < numShifts; i++)	
		queueAr[i] = value;

	nextEmpty = 0;
}

void RunningAvg::queue(int numQueues, double val)
{
	for (int i = 0; i < numQueues; i++)
	{
		currentAvg = currentAvg - queueAr[nextEmpty] + val;

		queueAr[nextEmpty] = val;
		nextEmpty = (nextEmpty + 1) % numShifts;
	}
}

void RunningAvg::next(int num)
{
	nextEmpty = (nextEmpty + num) % numShifts;
}

double RunningAvg::getAvg(void)
{
	return currentAvg;
}

/* end RunningAvg functions */

// Maybe this should be another class: Processing?

/* begin Process functions */

Process::Process(int bedNum, int maxWindow, int shift, int window, vector<Bed *> * &bedsByChr, vector<Wig *> &wigsByChr, ThreadInfo * threads, map<string, vector<int> *> * &chrs)
{
	bedNum = bedNum;
	maxWindow = maxWindow;
	shift = shift;
	window = window;
	chromsPerThread = chrs->size() / threads->chromThreads;
	bedsPerThread = bedNum / threads->bedThreads;
	divsPerThread = threads->divThreads;
	bedsByChr = bedsByChr;
	threads = threads;
	wigsByChr = wigsByChr;
	chrs = chrs;
	region = new MetaplotRegion*[bedNum];

	for (int i = 0; i < bedNum; i++)
		region[i] = new MetaplotRegion(maxWindow);

	makeChromVec();
}

void Process::makeChromVec(void)
{
	for (map<string, vector<int> * >::iterator iter = chrs->begin(); iter != chrs->end(); iter++)
		chromVec.push_back(iter->first);
}

MetaplotRegion * Process::mergeMetaplotRegions(MetaplotRegion ** array, int numRegions)
{
	MetaplotRegion * merge = new MetaplotRegion(maxWindow);	

	for (int i = 0; i < numRegions; i++)
	{
		for (int j = 0; j < maxWindow; j++)
		{
			merge->basePairs[j][0] += array[i]->basePairs[j][0];
			merge->basePairs[j][1] += array[i]->basePairs[j][1];
		}
	}
	
	return merge;
}

void Process::mergeMetaplotRegionsByBed(MetaplotRegion ** &oldArray, MetaplotRegion ** newArray)
{
	for (int i = 0; i < bedNum; i++)
	{
		for (int j = 0; j < maxWindow; j++)
		{
			oldArray[i]->basePairs[j][0] += newArray[i]->basePairs[j][0];
			oldArray[i]->basePairs[j][1] += newArray[i]->basePairs[j][1];
		}
	}
}

void Process::calculate(void)
{

//	if (threads->bedThreads == 0)
	if (threads->chromThreads == 0)
	{
		chromThread(chromVec.begin(), chromVec.end(), region); 
	}
	else
	{
		thread threadAr[threads->chromThreads];
//		MetaplotRegion ** regionAr = new MetaplotRegion*[threads->chromThreads];		

//		thread threadAr[threads->bedThreads];
//		MetaplotRegion ** regionAr = new MetaplotRegion*[bedNum];

//		for (int k = 0; k < bedNum; k++)
//			regionAr[k] = new MetaplotRegion(maxWindow);

		MetaplotRegion *** chromRegion = new MetaplotRegion**[threads->chromThreads];

		for (int i = 0; i < threads->chromThreads; i++)
		{
//			regionAr[i] = new MetaplotRegion();
			chromRegion[i] = new MetaplotRegion*[bedNum];
		
			for (int k = 0; k < bedNum; k++) // VERY EXPENSIVE init for Monte Carlo! Horrible! O(numChrom * numBed) so basically n^2!
				chromRegion[i][k] = new MetaplotRegion(maxWindow);


			threadAr[i] = thread(&Process::chromThread, this, chromVec.begin() + chromsPerThread * i, chromVec.begin() + chromsPerThread * (i + 1), ref(chromRegion[i]));
		}

		for (int j = 0; j < threads->chromThreads; j++)
			threadAr[j].join();

//		region = mergeMetaplotRegions(regionAr, threads->chromThreads);
		for (int m = 0; m < threads->chromThreads; m++)
			mergeMetaplotRegionsByBed(region, chromRegion[m]);				

	}
}

void Process::chromThread(vector<string>::iterator chromStart, vector<string>::iterator chromEnd, MetaplotRegion ** &regionMerge)
{
//	MetaplotRegion ** chromRegionMerge = new MetaplotRegion*[chromsPerThread];
	if (threads->bedThreads == 0)
	{
		for (vector<string>::iterator iter = chromStart; iter != chromEnd; iter++)
		{
//			chromRegionMerge[iter - chromStart] = new MetaplotRegion(maxWindow);

			MetaplotRegion ** bedRegion = new MetaplotRegion*[bedNum]; // for threadsafe

			int wigChrPos = *(chrs->find(*iter)->second->begin());
			wigsByChr[wigChrPos]->readPeaks();

			Bed ** bedAr = new Bed*[bedNum];

			for (vector<int>::iterator it = chrs->find(*iter)->second->begin() + 1; it != chrs->find(*iter)->second->end(); it++)
			{
				int index = it - chrs->find(*iter)->second->begin()-1;
				bedAr[index] = *(bedsByChr[index].begin() + *it);
			}

			bedThread(bedAr, wigsByChr.begin() + wigChrPos, bedRegion, 0);
			mergeMetaplotRegionsByBed(regionMerge, bedRegion);
		}

	}
	else
	{
//		for (int i = 0; i < bedNum; i++)
//			bedRegion[i] = new MetaplotRegion(maxWindow);


		for (vector<string>::iterator iter = chromStart; iter != chromEnd; iter++)
		{
//			chromRegionMerge[iter - chromStart] = new MetaplotRegion(maxWindow);
			

			MetaplotRegion ** bedRegion = new MetaplotRegion*[bedNum];
			int wigChrPos = *(chrs->find(*iter)->second->begin());
			wigsByChr[wigChrPos]->readPeaks();
			
			Bed *** bedAr = new Bed**[threads->bedThreads];
			 
			for (int i = 0; i < threads->bedThreads; i++)
			{
				bedAr[i] = new Bed*[bedsPerThread];
				for (vector<int>::iterator it = chrs->find(*iter)->second->begin() + (i + 1) * bedsPerThread; it != chrs->find(*iter)->second->begin() + (i + 1 + 1) * bedsPerThread; it++)
				{
					int index = it - chrs->find(*iter)->second->begin()-1;
					bedAr[i][index] = *(bedsByChr[index].begin() + *it);		
				}
			}
			

			thread threadAr[threads->bedThreads];
//			MetaplotRegion ** regionAr = new MetaplotRegion*[threads->bedThreads];
			for (int i = 0; i < threads->bedThreads; i++)
			{
//				regionAr[i] = new MetaplotRegion(maxWindow);
//				threadAr[i] = thread(&Process::bedThread, this, bedAr[i], wigsByChr.begin() + wigChrPos, ref(regionAr[i]));	
				threadAr[i] = thread(&Process::bedThread, this, bedAr[i], wigsByChr.begin() + wigChrPos, ref(bedRegion), i);	
			}

			for (int i = 0; i < threads->bedThreads; i++)
				threadAr[i].join();

//			chromRegionMerge[iter - chromStart] = mergeMetaplotRegions(regionAr, threads->bedThreads); 
/*
			if (iter != chromStart)
			{
				// merge this and the previous chr
				for (int i = 0; i < bedNum; i++)
				{
					mergeMetaplotRegionsByBed(bedRegion[0], bedRegion[1]);
				}

				delete bedRegion[1];
			}
*/
			mergeMetaplotRegionsByBed(regionMerge, bedRegion);
		}

		
	}

//	regionMerge = mergeMetaplotRegions(chromRegionMerge, chromsPerThread);
}

void Process::bedThread(Bed ** bedAr, vector<Wig *>::iterator wigIter, MetaplotRegion ** &regionMerge, int whichSlice)
{
	if (threads->divThreads == 0)
	{
//		MetaplotRegion ** bedRegionMerge = new MetaplotRegion*[bedsPerThread];

		for (int i = 0; i < bedsPerThread; i++)
		{
//			bedRegionMerge[i] = new MetaplotRegion(maxWindow);
			regionMerge[i + whichSlice] = new MetaplotRegion(maxWindow);

			bedAr[i]->readPeaks();
			
			vector<Peak>::iterator startBed, endBed, startWig, endWig;
			bedAr[i]->getPeakDiv(0, 0, startBed, endBed);
			(*wigIter)->getPeakDiv(startBed->start, endBed->start, startWig, endWig);
//			divThread(startBed, endBed, startWig, endWig, bedRegionMerge[i]);

			divThread(startBed, endBed, startWig, endWig, regionMerge[i + whichSlice]);
		}

//		regionMerge = mergeMetaplotRegions(bedRegionMerge, bedsPerThread);
	}
	else
	{
//		MetaplotRegion ** bedRegionMerge = new MetaplotRegion*[bedsPerThread];

		for (int i = 0; i < bedsPerThread; i++)
		{
			regionMerge[i + whichSlice] = new MetaplotRegion(maxWindow);

			MetaplotRegion ** divRegions = new MetaplotRegion*[threads->divThreads];
			bedAr[i]->readPeaks();
			thread threadAr[threads->divThreads];
			
			for (int j = 0; j < threads->divThreads; j++)
			{
				vector<Peak>::iterator startBed, endBed, startWig, endWig;
				bedAr[i]->getPeakDiv(threads->divThreads, j, startBed, endBed);
				(*wigIter)->getPeakDiv(startBed->start, endBed->start, startWig, endWig);

				threadAr[i] = thread(&Process::divThread, this, startBed, endBed, startWig, endWig, ref(divRegions[j])); 
			}

			for (int k = 0; k < threads->divThreads; k++)
				threadAr[k].join();

//			bedRegionMerge[i] = mergeMetaplotRegions(divRegions, threads->divThreads);
			regionMerge[i + whichSlice] = mergeMetaplotRegions(divRegions, threads->divThreads);
		}

//		regionMerge = mergeMetaplotRegions(bedRegionMerge, bedsPerThread);
	}
}


void Process::divThread(vector<Peak>::iterator bedPeakStart, vector<Peak>::iterator bedPeakEnd, vector<Peak>::iterator wigPeakStart, vector<Peak>::iterator wigPeakEnd, MetaplotRegion * &regionMerge)
{
	
	for (vector<Peak>::iterator iter = bedPeakStart; iter != bedPeakEnd; iter++)
	{
		// init to values that they could not possibly be in the wig ar
		vector<Peak>::iterator thisWigPeakStart = iter; 
		vector<Peak>::iterator thisWigPeakEnd = iter;

		vector<Peak>::iterator min = wigPeakStart;
		vector<Peak>::iterator max = wigPeakEnd;
		vector<Peak>::iterator mid = wigPeakEnd; // init only

		while (min < max) // note: modified from <= to <
		{
			mid = (max - min) / 2 + min;
	
			if (mid->start == iter->start)
				break;
			else if (mid->start < iter->start)
				min = mid + 1;
			else if (mid->start > iter->start)
				max = mid - 1;
		} // will almost never get an exact match

		if (mid->end < iter->start)
			mid++;
		if (mid != wigPeakStart && (mid - 1)->end > iter->start)
			mid--;

		thisWigPeakStart = mid;

		// bsearch #2 for end
		min = mid;
		max = wigPeakEnd;
		mid = wigPeakEnd;
	
		while (min < max)
		{
			mid = (max - min) / 2 + min;
	
			if (mid->end == iter->end)
				break;
			else if (mid->end < iter->end)
				min = mid + 1;
			else if (mid->end > iter->end)
				max = mid - 1;
		}

		if (mid->start > iter->end)
			mid--;
		if (mid != wigPeakEnd && mid->end < iter->end)
			mid++; 
	
		thisWigPeakEnd = mid;		
	
		mapWig(thisWigPeakStart, thisWigPeakEnd, iter->start, iter->strand, regionMerge);	
	} // for	
}

/* FIXME THIS ENTIRE ALGORITHM RESTS ON THE FACT THAT WINDOW % SHIFT == 0
 * EVERYTHING WILL BREAK IF THIS IS NOT TRUE */

void Process::mapWig(vector<Peak>::iterator wigPeakStart, vector<Peak>::iterator wigPeakEnd, int bedStart, char bedStrand, MetaplotRegion * &region)
{
	for (vector<Peak>::iterator iter = wigPeakStart; iter != wigPeakEnd; iter++)
	{	
		RunningAvg avg(window, shift);
		
		// set backIter and fwdIter
		vector<Peak>::iterator backIter = iter;
		vector<Peak>::iterator fwdIter = iter + 1;
		vector<Peak>::iterator fwdIterStop = iter;
		
		while (backIter != wigPeakStart && backIter->end > (iter->start - window))
			backIter--;

		if (backIter->end < iter->start - window)
			backIter++;

		while(fwdIterStop != wigPeakEnd && fwdIterStop->start < (iter->end + window))
			fwdIterStop++;

		if (fwdIterStop->start > iter->end + window)
			fwdIterStop--;

		int pos = iter->start - window;

		// init queue
		// can't addSignal until queue is fully init
		while (backIter != iter)
		{
			int backSignalLen = backIter->end - pos;
			
			// #1
			avg.queue(backSignalLen / shift, backIter->value * shift);
			
			// #2
			avg.queue((backSignalLen - (backSignalLen / shift) * shift + (shift-1) ) / shift, backIter->value * (backSignalLen % shift)); 
			
			int shift_end_num = (backSignalLen - 1) / shift + 1;	
			int shift_end_pos = pos + shift_end_num * shift;

			backIter++;	
			
			// #3
			avg.next((backIter->start - shift_end_pos) / shift);

			// #4
			// if we have hit peakStart, this will not add anything to queue b/c it is always a multiple of shift away from start - window.
			avg.queue( ( (backIter->start - shift_end_pos) + (shift - 1) ) / shift, backIter->value * (shift - ( (backIter->start - shift_end_pos) % shift) ) );
		} // while

		// slide forward until end of window == peakStart

		// start addSignal here
		region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);

		int numSlidingWindow = 0;
		pos += shift;

		// FIXME check this	
		int nextOverlap = iter->end - window;
 
		if (fwdIterStop != iter)
			nextOverlap = fwdIterStop->start - window;

		if (nextOverlap < iter->start)
		{ // #4a
			numSlidingWindow = (nextOverlap + window - iter->start) / shift;

			for (int i = 0; i < numSlidingWindow; i++)	
			{
				avg.queue(1, iter->value * shift);
				region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);
				pos += shift;
			}
		}
		else
		{
			// #4a
			numSlidingWindow = window / shift;
			
			for (int i = 0; i < numSlidingWindow; i++)
			{
				avg.queue(1, iter->value * shift);
				region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);
				pos += shift;
			}

			// #5
			int nextEnd = iter->end;

			if (fwdIterStop != iter)
				nextEnd = fwdIterStop->start - window;				

			int middle_len = (nextEnd - pos) / shift;

			region->addSignal(pos - bedStart, middle_len, avg.getAvg(), bedStrand);
	
			pos += middle_len * shift;
		}

		if (fwdIterStop == iter)
		{
			// there is no nextIter that overlaps
			int shift_start = numSlidingWindow * shift + iter->start;
			avg.queue( (iter->end - shift_start + (shift - 1) ) / shift, ((iter->end - shift_start) % shift ) * iter->value); // FIXME check this!
			region->addSignal( pos - bedStart, 1 * ((iter->end - shift_start + shift - 1) / shift), iter->value * ((iter -> end - shift_start) % shift), bedStrand);

			while (pos < shift_start) // ??
			{
				avg.queue(1, 0);
				region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);
				pos += shift;
			} 	

		}
/* Don't actually need all this as it is taken care of with the next peak ?? only if there IS a next peak*/
//		int shift_start = numSlidingWindow * shift + iter->start;
//		pos = iter->end; // FIXME reusing variable differently!
/*
		while(fwdIter != fwdIterStop + 1 && fwdIter->start < shift_start + window)	
		{
			// #6
			avg.queue( (pos - shift_start + (shift - 1) ) / shift, ((pos - shift_start) % shift ) * iter->value);	
			// FIXME addsignal!
			int shift_end = shift_start + shift;

			// #7
			// FIXME make sure fwdIter is always > iter
		//	avg.next( (fwdIter->start - shift_end ) / shift);

			//in the blocks that follow is shift_end being used correctly in addSignal???
			int numZeros = ( fwdIter->start - shift_end ) / shift;

			for (int i = 0; i < numZeros; i++)
			{
				avg.queue(1,0);

				region->addSignal( shift_end - window - bedStart, 1, avg.getAvg(), bedStrand); 
				shift_end += shift;
			}

			// #8	
			int conditional = ( (fwdIter->start - shift_end) * shift + (shift - 1) ) / shift;

			avg.queue( conditional, ((fwdIter->start - shift_end) % shift) * fwdIter->value);
			region->addSignal( shift_end - window - bedStart, conditional, avg.getAvg(), bedStrand);

			shift_end += shift * conditional; // should be same as end_pos
				

			// #9
			// this is wrong
			if (fwdIter != fwdIterStop + 1)
			{
				int numShifts = (fwdIter->end - shift_end) / shift;

				for (int i = 0; i < numShifts; i++)
				{
					avg.queue(1, fwdIter->value * shift);
					region->addSignal(shift_end - window - bedStart, 1, avg.getAvg(), bedStrand);
				}
				fwdIter++;
			}
			else
			{
				int numShifts = ( shift_start + window - shift_end ) / shift;

				for (int i = 0; i < numShifts; i++)
				{
					avg.queue(1, fwdIter->value * shift);
					region->addSignal(shift_end - window - bedStart, 1, avg.getAvg(), bedStrand);
					shift_end += shift;
				}
			}
				
		}

		
*/	
	} // for
}

string Process::printResults(string nameStr, string nameStrR)
{
	print("Printing out results");
	
	string outfileName = "metaplot_outfile.txt";
	ofstream file(outfileName.c_str());
	
	string header = "bp\t" + nameStr + "\n";
	file << header;

	int h = 0 - (maxWindow / 2);
	
	for (int i = 0; i < maxWindow; i++)
	{
		int addition = h + 1;
		stringstream hToString;
		hToString << addition;
		file << hToString.str() << "\t";

		for (int j = 0; i < bedNum; j++)
		{
		}
	} // for
}
/*
string printResults(MetaplotRegion * result, int bedNum, std::string nameStr, std::string nameStrR, int maxWindow)
{

}
*/
void monteCarlo(std::string outfile, int bedNum)
{

}

/* end Process functions */

/*
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
*/
/*
MetaplotRegion * map(int maxWindow, vector<Peak>::iterator wigStart, vector<Peak>::iterator wigEnd, int bedStart, char bedStrand, int window, int shift)
{
	MetaplotRegion * region;
	
	for(vector<Peak>::iterator iter = wigStart; iter != wigEnd; iter++)
	{
		RunningAvg avg(window, shift);
		int pos = iter->start - (window - 1);
		
		// this is wrong
//		avg.populate( (double)iter->value );
		
		vector<Peak> backIter = iter;

		if (iter != wigStart)
		{
			while( backIter != wigStart && backIter->end > pos )
				backIter--;

			if (backIter->end < pos)
				backIter++;
		}

		// initialize avg with first window
		while (backIter != iter && pos < iter->start)
		{
			int mod = shift;
			if (backIter->start < pos)
			{
				for (int i = 0; i < (backIter->end - pos) / shift; i++)
					avg.queue(backIter->value * shift);	
				
				// DO NOT remove /shift * shift until you understand why it is there! (hint: INTEGER MATH)
				pos += (backIter->end - pos) / shift * shift; 
	
				if ((mod = (backIter-> end - pos) % shift) != 0)
				{
					avg.queue(backIter->value * mod);
					pos += shift;
				}
			}
			
			backIter++;
	
			// THIS WORKS PROBABLY DO NOT MODIFY UNTIL YOU UNDERSTAND WHY
			// can make this into a single arithmetic operation
//for (int i = 0; i < (backIter->start - (backIter - 1)->end + (shift - mod)) / shift; i++)
//				avg.next();
			int weirdEq = backIter->start - (backIter - 1)->end + (shift - mod);

			avg.next(weirdEq / shift);

			pos += weirdEq / shift * shift;

			if (weirdEq % shift != 0)
			{
				avg.queue(backIter->value * (weirdEq % shift));
				pos += shift;
			}
		}
			
		// pos SHOULD be within [start, start + shift) here. 

		while (pos - (window-1) < iter->start && pos + shift <= iter->end)
		{
			avg->queue(iter->value * shift);
			region->addSignal(pos - bedStart, shift, (double) (avg->getAvg) / (double)window, bedStrand);
			pos += shift;
		}		

		// now calculate over peak
		if (pos + shift <= iter->end)
		{
				
		}
			
					
	} // for
} 
*/

/*
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
*/
/*
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
*/
