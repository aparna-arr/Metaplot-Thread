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

Process::Process(int bedNumArg, int maxWindowArg, int shiftArg, int windowArg, vector<Bed *> * &bedsByChrArg, vector<Wig *> &wigsByChrArg, ThreadInfo * threadsArg, map<string, vector<int> *> * &commonChrs)
{
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("Process(): start" + debug);

	bedNum = bedNumArg;
	maxWindow = maxWindowArg;
	shift = shiftArg;
	window = windowArg;

	stringstream iToS4;
	iToS4 << window;
	print("Process(): window is " + iToS4.str() + debug);

//	print("Process(): segfault?" + debug);

//	chromsPerThread = (commonChrs->size() + 1) / (threadsArg->chromThreads + 1) ; // otherwise div by 0
	if (threadsArg->chromThreads == 0)
		chromsPerThread = commonChrs->size();
	else
		chromsPerThread = commonChrs->size() / threadsArg->chromThreads;

//	bedsPerThread = (bedNum + 1) / (threadsArg->bedThreads + 1) ; // otherwise div by 0

	if (threadsArg->bedThreads == 0)
		bedsPerThread = bedNum;
	else
		bedsPerThread = bedNum / threadsArg->bedThreads;


	divsPerThread = threadsArg->divThreads;

//	print("Process(): segfault?" + debug);

	bedsByChr = bedsByChrArg;
	threads = new ThreadInfo;
	threads->chromThreads = threadsArg->chromThreads;
	threads->bedThreads = threadsArg->bedThreads;
	threads->divThreads = threadsArg->divThreads;

//	(*threads) = (*threads);
	wigsByChr = wigsByChrArg;

	stringstream iToS;
	iToS << wigsByChr.size();
	print("Process(): wigsByChr.size() " + iToS.str() + debug);
	print("Process(): wigsByChr[0]->getChr() " + wigsByChr[0]->getChr() + debug);

//	*chrs = *chrs;
	region = new MetaplotRegion*[bedNum];

	chrs = new map< string, vector<int> * >;

	print("Process(): after init members" + debug);

	for (int i = 0; i < bedNum; i++)
		region[i] = new MetaplotRegion(maxWindow);

	print("Process(): after init region" + debug);

	makeChromVec(commonChrs);
/*
	stringstream iToS2, iToS3;
	iToS2 << chrs->find("chr1")->second->size();
	iToS3 << commonChrs->find("chr1")->second->size();

	print("Process(): size of chrs vec (chr1) " + iToS2.str() + debug);
	print("Process(): size of commonChrs vec (chr1) " + iToS3.str() + debug);
*/
	print("Process(): end" + debug);
}

void Process::makeChromVec(map< string, vector<int> * > * commonChrs)
{
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;
	print("makeChromVec(): start" + debug);

	for (map<string, vector<int> * >::iterator iter = commonChrs->begin(); iter != commonChrs->end(); iter++)
	{
		print("makeChromVec(): on chr " + iter->first + debug);
		chromVec.push_back(iter->first);

		print("makeChromVec(): after push back" + debug);
		
		(*chrs)[iter->first] = new vector<int>;

		for (vector<int>::iterator it = (iter->second)->begin(); it != (iter->second)->end(); it++)
		{
			(*chrs)[iter->first]->push_back(*it);

			stringstream iToS;
			iToS << *it;

			print("\tmakeChromVec(): on chr " + iter->first + " with index " + iToS.str() + debug);
		}

		print("makeChromVec(): after loop" + debug);
	}

	print("makeChromVec(): end" + debug);
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
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("calculate(): start" + debug);

//	if (threads->bedThreads == 0)
	if (threads->chromThreads == 0)
	{
		print("calculate(): no chrom thread start" + debug);
		chromThread(chromVec.begin(), chromVec.end(), region); 
		print("calculate(): no chrom thread end" + debug);
	}
	else
	{
		print("calculate(): chrom thread start" + debug);
		thread threadAr[threads->chromThreads];

		print("calculate(): made threadAr" + debug);
//		MetaplotRegion ** regionAr = new MetaplotRegion*[threads->chromThreads];		

//		thread threadAr[threads->bedThreads];
//		MetaplotRegion ** regionAr = new MetaplotRegion*[bedNum];

//		for (int k = 0; k < bedNum; k++)
//			regionAr[k] = new MetaplotRegion(maxWindow);

		MetaplotRegion *** chromRegion = new MetaplotRegion**[threads->chromThreads];

		print("calculate(): made chromRegion " + debug);

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
		
		print("calculate(): chrom thread end" + debug);

	}

	print("calculate(): end" + debug);
}

void Process::chromThread(vector<string>::iterator chromStart, vector<string>::iterator chromEnd, MetaplotRegion ** &regionMerge)
{
//	MetaplotRegion ** chromRegionMerge = new MetaplotRegion*[chromsPerThread];
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;
	print("chromThread(): start" + debug);

	if (threads->bedThreads == 0)
	{

		print("chromThread(): no bed thread start" + debug);
		for (vector<string>::iterator iter = chromStart; iter != chromEnd; iter++)
		{
			print ("\tchromThread(): for start chr is " + *iter + debug);
//			chromRegionMerge[iter - chromStart] = new MetaplotRegion(maxWindow);

			MetaplotRegion ** bedRegion = new MetaplotRegion*[bedNum]; // for threadsafe

			print("\tchromThread(): after bed region init" + debug);

			int wigChrPos = *(chrs->find(*iter)->second->begin());
			

			stringstream iToS;	
			iToS << wigChrPos;	
			print("\tchromThread(): after wigChrPos init is " + iToS.str() + debug);
			
			wigsByChr[wigChrPos]->readPeaks();

			print("\tchromThread(): after wig read peaks" + debug);

			Bed ** bedAr = new Bed*[bedNum];

			print("\tchromThread(): after bed init" + debug);

			stringstream iToS2;
			iToS2 << chrs->find(*iter)->second->size();

			print("\tchromThread(): chrs size is " + iToS2.str() + debug);

			for (vector<int>::iterator it = chrs->find(*iter)->second->begin() + 1; it != chrs->find(*iter)->second->end(); it++)
			{
				print("\t\tchromThread(): in second for start" + debug);
				int index = it - chrs->find(*iter)->second->begin()-1;

//				bedAr[index] = new Bed();

				stringstream iToS3;
				iToS3 << index;
				print("\t\tchromThread(): index is " + iToS3.str() + debug);

				print("\t\tchromThread(): chr " + (*(bedsByChr[index].begin() + *it))->getChr() + debug);

				bedAr[index] = *(bedsByChr[index].begin() + *it);
				print("\t\tchromThread(): in second for end" + debug);
			}

			print("\tchromThread(): calling bedThread()" + debug);
			bedThread(bedAr, wigsByChr.begin() + wigChrPos, bedRegion, 0);

			print("\tchromThread(): after call to bedThread()" + debug);

			mergeMetaplotRegionsByBed(regionMerge, bedRegion);

			print("\tchromThread(): for end" + debug);
		}
		print("chromThread(): no bed thread end" + debug);

	}
	else
	{
//		for (int i = 0; i < bedNum; i++)
//			bedRegion[i] = new MetaplotRegion(maxWindow);

		print("chromThread(): bed thread start" + debug);

		for (vector<string>::iterator iter = chromStart; iter != chromEnd; iter++)
		{
//			chromRegionMerge[iter - chromStart] = new MetaplotRegion(maxWindow);
			
			print("\tchromThread(): in first for start" + debug);
			MetaplotRegion ** bedRegion = new MetaplotRegion*[bedNum];
			int wigChrPos = *(chrs->find(*iter)->second->begin());
			wigsByChr[wigChrPos]->readPeaks();
			
			print("\tchromThread(): after wig read peaks" + debug);

			Bed *** bedAr = new Bed**[threads->bedThreads];
			 
			for (int i = 0; i < threads->bedThreads; i++)
			{
				print("\t\tchromThread(): start inner for #1" + debug);
				bedAr[i] = new Bed*[bedsPerThread];

				/* debug */
				for (vector<int>::iterator deb = chrs->find("chr1")->second->begin(); deb != chrs->find("chr1")->second->end(); deb++)
				{
					stringstream d;
					d << deb - (chrs->find("chr1")->second->begin());
					print("\t\t\tchomThread(): debug index for chr1 is " + d.str() + debug);
				}


				for (vector<int>::iterator it = chrs->find(*iter)->second->begin() + (i * bedsPerThread) + 1; it != chrs->find(*iter)->second->begin() + (i + 1) * bedsPerThread + 1; it++)
				{

					print("\t\t\tchromThread(): start innermost for" + debug);
					int index = it - (chrs->find(*iter)->second->begin() + i * bedsPerThread + 1);

					print("\t\t\tchromThread(): found index for chr " + *iter + debug);
					stringstream indexToS, iToS, itDebug;
					indexToS << index;
					iToS << i;

					stringstream itFromStart;
					itFromStart << it - chrs->find(*iter)->second->begin();
				// FIXME HERE IT IS! MEMORY LEAK AGAIN
				// chromThread(): index is 1 i is 0 it 11053 My id: 47475403069440 
					itDebug << *it;
					print("\t\t\tchromThread(): index is " + indexToS.str() + " i is " + iToS.str() + " it " + itDebug.str() + " it from start " + itFromStart.str() + debug);

// FIXME the logic in this loop is wrong re: threading
// somehow
					stringstream debugSize;
					debugSize << bedsByChr[index + i].size();

					print("\t\t\tchromThread(): bedsbyChr size " + debugSize.str() + debug);

					bedAr[i][index] = *(bedsByChr[index + i].begin() + *it);

					stringstream iToStr, indexToStr;
					iToStr << i;
					indexToStr << index;

				// FIXME all threads get the same file!!! WHY -> made bedsByChr[index] to [index + i]
					print("\t\t\tchromThread(): for i is " + iToStr.str() + " index is " + indexToStr.str() + " file is " + (*(bedsByChr[index + i].begin()))->getFilename() + debug);	
//					print("\t\t\t bedsByChr[index].begin() + *it -> getChr() is " + (*(bedsByChr[index].begin() + *it))->getChr() + debug);	
					print("\t\t\tchromThread(): end innermost for" + debug);
				}

				print("\t\tchromThread(): end inner for #1" + debug);
			}
		
			print("\tchromThread(): after inner for #1" + debug);	

			thread threadAr[threads->bedThreads];
//			MetaplotRegion ** regionAr = new MetaplotRegion*[threads->bedThreads];

			for (int d = 0; d < bedsPerThread; d++)
			{
				print("\t\tchromThread(): debug chr is " + bedAr[0][d]->getChr());
			}



			print("\tchromThread(): before fill threadAr" + debug);
			for (int i = 0; i < threads->bedThreads; i++)
			{
//				regionAr[i] = new MetaplotRegion(maxWindow);
//				threadAr[i] = thread(&Process::bedThread, this, bedAr[i], wigsByChr.begin() + wigChrPos, ref(regionAr[i]));	
				threadAr[i] = thread(&Process::bedThread, this, bedAr[i], wigsByChr.begin() + wigChrPos, ref(bedRegion), i);	
			}

			print("\tchromThread(): before join()" + debug);

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
			print("\tchromThread(): before merge" + debug);
			mergeMetaplotRegionsByBed(regionMerge, bedRegion);
			print("\tchromThread(): in first for end" + debug);
		}

		print("chromThread(): bed thread end" + debug);
		
	}

	print("chromThread(): end" + debug);
//	regionMerge = mergeMetaplotRegions(chromRegionMerge, chromsPerThread);
}

void Process::bedThread(Bed ** bedAr, vector<Wig *>::iterator wigIter, MetaplotRegion ** &regionMerge, int whichSlice)
{
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("bedThread(): start" + debug);

	if (threads->divThreads == 0)
	{
		print("bedThread(): no div thread start" + debug);
//		MetaplotRegion ** bedRegionMerge = new MetaplotRegion*[bedsPerThread];

		for (int i = 0; i < bedsPerThread; i++)
		{
			print ("bedThread(): in for " + debug);
			Bed * bed;
			Wig * wig;

			stringstream iToStr, bpt, sliceToStr, startStr;
			iToStr << i;

			bpt << bedsPerThread;

			sliceToStr << whichSlice;

			
			print("\tbedThread(): in for start i is " + iToStr.str() + " beds per thread: " + bpt.str() + " whichSlice : " + sliceToStr.str() +  debug);
//			bedRegionMerge[i] = new MetaplotRegion(maxWindow);
			regionMerge[i + whichSlice] = new MetaplotRegion(maxWindow);

			print("\tbedThread(): after regionMerge init" + debug);

			print("\tbedThread(): bedAr[i]->getChr() is " + bedAr[i]->getChr() + debug);

			bedAr[i]->readPeaks();
			startStr << bedAr[i]->firstPeak()->start;

			print("\tbedThread(): after readPeaks() firstPeak start is " + startStr.str() + debug);
			
			vector<Peak>::iterator startBed, endBed, startWig, endWig;
			bedAr[i]->getPeakDiv(0, 0, startBed, endBed, bed);
//				stringstream bedstarttest, bedendtest, bedEndEnd;
//				bedstarttest << startBed->start;
//				bedendtest << endBed->start;
//				bedEndEnd << endBed->end;
//				print("\t\tbedThread(): startbed start is " + bedstarttest.str() + " endbed start is " + bedendtest.str() + " bedEnd end " + bedEndEnd.str() + debug);
//				print("\t\tbedThread(): startbed start is " + bedstarttest.str() + debug);
			print("\tbedThread(): after bedAr->getPeakDiv()" + debug);

// Issue:
// *** Error in `	bedThread(): after bedAr->getPeakDiv() My id: 47642893563648 
// 	bedThread(): endPeak() == peaks.end() My id: 47642893563648 
// 	bedThread(): endPeak() == peaks.end() My id: 47642893563648 
//  --threadBeds 2
//  appears that a thread is duplicated somehow?

			if (bed->endPeak() == bed->firstPeak() - bed->getPeakSize())
			{
				print("\tbedThread(): endPeak() == peaks.end()" + debug);
				(*wigIter)->getPeakDiv(bed->firstPeak()->start, (bed->endPeak() + 1)->end, wig);
			}
			else
			{
				print("\tbedThread(): endPeak() != peaks.end()" + debug);
				(*wigIter)->getPeakDiv(bed->firstPeak()->start, (bed->endPeak())->start, wig);
			}


//			(*wigIter)->getPeakDiv(startBed->start, endBed->start, wig);
			print("\tbedThread(): after wigIter->getPeakDiv()" + debug);
//			divThread(startBed, endBed, startWig, endWig, bedRegionMerge[i]);

			print("\tbedThread(): before call to divThread()" + debug);

			divThread(bed, wig, regionMerge[i + whichSlice]);
			print("\tbedThread(): in for end" + debug);
		}

//		regionMerge = mergeMetaplotRegions(bedRegionMerge, bedsPerThread);
		print("bedThread(): no div thread end" + debug);
	}
	else
	{
		print("bedThread(): else: div thread start" + debug);
//		MetaplotRegion ** bedRegionMerge = new MetaplotRegion*[bedsPerThread];

		for (int i = 0; i < bedsPerThread; i++)
		{
			stringstream iToS;
			iToS << i;

			print("bedThread(): start of bedsPerThread for i is " + iToS.str() + debug);
			regionMerge[i + whichSlice] = new MetaplotRegion(maxWindow);

			MetaplotRegion ** divRegions = new MetaplotRegion*[threads->divThreads];
			bedAr[i]->readPeaks();
			thread threadAr[threads->divThreads];
			
			for (int j = 0; j < threads->divThreads; j++)
			{
				print("bedThread(): start of threads->divThreads for" + debug);
				Bed * bed;
				Wig * wig;
				divRegions[j] = new MetaplotRegion(maxWindow);
				vector<Peak>::iterator startBed, endBed, startWig, endWig;
				
				// FIXME
				// ah, here is the problem
				// peakDiv divides # of peaks
				// but this bed file has only 1 peak
				// and we want 2 minimum to do peakDiv (--threadDivisions 2)
				// basically that last div, which may not be perfectly divisible by thread #, is a problem
				// for now, fixed the 0 case with return if bed->getPeakSize == 0
				// but any case < # peaks per div will cause array overflow in map()
	
				bedAr[i]->getPeakDiv(threads->divThreads, j, startBed, endBed, bed);
// FIXME debugging here
// I think the issue is that I didn't update getPeakDiv to handle bed peaks as a stack

//				stringstream bedstarttest, bedendtest, bedEndEnd;
//				bedstarttest << startBed->start;
//				bedendtest << endBed->start;
//				bedEndEnd << endBed->end;
//				print("\t\tbedThread(): startbed start is " + bedstarttest.str() + " endbed start is " + bedendtest.str() + " endBed end " + bedEndEnd.str()+ debug);
//			print("\t\tbedThread(): startbed start is " + bedstarttest.str() + debug);

				// FIXME EndBed is not always valid!! 
//				(*wigIter)->getPeakDiv(startBed->start, (endBed-1)->end, wig);
//				int start, end;
//				start = bed->firstPeak()->start;
//				end = bed->endPeak()->end;

				if (bed->getPeakSize() == 0)
					return;

				print("bedThread(): before wig->getPeakDiv" + debug);
				if (bed->endPeak() == bed->firstPeak() - bed->getPeakSize())
				{
					stringstream sizeToI;
					sizeToI << bed->getPeakSize();
					print("\t\tbedThread(): endPeak() == peaks.end()! Peak size is " + sizeToI.str() + " bed chr is " + bed->getChr() + debug);
					(*wigIter)->getPeakDiv(bed->firstPeak()->start, (bed->endPeak() + 1)->end, wig);
				}
				else
					(*wigIter)->getPeakDiv(bed->firstPeak()->start, (bed->endPeak())->start, wig);
			
				print("bedThread(): after wig->getPeakDiv" + debug);
				stringstream wigSize;
				wigSize << wig->getPeakSize();
				print("\tbedThread(): wig size is " + wigSize.str() + debug);


				print("\tbedThread(): Before thread()" + debug);
				threadAr[j] = thread(&Process::divThread, this, bed, wig, ref(divRegions[j])); 
				print("\tbedThread(): After thread()" + debug);
			}

			for (int k = 0; k < threads->divThreads; k++)
				threadAr[k].join();

			print("bedThread(): After join()" + debug);

//			bedRegionMerge[i] = mergeMetaplotRegions(divRegions, threads->divThreads);
			regionMerge[i + whichSlice] = mergeMetaplotRegions(divRegions, threads->divThreads);
		}

		print("bedThread(): div thread end" + debug);

//		regionMerge = mergeMetaplotRegions(bedRegionMerge, bedsPerThread);
	}
}


void Process::divThread(Bed * bed, Wig * wig, MetaplotRegion * &regionMerge)
{
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("divThread(): start" + debug);

	while(bed->isValid())
	{
		print("\tdivThread(): in while start" + debug);
		vector<Peak>::iterator iter = bed->getCurrPeak();

		Wig * smallWig;

		print("\tdivThread(): before wig getPeakDiv call" + debug);

		wig->getPeakDiv(iter->start, iter->end, smallWig);
		
		print("\tdivThread(): before call to map" + debug);
//		skip for now
//		mapWig(smallWig, iter->start, iter->strand, regionMerge);	

		bed->nextPeak();
		print("\tdivThread(): in while end" + debug);
	}

/*
	
	for (vector<Peak>::iterator iter = bedPeakStart; iter != bedPeakEnd; iter++)
	{
		print("\tdivThread(): in for start" + debug);
		// init to values that they could not possibly be in the wig ar
		vector<Peak>::iterator thisWigPeakStart = iter; 
		vector<Peak>::iterator thisWigPeakEnd = iter;

		vector<Peak>::iterator min = wigPeakStart;
		vector<Peak>::iterator max = wigPeakEnd;
		vector<Peak>::iterator mid = wigPeakEnd; // init only

		print("\tdivThread(): iter init finished" + debug);

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

		print("\tdivThread(): after first bsearch" + debug);

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
	
		print("\tdivThread(): before call to map" + debug);
		mapWig(thisWigPeakStart, thisWigPeakEnd, iter->start, iter->strand, regionMerge);	
		print("\tdivThread(): in for end" + debug);
	} // for	

*/
	print("divThread(): end" + debug);
}

/* FIXME THIS ENTIRE ALGORITHM RESTS ON THE FACT THAT WINDOW % SHIFT == 0
 * EVERYTHING WILL BREAK IF THIS IS NOT TRUE */

void Process::mapWig(Wig * wig, int bedStart, char bedStrand, MetaplotRegion * &region)
{
	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("mapWig(): start" + debug);

	while(wig->isValid())
	{
		vector<Peak>::iterator iter = wig->getCurrPeak();

		RunningAvg avg(window, shift);

		vector<Peak>::iterator backIter = iter; 	
		vector<Peak>::iterator fwdIter = iter; 

		while (backIter != wig->firstPeak() && backIter->end > (iter->start - window))	
			backIter++; // ++ is -- !!! because stack-like ops. 

		if (backIter->end < iter->start - window)
			backIter--;
	
		while(fwdIter != wig->endPeak() && fwdIter->start < (iter->end + window))
			fwdIter--;

		if (fwdIter->start > iter->end + window)
			fwdIter++;

		int pos = iter->start - window;

//		while((backIter == wig->firstPeak()) || backIter < iter)
		while(backIter < iter)
		{
			int backSignalLen = backIter->end - pos;
			
			avg.queue(backSignalLen / shift, backIter->value * shift);

			pos += backSignalLen / shift * shift;

			int cond = ( backSignalLen - backSignalLen / shift * shift + (shift - 1) ) / shift;
	
			avg.queue(cond, backIter->value * (backSignalLen % shift) );

			int shift_end = pos + shift - backSignalLen % shift + window;

			pos += backSignalLen % shift;

			// this will crash if backIter == iter : now we're past iter 
			backIter--;

			avg.next( (backIter->start - shift_end) / shift );

			cond = ( (backIter->start - shift_end) + (shift - 1) ) / shift;
			int rem = (shift - ( (backIter->start - shift_end) % shift) );

			avg.queue(cond, backIter->value * rem);

			pos += rem;
		}		

		region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);

// FIXME no idea what's going on after this point	

		int numSlidingWindow = 0;
		pos += shift;

		int nextOverlap = iter->end - window;	
		
		if (fwdIter != iter)
			nextOverlap = fwdIter->start - window;

		if (nextOverlap < iter->start)
		{
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
			numSlidingWindow = window / shift;

			for (int i = 0; i < numSlidingWindow; i++)
			{
				avg.queue(1, iter->value * shift);
				region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);
				pos += shift;
			}

			int nextEnd = iter->end;

			if (fwdIter != iter)
				nextEnd = fwdIter->start - window;

			int middle_len = (nextEnd - pos) / shift;
			
			region->addSignal(pos - bedStart, middle_len, avg.getAvg(), bedStrand);

			pos += middle_len * shift;
		}
	
		if (fwdIter == iter)
		{
			int shift_start = numSlidingWindow * shift + iter->start;

			int len = ( iter->end - shift_start + (shift - 1) ) / shift;
			int rem = ( iter->end - shift_start ) % shift;

			avg.queue(len, rem * iter->value);
		
			region->addSignal( pos - bedStart, len, rem * iter->value, bedStrand);

			while(pos < shift_start)
			{
				avg.queue(1, 0);
				region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);
				pos += shift;
			}

		}


		wig->nextPeak();
	}
/*	




	for (vector<Peak>::iterator iter = wigPeakStart; iter != wigPeakEnd; iter++)
	{
		print("\tmapWig(): in for start" + debug);	
		RunningAvg avg(window, shift);

		print("\tmapWig(): after init avg" + debug);		
		// set backIter and fwdIter
		vector<Peak>::iterator backIter = iter;
		vector<Peak>::iterator fwdIter = iter + 1;
		vector<Peak>::iterator fwdIterStop = iter;

		print("\tmapWig(): after init iters" + debug);		

		while (backIter != wigPeakStart && backIter->end > (iter->start - window))
			backIter--;
	

		if (backIter->end < iter->start - window)
			backIter++;

		print("\tmapWig(): after set backIter" + debug);

		while(fwdIterStop != wigPeakEnd && fwdIterStop->start < (iter->end + window))
			fwdIterStop++;

		if (fwdIterStop->start > iter->end + window)
			fwdIterStop--;

		print("\tmapWig(): after set fwdIterStop" + debug);

		stringstream iToS;
		iToS << window;
		print("\tmapWig(): window is " + iToS.str() + debug);

		int pos = iter->start - window;

		if (iter == wigPeakStart && wigPeakStart->start > bedStart)
			pos = iter->start;
		else if (iter == wigPeakStart && wigPeakStart->start < bedStart)
			pos = bedStart;

		print("\tmapWig(): after init pos" + debug);


		stringstream posTest, bedTest, wigTest, iterTest;
		posTest << pos;
		bedTest << bedStart;
		wigTest << wigPeakStart->start;
		iterTest << iter->start;

		print("\tmapWig(): bed start " + bedTest.str() + " pos " + posTest.str() + " wig " + wigTest.str() + " iter " + iterTest.str() + debug);
		// init queue
		// can't addSignal until queue is fully init
		// FIXME if this is the first peak we will not go into this loop, thus runningAvg will not be init (well, it will be filled with 0's so okay?? but POS will be very very wrong
//		while (backIter != iter)
		while ((backIter == wigPeakStart) || backIter < iter)
		{
			print("\t\tmapWig(): while start" + debug);
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

			print("\t\tmapWig(): while end" + debug);
		} // while

		print("\tmapWig(): after while" + debug);
		// slide forward until end of window == peakStart

		// start addSignal here
		region->addSignal(pos - bedStart, 1, avg.getAvg(), bedStrand);

		print("\tmapWig(): after addSignal" + debug);

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
*/
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
/*

		print("\tmapWig(): in for end" + debug);	
	} // for
	print("mapWig(): end" + debug);
*/
}

string Process::printResults(string nameStr, string nameStrR)
{
	print("Printing out results");
	
	string outfileName = "metaplot_outfile.txt";
	ofstream file(outfileName.c_str());
	
	string header = "bp\t" + nameStr + "\n";
	file << header;

	int h = 0 - (maxWindow / 2);
	
//	print("before for");
	for (int i = 0; i < maxWindow; i++)
	{
		h++;
		stringstream hToString;
		hToString << h;
		file << hToString.str() << "\t";

		for (int j = 0; j < bedNum; j++)
		{
//			print("start of inner for");
			if (region[j]->basePairs[i][1] <= 0)
				file << "NA\t";
			else
			{
				double avg = (double)region[j]->basePairs[i][0] / (double)region[j]->basePairs[i][1];
				stringstream avgToString;	
				avgToString << avg;
				file << avgToString.str() << "\t";
			}
//			print("end of inner for");
		}
		file << endl;
	} // for
//	print("after for");
	
	ofstream rstream("metaplot_outfile.R");
	rstream << "library(ggplot2)\nlibrary(reshape2)\n";
	rstream << "pdf(file=\"metaplot_outfile.pdf\", width=12, height=8)\n";
	rstream << "plot<-read.table(\"metaplot_outfile.txt\", header=T)\n";
	rstream << "plot.melt<-melt(plot[,c('bp', " << nameStrR << ")], id.vars=1)\n";
	rstream << "ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) +\n";
	rstream << "geom_line() +\n";
	rstream << "geom_smooth() +\n";
	rstream << "theme_bw() +\n";
	rstream << "ggtitle(\"Metaplot\") +\n";
	rstream << "theme(panel.grid.minor=element_blank()) +\n";
	rstream << "scale_colour_brewer(palette=\"Set1\", name=\"Bed\")\n";

	rstream.close(); 

	file.close();

	return outfileName;
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
