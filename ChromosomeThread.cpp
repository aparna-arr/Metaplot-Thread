#include"ChromosomeThread.h"
using namespace std;

void Chromosome::clearPeaks(void)
{
	peaks.clear();
	currPeakValid = false;
}

void Chromosome::setFilename(string filename)
{
	fileName = filename;
}

string Chromosome::getChr(void)
{
	return chrName;
}

void Chromosome::addPeak(void) 
{
	currPeakValid = true;

	currPeak = peaks.begin() + peaks.size() - 1;
}

void Chromosome::printPeaks(void)
{
	currPeakValid = false;
}

void Chromosome::readPeaks(void)
{

}

vector<Peak>::iterator Chromosome::getCurrPeak(void)
{
	if (peaks.empty() || currPeak == peaks.end())
	{
		currPeakValid = false;
		return peaks.end(); // DO NOT USE THIS VALUE
	}

	return currPeak;
}

void Chromosome::nextPeak(void)
{ // NOTE operating vector like a stack!
	if (!peaks.empty() && currPeak != peaks.begin())
		currPeak--;
	else
		currPeakValid = false; 
}

void Chromosome::generateFilename(void) { };

void Chromosome::reset(void)
{	
	if (!peaks.empty())
		currPeak = peaks.begin() + peaks.size() - 1;
	else
		currPeakValid = false;
}

bool Chromosome::isValid(void)
{
	return currPeakValid;
}

int Chromosome::getPeakSize(void)
{
	return peaks.size();
}

/* End Chromosome functions */

/* Wig functions */
/*
void Wig::unstack(void)
{
	unstacked.push(peaks.top());
	peaks.pop();
}

void Wig::restack(void)
{
	while(!unstacked.empty())
	{
		peaks.push(unstacked.top());
		unstacked.pop();
	}
}
*/

vector<Peak>::iterator Wig::firstPeak(void)
{
	if (!peaks.empty())
		return peaks.end() - 1;
	else
	{
		currPeakValid = false;
		return peaks.end();
	}
}

vector<Peak>::iterator Wig::endPeak(void)
{
	if (!peaks.empty())
		return peaks.begin() - 1;
	else
	{
		currPeakValid = false;
		return peaks.end();
	}
}

void Wig::printPeaks(void)
{
	ofstream fstream(fileName.c_str());

	while(!peaks.empty())
	{
		Peak tmp = *currPeak;
		fstream << chrName << "\t" << tmp.start << "\t" << tmp.end << "\t" << tmp.value << endl;
//		cerr << "printing " << chrName << "\t" << tmp.start << endl;

		nextPeak();
		peaks.pop_back(); // remove last element
	}

	fstream.close();	
	
	Chromosome::printPeaks(); // handles currPeak iter
}

void Wig::readPeaks(void)
{
	ifstream fstream(fileName.c_str());
	string line;

	while(getline(fstream, line))
	{
		stringstream linestream(line);
		int start, end;
		double value;
		string chrom;
		
		linestream >> chrom >> start >> end >> value;

		addPeak(start, end, value);
	}	
	
	fstream.close();
}

void Wig::addPeak(int start, int end, double value)
{
	Peak tmp;

	tmp.start = start;
	tmp.end = end;
	tmp.value = value;

	peaks.push_back(tmp);

//	cerr << "adding wig peak " << chrName << "\t" << start << endl;

	Chromosome::addPeak(); // handles currPeak
}

void Wig::generateFilename(void)
{
	fileName = chrName + "_wig" + ".bed";
}

void Wig::getPeakDiv(int startPos, int endPos, Wig * &div)
{
//	cerr << "in wig Get Peak Div" << endl;
	// do bsearch on start and end
	// return segment of wigpeaks between [ ) start and end

//	cerr << "in wig get peak div" << endl;
	div = new Wig(chrName);

	vector<Peak>::iterator max = peaks.begin();
	vector<Peak>::iterator min = peaks.end();
	vector<Peak>::iterator mid = (min - max) / 2 + max;

//	cerr << "mid - max is " << mid - max << endl;

//	if (min == max)
//		cerr << "begin() == end()" << endl;
//	else
//	{
//		cerr << "min - max is " << min - max << endl;
//		for (vector<Peak>::iterator iter = peaks.begin(); iter != peaks.end(); iter++)
//		{
//			cerr << "\tpeak is " << iter->start << endl;
//		}
//	}

	vector<Peak>::iterator startIter, endIter;


//	if (mid == min)
//		cerr << "mid == min!" <<endl;
//	else if (mid == max)
//		cerr << "mid == max!" << endl;
//	else
//		cerr << "mid != min | max" << endl;

//	cerr << "after init startpos is " <<startPos << " endpos is " << endPos << endl;
//	cerr << "mid start is " << mid->start << " end is " << mid->end << endl;
	
	while (max <= min)
	{
//		cerr << "in while" << endl;
		mid = (min - max) / 2 + max;
//	cerr << "mid start is " << mid->start << " end is " << mid->end << endl;
		
		if (mid->start == startPos)
			break;
		else if (mid->start < startPos)
			min = mid - 1;
		else if (mid->start > startPos)
			max = mid + 1;
	}
	
//	cerr << "After while" << endl;
//	cerr << "mid start is " << mid->start << " end is " << mid->end << endl;

//	cerr << "Mid start is " << mid->start << endl;
//	cerr << "before first if startPos is " << startPos << " mid->start is " << mid->start << " mid->end is " << mid->end << endl;
	if (mid != peaks.begin() && mid->end < startPos)
		mid--;


// FIXME modifying this caused the dreaded terminate called without an active exception
// error -> fixed and this had nothing to do with it!	
	if (mid < peaks.end() - 1 && mid->start > startPos && (mid+1)->end > startPos)
		mid++;

	startIter = mid;

//	cerr << "after 1st bsearch startIter is " << startIter->start << endl;

	max = peaks.begin();
	min = peaks.end();
	
	while(max <= min)
	{
		mid = (min - max) / 2 + max;
//	cerr << "mid start is " << mid->start << " end is " << mid->end << endl;

		if(mid->start == endPos)
			break;
		else if (mid->start < endPos)
			min = mid - 1;
		else if (mid->start > endPos)
			max = mid + 1;
	}	

//	cerr << "after 2nd bsearch before ifs mid is " << mid->start << endl;
	if (mid < peaks.end() - 1 && (mid-1)->start > endPos)
		mid++;
	if (mid != peaks.begin() && mid->end < endPos)
		mid--;

	endIter = mid;

//	cerr << "after 2nd bsearch endIter is " << endIter->start << endl;

//	if (startIter < endIter)
//		cerr << "startIter < endIter" << endl;

	// FIXME should this be iter <= endIter??? for the case where start and end are the same? Which should never happen unless NO wig peaks over bed ...

//	cerr << "startIter - endIter is " << startIter - endIter << endl;
	for (vector<Peak>::iterator iter = endIter; iter <= startIter; iter++)
		div->addPeak(iter->start, iter->end, iter->value);
	
//	cerr << "end of getPeakDiv" << endl;
}

/* End Wig functions */

/* Bed functions */

void Bed::printPeaks(void)
{
	ofstream fstream(fileName.c_str());
	while(!peaks.empty())
	{
		Peak tmp = *currPeak;
		fstream << chrName << "\t" << tmp.start << "\t" << tmp.end << "\t" << tmp.strand << endl;
		nextPeak();
		peaks.pop_back(); // remove last element
	}
	fstream.close();	

	Chromosome::printPeaks();	
}

void Bed::readPeaks(void)
{
	ifstream fstream(fileName.c_str());
	string line;

	while(getline(fstream, line))
	{
		stringstream linestream(line);
		int start, end;
		char strand;
		string chrom;
		
		linestream >> chrom >> start >> end >> strand;

		addPeak(start, end, strand);
	}	
	
	fstream.close();
}

void Bed::addPeak(int start, int end, char strand)
{
	Peak tmp;

	tmp.start = start;
	tmp.end = end;
	tmp.strand = strand;

	peaks.push_back(tmp);

	Chromosome::addPeak();
}

void Bed::generateFilename(int filenum)
{
	stringstream intToString;
	intToString << filenum;
	fileName = chrName + "_" + intToString.str() + ".bed";
}

void Bed::getPeakDiv (int numDivs, int iteration, vector<Peak>::iterator &startIter, vector<Peak>::iterator &endIter, Bed * &div)
{
	// divide size() by numDivs, return [ iteration * numDivs, iteration * numDivs + (iteration+1) * numDivs OR peaks.end() )

	div = new Bed(chrName, peakLen);

//	int sizeOfDiv = (peaks.size() + 1) / (numDivs + 1);

	int sizeOfDiv;
	if (numDivs == 0)
		sizeOfDiv = peaks.size();
	else
		sizeOfDiv = peaks.size() / numDivs;

//	cerr << "Bed():getPeakDiv(): sizeOfDiv is " << sizeOfDiv << endl;
//	cerr << "peaks size is " << peaks.size() << endl;

	startIter = peaks.begin() + sizeOfDiv * iteration;

	if (startIter + sizeOfDiv > peaks.end())
		endIter = peaks.end();
	else
		endIter = peaks.begin() + sizeOfDiv * (iteration + 1);

//	if (endIter == peaks.end())
//		cerr << "endIter == peaks.end()" << endl;

	// FIXME same question as Wig::getPeakDiv
	for (vector<Peak>::iterator iter = startIter; iter != endIter; iter++)
		div->addPeak(iter->start, iter->end, iter->strand);
}
	
vector<Peak>::iterator Bed::firstPeak(void)
{
	if (!peaks.empty())
		return peaks.end() - 1;
	else
	{
		currPeakValid = false;
		return peaks.end();
	}
}

vector<Peak>::iterator Bed::endPeak(void)
{
	if (!peaks.empty())
		return peaks.begin() - 1;
	else
	{
		currPeakValid = false;
		return peaks.end();
	}
}
/* End Bed functions */
