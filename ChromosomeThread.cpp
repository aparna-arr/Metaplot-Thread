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

	cerr << "in wig get peak div" << endl;
	div = new Wig(chrName);

	vector<Peak>::iterator min = peaks.begin();
	vector<Peak>::iterator max = peaks.end();
	vector<Peak>::iterator mid = (max - min) / 2 + min;

	vector<Peak>::iterator startIter, endIter;

	cerr << "after init" << endl;
	
	while (min < max)
	{
		mid = (max - min) / 2 + min;
		
		if (mid->start == startPos)
			break;
		else if (mid->start < startPos)
			max = mid - 1;
		else if (mid->start > startPos)
			min = mid + 1;
	}

	if (mid->end < startPos)
		mid++;
	
	if (mid != peaks.begin() && (mid - 1)->end > startPos)
		mid--;

	startIter = mid;

	cerr << "after 1st bsearch" << endl;

	min = mid;
	max = peaks.end();
	
	while(min < max)
	{
		mid = (max - min) / 2 + min;

		if(mid->end == endPos)
			break;
		else if (mid->end < endPos)
			max = mid - 1;
		else if (mid->end > endPos)
			min = mid + 1;
	}	

	if (mid != startIter && (mid-1)->start > endPos)
		mid--;
	if (mid != peaks.end() && mid->end < endPos)
		mid++;

	endIter = mid;

	cerr << "after 2nd bsearch" << endl;

	// FIXME should this be iter <= endIter??? for the case where start and end are the same? Which should never happen unless NO wig peaks over bed ...
	for (vector<Peak>::iterator iter = startIter; iter != endIter; iter++)
		div->addPeak(iter->start, iter->end, iter->value);
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

	int sizeOfDiv = (peaks.size() + 1) / (numDivs + 1);

	startIter = peaks.begin() + sizeOfDiv * iteration;

	if (startIter + sizeOfDiv > peaks.end())
		endIter = peaks.end();
	else
		endIter = peaks.begin() + sizeOfDiv * (iteration + 1);

	// FIXME same question as Wig::getPeakDiv
	for (vector<Peak>::iterator iter = startIter; iter != endIter; iter++)
		div->addPeak(iter->start, iter->end, iter->strand);
}
	
/* End Bed functions */
