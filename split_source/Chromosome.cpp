#include"Chromosome.h"
using namespace std;

/* Chromosome functions */


void Chromosome::clearPeaks(void)
{
	while(!peaks.empty())
	{
		peaks.pop();
	}
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

}

void Chromosome::printPeaks(void)
{

}

void Chromosome::readPeaks(void)
{

}

Peak * Chromosome::getCurrPeak(void)
{
	if (peaks.empty())
		return NULL;

	return &(peaks.top());
}

void Chromosome::nextPeak(void)
{
	if (!peaks.empty())
		peaks.pop();
}

void Chromosome::generateFilename(void)
{

}

/* End Chromosome functions */

/* Wig functions */

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

void Wig::printPeaks(void)
{
	ofstream fstream(fileName.c_str());
	while(!peaks.empty())
	{
		Peak tmp = peaks.top();
		fstream << chrName << "\t" << tmp.start << "\t" << tmp.end << "\t" << tmp.value << endl;
		peaks.pop();
	}
	fstream.close();		
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

	peaks.push(tmp);

}

void Wig::generateFilename(void)
{
	fileName = chrName + "_wig" + ".bed";
}

/* End Wig functions */

/* Bed functions */

void Bed::printPeaks(void)
{
	ofstream fstream(fileName.c_str());
	while(!peaks.empty())
	{
		Peak tmp = peaks.top();
		fstream << chrName << "\t" << tmp.start << "\t" << tmp.end << "\t" << tmp.strand << endl;
		peaks.pop();
	}
	fstream.close();		
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

	peaks.push(tmp);
}

void Bed::generateFilename(int filenum)
{
	stringstream intToString;
	intToString << filenum;
	fileName = chrName + "_" + intToString.str() + ".bed";
}

/* End Bed functions */
