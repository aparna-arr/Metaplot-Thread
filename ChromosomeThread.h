#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include<iostream>
#include<string>
#include<stack>
#include<vector>
#include<iterator>
#include<exception>
#include<sstream>
#include<fstream>

typedef struct Peak
{
	int start;
	int end;
	double value;
	char strand;
	std::string chr;
} Peak;

class Chromosome
{
	public:
	Chromosome() {};
	Chromosome(std::string chr) 
	{
		chrName = chr;
		currPeakValid = false;
	}	

	virtual void addPeak(void);
	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	
	std::string getChr(void);
	void nextPeak(void);
	std::vector<Peak>::iterator getCurrPeak(void);	
	int getPeakStart(void);
	int getPeakEnd(void);
	void setFilename(std::string filename);
	void clearPeaks(void);	
	void reset(void);	
	bool isValid(void);
	int getPeakSize(void);
	virtual void getPeakDiv(void) { } ;
	
	protected:
	virtual void generateFilename(void);

	std::string chrName;
	std::string fileName;
	std::vector<Peak> peaks;
	std::vector<Peak>::iterator currPeak;
	
	bool currPeakValid;
};

class Wig: public Chromosome
{
	public:
	Wig() {};
	Wig(std::string chr) : Chromosome(chr) { };
	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	virtual void addPeak(int start, int end, double value);
	virtual void generateFilename(void);
	virtual void getPeakDiv(int startPos, int endPos, std::vector<Peak>::iterator &startIter, std::vector<Peak>::iterator &endIter); // to do bsearch on
//	void unstack(void);
//	void restack(void);

	private:
	
//	std::stack<Peak> unstacked;
};

class Bed : public Chromosome
{
	public:
	Bed() {};
	Bed(std::string chr, int len) : Chromosome(chr)
	{
		peakLen = len;
	}

	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	virtual void addPeak(int start, int end, char strand);
	virtual void generateFilename(int filenum);
	virtual void getPeakDiv(int numDivs, int iteration, std::vector<Peak>::iterator &startIter, std::vector<Peak>::iterator &endIter);	

	private: 
	int peakLen;
};


#endif
