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
#include<cstdio>
#include<dirent.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

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
	}	

	virtual void addPeak(void);
	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	
	std::string getChr(void);
	void nextPeak(void);
	Peak * getCurrPeak(void);	
	int getPeakStart(void);
	int getPeakEnd(void);
	void setFilename(std::string filename);
	void clearPeaks(void);	

	protected:
	virtual void generateFilename(void);

	std::string chrName;
	std::string fileName;
	std::stack<Peak> peaks;
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

	void unstack(void);
	void restack(void);

	private:
	
	std::stack<Peak> unstacked;
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
	
	private: 
	int peakLen;
};
#endif
