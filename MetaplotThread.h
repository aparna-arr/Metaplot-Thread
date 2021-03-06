#ifndef METAPLOT_THREAD_H
#define METAPLOT_THREAD_H

#include<iostream>
#include<mutex>
#include<thread>
#include<vector>
#include<iterator>
#include<string>
#include<map>

#include"MetaplotRegionThread.h"
#include"ChromosomeThread.h"

static std::mutex printMutex;

typedef struct ThreadInfo
{
	int bedThreads;
	int chromThreads;
	int divThreads;
} ThreadInfo;

class RunningAvg
{
	public: 
	RunningAvg();
	RunningAvg(int window, int shift);
	void populate(double value);
	void queue(int numQueues, double val);
	void next(int num);
	double getAvg(void);	

	private:
	double * queueAr; // this should be circular
	int nextEmpty; // next empty block to use -- same as first???
	int numShifts;	
	int windowSize;
	double currentAvg; 		
};

class Process
{
	public:
	Process() { } ;
	Process(int bedNumArg, int maxWindowArg, int shiftArg, int windowArg, std::vector<Bed *> * &bedsByChrArg, std::vector<Wig *> &wigsByChrArg, ThreadInfo * threadsArg, std::map< std::string, std::vector<int> * > * &commonChrs);

	void calculate(void);
	std::string printResults(std::string nameStr, std::string nameStrR);

	private:
	void chromThread(std::vector<std::string>::iterator chromStart, std::vector<std::string>::iterator chromEnd, MetaplotRegion ** &regionMerge);
	void bedThread(Bed ** bedAr, int arraySize, std::vector<Wig *>::iterator wigIter, MetaplotRegion ** &regionMerge, int whichSlice);
	void divThread(Bed * bed, Wig * wig, MetaplotRegion * &regionMerge);
	void mapWig(Wig * wig, int bedStart, char bedStrand, MetaplotRegion * &regionMerge);
	void slidingWindow(Wig * wig, int bedStart, char bedStrand, MetaplotRegion * &regionMerge);
	
	void makeChromVec(std::map< std::string, std::vector<int> * > * commonChrs);

	MetaplotRegion * mergeMetaplotRegions(MetaplotRegion ** array, int numRegions);
	void mergeMetaplotRegionsByBed(MetaplotRegion ** &oldArray, MetaplotRegion ** newArray);

	int bedNum;	
	int maxWindow;
	int shift;
	int window;
	int chromsPerThread;
	int bedsPerThread;
	int divsPerThread;
	ThreadInfo * threads;
	std::vector<Bed *> * bedsByChr;
	std::vector<Wig *> wigsByChr; 
	std::map< std::string, std::vector<int> * > * chrs;
	std::vector<std::string> chromVec;
	MetaplotRegion ** region; // this needs to be an array ... one per bed file.
};

void print(std::string msg);

std::string convertI(int num);
std::string convertD(double num);

/*
MetaplotRegion * calculate(int bedNum, int maxWindow, std::vector<Bed *> * bedsByChr, std::vector<Wig *> wigsByChr, ThreadInfo * threads, std::map< std::string, std::vector<int>* > * chrs);

void chromThread(int bedNum, int maxWindow, std::vector<Bed *> * bedsByChr, std::vector<Wig *> wigsByChr, ThreadInfo * threads, std::map< std::string, std::vector<int> *> *chrs, std::vector<std::string> chromVec, int startChrom, int chromsPerThread);

void bedThread(int maxWindow, Bed ** bedAr, Wig * wig, int divThreads, int startBed, int bedsPerThread);

void divThread(int maxWindow, std::vector<Peak>::iterator bedStart, std::vector<Peak>::iterator bedEnd, std::vector<Peak>::iterator wigStart, std::vector<Peak>::iterator wigEnd);

std::string printResults(MetaplotRegion * result, int bedNum, std::string nameStr, std::string nameStrR, int maxWindow);
*/
void monteCarlo(std::string outfile, int bedNum);

#endif
