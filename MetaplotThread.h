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
	void queue(double val);
	void dequeue(void);
	double getAvg(void);	

	private:
	double * queue; // this should be circular
	int first; // thing to dequeue next
	int nextEmpty; // next empty block to use -- same as first???
	int numShifts;	

	double currentAvg; 		
};

void print(std::string msg);

MetaplotRegion * calculate(int bedNum, int maxWindow, std::vector<Bed *> * bedsByChr, std::vector<Wig *> wigsByChr, ThreadInfo * threads, std::map< std::string, std::vector<int>* > * chrs);

void chromThread(int bedNum, int maxWindow, std::vector<Bed *> * bedsByChr, std::vector<Wig *> wigsByChr, ThreadInfo * threads, std::map< std::string, std::vector<int> *> *chrs, std::vector<std::string> chromVec, int startChrom, int chromsPerThread);

void bedThread(int maxWindow, Bed ** bedAr, Wig * wig, int divThreads, int startBed, int bedsPerThread);

void divThread(int maxWindow, std::vector<Peak>::iterator bedStart, std::vector<Peak>::iterator bedEnd, std::vector<Peak>::iterator wigStart, std::vector<Peak>::iterator wigEnd);

std::string printResults(MetaplotRegion * result, int bedNum, std::string nameStr, std::string nameStrR, int maxWindow);

void monteCarlo(std::string outfile, int bedNum);

#endif
