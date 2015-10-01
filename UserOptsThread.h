#ifndef USEROPTS_H
#define USEROPTS_H

#include<iostream>
#include<string>
#include<exception>
#include<sstream>
#include<fstream>
#include<vector>
#include<iterator>
#include<stack>
#include<map>
#include<cstdio>
#include<dirent.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

#include"ChromosomeThread.h"
#include"MetaplotThread.h"

class UserOpts
{
	public:
	UserOpts(int argc, char *argv[]);

	void printUsage(void);
	void preprocess(void);
	void split(std::vector<Bed *> * &bedSplit, std::vector<Wig *> &wigSplit);
	
	int getBedNum(void);
	int getMaxWindow(void);
	std::string getNameStr(void);
	std::string getNameStrR(void);
	bool isMonteCarlo(void);
	ThreadInfo * getThreadInfo(void);
	std::map< std::string, std::vector<int> * > * getCommonChrs(void);

	private:
		
	std::stack<std::string> getFilesInDir(std::string dir);

	std::string handleOpts(int argc, char * argv[]);
	void preprocessBed(void);
	void preprocessBedThread(int bedsPerThread, int startBed);
	void preprocessWig(void);
	void calcSlidingWindow(std::vector<Peak> * wigBlocks);

	void splitBeds(std::vector<Bed *> * &arOfBedfiles, int startBed, int bedsPerThread);
	void splitWig(std::vector<Wig *> &wigSplit);
	void readBedSplit(std::vector<Bed *> * &arOfBedfiles);
	void readWigSplit(std::vector<Wig *> &wigSplit);
//	void splitSomeBeds(vector<Bed *> * &arOfBedfiles, int startBed, int bedsPerThread);

	void readBedsFromDir(std::vector<std::string> &bedsAndNames);

	void commonChrsWig(std::map<std::string,int> &allChrs);	
	void commonChrsBed(std::map<std::string,int> &allChrs, std::string filename);
	void analyzeAndPrintChrs(std::map<std::string,int> allChrs);

	std::string nameStr;
	std::string nameStrR;

	int maxWindow;
	int bedNum;
	int chromNum;	

	bool preprocessWigOpt;
	int step;
	int window;

	bool preprocessBedOpt;
	int mode;

	std::string wigSplitDir;
	std::string bedSplitDir;

	bool isMonteCarloOpt;

	std::string bedDir;

	bool onlyFindChrs;	
	std::string inputChrs;

	ThreadInfo threads;	
	
	std::string wigFile;
	std::vector<std::string> bedFiles;
	std::stack<std::string> names;

	std::map<std::string, std::vector<int> *> commonChrs;
};

#endif
