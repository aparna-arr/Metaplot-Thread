#ifndef USEROPTS_H
#define USEROPTS_H

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

#include"Chromosome.h"

typedef struct Freq
{
	std::string chr;
	int freq;
} Freq;

class UserOpts
{
	public:
	UserOpts(int argc, char *argv[]);

	void printUsage(void);
	int getBedNumber(void);
	int getMaxWindow(void);
	std::string getNextBedFile(void);
	void resetIter(void);
	std::vector<Bed *> * splitBedFiles(void);
	std::vector<Wig *> splitWigFile(void);
	void debugVectorSize(void);

	std::string getNameString(void);
	std::string getNameStringR(void);
	std::vector<std::string> commonChroms(void);

	bool monteCarlo(void);
	
	private:
	bool validateBeds(void);
	bool validateWig(void);
	bool validateNames(void);

	std::string handleOpts(int argc, char * argv[]);
	void preprocessBed(void);
	void preprocessWig(void);

	void calcSlidingWindow(std::vector<Peak> * wigBlocks);

	int numBeds;
	int maxWindow;

	int step;	
	int window;
	int mode;

	bool preprocessBedOpt;
	bool preprocessWigOpt;

	std::string splitWigDir;
	std::string splitBedDir;

	bool wigsSplit;	
	bool bedsSplit;

	bool isMonteCarlo;
	
	bool bedsFromDir;
	std::string bedDir;

	std::vector<Bed *> * readBedSplit(std::vector<Bed *> * &arOfBedfiles);
	std::vector<Wig *> readWigSplit(std::vector<Wig *> &wigSplit);
	void readBedsFromDir(std::vector<std::string> & tmpNamesAndBeds, std::string dir);

	std::stack<std::string> getFilesInDir(std::string dir);
	
	std::vector<Freq> chroms;	

	std::string wigFile;
	std::vector<std::string> bedFiles;
	std::stack<std::string> names;
	std::stack<std::string> namesR;
	// add more like debug mode later
};
#endif
