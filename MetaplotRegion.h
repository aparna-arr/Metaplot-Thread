#ifndef METAPLOTREGION_H
#define METAPLOTREGION_H

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

class MetaplotRegion
{
	public:
	MetaplotRegion() { };
	MetaplotRegion(int len);

	void resetIndicies(void);
	void addSignal(int offset, int len, double value, char strand);

	double ** basePairs;
	private:
	int length;
};

#endif
