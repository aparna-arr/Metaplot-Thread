#ifndef METAPLOTREGION_H
#define METAPLOTREGION_H

#include<iostream>
#include<string>
#include<exception>

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
