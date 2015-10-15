#include"MetaplotRegionThread.h"
using namespace std;

MetaplotRegion::MetaplotRegion(int len)
{
	length = len;

	basePairs = new double*[len];
	
	for (int i = 0; i < len; i++)
	{
		basePairs[i] = new double[2]; // 0 = height, 1 = count
		basePairs[i][0] = 0;
		basePairs[i][1] = 0;
	}
}

void MetaplotRegion::addSignal(int offset, int len, double value, char strand)
{
	cerr << "In addSignal" << endl;
	cerr << "offset is " << offset << " len is " << len << endl;
	
	if (offset < 0 || len <= 0)
		return;

	cerr << "going to add " << value << endl;

	if (strand == '+')
	{
		for (int i = offset; i < len + offset; i++)
		{
			basePairs[i][0] += value;
			basePairs[i][1]++;
		}		
	}
	else 
	{
		for (int i = length - offset - 1; i > length - (offset + len) - 1; i--)
		{
			basePairs[i][0] += value;
			basePairs[i][1]++;
		}
	}
//	cerr << "End of addSignal" << endl;
}
