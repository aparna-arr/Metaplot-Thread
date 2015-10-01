#include"MetaplotRegion.h"
using namespace std;
/* MetaplotRegion functions */

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
	if (strand == '+')
	{
		for (int i = offset; i < len + offset; i++)
		{
		//	cerr << "DEBUG adding " << value << " to " << i << endl;
			basePairs[i][0] += value;
			basePairs[i][1]++;
		//	cerr << "basePairs " << i << " height is now " << basePairs[i][1] << " signal is " << basePairs[i][0] <<  endl;
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
}

/* End MetaplotRegion functions */
