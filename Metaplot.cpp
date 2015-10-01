#include"Metaplot.h"
using namespace std;

void calculateMetaplot(int bedNum, vector<Bed *> * &bedByChrs, vector<Wig *> &wigByChrs, vector<string> commonChrs, MetaplotRegion * &region)
{
//	cerr << "DEBUG in calculateMetaplot" << endl;
	for (vector<string>::iterator iter = commonChrs.begin(); iter != commonChrs.end(); iter++)
	{
		cout << endl << "On chromosome " << *iter << endl;
		// read in wig
		// then read in each bed

		cout << "Reading in wig file ... " << endl;
		vector<Wig *>::iterator wigIter;

		for (wigIter = wigByChrs.begin(); wigIter != wigByChrs.end(); wigIter++)
		{
			if ((*wigIter)->getChr() == (*iter))
				break;
		}
		(*wigIter)->readPeaks();

		cout << "\tdone." << endl;
		for (int bedFileNum = 0; bedFileNum < bedNum; bedFileNum++)
		{
			cout << "\tOn bedfile " << bedFileNum << endl;

			cout << "\t\tReading in bedfile ... " << endl;
			vector<Bed *>::iterator bedIter;
			for (bedIter = bedByChrs[bedFileNum].begin(); bedIter !=bedByChrs[bedFileNum].end(); bedIter++)
			{
				if ((*bedIter)->getChr() == (*iter))
					break;
			}
			
			(*bedIter)->readPeaks();

			cout << "\t\t\tdone." << endl;

			Peak * currBedPeak = (*bedIter)->getCurrPeak();
			Peak * currWigPeak = (*wigIter)->getCurrPeak();
			
//			cerr << "DEBUG: segfault?" << endl;
			/* now begin actual processing */
			// assumes all peaks are in ASCENDING order (check the stacks!)
			// okay first hurdle: we need to go through the wig multiple times, but each bed only once. But they're all stacks.
			// unstack/restack wig
			 
			// call restack wig at the end of every bedfile processed
			// do NOT nextPeak() on wig file

			while(currBedPeak != NULL && currWigPeak != NULL)
			{
				// this loop runs through all bed peaks for this file/chr

				// find first wig peak that is not < bed peak
				while (currWigPeak != NULL && currBedPeak->start > currWigPeak->end)
				{
					(*wigIter)->unstack();
					currWigPeak = (*wigIter)->getCurrPeak();
				}
				// assumes wig peaks are smaller than bed peaks
				while(currWigPeak != NULL && currWigPeak->start < currBedPeak->end)
				{
					// match and add up
					int offset = 0;
					int len = 0;
					if (currWigPeak->start > currBedPeak->start)
					{
						offset = currWigPeak->start - currBedPeak->start;

						if (currWigPeak->end >= currBedPeak->end)
							len = currBedPeak->end - currWigPeak->start;
						else
							len = currWigPeak->end - currWigPeak->start;	
					}					
					else
					{
						offset = 0;
						
						if (currWigPeak->end > currBedPeak->end)
							len = currBedPeak->end - currBedPeak->start;
						else
							len = currWigPeak->end - currBedPeak->start;
					}
					region[bedFileNum].addSignal(offset, len, currWigPeak->value, currBedPeak->strand);
//					cerr << "DEBUG region[" << bedFileNum << "] has signal " << region[bedFileNum].basePairs[0][0] << " and count " << region[bedFileNum].basePairs[0][1] << endl;

					(*wigIter)->unstack();
					currWigPeak=(*wigIter)->getCurrPeak();
				}

				(*bedIter)->nextPeak();
				currBedPeak = (*bedIter)->getCurrPeak();
			}
			// now at end for this bedfile
			(*wigIter)->restack();
			(*bedIter)->clearPeaks();
		}
		(*wigIter)->clearPeaks();		
	}
//	cerr << "DEBUG end of calculateMetaplot"<<endl;
//	cerr << "DEBUG region[0][1] is " << region[0].basePairs[0][1] << endl;
}


void monteCarloMetaplot(string file, int bedNum)
{
 // avg horizontally
 // no need for int reps 	
//	cerr << "DEBUG: in MonteCarloMetaplot" << endl;
	ifstream infile(file.c_str());
	ofstream outfile("metaplot_outfile.tmp");	

	outfile << "bp\tsimulation" << endl;

	string line;

//	cerr << "DEBUG: before while" << endl;
	while (getline(infile, line))
	{
		// have to deal with NA's
		double num;
		stringstream linestream(line);

//		cerr << "linestream is " << linestream.str() << endl;		

		if (!(linestream >> num))
			continue; // like next?
		
		// num == bp# right now
		double bp = num;
		double avg = 0;
		string NA;

// want to treat NAs like 0s
		for (int i = 0; i < bedNum; i++)
		{
			linestream >> NA;
			stringstream toNum(NA);

			if ((toNum >> num))
				avg += num;
		} 

		avg /= (double)bedNum;
		
		outfile << bp << "\t" << avg << endl;
	}
//	cerr << "DEBUG: after while" << endl;
	
	outfile.close();
	infile.close();
	rename("metaplot_outfile.txt", "tmp");
	rename("metaplot_outfile.tmp", "metaplot_outfile.txt");
	rename("tmp", "metaplot_outfile.tmp");
//	cerr << "DEBUG: end of function" << endl;
}


void debug(MetaplotRegion * &region, int bedNumber, string nameStr, string nameStrR)
{
	cerr << "In DEBUG function" << endl;
}

// returns outfile name
string printResults(MetaplotRegion * &region, int bedNumber, string nameStr, string nameStrR, int maxWindow)
{
//	debug();
	cout << endl << "Printing out results" << endl;

	string outfileName = "metaplot_outfile.txt";

	ofstream fstream(outfileName.c_str()); 

	string header = "bp\t" + nameStr + "\n";
	fstream << header;

	int h = 0 - (maxWindow / 2);

	for (int i = 0; i < maxWindow; i++)
	{
		int addition = h + i;
		stringstream hToString;
		hToString << addition;
		fstream << hToString.str() << "\t";

		for (int j = 0; j < bedNumber; j++)
		{
			if (region[j].basePairs[i][1] <= 0)
				fstream << "NA\t";
			else
			{
				double avg = (double)region[j].basePairs[i][0] / (double)region[j].basePairs[i][1];
				stringstream avgToString;
				avgToString << avg;
				fstream << avgToString.str() << "\t";
			}	
		}
		fstream << endl;
	}

	fstream.close();

	ofstream rstream("metaplot_outfile.R");
	rstream << "library(ggplot2)\nlibrary(reshape2)\n";
	rstream << "pdf(file=\"metaplot_outfile.pdf\", width=12, height=8)\n";
	rstream << "plot<-read.table(\"metaplot_outfile.txt\", header=T)\n";
	rstream << "plot.melt<-melt(plot[,c('bp', " << nameStrR << ")], id.vars=1)\n";
	rstream << "ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) +\n";
	rstream << "geom_line() +\n";
	rstream << "geom_smooth() +\n";
	rstream << "theme_bw() +\n";
	rstream << "ggtitle(\"Metaplot\") +\n";
	rstream << "theme(panel.grid.minor=element_blank()) +\n";
	rstream << "scale_colour_brewer(palette=\"Set1\", name=\"Bed\")\n";

	rstream.close(); 

	return outfileName;
}
