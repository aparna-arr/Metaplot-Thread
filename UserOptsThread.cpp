#include"UserOptsThread.h"
using namespace std; 

UserOpts::UserOpts(int argc, char *argv[])
{
	if (argc < 5)
	{
		printUsage();
		throw 1;
	}

	/* INIT variables */

	nameStr = "";
	nameStrR = "";
	maxWindow = 0;
	bedNum = 0;
	chromNum = 0;
	preprocessWigOpt = false;
	step = 1;
	window = 100;
	preprocessBedOpt = false;
	mode = -1;
	wigSplitDir = "";
	bedSplitDir = "";
	isMonteCarloOpt = false;
	bedDir = "";
	onlyFindChrs = false;
	inputChrs = "";
	threads.bedThreads = 0;
	threads.chromThreads = 0;
	threads.divThreads = 0;
	
	/* Handle opts */
	string msg = "Handling opts";
	print(msg);

	string args = handleOpts(argc, argv);
	stringstream input(args);	
	
	if (!(input >> maxWindow))
	{
		string err = "Your maxWindow is not a valid number!";
		print (err);	
		throw 1;
	}	

	stringstream iToS;
	iToS << maxWindow;
	msg = "maxWindow: " + iToS.str();

	input >> wigFile;

	msg += "\nwigFile: " + wigFile + "\n";
	
	vector<string> bedsAndNames;
	string bedOrName;

	if (bedDir == "")
	{
		while((input >> bedOrName))
			bedsAndNames.push_back(bedOrName);
	}
	else
		readBedsFromDir(bedsAndNames);

	if (bedsAndNames.size() % 2 != 0)
	{
		string err = "Something is wrong with your bedfiles and names. Odd number returned.\n";
		
		iToS.str("");

		iToS << bedsAndNames.size();
		err += "size is " + iToS.str() + "\n";
		err += "Your input is:\n";
		for (vector<string>::iterator iter = bedsAndNames.begin(); iter != bedsAndNames.end(); iter++)	 
			err += *iter + "\n"; 

		print(err);	
		throw 1;
	}
	
	bedNum = bedsAndNames.size() / 2;

	if (threads.bedThreads > bedNum)
	{
		string err = "Your number of bedThreads is > bedNum!";
		print(err);	
		throw 1;
	}
	else if (threads.bedThreads == -1)
		threads.bedThreads = bedNum;


	msg += "Bedfiles:\n";

	for(vector<string>::iterator iter = bedsAndNames.begin(); iter != bedsAndNames.end() - bedNum; iter++)
	{
		msg += *iter + "\n";
		bedFiles.push_back(*iter);
	}
	
	msg += "Names:\n";

	for (vector<string>::iterator iter = bedsAndNames.begin() + bedNum; iter != bedsAndNames.end(); iter++)
	{
		msg += *iter + "\n";
		names.push(*iter);
	}

	/* FIXME later get rid of the name stack and generate strings in above for loop */

	nameStrR = "'" + names.top() + "'";
	nameStr = names.top() + "\t" + nameStr;
	names.pop();

	while(!names.empty())
	{
		nameStr = names.top() + "\t" + nameStr;
		nameStrR = "'" + names.top() + "', " + nameStr;
		names.pop();
	}

	print(msg);
}

void UserOpts::printUsage(void)
{
	string usage = "usage: metaplot-thread [opts] <max_window> <wigfile> <bedfile1 bedfile2 ... > <bedname1 bedname2 ... >\n";
	usage += "opts are:\n";
	usage += "\t--preprocessWig: preprocess wig (smooth)\n";
	usage += "\t--window <INT>: default 100. Window size to smooth wig file\n";
	usage += "\t--step <INT>: default 1. Step size to smooth wig file\n";
	usage += "\t--preprocessBed <INT>: preprocess bed. Modes are:\n\t\t0: TSS\n\t\t1: Centered\n";
	usage += "\t--readSplitWig <DIR>: read a chromosome-split wig file (already preprocessed) from DIR\n";
	usage += "\t--readSplitBed <DIR>: read chromosome-split bed files (already preprocessed) from DIR\n";
	usage += "\t--readAllBedsInDir <DIR>: use all bed files in DIR as input, instead of listing on CLI. Do NOT list bedfiles or bednames on CLI--names will be filenames\n";
	usage += "\t--findChrs: pre-find common chrs in all beds and wig, output file of common chrs and exit.";	
	usage += "\t--inputChrs <FILENAME>: input file produced by --findChrs so common chrs do not need to be found during processing\n";
	usage += "\t--monteCarlo: This run is a Monte Carlo simulation. Final outfile will be horizontally averaged.\n";
	usage += "\nTHREADING\n";
	usage += "\t--threadBeds <OPTIONAL INT>: For each and any processing step, instead of FOR $BEDS, thread bed files. If no argument, produce 1 thread per bedfile. If argument, then run only <INT> threads (must be < # bedfiles and # beds % INT MUST == 0!)\n";
	usage += "\t--threadChroms <OPTIONAL INT>: For each and any processing step, instead of FOR $CHROMS, thread chromosomes. If no argument, produce 1 thread per chromosome. If argument, then run only <INT> threads (must be < # chroms and # chroms % INT MUST == 0!)\n";
	usage += "\t--threadDivisions <INT>: Divide each bedfile into <INT> divisions for the map() and calculate() processing steps\n";
	print(usage);
}

void UserOpts::preprocess(void)
{
	/* get common chrs */
	
	map<string,int> allChrs;

	if (inputChrs != "")
	{
		// we already have common chrs in file
		// read and set
		ifstream in(inputChrs.c_str());
		string line;
	
		while(getline(in,line))
			commonChrs[line] = new vector<int>(bedNum + 1); 
		
		in.close();

		chromNum = commonChrs.size();		
	}
	else
	{
		// find common chrs
		commonChrsWig(allChrs);
		
		for (vector<string>::iterator iter = bedFiles.begin(); iter != bedFiles.end(); iter++)
			commonChrsBed(allChrs, *iter);

		analyzeAndPrintChrs(allChrs);

		chromNum = commonChrs.size();

		if (onlyFindChrs)
			throw 0; // exit cleanly
	}		

	if (threads.chromThreads > chromNum)
	{
		string err = "Error! ChromThreads > chromNum!";
		print(err);
		throw 1;
	}
	else if (threads.chromThreads == -1)
		threads.chromThreads = chromNum;

	/* preprocess bed and wig */

	// FIXME does wig need to be preprocessed? Since we don't do sliding window until after map. It's an empty function right now

	if (threads.bedThreads == 0) // do not thread bed and wig
	{
//		preprocessBed();
		preprocessBedThread(bedNum, 0); // works the same
		preprocessWig();
	}
	else
	{
		// thread bed and wig
//		int bedsPerThread = bedNum;
		
		print("In Preprocess, going to thread beds");
		thread threadAr[threads.bedThreads+1];

		threadAr[0] = thread(&UserOpts::preprocessWig, this);

		int bedsPerThread = bedNum / threads.bedThreads; // FIXME require bedNum to be a multiple of bedThreads?
		for (int i = 1; i <= threads.bedThreads; i++)
			threadAr[i] = thread(&UserOpts::preprocessBedThread, this, bedsPerThread, bedsPerThread * (i - 1));
				
		for (int j = 0; j <= threads.bedThreads; j++)
			threadAr[j].join();
	}

	print("Done processing beds and wig");
}

void UserOpts::split(vector<Bed *> * &bedSplit, vector<Wig *> &wigSplit)
{

	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;


	print("In split function" + debug);
	if (threads.bedThreads == 0) // do not thread bed and wig
	{
		print("Not threading" + debug);
		if (bedSplitDir != "")
		{
			print ("Beds already split" + debug);
			// bedfiles are already split
			// read bed files
			readBedSplit(bedSplit); // FIXME copy over this function -> done
		}
		else 	
			splitBeds(bedSplit, 0, bedNum);

		print("done splitting beds" + debug);

		splitWig(wigSplit);

		print("done splitting wigs" + debug);
	}
	else
	{
		print("split: going to thread beds" + debug);
		thread threadAr[threads.bedThreads+1];
		
		// thread Wig
		threadAr[0] = thread(&UserOpts::splitWig, this, ref(wigSplit));

		if (bedSplitDir != "")
		{
			// bedfiles are already split
			// read bed files
			// can't do this threaded
			readBedSplit(bedSplit); // FIXME copy over this function -> done	
			threadAr[0].join();
		}
		else
		{
			int bedsPerThread = bedNum / threads.bedThreads;

			

			for (int i = 1; i < threads.bedThreads + 1; i++)
			{
				print("making a thread");
/*
				stringstream debug;
				debug << this_thread::get_id();

				string d = "id: " + debug.str();

				print(d);				
*/
				// FIXME EVERYTHING PUSHES TO bedSplit--NOT THREAD SAFE
				threadAr[i] = thread(&UserOpts::splitBeds, this, ref(bedSplit), (i - 1) * bedsPerThread, bedsPerThread);
			}

			for (int j = 0; j < threads.bedThreads + 1; j++)
				threadAr[j].join();
		}

		print("Done threading");
	}
	print("Done splitting beds and wigs");
}
	
int UserOpts::getBedNum(void)
{
	return bedNum;
}

int UserOpts::getMaxWindow(void)
{
	return maxWindow;
}

int UserOpts::getShift(void)
{
	return step;
}

int UserOpts::getWindow(void)
{
	return window;
}

std::string UserOpts::getNameStr(void)
{
	return nameStr;	
}

std::string UserOpts::getNameStrR(void)
{
	return nameStrR;
}

bool UserOpts::isMonteCarlo(void)
{
	return isMonteCarloOpt;
}

ThreadInfo * UserOpts::getThreadInfo(void)
{
	return &threads;
}

map< string, vector<int> *> * UserOpts::getCommonChrs(void)
{
	return &commonChrs;
}

/* PRIVATE FUNCTIONS */

stack<string> UserOpts::getFilesInDir(string dir)
{
	stack<string> filenames;
	
	DIR *dp;
	struct dirent *dirp;
	struct stat filestat;
	string filepath;

	dp = opendir(dir.c_str());
	
	if (dp == NULL)
	{
		string err = "Error opening dir " + dir;
		print(err);
		throw 1;
	}

	while ((dirp = readdir(dp)))
	{
		filepath = dir + "/" + dirp->d_name;
		// skip dirs or weird files
		if (stat(filepath.c_str(), &filestat)) continue;
		if (S_ISDIR(filestat.st_mode)) continue;
		filenames.push(filepath);
	}

	return filenames;
}

string UserOpts::handleOpts(int argc, char * argv[])
{
	string currArg(argv[1]); // 0 is metaplot-thread
	int pos; // within currArg
	int i = 1; // which arg we're on

	while (i < argc && (pos = currArg.find("--")) != (signed) string::npos)
	{
		string opt = currArg.substr(pos+2, currArg.length());
	
		// Opts with an INT argument
		if (opt == "step" || opt == "window" || opt == "preprocessBed" || opt == "threadBeds" || opt == "threadChroms" || opt == "threadDivisions")
		{
			bool testNext = true;
			bool isNextInt = false;
			int num;
			// make sure we can test next arg without segfault
			if (i == argc -1)
			{
				if (opt != "threadBeds" && opt != "threadChroms")
					throw 2; // arg requires int
				else
					isNextInt = false;	
			}

			if (testNext) // can test next arg
			{
				stringstream test;
				test << argv[i+1];
				isNextInt = (test >> num);

				if (!isNextInt && (opt != "threadBeds" && opt != "threadChroms"))
					throw 2; // arg requires int
				else if (isNextInt)
					i++;
			}

			if (opt == "step")
				step = num; 
			else if (opt == "window")
				window = num;
			else if (opt == "threadDivisions")
				threads.divThreads = num;
			else if (opt == "preprocessBed")
			{
				preprocessBedOpt = true;
				mode = num;
			}
			else if (opt == "threadBeds")
			{
				if (!isNextInt)
					threads.bedThreads = -1; // thread all
				else
					threads.bedThreads = num; // FIXME check later in code if < #beds! -> done in UserOpts()
			}
			else if (opt == "threadChroms")
			{
				if (!isNextInt)
					threads.chromThreads = -1; // thread all
				else
					threads.chromThreads = num; // FIXME check later in code if < #chroms! -> done in preprocess()
			}
			else
			{
				string err = "Unrecognized option: " + opt;
				print(err);	
				throw 1;
			}
		}
		else if (opt == "readSplitWig" || opt == "readSplitBed" || opt == "readAllBedsInDir" || opt == "inputChrs")
		{ // Opts with a STRING argument
			if (i == argc -1)
			{
				throw 3; // arg requires string
			}
			
			string fileOrDir(argv[i+1]);
			
			if (opt == "readSplitWig")
				wigSplitDir = fileOrDir;
			else if (opt == "readSplitBed")
				bedSplitDir = fileOrDir;
			else if (opt == "readAllBedsInDir")
				bedDir = fileOrDir;
			else if (opt == "inputChrs")
				inputChrs = fileOrDir;	
			else 
			{
				string err = "Unrecognized option: " + opt;
				print(err);
				throw 1;
			}
		}
		else if (opt == "preprocessWig" || opt == "findChrs" || opt == "monteCarlo")
		{ // Opts with no argument
			if (opt == "preprocessWig")
				preprocessWigOpt = true;
			else if (opt == "findChrs")
				onlyFindChrs = true;
			else if (opt == "monteCarlo")
				isMonteCarloOpt = true;
			else
			{
				string err = "Unrecognized option: " + opt;
				print(err);
				throw 1;
			}
		}
		else
		{
			string err = "Unrecognized opt: " + opt;
			print(err);
			throw 1; 
		}

		i++;

		if (i < argc)
			currArg = string(argv[i]);
	} // while
	
	string rem = "";

	for (int j = i; j < argc; j++)
		rem += string(argv[j]) + " ";

	return rem;
}

// PROBABLY OBSOLETE
void UserOpts::preprocessBed(void)
{

}

void UserOpts::preprocessWig(void)
{

}

void UserOpts::preprocessBedThread(int bedsPerThread, int startBed)
{
	print("In PreProcessBedThread");

	if (mode == -1)
		return;

	for (vector<string>::iterator iter = bedFiles.begin() + startBed; iter != bedFiles.begin() + startBed + bedsPerThread; iter++)
	{
		ifstream file((*iter).c_str());
		string line;
		(*iter) = (*iter) + ".tmp";
		ofstream tmpfile((*iter).c_str());
		
		while (getline(file, line))
		{
			int start, end;
			char strand;
			string chr;
			stringstream linestream(line);
			linestream >> chr >> start >> end >> strand;	
			
			int more, less;

			if (mode == 0) // TSS
			{
				if (strand == '+')
				{
					more = start + (maxWindow / 2);
					less = start - (maxWindow / 2);
				}
				else	
				{
					more = end + (maxWindow / 2);	
					less = end - (maxWindow / 2);
				}	
			}
			else if (mode == 1)
			{
				int middle = (end - start) / 2 + start;
				more = middle + (maxWindow / 2);
				less = middle - (maxWindow / 2);
			}

			tmpfile << chr << "\t" << less << "\t" << more << "\t" << strand << endl;
		} // while
		file.close();
		tmpfile.close();
	} // for
	print("End preprocess bed");
}

void UserOpts::calcSlidingWindow(std::vector<Peak> * wigBlocks)
{

}

void UserOpts::splitBeds(vector<Bed *> * &arOfBedfiles, int startBed, int bedsPerThread)
{


	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("In splitBeds " + debug);

	stringstream startBedToI, bedsPerThreadToI;
	startBedToI << startBed;
	bedsPerThreadToI << bedsPerThread;

	print("splitBeds(): startBed is " + startBedToI.str() + " bedsPerThread is " + bedsPerThreadToI.str() + debug);


	for (vector<string>::iterator iter = bedFiles.begin() + startBed; iter != bedFiles.begin() + startBed + bedsPerThread; iter++)
	{
		print("\tIn for loop" + debug);
		ifstream file((*iter).c_str());
//		print("\tAfter iter call" + debug);

		string line;
		string prevChr = "INIT";

		Bed * currBed = new Bed();
		
		stringstream iToS;
		iToS << startBed;
/*		
		string start = "startBed is " + iToS.str();
		//string debug = "File is " + (*iter);
		print(start + debug);
		print("\tbefore while loop" + debug);
*/
		while (getline(file, line))
		{
//			print("\tin while loop"+ debug);

			print("\ton line " +line + debug);
			stringstream linestream(line);	
			string chr; 
			int start, end;
			char strand;
			linestream >> chr >> start >> end >> strand;
				
//			print("\tafter get line and split line" + debug);

			if (commonChrs.find(chr) == commonChrs.end())
				continue;

			if (chr != prevChr)
			{
//				print("\t\tin first if" + debug);
				if (prevChr != "INIT")
				{
//					print("\t\t\tin if" + debug);
					currBed->printPeaks();
					Bed *  tmp = new Bed();
					*tmp = *currBed;
		
//					print("\t\t\tafter assigning tmp" + debug);

					arOfBedfiles[iter - bedFiles.begin()].push_back(tmp);
//					*(commonChrs[prevChr]->begin() + (iter - bedFiles.begin() + 1)) = arOfBedfiles[iter - bedFiles.begin()].size() - 1;
					print("\t\t\tafter arOfBedfiles push" + debug);
/*

					// FIXME NEED A MUTEX ON THIS, OR COPY OVER!	
					// well each should be modifying--inserting--a unique place ...
					stringstream iToS;
					iToS << iter - bedFiles.begin();
					print("\t\t\tcurrBed is " + iToS.str() + debug);					
*/						

					*(commonChrs[prevChr]->begin() + (iter - bedFiles.begin() + 1)) = arOfBedfiles[iter - bedFiles.begin()].size() - 1;
//					commonChrs[prevChr]->insert(commonChrs[prevChr]->begin() + (iter - bedFiles.begin() + 1), arOfBedfiles[iter-bedFiles.begin()].size() - 1);
//					print("\t\t\tend of if" + debug);
					
					stringstream bedPosS;	
					bedPosS << arOfBedfiles[iter - bedFiles.begin()].size() - 1;
					print("Adding to chr " + prevChr + " the pos " + bedPosS.str() + debug);

				}
			
				*currBed = Bed(chr, maxWindow);
				currBed->generateFilename(iter - bedFiles.begin());
				prevChr = chr;
			}
			
			currBed->addPeak(start, end, strand);
//			print("\t\tend of first if" + debug);
		}
		print("After while");
		currBed->printPeaks();
		Bed * tmp = new Bed();
		*tmp = *currBed;
		arOfBedfiles[iter - bedFiles.begin()].push_back(tmp);
//		*(commonChrs[prevChr]->begin()+ (iter - bedFiles.begin() + 1)) = arOfBedfiles[iter - bedFiles.begin()].size() - 1;
		*(commonChrs[prevChr]->begin() + (iter - bedFiles.begin() + 1)) = arOfBedfiles[iter - bedFiles.begin()].size() - 1;
//		commonChrs[prevChr]->insert(commonChrs[prevChr]->begin() + (iter - bedFiles.begin() + 1), arOfBedfiles[iter-bedFiles.begin()].size() - 1);

		print("splitBeds(): before file close" + debug);

		file.close();
	}

	print("splitBeds(): done with split beds" + debug);
}

void UserOpts::splitWig(vector<Wig *> &wigSplit)
{
	if (wigSplitDir != "")
		return readWigSplit(wigSplit); // FIXME copy over this function -> done

	// can't split this function to thread it
	// FIXME copy over rest of this function -> done

	stringstream id;
	id << this_thread::get_id();

	string debug = " My id: " + id.str() + " " ;

	print("splitWig(): start" + debug);

	ifstream fstream(wigFile.c_str());
	string line;

	Wig * currWig = new Wig();
	int currSpan = 0;
	string currChr = "INIT";

	while(getline(fstream, line))
	{

//		print("splitWig(): in while" + debug);

		print("splitWig(): while start line is " + line + debug);
		stringstream linestream(line);
		stringstream test(line);
		int pos;
		double val;
		string substring;
		if (!(test >> pos))
		{
			// on a non-integer line
			// FIXME do a check here for variableStep, track lines, and #'s
			linestream >> substring; // variableStep
			linestream >> substring; // chrom=
	
//			print("splitWig(): substring is " + substring + debug);
			string chr = substring.substr(substring.find("=") + 1, substring.length());

// DO NOT UNCOMMENT
// AT ALL
// OR ELSE CRAZY STUFF THAT IS HARD TO DEBUG HAPPENS
//
//			if (commonChrs.find(chr) == commonChrs.end())
//				continue;

			linestream >> substring; // span=
//			print("splitWig(): substring is " + substring + debug);

			string span = substring.substr(substring.find("=") + 1, substring.length());

			if (chr != currChr)
			{
				if (currChr != "INIT" && commonChrs.find(currChr) != commonChrs.end())
				{
					// print prev Chr peaks
					currWig->printPeaks();
					Wig * tmp = new Wig();
					*tmp = *currWig;
					wigSplit.push_back(tmp);	

					*(commonChrs[currChr]->begin()) = wigSplit.size() - 1;
					
//					commonChrs[currChr]->insert(commonChrs[currChr]->begin(), wigSplit.size() - 1);
				}
				else if (currChr != "INIT" && commonChrs.find(currChr) == commonChrs.end())
				{
					currWig->clearPeaks();
				}
				
				currChr = chr;
				*currWig = Wig(chr);
				currWig->generateFilename();
				
			}
			stringstream itoS(span);
			itoS >> currSpan;
//			print("splitwig(): currSpan is " + itoS.str() + debug);
		}
		else
		{
			linestream >> pos >> val;
		
			print ("splitwig(): before addPeak linestream is " + linestream.str() + debug);	
//			stringstream spanStr;
//			spanStr << currSpan;
//			print("In Else currSpan is " + spanStr.str() + debug);
			currWig->addPeak(pos, pos+currSpan, val);
		}

//		print("splitWig(): end of while" + debug);
	}


	print("splitWig(): after while" + debug);


	if (commonChrs.find(currChr) != commonChrs.end())
	{

		currWig->printPeaks();

		print("splitWig(): After print peaks" + debug);
	
		Wig * tmp = new Wig();

		print("splitWig(): after new" + debug);

		*tmp = *currWig;
		wigSplit.push_back(tmp);	

		print ("splitWig(): after handle tmp" + debug);

		print("splitWig(): currChr is " + currChr + debug);

		*(commonChrs[currChr]->begin()) = wigSplit.size() - 1;
//		commonChrs[currChr]->insert(commonChrs[currChr]->begin(), wigSplit.size() - 1);

		print("splitWig(): after assign commmonChr" + debug);
	}


	fstream.close();
	print("splitWig(): Done" + debug);
}

void UserOpts::readBedSplit(std::vector<Bed *> * &arOfBedfiles)
{
	stack<string> splitFiles = getFilesInDir(bedSplitDir);
	int pos;
	while(!splitFiles.empty())
	{
		if ( (pos = (splitFiles.top()).rfind("_")) != (signed) string::npos )
		{
			string chr = (splitFiles.top()).substr(0,pos);
			chr = chr.substr(chr.rfind("/")+1);
		
			if (commonChrs.find(chr) == commonChrs.end())
				continue;

			stringstream findBedNum((splitFiles.top()).substr(pos+1));
			int loc_bedNum;

			if (!(findBedNum >> loc_bedNum))
			{
				string err = "There is no bednum after _ in file " + splitFiles.top();
				throw 1;

			}

			Bed * tmp = new Bed(chr, maxWindow);
			tmp->setFilename(splitFiles.top());
			
			arOfBedfiles[loc_bedNum].push_back(tmp);	
			//*(commonChrs[chr]->begin() + loc_bedNum + 1) = arOfBedfiles[loc_bedNum].size() - 1;
			commonChrs[chr]->insert(commonChrs[chr]->begin() + (loc_bedNum + 1), arOfBedfiles[loc_bedNum].size() - 1);
		}
		else
		{
			string err = "Split file " + splitFiles.top() + " has a weird filename (not *_*)!";	
			print(err);
			throw 1;

		}
		splitFiles.pop();
	}
}

void UserOpts::readWigSplit(std::vector<Wig *> &wigSplit)
{
	stack<string> splitFiles = getFilesInDir(wigSplitDir);
	int pos;
	while(!splitFiles.empty())
	{
		if ( (pos = (splitFiles.top()).find("_wig.bed")) != (signed) string::npos )
		{
			string chr = (splitFiles.top()).substr(0,pos);
			chr = chr.substr(chr.rfind("/")+1);

			if (commonChrs.find(chr) != commonChrs.end())
			{
				Wig * tmp = new Wig(chr);
				tmp->setFilename(splitFiles.top());
		
				wigSplit.push_back(tmp);
//				*(commonChrs[chr]->begin()) = wigSplit.size() - 1;
				commonChrs[chr]->insert(commonChrs[chr]->begin(), wigSplit.size() - 1);

			}
		}
		else
		{
			string err = "Split file " + splitFiles.top() + " has a weird filename (not *_wig.bed)!";	
			print(err);
			throw 1;
		}
		
		splitFiles.pop();
	}
}

void UserOpts::readBedsFromDir(std::vector<std::string> &bedsAndNames)
{
	stack<string> files = getFilesInDir(bedDir);

	vector<string> names;

	while(!files.empty())
	{
		bedsAndNames.push_back(files.top());
		names.push_back(files.top());
		files.pop();			
	}

	bedsAndNames.insert(bedsAndNames.end(), names.begin(), names.end());
}

void UserOpts::commonChrsWig(map<string,int> &allChrs)
{
	ifstream file(wigFile.c_str());
	string line;

	map<string,int> myChrs;

	while(getline(file,line))
	{
		stringstream linestream(line);
		string word;
		linestream >> word;

		if (word != "variableStep")
			continue;
	
		linestream >> word;
		string chr = word.substr(word.find("=") + 1, word.length());

		if (myChrs.find(chr) == myChrs.end())
			allChrs[chr]++;	// will create element & init it to 0 if not exist
		// http://stackoverflow.com/questions/5616421/increment-mapstring-int-using-operator	

		myChrs[chr]++;
	}	

	file.close();
}
// FIXME: multiple chrs in one file i.e:
// chr1 10 20
// chr1 30 40
// breaks these commonchrs functions!
// FIXED
void UserOpts::commonChrsBed(map<string,int> &allChrs, string filename)
{
	ifstream file(filename);
	string line;

	map<string,int> myChrs;	

	while(getline(file,line))
	{
		stringstream linestream(line);
		string chr;
		linestream >> chr;

		// only increment chrs that are already there
		if (allChrs.find(chr) != allChrs.end() && myChrs.find(chr) == myChrs.end())
			allChrs[chr]++;

		myChrs[chr]++;
	}
	file.close();
}

void UserOpts::analyzeAndPrintChrs(map<string,int> allChrs)
{
	ofstream out("common_chrs.txt");

	for (map<string,int>::iterator iter = allChrs.begin(); iter != allChrs.end(); iter++)
	{
		if (iter->second == bedNum + 1)
		{
			commonChrs[iter->first] = new vector<int>(bedNum + 1);
			out << iter->first << endl;
		}
	}

	out.close();
}
