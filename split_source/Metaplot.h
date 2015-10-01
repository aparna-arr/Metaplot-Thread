#ifndef METAPLOT_H
#define METAPLOT_H

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

#include"MetaplotRegion.h"
#include"Chromosome.h"

//bool calculateMetaplot(UserOpts &input, std::vector<Bed *> * &bedByChrs, std::vector<Wig *> &wigByChrs, std::vector<std::string> commonChrs, MetaplotRegion * &region);
void calculateMetaplot(int bedNum, std::vector<Bed *> * &bedByChrs, std::vector<Wig *> &wigByChrs, std::vector<std::string> commonChrs, MetaplotRegion * &region);

void monteCarloMetaplot(std::string file, int bedNum);


void debug(MetaplotRegion * &region, int bedNumber, std::string nameStr, std::string nameStrR);

std::string printResults(MetaplotRegion * &region, int bedNumber, std::string nameStr, std::string nameStrR, int maxWindow);

#endif
