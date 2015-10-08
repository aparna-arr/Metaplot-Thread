CXX := g++
CXXFLAGS := -Wall -std=c++11 -g -pthread
OBJECTS := Main.o MetaplotThread.o ChromosomeThread.o MetaplotRegionThread.o UserOptsThread.o

metaplot-thread: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

Main.o: Main.cpp MetaplotThread.h ChromosomeThread.h MetaplotRegionThread.h UserOptsThread.h
	$(CXX) $(CXXFLAGS) -c Main.cpp

MetaplotThread.o: MetaplotThread.cpp MetaplotThread.h ChromosomeThread.h MetaplotRegionThread.h
	$(CXX) $(CXXFLAGS) -c MetaplotThread.cpp

ChromosomeThread.o: ChromosomeThread.cpp ChromosomeThread.h 
	$(CXX) $(CXXFLAGS) -c ChromosomeThread.cpp

MetaplotRegionThread.o: MetaplotRegionThread.cpp MetaplotRegionThread.h 
	$(CXX) $(CXXFLAGS) -c MetaplotRegionThread.cpp

UserOptsThread.o: UserOptsThread.cpp UserOptsThread.h Chromosome.h
	$(CXX) $(CXXFLAGS) -c UserOptsThread.cpp

clean-thread:
	rm -f $(OBJECTS) metaplot-thread

metaplot: RunMetaplot.o Metaplot.o Chromosome.o MetaplotRegion.o UserOpts.o
	g++ -Wall -o $@ RunMetaplot.o Metaplot.o Chromosome.o MetaplotRegion.o UserOpts.o
RunMetaplot.o: RunMetaplot.cpp Metaplot.h Chromosome.h MetaplotRegion.h UserOpts.h
	g++ -Wall -c RunMetaplot.cpp
Metaplot.o: Metaplot.cpp Metaplot.h MetaplotRegion.h
	g++ -Wall -c Metaplot.cpp
Chromosome.o: Chromosome.cpp Chromosome.h
	g++ -Wall -c Chromosome.cpp
MetaplotRegion.o: MetaplotRegion.cpp MetaplotRegion.h
	g++ -Wall -c MetaplotRegion.cpp
UserOpts.o: UserOpts.cpp UserOpts.h Chromosome.h
	g++ -Wall -c UserOpts.cpp
clean: 
	rm -f Metaplot.o RunMetaplot.o Chromosome.o MetaplotRegion.o UserOpts.o metaplot
