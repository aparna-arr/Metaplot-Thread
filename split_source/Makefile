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
