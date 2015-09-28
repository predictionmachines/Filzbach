CPPFLAGS = -I../include -DFILZBACH_LIB -std=c++11
objects = catherine_loader.o metropolis.o params.o stats.o table.o
libfilzbach.a: $(objects)
	ar rs libfilzbach.a $(objects)
# Backup.h
catherine_loader.o: catherine_loader.cpp catherine_loader.h stdafx.h
metropolis.o: metropolis.cpp dwp_metropolis.h stdafx.h
params.o: params.cpp params.h stdafx.h
stats.o: stats.cpp dwp_stats.h stdafx.h
table.o: table.cpp dwp_table.h stdafx.h
# targetver.h


clean:
	rm -f libfilzbach.a $(objects)
