#ifndef BATCH_TEST_RUN
	#include "examples.h"
#ifdef WIN32
	#define PAUSE system("pause");
#else
	#define PAUSE system("read -p \"Press ENTER key to continue...\" __dummy__");
#endif
#else
	#define PAUSE
#endif
#ifdef _OPENMP
	#include "omp.h"
#endif
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
