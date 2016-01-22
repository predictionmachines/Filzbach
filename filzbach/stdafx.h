// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#pragma once

#ifndef STDAFX_H
#define STDAFX_H

// #define DEBUG
#define _USE_MATH_DEFINES
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <limits>
const double MISSING_VALUE = -999.0;
const char ws_dir[] = "./workspace";


#ifdef _DEBUG
#define ASSERT(cond,msg) if (!(cond)) {throw (msg);}
#define CHECK(cond,msg) ASSERT(cond,msg)
#define HEAP_CHECK ASSERT(_heapchk()==_HEAPOK, "Heap detroyed!")
#else
#define ASSERT(cond,msg)
#ifdef WIN32
#define CHECK(cond,msg) if (!(cond)) {printf("\nERROR: *** %s ***\nIn %s\n(%s:%d)\n", msg, __FUNCSIG__,  __FILE__, __LINE__);exit(2);}
#else
#define CHECK(cond,msg) if (!(cond)) {printf("\nERROR: *** %s ***\nIn %s:%d\n", msg,  __FILE__, __LINE__);exit(2);}
#endif
#define HEAP_CHECK
#endif


#ifndef WIN32
#define __max(a,b) (((a) > (b)) ? (a) : (b))
#define __min(a,b) (((a) < (b)) ? (a) : (b))
#define strcpy_s(dest,bufsize,source) strncpy(dest,source,bufsize)
#define strcat_s(dest,bufsize,source) strncat(dest,source,bufsize)
#include <cstdarg>
template <size_t size>
int sprintf_s(
   char (&buffer)[size],
   const char *format, ... 
	) {
	va_list vargs; 
	va_start(vargs, format);  
	return vsnprintf(buffer, size, format, vargs);
}
#endif

#endif // TDAFX_H_H 
