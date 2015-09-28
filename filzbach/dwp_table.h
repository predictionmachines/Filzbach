#ifndef DWP_TABLE_H
#define DWP_TABLE_H
typedef char name_t[256];
struct table
{
	name_t name;
	name_t sourcename;
	int numcolumns;
	int numrows;
	name_t *columnname;
	double **data;
};

#include "filzbach.h"

void table_showinfo(char[]);

#endif // DWP_TABLE_H