#include "stdafx.h"
#include "dwp_table.h"

const int SHRINK = 0;
const int GROW = 1;

int num_tables = 0; //current number of tables
int maxrows = 100000;
table *mytable;


/**************************/
/**  Local Declarations  **/
/**************************/

/*
 * Allocates memory for a new table and returns tid of the new table
 */
int talloc()
{
	// try find a deleted table
	for (int t = 0; t < num_tables; t++)
	{
		if (mytable[t].numcolumns == 0
			&& mytable[t].numrows == 0
			&& mytable[t].data == NULL
			&& mytable[t].columnname == NULL
			&& mytable[t].name[0] == '\0'
			&& mytable[t].sourcename == '\0')
			return t;
	}
	table *temptable = new table[num_tables + 1];
	CHECK(temptable != NULL, "Memory allocation failed");

	for (int i = 0; i < num_tables; i++)
		temptable[i] = mytable[i];
	temptable[num_tables].numcolumns = 0;
	temptable[num_tables].numrows = 0;
	temptable[num_tables].columnname = NULL;
	temptable[num_tables].data = NULL;
	temptable[num_tables].name[0] = '\0';
	temptable[num_tables].sourcename[0] = '\0';

	delete[] mytable;
	mytable = temptable;

	HEAP_CHECK;
	return num_tables++;
}

/*
 * Frees the memory for all columns of a table with the given id.
 */
void tfreecolumns(int tid)
{
	// delete old column names and column data
	ASSERT((tid >= 0) && (tid < num_tables), "tid out of range");
	for (int i = 0; i < mytable[tid].numcolumns; i++)
	{
		//delete[] mytable[tid].columnname[i];
		delete[] mytable[tid].data[i];
	}
	delete[] mytable[tid].columnname;
	delete[] mytable[tid].data;
	HEAP_CHECK;
}

/*
 * Frees the memory used by the specified table.
 *
 * index: The table id for the table to free.
 */
void tfree(int tid)
{
	// cannot shrink mytable[] array as this will change tid's of existing tables
	// free memory occupied by internal table structures
	ASSERT((tid >= 0) && (tid < num_tables), "tid out of range");
	tfreecolumns(tid);
	mytable[tid].columnname = 0;
	mytable[tid].data = 0;
	mytable[tid].numcolumns = 0;
	mytable[tid].numrows = 0;
	mytable[tid].name[0] = '\0';
	mytable[tid].sourcename[0] = '\0';
}

/*
 * Shrinks or growths a table.
 * mode: SHRINK or GROW
 * count: the number to change the table.
 */
int tchange(int tid, int mode, int count)
{
	// get current number of columns and rows for this table
	ASSERT((tid >= 0) && (tid < num_tables), "tid out of range");
	ASSERT((mode == GROW) || (count <= mytable[tid].numrows), "invalid combination of mode and count");
	int ncol = mytable[tid].numcolumns;
	int nrow = mytable[tid].numrows;

	// nothing to change
	if (count == 0)	return nrow;

	int newrowscount = mode == GROW ? nrow + count : nrow - count;

	// allocate temporary memory for column names and new column data 
	name_t *tempcolname = new name_t[ncol];
	double **tempcoldata = new double*[ncol];
	for (int i = 0; i < ncol; i++)
	{
		tempcoldata[i] = new double[newrowscount];
	}

	// copy old column names and column data into this new space
	for (int i = 0; i < ncol; i++)
	{
		strcpy(tempcolname[i], mytable[tid].columnname[i]);

		// do we need copy everything or just a part?
		int max = mode == GROW ? nrow : newrowscount;

		for (int j = 0; j < max; j++) // copy existing data
			tempcoldata[i][j] = mytable[tid].data[i][j];

		if (mode == GROW)
		{
			for (int j = nrow; j < newrowscount; j++) // and fill with missing value
				tempcoldata[i][j] = MISSING_VALUE;
		}
	}

	// delete old columns
	tfreecolumns(tid);

	// assign new column names and data
	mytable[tid].columnname = tempcolname;
	mytable[tid].data = tempcoldata;
	mytable[tid].numrows = newrowscount;
	HEAP_CHECK;
	return newrowscount;
}

/***************************************************************/
void table_allocate_memory_all(int tid, int ncol, int numrows)
{
	ASSERT((tid >= 0) && (tid < num_tables), "tid out of range");
	tfreecolumns(tid);
	mytable[tid].columnname = new name_t[ncol + 1];
	mytable[tid].data = new double*[ncol + 1];
	mytable[tid].numcolumns = ncol;
	for (int i = 0; i < ncol + 1; i++)
	{
		//mytable[tid].columnname[i] = new char[256];
		mytable[tid].data[i] = new double[numrows];
	}
	mytable[tid].numrows = numrows;
	HEAP_CHECK;
}

/***************************************************************/
void table_allocate_memory_column(int tid, int cid, int numrows)
{
	ASSERT((tid >= 0) && (tid < num_tables), "tid out of range");
	ASSERT((cid >= 0) && (cid < mytable[tid].numcolumns), "cid out of range");
	//delete[] mytable[tid].columnname[cid];
	//mytable[tid].columnname[cid] = new char[256];
	//delete[] mytable[tid].data[cid];
	mytable[tid].data[cid] = new double[numrows];
	HEAP_CHECK;
}

/***************************/
/** External Declarations **/
/***************************/


/*
 * Returns the ID of a local table.
 * If no or more than one table with this name are loaded, -1 is returned.
 *
 * TODO: It should not be possible to load a table with a name that is already taken.
 */
int table_getID(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = -1;
	int hits = 0;

	for (int i = 0; i < num_tables; i++)
	{
		if (strcmp(mytable[i].name, local_name) == 0)
		{
			if (++hits > 1)
				break;
			tid = i;
		}
	}
	switch (hits)
	{
	case 1:
		return tid;
	case 0:
		printf("\n Error: Looking up table ID by local name %s, but no table with that local name exist.\n \n", local_name);
		return -1;
	default:
		printf("\n Error: Looking up table ID by local name %s, but more than one table with that local name exists \n \n", local_name);
		return -1;
	}
}

/*
 * Returns the ID of a column of a local table.
 * If no or more than one column with this name exist, -1 is returned.
 */
int table_getcolumnID(int tid, const char cname[])
{
	CHECK((tid >= 0) && (tid < num_tables), "tid out of range");
	CHECK(cname != NULL, "Column name cannot be null");
	int cid = -1;
	int hits = 0;

	for (int i = 0; i < mytable[tid].numcolumns; i++)
	{
		if (strcmp(mytable[tid].columnname[i], cname) == 0)
		{
			if (++hits > 1)
				return -1;
			cid = i;
		}
	}

	switch (hits)
	{
	case 1:
		return cid;
	case 0:
		printf("\n Error: Looking up column ID, but no column named %s within table %s \n \n", cname, mytable[tid].name);
		return -1;
	default:
		printf("\n Error: Looking up column ID, but more than one column named %s within table %s \n \n", cname, mytable[tid].name);
		return -1;
	}
}

/*
 * Prints information of a table.
 */
void table_showinfo(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int id = table_getID(local_name);

	if (id == -1)
	{
		printf("\nError: No table with local name %s.", local_name);
		exit(-1);
		return;
	}

	printf("\n Info for table with local name: %s", local_name);
	printf("\n Source: %s", mytable[id].sourcename);
	printf("\n Number of columns: %d and rows: %d ", mytable[id].numcolumns, mytable[id].numrows);
	printf("\n Column headings:");

	for (int i = 0; i < mytable[id].numcolumns; i++)
		printf("%s ", mytable[id].columnname[i]);

	printf("\n");
}

/*
 * Puts the table with the given name into a file.
 */
void table_output(const char local_name[], const char file_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	FILE *ofile;

	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\nError: Couldn't output table with local name %s \n \n", local_name);
		return;
	}

#pragma warning(push)
#pragma warning(disable:4996)
	ofile = fopen(file_name, "w");
#pragma warning(pop)
	if (!ofile)
		printf("\nError: couldn't open output file %s \n \n", file_name);
	else
		printf("\nOK: opened output file %s \n \n", file_name);

	// output header
	for (int c = 0; c < mytable[tid].numcolumns; c++)
	{
		fprintf(ofile, "%s\t", mytable[tid].columnname[c]);
	}
	fprintf(ofile, "\n");

	// output data
	for (int r = 0; r < mytable[tid].numrows; r++)
	{
		for (int c = 0; c < mytable[tid].numcolumns; c++)
		{
			fprintf(ofile, "%lf\t", mytable[tid].data[c][r]);
		}
		fprintf(ofile, "\n");
	}

	fclose(ofile);

	return;
}

/*
 * Creates a table with the given name based on the given input file.
 *
 * local_name: The name for the table to create.
 * file_name:  The file to read the data from.
 * numcolumns: The number of columns to read teh data into.
 */
void table_read(const char local_name[], const char file_name[], int numcolumns)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	FILE *ifile;

	int tableindex = talloc();

#pragma warning(push)
#pragma warning(disable:4996)
	ifile = fopen(file_name, "r");
#pragma warning(pop)
	if (!ifile)
	{
		printf("\nError: Can't open input file named %s.\n\n", file_name);
		exit(1);
	}
	else
		printf("\nOK: Opened input file named %s.\n\n", file_name);

	// determine number of rows - numcolumns already provided by user
	int cr_count = 0;
	int c;
	do {
		c = getc(ifile);
		if (c == 10)
			cr_count++;
	} while (c != EOF);
	rewind(ifile);

	// we don't count the header
	int numrows = cr_count - 1;

	// create table
	strcpy(mytable[tableindex].name, local_name);
	strcpy(mytable[tableindex].sourcename, file_name);

	// allocate memory for this table
	table_allocate_memory_all(tableindex, numcolumns, numrows);

	// read header
	char header[256];
	for (int col = 0; col < numcolumns; col++)
	{
#ifdef WIN32
		fscanf_s(ifile, "%255s", header, 256);
#else
		fscanf(ifile, "%255s", header);
#endif
		strcpy(mytable[tableindex].columnname[col], header);
	}

	// read values
	for (int r = 0; r < numrows; r++) // foreach row, skipping first row as already pointing to row 1. 
	{
		for (int col = 0; col < numcolumns; col++) // and each column
		{
			if (!feof(ifile)) // if there is still file left, read in tha values
			{
#ifdef WIN32
				fscanf_s(ifile, "%lf", &mytable[tableindex].data[col][r]);
#else
				fscanf(ifile, "%lf", &mytable[tableindex].data[col][r]);
#endif
			}
			else // fill the rest of the table  with missing values 
			{
				mytable[tableindex].data[col][r] = MISSING_VALUE;
			}
		}
	}

	// output to screen
	table_showinfo(local_name);

	fclose(ifile);

	return;
}

/*****************************************************************/
double table_getvalue(const char local_name[], const char cname[], int row)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	CHECK(cname != NULL, "Column name cannot be null");
	int tid = table_getID(local_name);

	if (tid > -1)
	{
		int column = table_getcolumnID(tid, cname);

		return table_getvalue(local_name, column, row);
	}
	else
	{
		printf("\n Error: Couldn't find a table by name.\n\n");
		return MISSING_VALUE;
	}
}

double table_getvalue(const char local_name[], int cid, int row)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (cid > -1)
	{
		return mytable[tid].data[cid][row];
	}
	else
	{
		printf("\n ERROR: Couldn't get value.\n\n");
		return MISSING_VALUE;
	}
}


/*
 * Returns the number of rows for a table with the given name.
 * If the table is not found, -1 is returned.
 */
int table_numrows(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid;
	tid = table_getID(local_name);

	if (tid > -1)
		return mytable[tid].numrows;

	printf("\n Error: couldn't find table %s\n \n", local_name);
	return -1;

}

/*
 * Returns the number of rows for a table with the given id.
 * If the table is not found, -1 is returned.
 */
int table_numrows(int tid)
{
	if (tid >= 0 && tid < num_tables)
		return mytable[tid].numrows;

	printf("\n Error: couldn't find table with ID %i\n \n", tid);
	return -1;
}

/*
 * Returns the number of columns for a table with the given name.
 * If the table is not found, -1 is returned.
 */
int table_numcolumns(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid;
	tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\n Error: couldn't find table %s\n \n", local_name);
		return -1;
	}

	return table_numcolumns(tid);
}

/*
 * Returns the number of columns for a table with the given ID.
 * If the table is not found, -1 is returned.
 */
int table_numcolumns(int tid)
{
	if (tid > -1)
		return mytable[tid].numcolumns;

	printf("\n Error: couldn't find table with ID %d\n \n", tid);
	return -1;
}


/*
 * Creates a new table with the given name.
 * Return the ID of the table.
 */
int table_create(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = talloc();

	strcpy(mytable[tid].name, local_name);

	mytable[tid].numcolumns = 0;
	mytable[tid].numrows = 0;

	strcpy(mytable[tid].sourcename, "internally_generated");

	return tid;
}

/***************************************************************/
int table_addcolumn(const char local_name[], const char cname[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	CHECK(cname != NULL, "Column name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("Error. Cannot add column %s to table %s. Table cannot be found. \n", cname, local_name);
		return -1;
	}

	// get current number of columns for this table
	int ncol = mytable[tid].numcolumns;

	// allocate temporary memory for column names and column data 
	name_t *tempcolname = new name_t[ncol + 1];
	double **tempcoldata = new double*[ncol + 1];
	for (int i = 0; i < ncol + 1; i++)
	{
		//tempcolname[i] = new char[256];
		tempcoldata[i] = new double[mytable[tid].numrows];
	}

	// copy old column names and column data into this new space
	for (int i = 0; i < ncol; i++)
	{
		strcpy(tempcolname[i], mytable[tid].columnname[i]);

		for (int j = 0; j < mytable[tid].numrows; j++)
			tempcoldata[i][j] = mytable[tid].data[i][j];
	}

	// add new column name
	strcpy(tempcolname[ncol], cname); // if old table had 6 columns, this would write to [6], which is the 7th column

	// fill new columndata with missing values
	for (int i = 0; i < mytable[tid].numrows; i++)
		tempcoldata[ncol][i] = MISSING_VALUE;

	// delete old columns
	tfreecolumns(tid);

	// assign new column names and data
	mytable[tid].columnname = tempcolname;
	mytable[tid].data = tempcoldata;
	mytable[tid].numcolumns++;

	// TODO: We still have to increase the allocated memory for the table

	//int cc,yy;
	//if(mytable[tid].numcolumns<20)
	//{
	//	mytable[tid].numcolumns++;
	//	cc = mytable[tid].numcolumns-1;
	//						
	//	//printf(mytable[tt].columnname[cc]);
	//	exit(-1);

	//	strcpy(mytable[tid].columnname[cc], "bar"/*cname*/);
	//	/* write missing values into cells, if nrows >0 */
	//	if(mytable[tid].numrows>0)
	//	{
	//		for(yy=1 ; yy<=mytable[tid].numrows ; yy++)
	//		{
	//			mytable[tid].data[cc][yy]=-888.0;
	//		}
	//	}
	//	return cc;
	//	}
	//	else
	//	{
	//		printf("\n dwp_table error: can't add column 'cos table has max no. columns already \n");
	//		return -1;
	//	}	
	return 0;
}

/*
 * Adds additional rows to the table with the given name.
 * Returns the number of all rows after adding, -1 if the table cannot be found.
 */
int table_addrows(const char local_name[], int count)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("Error: Cannot add columns to table %s. Table cannot be found.\n", local_name);
		return -1;
	}

	return tchange(tid, GROW, count);
}

/*
 * Removes corws from the table wit hte given name.
 * Returns the number of rows after removing, -1 if the table cannot be found
 * or more rows were tried to be removed than are available.
 * local_name: The name of the table to modify.
 * count: the number of rows to remove.
 */
int table_removerows(const char local_name[], int count)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("Error: Cannot add rows from table %s. Table cannot be found.\n", local_name);
		return -1;
	}

	if (count > table_numrows(local_name))
	{
		printf("Error: Cannot remove rows from table %s. Table has that many rows.\n", local_name);
		return -1;
	}

	return tchange(tid, SHRINK, count);
}

/*
 * Sets the number of rows for a table.
 * If the number is larger than the current nuber of rows, additional
 * rows are added, otherwise the rows are deleted.
 */
void table_setnumrows(const char local_name[], int count)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	if (count < 0)
	{
		printf("\nError: Cannot set number of rows for table %s to negative value %d.\n", local_name, count);
		return;
	}

	if (count < 0)
	{
		printf("\nError: Cannot set number of rows for table %s. Table cannot be found.\n", local_name);
		return;
	}

	if (count > table_numcolumns(local_name))
	{
		int diff = count - table_numcolumns(local_name);
		table_addrows(local_name, diff);
		return;
	}
	else
	{
		int diff = table_numcolumns(local_name) - count;
		table_removerows(local_name, diff);
	}
}

/*
 * Writes the given value into the column with the specified name to the specified row.
 * If the value is written into a row exeeding the current size of the table, the
 * numer of rows is extended accordingly.
 *
 * local_name: The name of the table to write the value in.
 * column_name: The name of the column to write the value in.
 * row: The row to write the value in.
 */
void table_writevalue_multichain(const char local_name[], const char cname[], int row, double value)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	CHECK(cname != NULL, "Column name cannot be null");
	CHECK(row >= 0, "Row number cannot be negative");
	int tid = table_getID(local_name);
	if (tid < 0)
	{
		printf("\n Error: Cannot write value, can't find table %s\n", local_name);
		return;
	}

	int cid = table_getcolumnID(tid, cname);
	if (cid < 0)
	{
		printf("\n Error: Cannot write value, can't find column %s in table %s \n", cname, local_name);
		return;
	}

	if (row < 0)
	{
		printf("\n Error: Cannot write value, row number outside of bounds.\n");
		return;
	}

	// add rows if required
	if (row >= mytable[tid].numrows)
		table_addrows(local_name, row - mytable[tid].numrows + 1);

	mytable[tid].data[cid][row] = value;
}

/*
 * Writes the given value into the column with the specified name to the specified row.
 * If the value is written into a row exeeding the current size of the table, the
 * numer of rows is extended accordingly.
 *
 * local_name: The name of the table to write the value in.
 * column_name: The name of the column to write the value in.
 * row: The row to write the value in.
 */
void table_writevalue_multichain(const char local_name[], int column, int row, double value)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);
	if (tid < 0)
	{
		printf("\n Error: Cannot write value, can't find table %s\n", local_name);
		return;
	}

	if (column < 0 || column > table_numcolumns(local_name))
	{

		printf("\n Error: Table %s has no column with index %i. ", local_name, column);
		return;
	}

	// add rows if required
	if (row >= mytable[tid].numrows)
		table_addrows(local_name, row - mytable[tid].numrows + 1);

	mytable[tid].data[column][row] = value;
}

extern int get_currentchain();

void table_writevalue(const char local_name[], const char cname[], int row, double value)
{
	if (get_currentchain() > 0) return;
	table_writevalue_multichain(local_name, cname, row, value);
}

void table_writevalue(const char local_name[], int column, int row, double value)
{
	if (get_currentchain() > 0) return;
	table_writevalue_multichain(local_name, column, row, value);
}

/*
 * Prints the table with the given name on the screen.
 */
void table_print(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\nError: No table with local name %s.", local_name);
		return;
	}

	for (int c = 0; c < mytable[tid].numcolumns; c++)
	{
		printf("%s ", mytable[tid].columnname[c]);
	}
	printf("\n");


	for (int r = 0; r < mytable[tid].numrows; r++)
	{
		for (int c = 0; c < mytable[tid].numcolumns; c++)
		{
			printf("%lf ", mytable[tid].data[c][r]);
		}
		printf("\n");
	}
}

/*
 * Delets a table with the given name.
 * Retursn the former ID of the deleted table,
 * -1 if the table does not exist.
 */
int table_delete(const char local_name[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	// table not found
	if (tid < 0)
	{
		printf("\n Error: Cannot delete table. Table %s does not exist.", local_name);
		return -1;
	}

	return table_delete(tid);
}


/*
 * Delets a table with the given ID.
 * Retursn the former ID of the deleted table,
 * -1 if the table does not exist.
 */
int table_delete(int tid)
{
	if (tid >= 0 && tid < num_tables) {
		mytable[tid].numcolumns = 0;
		mytable[tid].numrows = 0;
		mytable[tid].columnname = NULL;
		mytable[tid].data = NULL;
		mytable[tid].name[0] = '\0';
		mytable[tid].sourcename[0] = '\0';
		return tid;
	}
	else {
		printf("\n Error: Cannot delete table. Table with ID % does not exist.", tid);
		return -1;
	}

}


/*
 * Deletes a column with the given name out of the table with the given name.
 */
int table_deletecolumn(const char local_name[], const char cname[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	CHECK(cname != NULL, "Column name cannot be null");
	int tid = table_getID(local_name);
	if (tid < 0)
	{
		printf("\n Error: Cannot delete column %s. Table with name %s does not exist.", cname, local_name);
		return -1;
	}

	int cid = table_getcolumnID(tid, cname);
	if (cid < 0)
	{
		printf("\n Error: Cannot delete column %s from table %s. Column does not exist.", cname, local_name);
		return -1;
	}

	return table_deletecolumn(cid, tid);
}

/*
 * Deletes a column with the given name out of the table with the given ID.
 */
int table_deletecolumn(int tid, const char cname[])
{
	CHECK(cname != NULL, "Column name cannot be null");
	int cid = table_getcolumnID(tid, cname);
	if (cid < 0)
	{
		printf("\n Error: Cannot delete column %s from table with ID %i. Table does not exist.", cname, tid);
		return -1;
	}

	return table_deletecolumn(cid, tid);
}

/*
 * Deletes a column with the given ID out of the table with the given name.
 */
int table_deletecolumn(const char local_name[], int cid)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);
	if (tid < 0)
	{
		printf("\n Error: Cannot delete column from table %s. Table does not exist.", local_name);
		return -1;
	}

	return table_deletecolumn(cid, tid);
}

/*
 * Deletes a column with the given ID out of the table with the given ID.
 */
int table_deletecolumn(int tid, int cid)
{
	printf("\n Error: Cannot delete column with ID %i from table with ID %i: the operation not implemented.", cid, tid);
	return -1;
	////TODO!!	if (table_get(tid) == NULL)
	//{
	//	printf("\n Error: Cannot delete column with ID %i from table with ID %i. Table does not exist.", cid, tid);
	//	return -1;
	//}

	//if (cid > table_numcolumns(tid) - 1)
	//{
	//	printf("\n Error: Cannot delete column with ID %i from table with ID %i. Column does not exist.", cid, tid);
	//	return -1;
	//}

	//// get current number of columns for this table
	//int ncol = mytable[tid].numcolumns;

	//// allocate temporary memory for column names and column data 
	//name_t *tempcolname = new name_t[ncol - 1];
	//double **tempcoldata = new double*[ncol - 1];
	//for (int i = 0; i < ncol - 1; i++)
	//{
	//	//tempcolname[i] = new char[256];
	//	tempcoldata[i] = new double[mytable[tid].numrows];
	//}

	//// copy old column names and column data into this new space, skipping the column to delete
	//int newindex = 0;
	//for (int i = 0; i < ncol; i++)
	//{
	//	if (i != cid)
	//	{
	//		strcpy(tempcolname[newindex], mytable[tid].columnname[i]);

	//		for (int j = 0; j < mytable[tid].numrows; j++)
	//			tempcoldata[newindex][j] = mytable[tid].data[i][j];

	//		newindex++;
	//	}
	//}

	//// delete old columns
	//tfreecolumns(tid);

	//// assign new column names and data
	//mytable[tid].columnname = tempcolname;
	//mytable[tid].data = tempcoldata;
	//mytable[tid].numcolumns--;
	//return 0;
}

/*
 * Returns the maximum value (!= missing value)
 * from the column with the given name within the
 * table with the given name.
 * If the table or column does not exist, the missing value is returned.
 *
 * local_name: The name of the table to look at.
 * cname: The column name to look at.
 */
double table_getcolumnmax(const char local_name[], const char cname[])
{
	CHECK(cname != NULL, "Column name cannot be null");
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);
	if (tid < 0)
	{
		printf("\n Error: Cannot get max value from column %s within table %s. Table does not exist.", cname, local_name);
		return MISSING_VALUE;
	}

	int cid = table_getcolumnID(tid, cname);
	if (cid > table_numcolumns(tid) - 1)
	{
		printf("\n Error: Cannot get max value from column %s within table %s. Column does not exist.", cname, local_name);
		return MISSING_VALUE;
	}

	return table_getcolumnmax(local_name, cid);
}

/*
 * Returns the maximum value (!= missing value)
 * from the column with the given ID within the
 * table with the given name.
 * If the table or column does not exist, the missing value is returned.
 *
 * local_name: The name of the table to look at.
 * cid: The ID of the column to look at.
 */
double table_getcolumnmax(const char local_name[], int cid)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\n Error: Cannot get max value from column with ID %i from table %s. Table does not exist.", cid, local_name);
		return MISSING_VALUE;
	}

	if (cid < 0 || cid > table_numcolumns(tid) - 1)
	{
		printf("\n Error: Cannot get max value from column with ID %i from table %s. Column does not exist.", cid, local_name);
		return MISSING_VALUE;
	}

	if (table_numrows(tid) < 1)
	{
		printf("\n Error: Cannot get max value column with ID %i from table %s. Column does not exist.", cid, local_name);
		return MISSING_VALUE;
	}

	double* data = mytable[tid].data[cid];
	double retval = data[0];

	int nrows = table_numrows(tid);
	if (nrows == 1)
		return retval;

	for (int i = 1; i < nrows; i++)
	{
		double curval = data[i];

		// do we get something better than MISSING_VALUE
		if (retval == MISSING_VALUE && curval != MISSING_VALUE)
			retval = curval;

		// if its not the MISSING_VALUE we take the greater value
		if (curval != MISSING_VALUE && curval > retval)
			retval = curval;
	}

	return retval;
}

/*
 * Returns the minimum value (!= missing value)
 * from the column with the given name within the
 * table with the given name.
 * If the table or column does not exist, the missing value is returned.
 *
 * local_name: The name of the table to look at.
 * cname: The column name to look at.
 */
double table_getcolumnmin(const char local_name[], const char cname[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	CHECK(cname != NULL, "Column name cannot be null");
	int tid = table_getID(local_name);
	if (tid < 0)
	{
		printf("\n Error: Cannot get min value from column %s from table %s. Table does not exist.", cname, local_name);
		return MISSING_VALUE;
	}

	int cid = table_getcolumnID(tid, cname);
	if (cid > table_numcolumns(tid) - 1)
	{
		printf("\n Error: Cannot get min value from column %s from table %s. Column does not exist.", cname, local_name);
		return MISSING_VALUE;
	}

	return table_getcolumnmin(local_name, cid);
}


/*
 * Returns the maximum value (!= missing value)
 * from the column with the given ID within the
 * table with the given name.
 * If the table or column does not exist, the missing value is returned.
 *
 * local_name: The name of the table to look at.
 * cid: The ID of the column to look at.
 */
double table_getcolumnmin(const char local_name[], int cid)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\n Error: Cannot get min value from column with ID %i from table %s. Table does not exist.", cid, local_name);
		return MISSING_VALUE;
	}

	if (cid < 0 || cid > table_numcolumns(tid) - 1)
	{
		printf("\n Error: Cannot get min value from column with ID %i from table %s. Column does not exist.", cid, local_name);
		return MISSING_VALUE;
	}

	if (table_numrows(tid) < 1)
	{
		printf("\n Error: Cannot get min value frmo column with ID %i from table %s. Column does not exist.", cid, local_name);
		return MISSING_VALUE;
	}

	double* data = mytable[tid].data[cid];
	double retval = data[0];

	int nrows = table_numrows(tid);
	if (nrows == 1)
		return retval;

	for (int i = 1; i < nrows; i++)
	{
		double curval = data[i];

		// do we get something better than MISSING_VALUE
		if (retval == MISSING_VALUE && curval != MISSING_VALUE)
			retval = curval;

		// if its not the MISSING_VALUE we take the smaller value
		if (curval != MISSING_VALUE && curval < retval)
			retval = curval;
	}

	return retval;
}

/*
 * Returns the number of values within the given column.
 *
 * cname: The name of the column to count the values.
 * local_name The name of the table the column is within.
 */
int table_getvaluecount(const char local_name[], const char cname[])
{
	CHECK(local_name != NULL, "Table name cannot be null");
	CHECK(cname != NULL, "Column name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\n Error: Cannot get  column %s from table %s. Table does not exist.", cname, local_name);
		return (int)MISSING_VALUE;
	}

	int cid = table_getcolumnID(tid, cname);
	if (cid > table_numcolumns(tid) - 1)
	{
		printf("\n Error: Cannot delete column %s from table %s. Column does not exist.", cname, local_name);
		return (int)MISSING_VALUE;
	}

	return table_getvaluecount(local_name, cid);
}

/*
 * Returns the number of values within the given column.
 *
 * local_name The name of the table the column is within.
 * cid: The id of the column to count the values.
 */
int table_getvaluecount(const char local_name[], int cid)
{
	CHECK(local_name != NULL, "Table name cannot be null");
	int tid = table_getID(local_name);

	if (tid < 0)
	{
		printf("\n Error: Cannot get value vount from column with ID %i from table %s. Table does not exist.", cid, local_name);
		return (int)MISSING_VALUE;
	}

	if (cid < 0 || cid > table_numcolumns(tid) - 1)
	{
		printf("\n Error: Cannot get value count from column with ID %i from table %s. Column does not exist.", cid, local_name);
		return (int)MISSING_VALUE;
	}


	int count = 0;

	for (int i = 0; i < table_numrows(tid); i++)
	{
		if (table_getvalue(local_name, cid, i) != MISSING_VALUE)
			count++;
	}

	return count;
}
