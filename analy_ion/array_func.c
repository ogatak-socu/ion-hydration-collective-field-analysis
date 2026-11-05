#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cal_vector.h"

char *get_init_char_array(int len,char init_c)
{
	char *buf;
	int i;

	buf = (char *)calloc(len+1,sizeof(char ));
	for(i = 0; i < len; i++)
		*(buf+i) = init_c;
	*(buf+i) = '\0';
	return buf;
}

void set_init_char_array(char *buf,int len,char init_c)
{
	int i;

	for(i = 0; i < len; i++)
		*(buf+i) = init_c;
	*(buf+i) = '\0';
}

char **get_char_array_2d(int col,int line)
{
	char **buf;
	int i;

	buf = (char **)calloc(col,sizeof(char *));
	for(i = 0; i < col; i++)
		*(buf+i) = (char *)calloc(line+1,sizeof(char ));
	return buf;
}

char **get_init_char_array_2d(int col,int line,char init_c)
{
	char **buf;
	int i,j;

	buf = (char **)calloc(col,sizeof(char *));
	for(i = 0; i < col; i++){
		*(buf+i) = (char *)calloc(line+1,sizeof(char ));
		for(j = 0; j < line; j++)
			*(*(buf+i)+j) = init_c;
		*(*(buf+i)+j) = '\0';
	}
	return buf;
}

void set_init_char_array_2d(char **buf,int col,int line,char init_c)
{
	int i,j;

	for(i = 0; i < col; i++){
		for(j = 0; j < line; j++)
			*(*(buf+i)+j) = init_c;
		*(*(buf+i)+j) = '\0';
	}
}

void cp_char_array(char *a,char *b,int n)
{
	int i;

	for(i = 0; i < n; i++)
		*(a+i) = *(b+i);
}

void free_char_array_2d(char **buf,int col)
{
	int i;

	for(i = 0; i < col; i++)
		free(*(buf+i));
	free(buf);
}

int *get_init_int_array(int len,int init_c)
{
	int i,*v;

	v = (int *)calloc(len,sizeof(int ));
	for(i = 0; i < len; i++)
		*(v+i) = init_c;
	return v;
}

void set_init_int_array(int *v,int len,int init_c)
{
	int i;

	for(i = 0; i < len; i++)
		*(v+i) = init_c;
}

int **get_int_array_2d(int col,int line)
{
	int i,**v;

	v = (int **)calloc(col,sizeof(int *));
	for(i = 0; i < col; i++)
		*(v+i) = (int *)calloc(line,sizeof(int ));
	return v;
}

int **get_init_int_array_2d(int col,int line,int init_c)
{
	int i,j,**v;

	v = (int **)calloc(col,sizeof(int *));
	for(i = 0; i < col; i++){
		*(v+i) = (int *)calloc(line,sizeof(int ));
		for(j = 0; j < line; j++)
			*(*(v+i)+j) = init_c;
	}
	return v;
}

void set_init_int_array_2d(int **v,int col,int line,int init_c)
{
	int i,j;

	for(i = 0; i < col; i++){
		for(j = 0; j < line; j++)
			*(*(v+i)+j) = init_c;
	}
}

void cp_int_array(int *a,int *b,int n)
{
	int i;

	for(i = 0; i < n; i++)
		*(a+i) = *(b+i);
}

void free_int_array_2d(int **v,int col)
{
	int i;

	for(i = 0; i < col; i++)
		free(*(v+i));
	free(v);
}

unsigned long int *get_init_unsigned_long_int_array(int len,unsigned long int init_c)
{
	int i;
	unsigned long int *v;

	v = (unsigned long int *)calloc(len,sizeof(unsigned long int ));
	for(i = 0; i < len; i++)
		*(v+i) = init_c;
	return v;
}

void set_init_unsigned_long_int_array(unsigned long int *v,int len,unsigned long int init_c)
{
	int i;

	for(i = 0; i < len; i++)
		*(v+i) = init_c;
}

unsigned long int **get_unsigned_long_int_array_2d(int col,int line)
{
	int i;
	unsigned long int **v;

	v = (unsigned long int **)calloc(col,sizeof(unsigned long int *));
	for(i = 0; i < col; i++)
		*(v+i) = (unsigned long int *)calloc(line,sizeof(unsigned long int ));
	return v;
}

unsigned long int **get_init_unsigned_long_int_array_2d(int col,int line,unsigned long int init_c)
{
	int i,j;
	unsigned long int **v;

	v = (unsigned long int **)calloc(col,sizeof(unsigned long int *));
	for(i = 0; i < col; i++){
		*(v+i) = (unsigned long int *)calloc(line,sizeof(unsigned long int ));
		for(j = 0; j < line; j++)
			*(*(v+i)+j) = init_c;
	}
	return v;
}

void set_init_unsigned_long_int_array_2d(unsigned long int **v,int col,int line,unsigned long int init_c)
{
	int i,j;

	for(i = 0; i < col; i++){
		for(j = 0; j < line; j++)
			*(*(v+i)+j) = init_c;
	}
}

void cp_unsigned_long_int_array(unsigned long int *a,unsigned long int *b,int n)
{
	int i;

	for(i = 0; i < n; i++)
		*(a+i) = *(b+i);
}

void free_unsigned_long_int_array_2d(unsigned long int **v,int col)
{
	int i;

	for(i = 0; i < col; i++)
		free(*(v+i));
	free(v);
}


double *get_init_double_array(int len,double init_c)
{
	int i;
	double *v;

	v = (double *)calloc(len,sizeof(double ));
	for(i = 0; i < len; i++)
		*(v+i) = init_c;
	return v;
}

void set_init_double_array(double *v,int len,double init_c)
{
	int i;

	for(i = 0; i < len; i++)
		*(v+i) = init_c;
}

double **get_double_array_2d(int col,int line)
{
	int i;
	double **v;

	v = (double **)calloc(col,sizeof(double *));
	for(i = 0; i < col; i++)
		*(v+i) = (double *)calloc(line,sizeof(double ));
	return v;
}

double **get_init_double_array_2d(int col,int line,double init_c)
{
	int i,j;
	double **v;

	v = (double **)calloc(col,sizeof(double *));
	for(i = 0; i < col; i++){
		*(v+i) = (double *)calloc(line,sizeof(double ));
		for(j = 0; j < line; j++)
			*(*(v+i)+j) = init_c;
	}
	return v;
}

void set_init_double_array_2d(double **v,int col,int line,double init_c)
{
	int i,j;

	for(i = 0; i < col; i++){
		for(j = 0; j < line; j++)
			*(*(v+i)+j) = init_c;
	}
}

void free_double_array_2d(double **v,int col)
{
	int i;

	for(i = 0; i < col; i++)
		free(*(v+i));
	free(v);
}

void cp_double_array(double *a,double *b,int n)
{
	int i;

	for(i = 0; i < n; i++)
		*(a+i) = *(b+i);
}
