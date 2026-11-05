#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cal_vector.h"

#define GZCAT_COM "gunzip -c "
#define BUFFER_SIZ   1024

FILE *rd_file_open(char *fname)
{
	char buf[BUFFER_SIZ];
	FILE *fp;

	if((fp = fopen(fname,"r")) == NULL){
		(void)fprintf(stderr,"Cannot open file %s\n",fname);
		return NULL;
	}
	else
		fclose(fp);
        if(strstr(fname,".gz") == NULL && strstr(fname,".Z") == NULL){
		if((fp = fopen(fname,"r")) == NULL)
			return NULL;
        }
        else{
                sprintf(buf,"%s %s",GZCAT_COM,fname);
		if((fp = popen(buf,"r")) == NULL)
			return NULL;
        }
        return fp;
}

FILE *rd_binary_file_open(char *fname)
{
	FILE *fp;

	if((fp = fopen(fname,"r")) == NULL){
		(void)fprintf(stderr,"Cannot open file %s\n",fname);
		return NULL;
	}
        return fp;
}


FILE *wt_file_open(char *fname)
{
	FILE *fp;

	if((fp = fopen(fname,"w")) == NULL){
		fprintf(stderr,"Cannot write file %s\n",fname);
		exit(0);
	}
	return fp;
}

FILE *wt_binary_file_open(char *fname)
{
	FILE *fp;

	if((fp = fopen(fname,"wb")) == NULL){
		fprintf(stderr,"Cannot write file %s\n",fname);
		exit(0);
	}
	return fp;
}

FILE *append_file_open(char *fname)
{
	FILE *fp;

	if((fp = fopen(fname,"a")) == NULL){
		(void)fprintf(stderr,"Cannot open file %s\n",fname);
		return NULL;
	}
	return fp;
}

void rd_file_close(char *fname,FILE *fp)
{
        if(strstr(fname,".gz") == NULL && strstr(fname,".Z") == NULL)
                pclose(fp);
        else
                fclose(fp);
}

void wt_file_close(char *fname,FILE *fp)
{
	fclose(fp);
}
