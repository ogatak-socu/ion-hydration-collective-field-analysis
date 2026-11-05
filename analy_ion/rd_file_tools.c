#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cal_vector.h"

#define MAX_RD_LINE 1000000

char **rd_text_file(char *fname,int *l_cnt)
{
	char **buf_t,*buf;
	int i,len;
	FILE *fp;

	if((*l_cnt = get_text_line_cnt(fname)) == 0)
		return NULL;
	fp = rd_file_open(fname);
	buf_t = (char **)calloc(*l_cnt,sizeof(char *));
	buf = (char *)calloc(MAX_RD_LINE,sizeof(char ));
	for(i = 0; i < *l_cnt; i++){
		fgets(buf,MAX_RD_LINE,fp);
		len = strlen(buf);
		*(buf_t+i) = (char *)calloc(len+1,sizeof(char ));
		strcpy(*(buf_t+i),buf);
	}
	rd_file_close(fname,fp);
	free(buf);
	return buf_t;
}

char **rd_text_file_without_cr(char *fname,int *l_cnt)
{
	char **buf_t,*buf;
	int i,len;
	FILE *fp;

	if((*l_cnt = get_text_line_cnt(fname)) == 0)
		return NULL;
	fp = rd_file_open(fname);
	buf_t = (char **)calloc(*l_cnt,sizeof(char *));
	buf = (char *)calloc(MAX_RD_LINE,sizeof(char ));
	for(i = 0; i < *l_cnt; i++){
		fgets(buf,MAX_RD_LINE,fp);
		len = strlen(buf);
		*(buf_t+i) = (char *)calloc(len+1,sizeof(char ));
		*(buf+len-1) = '\0';
		strcpy(*(buf_t+i),buf);
	}
	rd_file_close(fname,fp);
	free(buf);
	return buf_t;
}

int get_text_line_cnt(char *fname)
{
	char *buf;
	int n;
	FILE *fp;

	if((fp = rd_file_open(fname)) == NULL)
		return 0;
	n = 0;
	buf = (char *)calloc(MAX_RD_LINE,sizeof(char ));
	while(fgets(buf,MAX_RD_LINE,fp) != NULL)
		n++;
	rd_file_close(fname,fp);
	free(buf);
	return n;
}

char **get_text_area_fix(int l_cnt)
{
        char **buf_t;
        int i;

        if((buf_t = (char **)calloc(l_cnt,sizeof(char *))) == NULL){
                fprintf(stderr,"Not enough memory\n");
                exit(0);
        }
        for(i = 0; i < l_cnt; i++){
                if((*(buf_t+i) = (char *)calloc(MAXLINE,sizeof(char ))) == NULL){
                        fprintf(stderr,"Not enough memory\n");
                        exit(0);
                }
        }
        return buf_t;
}

char **get_text_area_flexible(int l_cnt,int c_cnt)
{
        char **buf_t;
        int i;

        if((buf_t = (char **)calloc(l_cnt,sizeof(char *))) == NULL){
                fprintf(stderr,"Not enough memory\n");
                exit(0);
        }
        for(i = 0; i < l_cnt; i++){
                if((*(buf_t+i) = (char *)calloc(c_cnt,sizeof(char ))) == NULL){
                        fprintf(stderr,"Not enough memory\n");
                        exit(0);
                }
        }
        return buf_t;
}


void free_text_area(char **buf_t,int l_cnt)
{
        int i;

        for(i = 0; i < l_cnt; i++)
                free(*(buf_t+i));
        free(buf_t);
}
