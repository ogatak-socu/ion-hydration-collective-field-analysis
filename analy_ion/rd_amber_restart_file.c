/*                                              */
/*      Copyright : Koji Ogata                  */
/*                                              */
/*      Date      : 2005.08.21                  */
/*                                              */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cal_vector.h"

static char *get_title(char *);
static int get_line_cnt(FILE *);

RST *rd_amber_restart_file(char *fname)
{
	char buf[MAXLINE],buf2[MAXLINE],*flg;
	int i,line_n,k;
	FILE *fp;
	RST *rst;

	rst = (RST *)calloc(1,sizeof(RST ));
	line_n = get_text_line_cnt(fname);
	if((fp = fopen(fname,"r")) == NULL){
		(void)fprintf(stderr,"Cannot open file %s\n",fname);
		exit(0);
	}
	fgets(buf,MAXLINE,fp);
	rst->title = get_title(buf);
	fgets(buf,MAXLINE,fp);
	k = 2;
	rst->natom = atoi(buf);
	if(strlen(buf) >= 5)
		rst->t = atof(buf+5);
	else
		rst->t = 0.0;
	rst->natom_s = get_title(buf);
	rst->r = (XYZ_T *)calloc(rst->natom,sizeof(XYZ_T ));
	rst->v = (XYZ_T *)calloc(rst->natom,sizeof(XYZ_T ));
	for(i = 0; i < rst->natom; i++){
                if(i % 2 == 0){
                        fgets(buf,MAXLINE,fp);
                         (rst->r+i)->x = atof(buf);
                         (rst->r+i)->y = atof(buf+12);
                         (rst->r+i)->z = atof(buf+24);
                }
                else{
                         (rst->r+i)->x = atof(buf+36);
                         (rst->r+i)->y = atof(buf+48);
                         (rst->r+i)->z = atof(buf+60);
                }
	}
	if(line_n  != k){
		memset(buf2,'\0',sizeof(buf2));
		for(i = 0; i < rst->natom; i++){
			if(i % 2 == 0){
				if((flg = fgets(buf,MAXLINE,fp)) == NULL)
					break;
				strcpy(buf2,buf);
				(rst->v+i)->x = atof(buf);
				(rst->v+i)->y = atof(buf+12);
				(rst->v+i)->z = atof(buf+24);
			}
			else{
				(rst->v+i)->x = atof(buf+36);
				(rst->v+i)->y = atof(buf+48);
				(rst->v+i)->z = atof(buf+60);
			}
		}
	}
	if(i != rst->natom){
		if(flg == NULL && strlen(buf2) != 0){
			rst->box.x = atof(buf2);
			rst->box.y = atof(buf2+12);
			rst->box.z = atof(buf2+24);
			rst->box_a.x = atof(buf2+36);
			rst->box_a.y = atof(buf2+48);
			rst->box_a.z = atof(buf2+60);
		}
		else{
			rst->box.x = -1.0;
			rst->box.y = -1.0;
			rst->box.z = -1.0;
			rst->box_a.x = -1.0;
			rst->box_a.y = -1.0;
			rst->box_a.z = -1.0;
		}
	}
	else{
		if(fgets(buf,MAXLINE,fp) != NULL){
			rst->box.x = atof(buf);
			rst->box.y = atof(buf+12);
			rst->box.z = atof(buf+24);
			rst->box_a.x = atof(buf+36);
			rst->box_a.y = atof(buf+48);
			rst->box_a.z = atof(buf+60);
		}
		else{
			rst->box.x = -1.0;
			rst->box.y = -1.0;
			rst->box.z = -1.0;
			rst->box_a.x = -1.0;
			rst->box_a.y = -1.0;
			rst->box_a.z = -1.0;
		}
	}
	fclose(fp);
	return rst;
}

void free_amber_rst(RST *rst)
{
	free(rst->natom_s);
	free(rst->r);
	free(rst->v);
	free(rst->title);
	free(rst);
}

static char *get_title(char *buf)
{
	char *title;
	int len;

	len = strlen(buf);
	title = (char *)calloc(len+1,sizeof(char ));
	strcpy(title,buf);
	*(title+len) = '\0';
	return title;
}

static int get_line_cnt(FILE *fp)
{
	char buf[MAXLINE];
	int n;

	n = 0;
	while(fgets(buf,MAXLINE,fp) != NULL)
		n++;
	return n;
}
