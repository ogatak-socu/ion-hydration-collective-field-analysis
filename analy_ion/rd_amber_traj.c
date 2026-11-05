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

static void get_start_coord(char *,int ,int ,FILE *);
static double *get_coord_value(char *,int ,FILE *);
static XYZ_T *trans_coord_value(double *,int );
static XYZ_T *get_box_value(char *,FILE *,int *);

TRAJ *rd_amber_traj(char *fname,int natom,int sp,int num_s,int *t_cnt)
{
	char buf[MAXLINE];
	int i,flg;
	double *xyz;
	TRAJ *traj;
	FILE *fp;

	traj = (TRAJ *)calloc(num_s,sizeof(TRAJ ));
	xyz = (double *)calloc(natom*3,sizeof(double ));
	if((fp = rd_file_open(fname)) == NULL){
		(void)fprintf(stderr,"Cannot open mdcrd file %s\n",fname);
		exit(0);
	}
	get_start_coord(buf,natom,sp,fp);
	for(i = *t_cnt = 0; i < num_s; i++){
		xyz = get_coord_value(buf,natom,fp);
		(traj+i)->mdcrd = trans_coord_value(xyz,natom);
		free(xyz);
		(traj+i)->box = get_box_value(buf,fp,&flg);
		(*t_cnt)++;
		if(flg == -1)
			break;
	}
	rd_file_close(fname,fp);
	return traj;
}

TRAJ *rd_amber_traj_other(char *fname,int natom,int sp,int num_s,int *t_cnt)
{
	char buf[MAXLINE];
	int i,flg;
	double *xyz;
	TRAJ *traj;
	FILE *fp;

	traj = (TRAJ *)calloc(num_s,sizeof(TRAJ ));
	xyz = (double *)calloc(natom*3,sizeof(double ));
	if((fp = rd_file_open(fname)) == NULL){
		(void)fprintf(stderr,"Cannot open mdcrd file %s\n",fname);
		exit(0);
	}
	get_start_coord(buf,natom,sp,fp);
	for(i = *t_cnt = 0; i < num_s; i++){
		xyz = get_coord_value(buf,natom,fp);
		(traj+i)->mdcrd = trans_coord_value(xyz,natom);
		free(xyz);
		(*t_cnt)++;
		if(fgets(buf,MAXLINE,fp) == NULL)
			break;
	}
	rd_file_close(fname,fp);
	return traj;
}



void free_traj_area(TRAJ *traj,int t_cnt)
{
	int i;

	for(i = 0; i < t_cnt; i++){
		free((traj+i)->mdcrd);
		free((traj+i)->box);
	}
	free(traj);
}

void free_traj_area_other(TRAJ *traj,int t_cnt)
{
	int i;

	for(i = 0; i < t_cnt; i++)
		free((traj+i)->mdcrd);
	free(traj);
}

TRAJ *get_init_traj_area(int natom)
{
	TRAJ *traj;

	traj = (TRAJ *)calloc(1,sizeof(TRAJ ));
	traj->mdcrd = get_init_xyz_t(natom);
	traj->box = get_init_xyz_t(1);
	return traj;
}

static void get_start_coord(char *buf,int natom,int sp,FILE *fp)
{
	int m_cnt,k;

	m_cnt = 0;
	fgets(buf,MAXLINE,fp);
	while(sp-1 != m_cnt){
		k = 0;
		while(k < natom*3 && fgets(buf,MAXLINE,fp) != NULL){
			k += 10;
		}
		if((natom*3) % 10 != 0)
			fgets(buf,MAXLINE,fp);
		m_cnt++;
	}
	fgets(buf,MAXLINE,fp);
}

static double *get_coord_value(char *buf,int natom,FILE *fp)
{
	int i,k;
	double *xyz;

	xyz = (double *)calloc(natom*3,sizeof(double ));
	for(i = k = 0; i < natom*3; i++){
		*(xyz+i) = atof(buf+k*8);
		k++;
		if(i != natom*3-1 && (i+1) % 10 == 0){
			fgets(buf,MAXLINE,fp);
			k = 0;
		}
	}
	return xyz;
}

static XYZ_T *trans_coord_value(double *xyz,int natom)
{
	int i;
	XYZ_T *mdcrd;

	mdcrd = (XYZ_T *)calloc(natom,sizeof(XYZ_T ));
	for(i = 0; i < natom; i++){
		(mdcrd+i)->x = *(xyz+i*3);
		(mdcrd+i)->y = *(xyz+i*3+1);
		(mdcrd+i)->z = *(xyz+i*3+2);
	}
	return mdcrd;
}

static XYZ_T *get_box_value(char *buf,FILE *fp,int *flg)
{
	XYZ_T *box;

	box = (XYZ_T *)calloc(1,sizeof(XYZ_T ));
	*flg = 1;
	if(fgets(buf,MAXLINE,fp) == NULL)
		*flg = -1;
	if(strlen(buf) >= 30)
		return box;
	box->x = atof(buf);
	box->y = atof(buf+8);
	box->z = atof(buf+16);
	if(fgets(buf,MAXLINE,fp) == NULL)
		*flg = -1;
	return box;
}

void cp_traj_mdcrd(XYZ_T *w_mdcrd,XYZ_T *mdcrd,int natom)
{
        int i;

        for(i = 0; i < natom; i++)
                cp_vector_xyz(w_mdcrd+i,mdcrd+i);
}

