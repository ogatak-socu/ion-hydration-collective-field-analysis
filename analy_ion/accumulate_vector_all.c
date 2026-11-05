#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>
#include <limits.h>

#include "cal_vector.h"

static int max_layer;

static int rd_poly_a_values(POLY_A *,SW *,int );
static void rd_grid_values(GRID_T *,char *);
static void wt_grid_values_angle_all(POLY_A *,SW *);
static void wt_grid_values_angle(GRID_T *,char *,SW *);
static double cal_point_vector_angle(XYZ_T *,XYZ_T *);
static void wt_grid_values_angle_radical(double ,double ,int ,FILE *);

void accumulate_vector_all(SW *sw)
{
	char fname[MAXLINE];
        int k,max_g,n;
	POLY_A *poly_a;

        max_layer = (int )floor((double )sw->rn*(1.0/sw->delta_d)) + 1;

	poly_a = get_grid_system(max_layer,sw);

	fprintf(stderr,"Number of vectors: %d\n",sw->number_of_vectors);
        for(k = 0; k < sw->number_of_vectors; k++)
		rd_poly_a_values(poly_a,sw,k+1);

	wt_grid_values_angle_all(poly_a,sw);
}


static int rd_poly_a_values(POLY_A *poly_a,SW *sw,int no)
{
	char fname[MAXLINE];

	sprintf(fname,"%s.%d.out",sw->dip_vector_prefix,no);
	rd_grid_values(poly_a->grid_d,fname);

	sprintf(fname,"%s.%d.out",sw->frc_vector_prefix,no);
	rd_grid_values(poly_a->grid_f,fname);

	sprintf(fname,"%s.%d.out",sw->bck_vector_prefix,no);
	rd_grid_values(poly_a->grid_b,fname);

	sprintf(fname,"%s.%d.out",sw->fwd_vector_prefix,no);
	rd_grid_values(poly_a->grid_m,fname);

	return 1;
}

static void rd_grid_values(GRID_T *grid_t,char *fname)
{
	int n,layer,idx;
	FILE *fp;
	WGRID_T wgrid_t;

	fprintf(stderr,"Reading %s\n",fname);
	if((fp = rd_binary_file_open(fname)) != NULL){
		while((n = fread(&wgrid_t,sizeof(WGRID_T),1,fp)) == 1){
			layer = wgrid_t.layer;
			if(layer >= max_layer)
				continue;
			idx = wgrid_t.idx;
			((grid_t+layer)->dxyz_t+idx)->x += wgrid_t.x;
			((grid_t+layer)->dxyz_t+idx)->y += wgrid_t.y;
			((grid_t+layer)->dxyz_t+idx)->z += wgrid_t.z;
			if((*((grid_t+layer)->n+idx) += wgrid_t.n) == 0){
				fprintf(stderr,"Stack overflow in grid_t:n\n");
				exit(0);
			}
			*((grid_t+layer)->weight+idx) += wgrid_t.weight;
		}
		rd_file_close(fname,fp);
	}
}

static void wt_grid_values_angle_all(POLY_A *poly_a,SW *sw)
{
	char fname[MAXLINE];

	sprintf(fname,"%s.data.out",sw->dip_vector_prefix);
	wt_grid_values_angle(poly_a->grid_d,fname,sw);

	sprintf(fname,"%s.data.out",sw->frc_vector_prefix);
	wt_grid_values_angle(poly_a->grid_f,fname,sw);

	sprintf(fname,"%s.data.out",sw->bck_vector_prefix);
	wt_grid_values_angle(poly_a->grid_b,fname,sw);

	sprintf(fname,"%s.data.out",sw->fwd_vector_prefix);
	wt_grid_values_angle(poly_a->grid_m,fname,sw);
}

static void wt_grid_values_angle(GRID_T *grid_t,char *fname,SW *sw)
{
        int k,i,n_n;
        double ang,t_ang,len,t_len;
        XYZ_T w;
        FILE *fp;

        fprintf(stderr,"Writing %s\n",fname);
        fp = wt_file_open(fname);
        for(k = 0; k < max_layer; k++){
                fprintf(fp,"%6.2f ",(grid_t+k)->r);
                n_n = 0;
                t_ang = t_len = 0.0;
                for(i = 0; i < (grid_t+k)->npoint; i++){
                        if(*((grid_t+k)->n+i) != 0UL && fabs(distance0_xyz((grid_t+k)->p+i)) >= 1.0E-06){
                                divid_scalar2vector_xyz(&w,(grid_t+k)->dxyz_t+i,*((grid_t+k)->weight+i));
                                ang = cal_point_vector_angle((grid_t+k)->p+i,&w);
                                len = distance0_xyz(&w);
                                t_ang += ang;
                                t_len += len;
                                n_n++;
                        }
                }
		if(n_n != 0){
			t_ang /= (double )n_n;
			t_len /= (double )n_n;
		}
		wt_grid_values_angle_radical(t_ang,t_len,n_n,fp);
        }
        wt_file_close(fname,fp);
}

static double cal_point_vector_angle(XYZ_T *p0,XYZ_T *p1)
{
	double ang;
	XYZ_T v0,v1;

	if(distance0_xyz(p1) <= 1.0E-06)
		return -1.0;
	unit_vector_xyz(&v0,p0);
	scale_vector_xyz(&v0,-1.0);
	unit_vector_xyz(&v1,p1);
	ang = vec_angle0_xyz(&v0,&v1)*180.0/M_PI;

	return ang;
}


static void wt_grid_values_angle_radical(double t_ang,double t_len,int n,FILE *fp)
{
        if(n != 0)
		fprintf(fp,"%9.3f %9.5f %10d\n",t_ang,t_len,n);
        else
		fprintf(fp,"%9.3f %9.5f %10d\n",-1.0,-1.0,-1.0,-1.0,0);
}

