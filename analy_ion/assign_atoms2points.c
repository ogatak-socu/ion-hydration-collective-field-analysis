#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static double cal_grid_distance(int ,int ,int ,double );
static double distance_xyz2(XYZ_T *,XYZ_T *);
static double cal_distance_weight2(double ,double );

RANGE_P *get_init_range_p(int max_layer)
{
        int i;
        RANGE_P *range_p;

        range_p = (RANGE_P *)calloc(max_layer+1,sizeof(RANGE_P ));
        for(i = 0; i <= max_layer; i++){
                (range_p+i)->weight = -1.0;
                (range_p+i)->layer = -1;
                (range_p+i)->idx = -1;
        }
        return range_p;
}

void set_init_range_p(RANGE_P *range_p,int max_layer)
{
        int i;

        for(i = 0; i <= max_layer; i++){
                (range_p+i)->weight = -1.0;
                (range_p+i)->layer = -1;
                (range_p+i)->idx = -1;
        }
}

void get_range_point_polygon_grid_normal(RANGE_P *range_p,GRID_T *grid_t,XYZ_T *xyz_t,POLY_P *poly_p,XYZ_I *box_size,double grid_size,int delta_idx,double rc2,double alpha)
{
	int gi,gj,gk,i,j,k,ii,n,ip,il,idx_r,nn,dx,dy,dz;
	double d;
	XYZ_T g;

	nn = 0;
	zero_vector_xyz(&g);
	g_xyz2ijk_mass(xyz_t,&gi,&gj,&gk,&g,grid_size);
	for(i = gi - delta_idx; i <= gi + delta_idx; i++){
		dx = abs(i - gi);
		for(j = gj - delta_idx; j <= gj + delta_idx; j++){
			dy = abs(j - gj);
			for(k = gk - delta_idx; k <= gk + delta_idx; k++){
				dz = abs(k - gk);
				if(g_box_range(i,j,k,box_size) == -1)
					continue;
				if(cal_grid_distance(dx,dy,dz,grid_size) > rc2)
					continue;

				ii = g_ijk2idx(box_size,i,j,k);
				for(n = 0; n < (poly_p+ii)->n; n++){
					il = *((poly_p+ii)->layer+n);
					ip = *((poly_p+ii)->idx+n);
					d = distance_xyz2((grid_t+il)->p+ip,xyz_t);
					if(d < rc2){
						(range_p+nn)->weight = cal_distance_weight2(d,alpha);
						(range_p+nn)->layer = il;
						(range_p+nn)->idx = ip;

						nn++;
					}
				}
			}
		}
	}
}

static double cal_grid_distance(int i,int j,int k,double grid_size)
{
	double dx,dy,dz;

	if(i == 0)
		dx = ((double )i + 0.5)*grid_size;
	else
		dx = ((double )(i-1) + 0.5)*grid_size;
	if(j == 0)
		dy = ((double )j + 0.5)*grid_size;
	else
		dy = ((double )(j-1) + 0.5)*grid_size;
	if(j == 0)
		dz = ((double )k + 0.5)*grid_size;
	else
		dz = ((double )(k-1) + 0.5)*grid_size;
	return dx*dx + dy*dy + dz*dz;
}

static double distance_xyz2(XYZ_T *xyz_t0,XYZ_T *xyz_t1)
{
	double dx,dy,dz;

	dx = (xyz_t0->x - xyz_t1->x);
	dy = (xyz_t0->y - xyz_t1->y);
	dz = (xyz_t0->z - xyz_t1->z);
	return dx*dx + dy*dy + dz*dz;
}

static double cal_distance_weight2(double d,double sigma)
{
	if(sigma < 0.0)
		return 1.0;
	else
		return exp(-d*sigma);
}
