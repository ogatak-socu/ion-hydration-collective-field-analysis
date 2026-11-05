#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static POLY_P *get_init_poly_p(int );
static int get_init_number_layer_index(GRID_T *,POLY_P *,XYZ_I *,double ,int ,int );
static void alloc_layer_index_area(POLY_P *,int );
static int sign(double );

POLY_P *assign_spherical_system2cube_system(GRID_T *grid_t,XYZ_I *box_size,double grid_size,int max_g,int max_layer,int *max_grid_index)
{
        int i,k,n,idx,nn;
        XYZ_T g;
        POLY_P *poly_p;

	fprintf(stderr,"Assign all the points in spherical system to the points in cubic system\n");
        poly_p = get_init_poly_p(max_g);
        *max_grid_index = get_init_number_layer_index(grid_t,poly_p,box_size,grid_size,max_g,max_layer);
        alloc_layer_index_area(poly_p,max_g);
        zero_vector_xyz(&g);
        for(k = 0; k <= max_layer; k++){
                for(i = 0; i < (grid_t+k)->npoint; i++){
                        if((idx = g_xyz2idx((grid_t+k)->p+i,box_size,&g,grid_size)) == -1)
                                continue;
                        n = (poly_p+idx)->n;
                        *((poly_p+idx)->layer+n) = k;
                        *((poly_p+idx)->idx+n) = i;
                        (poly_p+idx)->n++;
                }
        }
        return poly_p;
}

static POLY_P *get_init_poly_p(int max_g)
{
	int i;
	POLY_P *poly_p;

	poly_p = (POLY_P *)calloc(max_g,sizeof(POLY_P ));
	for(i = 0; i < max_g; i++)
		(poly_p+i)->n = 0;
	return poly_p;
}

static int get_init_number_layer_index(GRID_T *grid_t,POLY_P *poly_p,XYZ_I *box_size,double grid_size,int max_g,int max_layer)
{
	int i,k,n,idx,max_i,max_k,gi,gj,gk;
	XYZ_T g,pos;

	zero_vector_xyz(&g);

	for(k = 0; k <= max_layer; k++){
		for(i = 0; i < (grid_t+k)->npoint; i++){
			if((idx = g_xyz2idx((grid_t+k)->p+i,box_size,&g,grid_size)) == -1)
				continue;
			(poly_p+idx)->n++;
		}
	}
	max_i = -100000;
	max_k = 0;
	for(k = 0; k < max_g; k++){
		if(max_i < (poly_p+k)->n){
			max_i = (poly_p+k)->n;
			max_k = k;
		}
	}
	return max_i;
}

static void alloc_layer_index_area(POLY_P *poly_p,int max_g)
{
	int i;

	for(i = 0; i < max_g; i++){
		(poly_p+i)->layer = get_init_int_array((poly_p+i)->n,-1);
		(poly_p+i)->idx = get_init_int_array((poly_p+i)->n,-1);
		(poly_p+i)->n = 0;
	}
}


void g_xyz2ijk_mass(XYZ_T *pos,int *i,int *j,int *k,XYZ_T *g,double grid_size)
{
	int sign_x,sign_y,sign_z;
	double x,y,z;

	x = pos->x - g->x;
	sign_x = sign(x);
	x = fabs(x);

	y = pos->y - g->y;
	sign_y = sign(y);
	y = fabs(y);

	z = pos->z - g->z;
	sign_z = sign(z);
	z = fabs(z);

	if(fabs(x - grid_size*floor(x/grid_size)) >= grid_size/2.0)
		*i = (int )ceil(fabs(x)/grid_size);
	else
		*i = (int )floor(fabs(x)/grid_size);

	if(fabs(y - grid_size*floor(y/grid_size)) >= grid_size/2.0)
		*j = (int )ceil(fabs(y)/grid_size);
	else
		*j = (int )floor(fabs(y)/grid_size);

	if(fabs(z - grid_size*floor(z/grid_size)) >= grid_size/2.0)
		*k = (int )ceil(fabs(z)/grid_size);
	else
		*k = (int )floor(fabs(z)/grid_size);

	*i *= sign_x;
	*j *= sign_y;
	*k *= sign_z;
}

static int sign(double x)
{
	if(x < 0.0)
		return -1;
	else
		return 1;
}

int g_xyz2idx(XYZ_T *pos,XYZ_I *box_size,XYZ_T *g,double grid_size)
{
	int gi,gj,gk;

	g_xyz2ijk_mass(pos,&gi,&gj,&gk,g,grid_size);
	if(g_box_range(gi,gj,gk,box_size) == -1)
		return -1;
	return g_ijk2idx(box_size,gi,gj,gk);
}

int g_ijk2idx(XYZ_I *box_size,int i,int j,int k)
{
        int nx,ny;

        nx = 2*box_size->x + 1;
        ny = 2*box_size->y + 1;
        return (k + box_size->z)*nx*ny + nx*(j + box_size->y) + (i + box_size->x);
}


int g_box_range(int x,int y,int z,XYZ_I *box_size)
{
        if(x < -box_size->x || x > box_size->x)
                return -1;
        if(y < -box_size->y || y > box_size->y)
                return -1;
        if(z < -box_size->z || z > box_size->z)
                return -1;
        return 1;
}

int cal_total_grid_index(SW *sw,int max_grid_index)
{
        int delta_idx,length;

        delta_idx = (int )floor(sw->rc / sw->grid_size) + 1;
        length = delta_idx*2 + 1;
        return max_grid_index * (length * length * length);
}

