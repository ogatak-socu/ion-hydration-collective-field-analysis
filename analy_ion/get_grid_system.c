#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static GRID_T *get_grid_t_area(int ,double );
static int get_total_point(int );
static void rd_grid_points(GRID_T *,int ,char *);
static void cp_xyz_t_array_scale(XYZ_T *,XYZ_T *,int ,double );

POLY_A *get_grid_system(int max_layer,SW *sw)
{
        int i,layer_n;
        POLY_A *poly_a;
        GRID_T *grid_t;

        poly_a = (POLY_A *)calloc(1,sizeof(POLY_A ));
	fprintf(stderr,"Reading grid points file %s\n",sw->spherical_grid_points_file);
	grid_t = get_init_grid_system_file(max_layer,sw->delta_d,sw->spherical_grid_points_file);

        poly_a->grid_d = cp_init_grid_t(grid_t,max_layer,sw,1);

        poly_a->grid_f = cp_init_grid_t(grid_t,max_layer,sw,0);
        mk_link_grid_coord(poly_a->grid_f,poly_a->grid_d,max_layer);

        poly_a->grid_m = cp_init_grid_t(grid_t,max_layer,sw,0);
        mk_link_grid_coord(poly_a->grid_m,poly_a->grid_d,max_layer);

        poly_a->grid_b = cp_init_grid_t(grid_t,max_layer,sw,0);
        mk_link_grid_coord(poly_a->grid_b,poly_a->grid_d,max_layer);

        free_grid_t_area(grid_t,max_layer);
        return poly_a;
}

GRID_T *get_init_grid_system_file(int max_layer,double scale,char *fname)
{
        int i,j,n,npoint;
        GRID_T *grid_t;

        grid_t = get_grid_t_area(max_layer,scale);
        rd_grid_points(grid_t,max_layer,fname);
        for(i = 0; i <= max_layer; i++)
                scaling_xyz((grid_t+i)->p,(grid_t+i)->npoint,(double )i*scale);
        return grid_t;
}

static GRID_T *get_grid_t_area(int max_layer,double scale)
{
        int i,j,n,npoint,base_n,line_n,inter_n;
        double r;
        GRID_T *grid_t;

        grid_t = (GRID_T *)calloc(max_layer+1,sizeof(GRID_T ));
        for(i = 0; i <= max_layer; i++){
                (grid_t+i)->npoint = npoint = get_total_point(i);
		(grid_t+i)->p = get_init_xyz_t((grid_t+i)->npoint);
                (grid_t+i)->dxyz_t = get_init_xyz_t(npoint);
                (grid_t+i)->n = get_init_unsigned_long_int_array(npoint,0);
                (grid_t+i)->weight = get_init_double_array(npoint,0.0);
        }
        return grid_t;
}

static int get_total_point(int layer_n)
{
	int base_n,line_n,inter_n;

        if(layer_n == 0)
                return 1;
        base_n = 12;
        line_n = 30*(layer_n - 1);
        if(layer_n <= 2)
                inter_n = 0;
        else
                inter_n = (layer_n - 2)*(layer_n - 1)/2*20;

        return base_n + line_n + inter_n;
}



static void rd_grid_points(GRID_T *grid_t,int max_layer,char *fname)
{
        char **buf_t;
        int i,l_cnt,layer,idx;
        double x,y,z;

        buf_t = rd_text_file_without_cr(fname,&l_cnt);
        if(l_cnt == 0){
                free(buf_t);
                exit(0);
        }
        else{
                for(i = 0; i < l_cnt; i++){
                        sscanf(*(buf_t+i),"%d %d %lf %lf %lf",&layer,&idx,&x,&y,&z);
                        if(max_layer < layer)
                                break;
                        ((grid_t+layer)->p+idx)->x = x;
                        ((grid_t+layer)->p+idx)->y = y;
                        ((grid_t+layer)->p+idx)->z = z;
                }
        }
        free_text_area(buf_t,l_cnt);
}

GRID_T *cp_init_grid_t(GRID_T *grid_t0,int max_layer,SW *sw,int flg)
{
        int i,npoint,layer_n;
        double radius;
        GRID_T *grid_t;

	grid_t = (GRID_T *)calloc(max_layer+2,sizeof(GRID_T ));
	for(i = 0; i <= max_layer; i++){
		radius = (double )i*sw->delta_d;
		(grid_t+i)->r = radius;
                (grid_t+i)->npoint = npoint = (grid_t0+i)->npoint;
		if(flg == 1){
			(grid_t+i)->p = get_init_xyz_t(npoint);
			cp_xyz_t_array_scale((grid_t+i)->p,(grid_t0+i)->p,npoint,radius);
		}
                (grid_t+i)->dxyz_t = get_init_xyz_t(npoint);
                (grid_t+i)->weight = get_init_double_array(npoint,0.0);
                (grid_t+i)->n = get_init_unsigned_long_int_array(npoint,0);
        }
        return grid_t;
}

static void cp_xyz_t_array_scale(XYZ_T *xyz_t,XYZ_T *xyz_t0,int natom,double scale)
{
        int i;

        for(i = 0; i < natom; i++){
                (xyz_t+i)->x = (xyz_t0+i)->x;
                (xyz_t+i)->y = (xyz_t0+i)->y;
                (xyz_t+i)->z = (xyz_t0+i)->z;
                unit_vector_self_xyz(xyz_t+i);
                scale_vector_xyz(xyz_t+i,scale);
        }
}

void mk_link_grid_coord(GRID_T *grid_t,GRID_T *grid_t0,int max_layer)
{
        int i;

        for(i = 0; i <= max_layer; i++)
                (grid_t+i)->p = (grid_t0+i)->p;
}

void free_grid_t_area(GRID_T *grid_t,int max_layer)
{
        int i;

        for(i = 0; i <= max_layer; i++){
                free((grid_t+i)->p);
                free((grid_t+i)->dxyz_t);
                free((grid_t+i)->weight);
                free((grid_t+i)->n);
        }
        free(grid_t);
}
