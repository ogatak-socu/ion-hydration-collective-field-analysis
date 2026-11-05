#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static void get_charge_array_area(GRID_T *,int );

CHARGE_T *get_init_charge_t(int max_layer,SW *sw)
{
	int i;
	CHARGE_T *charge_t;
	GRID_T *grid_t;

	charge_t = (CHARGE_T *)calloc(1,sizeof(CHARGE_T ));
	fprintf(stderr,"Reading grid points file %s\n",sw->spherical_grid_points_file);
	grid_t = get_init_grid_system_file(max_layer,sw->delta_d,sw->spherical_grid_points_file);

	charge_t->grid_p = cp_init_grid_t(grid_t,max_layer,sw,1);
	get_charge_array_area(charge_t->grid_p,max_layer);
	
	charge_t->grid_n = cp_init_grid_t(grid_t,max_layer,sw,0);
	mk_link_grid_coord(charge_t->grid_n,charge_t->grid_p,max_layer);
	get_charge_array_area(charge_t->grid_n,max_layer);

	free_grid_t_area(grid_t,max_layer);

        return charge_t;
}

static void get_charge_array_area(GRID_T *grid_t,int max_layer)
{
	int i;

	for(i = 0; i <= max_layer; i++)
		(grid_t+i)->charge = get_init_double_array((grid_t+i)->npoint,0.0);
}
