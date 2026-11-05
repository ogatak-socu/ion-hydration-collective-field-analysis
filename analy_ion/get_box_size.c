#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static double grid_size;
static XYZ_I *get_min_box_size(SW *);

XYZ_I *get_box_size(SW *sw,int *max_g)
{
        XYZ_I *box_size;

	box_size = get_min_box_size(sw);
        fprintf(stderr,"BOX size %d x %d x %d\n",2*box_size->x+1,2*box_size->y+1,2*box_size->z+1);
        *max_g = (2*box_size->x+1)*(2*box_size->y+1)*(2*box_size->z+1);
        fprintf(stderr,"The number of points %d\n",*max_g);
        return box_size;
}

static XYZ_I *get_min_box_size(SW *sw)
{
        char **buf_t;
        int no,l_cnt,i,w;
        double x,y,z;
        XYZ_I *box_size;

        box_size = (XYZ_I *)calloc(1,sizeof(XYZ_I));
        box_size->x = box_size->y = box_size->z = 100000;
        box_size->dx = box_size->dy = box_size->dz = 1.0E+10;
	if(sw->system_size != NULL)
		buf_t = rd_text_file_without_cr(sw->system_size,&l_cnt);
	else
		l_cnt = 0;
	if(l_cnt == 0)
		return box_size;

        for(i = 0; i < l_cnt; i++){
		if(*(*(buf_t+i)) == '#')
			continue;
                sscanf(*(buf_t+i),"%d %lf %lf %lf",&no,&x,&y,&z);
		x -= 1.0;
		y -= 1.0;
		z -= 1.0;

                if(x < box_size->dx)
                        box_size->dx = x;
                if(y < box_size->dy)
                        box_size->dy = y;
                if(z < box_size->dz)
                        box_size->dz = z;

                w = (int )floor((x/2.0)/sw->grid_size);
                if(w < box_size->x)
                        box_size->x = w;
                w = (int )floor((y/2.0)/sw->grid_size);
                if(w < box_size->y)
                        box_size->y = w;
                w = (int )floor((z/2.0)/sw->grid_size);
                if(w < box_size->z)
                        box_size->z = w;
        }
        free_text_area(buf_t,l_cnt);
	return box_size;
}
