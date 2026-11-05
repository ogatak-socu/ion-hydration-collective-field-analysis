#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static double   grid_size;
static double   delta_d;
static int      max_g;
static double   alpha;
static double   rc;
static double   rc2;
static int      target_idx;
static int      max_layer;
static int      total_grid_index;
static int      delta_idx;

static void cal_grid_charge_values(XYZ_T *,PRMTOP *,CHARGE_T *,POLY_P *,XYZ_I *);
static void store_grid_t_values4charge(GRID_T *,int ,double ,double );
static void wt_charge_values(CHARGE_T *,int ,SW *);

void cal_charge_mdcrd(PRMTOP *prmtop,SW *sw)
{
	int i,j,t_cnt,max_g,max_grid_index;
	FILE *fp;
	XYZ_T p;
	XYZ_I *box_size;
	POLY_P *poly_p;
	CHARGE_T *charge_t;
	TRAJ *traj;

        grid_size = sw->grid_size;
        box_size = get_box_size(sw,&max_g);
        target_idx = sw->target_atom_index-1;

        fprintf(stderr,"Target atom %s\n",*(prmtop->igraph+target_idx));

        rc = (double )sw->rc; //scan_dis_area
        rc2 = rc * rc;
        alpha = sw->alpha;   //g_sigma
        delta_d = sw->delta_d; // poly_radius
        delta_idx = floor(rc/grid_size) + 1;

        max_layer = (int )floor((double )sw->rn*(1.0/sw->delta_d)) + 1;
        fprintf(stderr,"Total %6d layers in %10d A\n",max_layer,sw->rn);

	charge_t = get_init_charge_t(max_layer,sw);
	poly_p = assign_spherical_system2cube_system(charge_t->grid_p,box_size,grid_size,max_g,max_layer,&max_grid_index);
	total_grid_index = cal_total_grid_index(sw,max_grid_index);

	for(i = 0; i < sw->tr_cnt; i++){
		fprintf(stderr,"%s ",*(sw->mdcrd+i));
		traj = rd_amber_traj(*(sw->mdcrd+i),prmtop->natom,1,2000,&t_cnt);
		fprintf(stderr,"%6d\n",t_cnt);

                for(j = 0; j < t_cnt; j++){

                        fprintf(stderr,"%4d",j);
                        if((j+1) % 25 == 0)
                                fprintf(stderr,"\n");
			cp_vector_xyz(&p,(traj+j)->mdcrd+target_idx);
			shift_traj_center(traj+j,&p,prmtop->natom);
			check_water_hydrogen_all((traj+j)->mdcrd,(traj+j)->box,prmtop);
			cal_grid_charge_values((traj+j)->mdcrd,prmtop,charge_t,poly_p,box_size);
                }
		free_traj_area(traj,t_cnt);
	}
	wt_charge_values(charge_t,max_layer,sw);
}

static void cal_grid_charge_values(XYZ_T *mdcrd,PRMTOP *prmtop,CHARGE_T *charge_t,POLY_P *poly_p,XYZ_I *box_size)
{
        int k,i,j,idx,sp,ep,layer;
	double r,w,c;
	RANGE_P *range_p;

	range_p = get_init_range_p(total_grid_index);
	for(k = 0; k < prmtop->nres; k++){
		if(target_idx == k || strcmp(*(prmtop->labres+k),"WAT") != 0)
			continue;
		get_prmtop_sp_ep_idx(prmtop,k+1,k+1,&sp,&ep);
		for(i = sp; i <= ep; i++){
			c = *(prmtop->chrg+i)/sqrt(EU);
			r = distance0_xyz(mdcrd+i);

			if((int )round(r/delta_d) > max_layer)
				continue;

			set_init_range_p(range_p,total_grid_index);
                        get_range_point_polygon_grid_normal(range_p,charge_t->grid_p,mdcrd+i,poly_p,box_size,grid_size,delta_idx,rc2,alpha);

                        for(j = 0; j <= total_grid_index; j++){
                                if((idx = (range_p+j)->idx) == -1)
                                        break;
                                layer = (range_p+j)->layer;
				w = (range_p+j)->weight;

				if(c >= 0.0)
                                	store_grid_t_values4charge(charge_t->grid_p+layer,idx,c,w);
				else
                                	store_grid_t_values4charge(charge_t->grid_n+layer,idx,c,w);
                        }
                }
        }
	free(range_p);
}

static void store_grid_t_values4charge(GRID_T *grid_t,int idx,double c,double weight)
{
        (*(grid_t->charge+idx)) += weight*c;
        (*(grid_t->weight+idx)) += weight;
        (*(grid_t->n+idx))++;
}

static void wt_charge_values(CHARGE_T *charge_t,int max_layer,SW *sw)
{
	char fname[MAXLINE];
	int k,i;
	FILE *fp;
	WCHARGE_T wcharge_t;

	sprintf(fname,"%s.%d.out",sw->charge_prefix,sw->file_no);
        fp = wt_binary_file_open(fname);
        for(k = 0; k < max_layer; k++){
                for(i = 0; i < (charge_t->grid_p+k)->npoint; i++){
                        if(*((charge_t->grid_p+k)->n+i) + *((charge_t->grid_n+k)->n+i) == 0)
                                continue;
                        wcharge_t.layer = k;
                        wcharge_t.idx = i;
                        wcharge_t.charge_p = *((charge_t->grid_p+k)->charge+i);
                        wcharge_t.charge_n = *((charge_t->grid_n+k)->charge+i);
                        wcharge_t.weight_p = *((charge_t->grid_p+k)->weight+i);
                        wcharge_t.weight_n = *((charge_t->grid_n+k)->weight+i);
                        wcharge_t.n_p = *((charge_t->grid_p+k)->n+i);
                        wcharge_t.n_n = *((charge_t->grid_n+k)->n+i);
                        fwrite(&wcharge_t,sizeof(WCHARGE_T ),1,fp);
                }
        }
        wt_file_close(fname,fp);
}
