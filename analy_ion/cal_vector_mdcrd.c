#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

#define MAX_EXP_V   900000
#define SCALE_INDEX 100000.0
#define SCALE_DIS   0.00001

static double	grid_size;
static double	delta_d;
static int	max_g;
static double	alpha;
static double	rc;
static double	rc2;
static int	target_idx;
static int	max_layer;
static int	total_grid_index;
static int	delta_idx;

static int case_value(int ,int ,int );
static void cal_dipole_vector(XYZ_T *,PRMTOP *,POLY_A *,POLY_P *,XYZ_I *);
static void cal_frc_diplacement_vector(XYZ_T *,XYZ_T *,XYZ_T *,XYZ_T *,PRMTOP *,POLY_A *,POLY_P *,XYZ_I *);
static void get_mass_center_water(XYZ_T *,XYZ_T *,int ,int );
static void get_dipole_water(XYZ_T *,PRMTOP *,int ,int ,XYZ_T *);
static void store_grid_values(GRID_T *,XYZ_T *,int ,double );
static void wt_poly_a_values(POLY_A *,int ,SW *);
static void wt_grid_t_values(GRID_T *,int ,char *);
static int num_grid_t_values(GRID_T *);

void cal_vector_mdcrd(PRMTOP *prmtop,SW *sw)
{
	int i,j,k,no,max_g,t_cnt,f_cnt,n,max_grid_index;
	FILE *fp;
	XYZ_I *box_size;
	XYZ_T *dxyz_f,*dxyz_b,p;
	TRAJ *traj,*frc,*wtraj,*wfrc;
	POLY_A *poly_a;
	POLY_P *poly_p;

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
	poly_a = get_grid_system(max_layer,sw);

	poly_p = assign_spherical_system2cube_system(poly_a->grid_m,box_size,grid_size,max_g,max_layer,&max_grid_index);

	total_grid_index = cal_total_grid_index(sw,max_grid_index);


        wtraj = get_init_traj_area(prmtop->natom);
        wfrc = get_init_traj_area(prmtop->natom);
        dxyz_b = get_init_xyz_t(prmtop->natom);
        dxyz_f = get_init_xyz_t(prmtop->natom);
	for(i = k = 0; i < sw->tr_cnt; i++){
                fprintf(stderr,"%s ",*(sw->mdcrd+i));
                traj = rd_amber_traj(*(sw->mdcrd+i),prmtop->natom,1,2000,&t_cnt);
                fprintf(stderr,"%6d\n",t_cnt);
                fprintf(stderr,"%s ",*(sw->mdfrc+i));
                frc = rd_amber_traj_other(*(sw->mdfrc+i),prmtop->natom,1,2000,&f_cnt);
                fprintf(stderr,"%6d\n",f_cnt);
                if(t_cnt != f_cnt){
                        fprintf(stderr,"Mismatch the number of frames   traj: %d  frc: %d\n",t_cnt,f_cnt);
                        exit(0);
                }

                for(j = 0; j < t_cnt; j++){
			fprintf(stderr,"%4d",j);
			if((j+1) % 25 == 0)
				fprintf(stderr,"\n");

			cp_vector_xyz(&p,(traj+j)->mdcrd+target_idx);
			shift_traj_center(traj+j,&p,prmtop->natom);
                        check_water_hydrogen_all((traj+j)->mdcrd,(traj+j)->box,prmtop);
			cal_dipole_vector((traj+j)->mdcrd,prmtop,poly_a,poly_p,box_size);
			switch (case_value(k,i == sw->tr_cnt-1, j == t_cnt-1)){
				case 0:
                                	get_delta_xyz_values_array(dxyz_f,(traj+j+1)->mdcrd,(traj+j)->mdcrd,(traj+j)->box,prmtop);
					cal_frc_diplacement_vector((traj+j)->mdcrd,dxyz_f,NULL,(frc+j)->mdcrd,prmtop,poly_a,poly_p,box_size);
					break;
				case 1:
					get_delta_xyz_values_array(dxyz_f,(traj+j+1)->mdcrd,(traj+j)->mdcrd,(traj+j)->box,prmtop);
					get_delta_xyz_values_array(dxyz_b,(traj+j)->mdcrd,(traj+j-1)->mdcrd,(traj+j)->box,prmtop);
					cal_frc_diplacement_vector((traj+j)->mdcrd,dxyz_f,dxyz_b,(frc+j)->mdcrd,prmtop,poly_a,poly_p,box_size);
					break;
				case 2:
					get_delta_xyz_values_array(dxyz_b,(traj+j)->mdcrd,(traj+j-1)->mdcrd,(traj+j)->box,prmtop);
					cal_frc_diplacement_vector((traj+j)->mdcrd,NULL,dxyz_b,(frc+j)->mdcrd,prmtop,poly_a,poly_p,box_size);
					break;
				case 3:
					get_delta_xyz_values_array(dxyz_f,(traj+j)->mdcrd,wtraj->mdcrd,wtraj->box,prmtop);
					cal_frc_diplacement_vector(wtraj->mdcrd,dxyz_f,dxyz_b,wfrc->mdcrd,prmtop,poly_a,poly_p,box_size);
					break;
				default:
					break;
			}
			cp_traj_mdcrd(dxyz_b,dxyz_f,prmtop->natom);
			cp_traj_mdcrd(wtraj->mdcrd,(traj+j)->mdcrd,prmtop->natom);
			cp_traj_mdcrd(wtraj->box,(traj+j)->box,1);
			cp_traj_mdcrd(wfrc->mdcrd,(frc+j)->mdcrd,prmtop->natom);
			k++;
		}
		free_traj_area(traj,t_cnt);
		free_traj_area_other(frc,f_cnt);
	}
	wt_poly_a_values(poly_a,max_layer,sw);
}

static int case_value(int k,int cond1,int cond2)
{
	if(k == 0)
		return 0;
	if(k == 1)
		return 1;
	if(cond1 == 1 && cond2 == 1)
		return 2;
	return 3;
}

static void cal_dipole_vector(XYZ_T *mdcrd,PRMTOP *prmtop,POLY_A *poly_a,POLY_P *poly_p,XYZ_I *box_size)
{
	int k,i,sp,ep,idx,layer;
	double r,weight;
        XYZ_T di,w;
        RANGE_P *range_p;

	range_p = get_init_range_p(total_grid_index);
	for(k = 0; k < prmtop->nres; k++){
		if(target_idx == k || strcmp(*(prmtop->labres+k),"WAT") != 0)
			continue;
		get_prmtop_sp_ep_idx(prmtop,k+1,k+1,&sp,&ep);
		get_dipole_water(mdcrd,prmtop,sp,ep,&di);
		get_mass_center_water(&w,mdcrd,sp,ep);

		r = distance0_xyz(&w);
		if((int )round(r/delta_d) > max_layer)
			continue;
		set_init_range_p(range_p,total_grid_index);
		get_range_point_polygon_grid_normal(range_p,poly_a->grid_m,&w,poly_p,box_size,grid_size,delta_idx,rc2,alpha);
		for(i = 0; i <= total_grid_index; i++){
			if((idx = (range_p+i)->idx) == -1)
				break;
			layer = (range_p+i)->layer;
			weight = (range_p+i)->weight;
			store_grid_values(poly_a->grid_m+layer,&di,idx,weight);
                }
	}
	free(range_p);
}

static void cal_frc_diplacement_vector(XYZ_T *mdcrd,XYZ_T *dxyz_f,XYZ_T *dxyz_b,XYZ_T *frc,PRMTOP *prmtop,POLY_A *poly_a,POLY_P *poly_p,XYZ_I *box_size)
{
	int k,i,j,sp,ep,layer,idx;
	double r,weight;
        XYZ_T mf,mb,fx;
        RANGE_P *range_p;

	range_p = get_init_range_p(total_grid_index);
	for(k = 0; k < prmtop->nres; k++){
		if(target_idx == k || strcmp(*(prmtop->labres+k),"WAT") != 0)
			continue;
		get_prmtop_sp_ep_idx(prmtop,k+1,k+1,&sp,&ep);
                for(i = sp; i <= ep; i++){
			if(strcmp(*(prmtop->igraph+i),"EPW") == 0 || strcmp(*(prmtop->igraph+i),"EP1") == 0 || strcmp(*(prmtop->igraph+i),"EP2") == 0)
				continue;
			r = distance0_xyz(mdcrd+i);
			if((int )round(r/delta_d) > max_layer)
				continue;

			if(dxyz_f != NULL)
				cp_vector_xyz(&mf,dxyz_f+i);
			if(dxyz_b != NULL)
				cp_vector_xyz(&mb,dxyz_b+i);
			cp_vector_xyz(&fx,frc+i);
			set_init_range_p(range_p,total_grid_index);
			get_range_point_polygon_grid_normal(range_p,poly_a->grid_m,mdcrd+i,poly_p,box_size,grid_size,delta_idx,rc2,alpha);
			for(j = 0; j <= total_grid_index; j++){
				if((idx = (range_p+j)->idx) == -1)
					break;
				layer = (range_p+j)->layer;
				weight = (range_p+j)->weight;

				if(dxyz_f != NULL)
					store_grid_values(poly_a->grid_m+layer,&mf,idx,weight);
				if(dxyz_b != NULL)
					store_grid_values(poly_a->grid_b+layer,&mb,idx,weight);
				store_grid_values(poly_a->grid_f+layer,&fx,idx,weight);
                        }
                }
	}
	free(range_p);
}

static void get_mass_center_water(XYZ_T *w,XYZ_T *mdcrd,int sp,int ep)
{
        int i;

        zero_vector_xyz(w);
        for(i = sp; i <= ep; i++)
                add_vector_xyz(w,mdcrd+i);
        divid_vector_xyz(w,(double )(ep-sp+1));
}

static void get_dipole_water(XYZ_T *mdcrd,PRMTOP *prmtop,int sp,int ep,XYZ_T *di)
{
	int i;
	double q;

	zero_vector_xyz(di);
	for(i = sp; i <= ep; i++){
		q = *(prmtop->chrg+i)/sqrt(EU);
		di->x += q*(mdcrd+i)->x*C2DEBYE;
		di->y += q*(mdcrd+i)->y*C2DEBYE;
		di->z += q*(mdcrd+i)->z*C2DEBYE;
	}
}

static void store_grid_values(GRID_T *grid_t,XYZ_T *vec,int idx,double weight)
{
	double len;

	if((len = distance0_xyz(vec)) > 1.0E-07){
		(grid_t->dxyz_t+idx)->x += weight * vec->x;
		(grid_t->dxyz_t+idx)->y += weight * vec->y;
		(grid_t->dxyz_t+idx)->z += weight * vec->z;
		(*(grid_t->weight+idx)) += weight;
		(*(grid_t->n+idx))++;
		if(*(grid_t->n+idx) < 0){
			fprintf(stderr,"Stack overflow grid_t->n\n");
			exit(0);
		}
	}
}

static void wt_poly_a_values(POLY_A *poly_a,int max_layer,SW *sw)
{
	char fname[MAXLINE];

	sprintf(fname,"%s.%d.out",sw->fwd_vector_prefix,sw->file_no);
	wt_grid_t_values(poly_a->grid_m,max_layer,fname);

	sprintf(fname,"%s.%d.out",sw->bck_vector_prefix,sw->file_no);
	wt_grid_t_values(poly_a->grid_b,max_layer,fname);

	sprintf(fname,"%s.%d.out",sw->frc_vector_prefix,sw->file_no);
	wt_grid_t_values(poly_a->grid_f,max_layer,fname);

	sprintf(fname,"%s.%d.out",sw->dip_vector_prefix,sw->file_no);
	wt_grid_t_values(poly_a->grid_d,max_layer,fname);
}

static void wt_grid_t_values(GRID_T *grid_t,int max_layer,char *fname)
{
	int k,i,n;
	FILE *fp;
	WGRID_T wgrid_t;

	for(k = n = 0; k < max_layer; k++)
		n += num_grid_t_values(grid_t+k);
	if(n != 0){
		fp = wt_binary_file_open(fname);
		for(k = 0; k < max_layer; k++){
			for(i = 0; i < (grid_t+k)->npoint; i++){
				if(*((grid_t+k)->n+i) == 0)
					continue;
				wgrid_t.layer = k;
				wgrid_t.idx = i;
				wgrid_t.x = ((grid_t+k)->dxyz_t+i)->x;
				wgrid_t.y = ((grid_t+k)->dxyz_t+i)->y;
				wgrid_t.z = ((grid_t+k)->dxyz_t+i)->z;
				wgrid_t.n = *((grid_t+k)->n+i);
				wgrid_t.weight = *((grid_t+k)->weight+i);
				fwrite(&wgrid_t,sizeof(WGRID_T ),1,fp);
			}
		}
		wt_file_close(fname,fp);
	}
}

static int num_grid_t_values(GRID_T *grid_t)
{
	int i,n;

	for(i = n = 0; i < grid_t->npoint; i++)
		n += *(grid_t->n+i);
	return n;
}

