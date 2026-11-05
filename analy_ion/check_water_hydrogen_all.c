#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static void check_water_hydrogen_atoms(XYZ_T *,int ,XYZ_T *,PRMTOP *);
static int get_oxygen_idx(PRMTOP *,int );
static double get_moving_atoms_in_box(double ,double ,double );

void check_water_hydrogen_all(XYZ_T *mdcrd,XYZ_T *box,PRMTOP *prmtop)
{
	int i,sp,ep;

	for(i = 0; i < prmtop->nres; i++){
		if(strcmp(*(prmtop->labres+i),"WAT") != 0)
			continue;
		check_water_hydrogen_atoms(mdcrd,i,box,prmtop);
	}
}

static void check_water_hydrogen_atoms(XYZ_T *mdcrd,int k,XYZ_T *box,PRMTOP *prmtop)
{
	int i,idx,sp,ep;
	XYZ_T w;

	idx = get_oxygen_idx(prmtop,k);
	get_prmtop_sp_ep_idx(prmtop,k+1,k+1,&sp,&ep);
	for(i = sp; i <= ep; i++){
		if(i == idx)
			continue;
		w.x = get_moving_atoms_in_box((mdcrd+idx)->x,(mdcrd+i)->x,box->x);
		w.y = get_moving_atoms_in_box((mdcrd+idx)->y,(mdcrd+i)->y,box->y);
		w.z = get_moving_atoms_in_box((mdcrd+idx)->z,(mdcrd+i)->z,box->z);
		cp_vector_xyz(mdcrd+i,&w);
	}
}

static int get_oxygen_idx(PRMTOP *prmtop,int target_idx)
{
        int i,sp,ep;

        get_prmtop_sp_ep_idx(prmtop,target_idx+1,target_idx+1,&sp,&ep);
        for(i = sp; i <= ep; i++){
                if(strcmp(*(prmtop->igraph+i),"O") == 0)
                        return i;
        }
        return 0;
}

static double get_moving_atoms_in_box(double x1,double x0,double box)
{
        double d;

        d = fabs(x1 - x0);
        if(d > box/2.0){
                if(x0 < 0.0)
                        return x0 + box;
                else
                        return x0 - box;
        }
        return x0;
}
