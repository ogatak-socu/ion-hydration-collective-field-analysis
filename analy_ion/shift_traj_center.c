#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static int check_center_in_box(XYZ_T *,XYZ_T *);
static void shift_mol_in_box_x(XYZ_T *,XYZ_T *,int ,int ,XYZ_T *);
static void shift_mol_in_box_y(XYZ_T *,XYZ_T *,int ,int ,XYZ_T *);
static void shift_mol_in_box_z(XYZ_T *,XYZ_T *,int ,int ,XYZ_T *);
static void shift_element_x(XYZ_T *,int ,int ,double );
static void shift_element_y(XYZ_T *,int ,int ,double );
static void shift_element_z(XYZ_T *,int ,int ,double );
static double atom_inbox_element(double ,double );

void shift_traj_center(TRAJ *traj,XYZ_T *p0,int natom)
{
	int i;
	double x,y,z;

	for(i = 0; i < natom; i++){
		x = (traj->mdcrd+i)->x - p0->x;
		y = (traj->mdcrd+i)->y - p0->y;
		z = (traj->mdcrd+i)->z - p0->z;
		(traj->mdcrd+i)->x = atom_inbox_element(x,(traj->box)->x);
		(traj->mdcrd+i)->y = atom_inbox_element(y,(traj->box)->y);
		(traj->mdcrd+i)->z = atom_inbox_element(z,(traj->box)->z);
	}
}

static double atom_inbox_element(double x,double box)
{
	if(x > box/2.0)
		return x - box;
	else{
		if(x < -box/2.0)
			return x + box;
		else
			return x;
	}
	fprintf(stderr,"Error  atom_inbox_element\n");
}

