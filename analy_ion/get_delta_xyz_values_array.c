#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

static double get_moving_atoms_element(double ,double ,double );

void get_delta_xyz_values_array(XYZ_T *dxyz_t,XYZ_T *mdcrd1,XYZ_T *mdcrd0,XYZ_T *box,PRMTOP *prmtop)
{
        int i;

        for(i = 0; i< prmtop->natom; i++){
                (dxyz_t+i)->x = get_moving_atoms_element((mdcrd1+i)->x,(mdcrd0+i)->x,box->x);
                (dxyz_t+i)->y = get_moving_atoms_element((mdcrd1+i)->y,(mdcrd0+i)->y,box->y);
                (dxyz_t+i)->z = get_moving_atoms_element((mdcrd1+i)->z,(mdcrd0+i)->z,box->z);
        }
}

static double get_moving_atoms_element(double x1,double x0,double box)
{
        double d;

        d = fabs(x1 - x0);
        if(d > box/2.0){
                if(x1 > x0){
                        return x1 - box - x0;
                }
                else{
                        return x1 + box - x0;
                }
        }
        return x1 - x0;
}
