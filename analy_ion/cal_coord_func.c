#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

XYZ_T *get_init_xyz_t(int max_g)
{
        int i;
        XYZ_T *xyz_t;

        xyz_t = (XYZ_T *)calloc(max_g,sizeof(XYZ_T ));
        for(i = 0; i < max_g; i++){
                (xyz_t+i)->x = 0.0;
                (xyz_t+i)->y = 0.0;
                (xyz_t+i)->z = 0.0;
        }
        return xyz_t;
}

void scaling_xyz(XYZ_T *p,int p_cnt,double r)
{
        int i;
        double dis;
        XYZ_T zero;

        zero_vector_xyz(&zero);
        for(i = 0; i < p_cnt; i++){
                unit_vector_self_xyz(p+i);
                (p+i)->x *= r;
                (p+i)->y *= r;
                (p+i)->z *= r;
        }
}

double distance_xyz(XYZ_T *a0,XYZ_T *a1)
{
	double xx,yy,zz;

	xx = (a0->x - a1->x)*(a0->x - a1->x);
	yy = (a0->y - a1->y)*(a0->y - a1->y);
	zz = (a0->z - a1->z)*(a0->z - a1->z);
	return sqrt(xx+yy+zz);
}

double distance_2_xyz(XYZ_T *a0,XYZ_T *a1)
{
	double xx,yy,zz;

	xx = (a0->x - a1->x)*(a0->x - a1->x);
	yy = (a0->y - a1->y)*(a0->y - a1->y);
	zz = (a0->z - a1->z)*(a0->z - a1->z);
	return xx+yy+zz;
}

double distance0_xyz(XYZ_T *a)
{
	double xx,yy,zz;

	xx = a->x*a->x;
	yy = a->y*a->y;
	zz = a->z*a->z;
	return sqrt(xx+yy+zz);
}

double distance0_2_xyz(XYZ_T *a)
{
	double xx,yy,zz;

	xx = a->x*a->x;
	yy = a->y*a->y;
	zz = a->z*a->z;
	return xx+yy+zz;
}

double vec_angle_xyz(XYZ_T *a0,XYZ_T *a1,XYZ_T *a2)
{
	XYZ_T p0,p1;
	double r,r0,r1,ang;

	pos_vector_xyz(&p0,a0,a1);
	pos_vector_xyz(&p1,a2,a1);
	r0 = distance_xyz(a0,a1);
	r1 = distance_xyz(a2,a1);
	r = inner_product_xyz(&p0,&p1)/(r0*r1);
	if(r > 1.0)
		r = 1.0;
	if(r < -1.0)
		r = -1.0;
	ang = acos(r);
	return ang;
}

double vec_angle0_xyz(XYZ_T *a0,XYZ_T *a1)
{
	XYZ_T p0,p1,zero;
	double r,r0,r1,ang;

	zero_vector_xyz(&zero);

	pos_vector_xyz(&p0,a0,&zero);
	pos_vector_xyz(&p1,a1,&zero);

	r0 = distance0_xyz(a0);
	r1 = distance0_xyz(a1);
	r = inner_product_xyz(&p0,&p1)/(r0*r1);
	if(r > 1.0)
		r = 1.0;
	if(r < -1.0)
		r = -1.0;
	ang = acos(r);
	return ang;
}

void cal_center_xyz(XYZ_T *xyz_t,int a_cnt,XYZ_T *g)
{
	int i;
	
	zero_vector_xyz(g);
	for(i = 0; i < a_cnt; i++){
		g->x += (xyz_t+i)->x;
		g->y += (xyz_t+i)->y;
		g->z += (xyz_t+i)->z;
	}
	g->x /= (double )a_cnt;
	g->y /= (double )a_cnt;
	g->z /= (double )a_cnt;
}

void zero_vector_xyz(XYZ_T *p)
{
	p->x = p->y = p->z = 0.0;
}

void x_axis_vector_xyz(XYZ_T *p)
{
	p->x = 1.0;
	p->y = p->z = 0.0;
}

void y_axis_vector_xyz(XYZ_T *p)
{
	p->y = 1.0;
	p->x = p->z = 0.0;
}

void z_axis_vector_xyz(XYZ_T *p)
{
	p->z = 1.0;
	p->x = p->y = 0.0;
}

void unit_vector_xyz(XYZ_T *u,XYZ_T *p)
{
	double r;

	r = distance0_xyz(p);
	if(r > 1.0E-38){
		u->x = p->x/r;
		u->y = p->y/r;
		u->z = p->z/r;
	}
	else
		zero_vector_xyz(u);
}

void unit_vector_self_xyz(XYZ_T *u)
{
	double r;

	r = distance0_xyz(u);
	if(r > 1.0E-38){
		u->x /= r;
		u->y /= r;
		u->z /= r;
	}
}


void pos_vector_xyz(XYZ_T *p,XYZ_T *a1,XYZ_T *a0)
{
	p->x = a1->x - a0->x;
	p->y = a1->y - a0->y;
	p->z = a1->z - a0->z;
}

void pos_unit_vector_xyz(XYZ_T *p,XYZ_T *a1,XYZ_T *a0)
{
	double dis;

	p->x = a1->x - a0->x;
	p->y = a1->y - a0->y;
	p->z = a1->z - a0->z;
	dis = distance0_xyz(p);
	divid_vector_xyz(p,dis);
}

void cp_vector_xyz(XYZ_T *a1,XYZ_T *a0)
{
	a1->x = a0->x;
	a1->y = a0->y;
	a1->z = a0->z;
}

void add_vector_xyz(XYZ_T *a1,XYZ_T *a0)
{
	a1->x += a0->x;
	a1->y += a0->y;
	a1->z += a0->z;
}

void add_vector2_xyz(XYZ_T *p,XYZ_T *a1,XYZ_T *a0)
{
	p->x = a1->x + a0->x;
	p->y = a1->y + a0->y;
	p->z = a1->z + a0->z;
}

void del_vector_xyz(XYZ_T *a1,XYZ_T *a0)
{
	a1->x -= a0->x;
	a1->y -= a0->y;
	a1->z -= a0->z;
}

void del_vector2_xyz(XYZ_T *p,XYZ_T *a1,XYZ_T *a0)
{
	p->x = a1->x - a0->x;
	p->y = a1->y - a0->y;
	p->z = a1->z - a0->z;
}

void center_vector_xyz(XYZ_T *p,XYZ_T *a0,XYZ_T *a1)
{
	p->x = (a0->x + a1->x)/2.0;
	p->y = (a0->y + a1->y)/2.0;
	p->z = (a0->z + a1->z)/2.0;
}

void scale_vector_xyz(XYZ_T *a,double r)
{
	a->x *= r;
	a->y *= r;
	a->z *= r;
}

void multi_scalar2vector_xyz(XYZ_T *p,XYZ_T *a,double r)
{
	p->x = a->x * r;
	p->y = a->y * r;
	p->z = a->z * r;
}

void divid_vector_xyz(XYZ_T *a,double r)
{
	a->x /= r;
	a->y /= r;
	a->z /= r;
}

void divid_scalar2vector_xyz(XYZ_T *p,XYZ_T *a,double r)
{
	p->x = a->x / r;
	p->y = a->y / r;
	p->z = a->z / r;
}

double inner_product_xyz(XYZ_T *a0,XYZ_T *a1)
{
	return a0->x*a1->x + a0->y*a1->y + a0->z*a1->z;
}

void outer_product_xyz(XYZ_T *p,XYZ_T *a0,XYZ_T *a1)
{
	double r;

	p->x = a0->y*a1->z - a0->z*a1->y;
	p->y = a0->z*a1->x - a0->x*a1->z;
	p->z = a0->x*a1->y - a0->y*a1->x;
	r = distance0_xyz(p);
	divid_vector_xyz(p,r);
}

