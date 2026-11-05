#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>
#include <limits.h>

#include "cal_vector.h"

static int max_layer;

static void rd_grid_charge_values(CHARGE_T *,SW *,int );
static void wt_grid_charge_values(CHARGE_T *,SW *);

void accumulate_charge_all(SW *sw)
{
	char fname[MAXLINE];
        int k,max_g,n;
	CHARGE_T *charge_t;

        max_layer = (int )floor((double )sw->rn*(1.0/sw->delta_d)) + 1;

	charge_t = get_init_charge_t(max_layer,sw);
	fprintf(stderr,"Number of vectors: %d\n",sw->number_of_vectors);
        for(k = 0; k < sw->number_of_vectors; k++)
		rd_grid_charge_values(charge_t,sw,k+1);

	wt_grid_charge_values(charge_t,sw);
}


static void rd_grid_charge_values(CHARGE_T *charge_t,SW *sw,int no)
{
	char fname[MAXLINE];
	int n,layer,idx;
	FILE *fp;
	WCHARGE_T wcharge_t;

	sprintf(fname,"%s.%d.out",sw->charge_prefix,no);
	fprintf(stderr,"Reading %s\n",fname);
	if((fp = rd_binary_file_open(fname)) != NULL){
		while((n = fread(&wcharge_t,sizeof(WCHARGE_T),1,fp)) == 1){
			layer = wcharge_t.layer;
			if(layer >= max_layer)
				continue;
			idx = wcharge_t.idx;
			*((charge_t->grid_p+layer)->charge+idx) += wcharge_t.charge_p;
			*((charge_t->grid_n+layer)->charge+idx) += wcharge_t.charge_n;
			*((charge_t->grid_p+layer)->weight+idx) += wcharge_t.weight_p;
			*((charge_t->grid_n+layer)->weight+idx) += wcharge_t.weight_n;
			*((charge_t->grid_p+layer)->n+idx) += wcharge_t.n_p;
			*((charge_t->grid_n+layer)->n+idx) += wcharge_t.n_n;
		}
		rd_file_close(fname,fp);
	}
}

static void wt_grid_charge_values(CHARGE_T *charge_t,SW *sw)
{
        char fname[MAXLINE];
        int k,i,n;
        double w,c,tc;
        FILE *fp;

        sprintf(fname,"%s.data.out",sw->charge_prefix);
        fprintf(stderr,"Writing %s\n",fname);

        fp = wt_file_open(fname);
        for(k = 0; k < max_layer; k++){
                fprintf(fp,"%6.2f ",(charge_t->grid_p+k)->r);
                n = 0;
                c = 0.0;
                for(i = 0; i < (charge_t->grid_p+k)->npoint; i++){
                        w = *((charge_t->grid_p+k)->weight+i) + *((charge_t->grid_n+k)->weight+i);
                        tc = *((charge_t->grid_p+k)->charge+i) + *((charge_t->grid_n+k)->charge+i);
                        c += tc/w;
                        n++;
                }
                if(n == 0)
                        fprintf(fp,"%9.5f\n",0.0);
                else
                        fprintf(fp,"%9.5f\n",c/(double )n);
        }
        wt_file_close(fname,fp);
}

