#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

int main(int argc,char **argv)
{
        SW *sw;
        PRMTOP *prmtop;
        RST *rst;

        if((sw = rd_options(argv,argc)) == NULL){
                exit(0);
        }

        prmtop = rd_amber_prmtop(sw->prmtop);
        fprintf(stderr,"PRMTOP file         %s\n",sw->prmtop);
//        rst = rd_amber_restart_file(sw->inpcrd);
//        fprintf(stderr,"Initial rstart file %s\n",sw->inpcrd);
        fprintf(stderr,"Target atom: %s\n",*(prmtop->igraph+sw->target_atom_index-1));
	disp_parameter(sw);

	cal_vector_mdcrd(prmtop,sw);

	return 0;
}
