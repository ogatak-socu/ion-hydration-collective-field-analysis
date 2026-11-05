#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>

#include "cal_vector.h"

int main(int argc,char **argv)
{
        SW *sw;

        if((sw = rd_options(argv,argc)) == NULL){
                exit(0);
        }
	accumulate_vector_all(sw);
	return 0;
}
