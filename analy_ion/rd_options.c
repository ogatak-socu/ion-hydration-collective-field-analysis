#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "cal_vector.h"

static void init_sw(SW *);
static void rd_option_file(SW *,char *);
static char *check_option_file(char **,int );
static int check_option_string(char **,int ,char *);
static char *get_option_single_string(char **,int ,char *);
static char **get_option_multi_string(char **,int ,int *,char *);
static int get_option_int_value(char **,int ,char *);
static double get_option_double_value(char **,int ,char *);
static int get_option_multi_tag_num(char **,int ,char *);
static void get_tag_val(char *,char *,char *);

SW *rd_options(char **buf,int l_cnt)
{
	char *options;
	SW *sw;

	sw = (SW *)calloc(1,sizeof(SW));
	init_sw(sw);
	if((options = check_option_file(buf,l_cnt)) != NULL){
		fprintf(stderr,"Input file :%s\n",options);
		rd_option_file(sw,options);
	}
	else{
		fprintf(stderr,"  Error : No input file\n");
		fprintf(stderr,"    Usage:\n\n");
		fprintf(stderr,"    executable-file -op input_file\n");
		fprintf(stderr,"    ex.\n\n");
		fprintf(stderr,"    cal_vector -op input_file\n");
	}
	free(options);
	return sw;
}

static void init_sw(SW *sw)
{
	sw->prmtop = NULL;
	sw->inpcrd = NULL;

	sw->system_size = NULL;
	sw->spherical_grid_points_file = NULL;
	sw->target_atom_index = 0;

	sw->tr_cnt = 0;
	sw->mdcrd = NULL;
	sw->f_cnt = 0;
	sw->mdfrc = NULL;

	sw->rn = 15;
	sw->delta_d = 0.2;


	sw->grid_size = 0.5;
	sw->dip_vector_prefix = NULL;
	sw->frc_vector_prefix = NULL;
	sw->bck_vector_prefix = NULL;
	sw->fwd_vector_prefix = NULL;
	sw->charge_prefix = NULL;
	sw->file_no = 0;
	sw->number_of_vectors = 1;
}

static void rd_option_file(SW *sw,char *options)
{
	char **buf_t;
	int l_cnt;

	buf_t = rd_text_file_without_cr(options,&l_cnt);
	if(check_option_string(buf_t,l_cnt,"prmtop") == 1)
		sw->prmtop		= 	get_option_single_string(buf_t,l_cnt,"prmtop");

	if(check_option_string(buf_t,l_cnt,"inpcrd") == 1)
		sw->inpcrd		= 	get_option_single_string(buf_t,l_cnt,"inpcrd");

	if(check_option_string(buf_t,l_cnt,"system_size") == 1)
		sw->system_size		= 	get_option_single_string(buf_t,l_cnt,"system_size");

	if(check_option_string(buf_t,l_cnt,"spherical_grid_points_file") == 1)
		sw->spherical_grid_points_file	= 	get_option_single_string(buf_t,l_cnt,"spherical_grid_points_file");

	if(check_option_string(buf_t,l_cnt,"target_atom_index") == 1)
		sw->target_atom_index	=	get_option_int_value(buf_t,l_cnt,"target_atom_index");

        if(check_option_string(buf_t,l_cnt,"mdcrd") == 1)
		sw->mdcrd 		=	get_option_multi_string(buf_t,l_cnt,&sw->tr_cnt,"mdcrd");

        if(check_option_string(buf_t,l_cnt,"mdfrc") == 1)
		sw->mdfrc		=       get_option_multi_string(buf_t,l_cnt,&sw->f_cnt,"mdfrc");

	if(check_option_string(buf_t,l_cnt,"rn") == 1)
		sw->rn			=	get_option_int_value(buf_t,l_cnt,"rn");

	if(check_option_string(buf_t,l_cnt,"delta_d") == 1)
		sw->delta_d		=	get_option_double_value(buf_t,l_cnt,"delta_d");

	if(check_option_string(buf_t,l_cnt,"rc") == 1)
		sw->rc			=	get_option_double_value(buf_t,l_cnt,"rc");

	if(check_option_string(buf_t,l_cnt,"alpha") == 1)
		sw->alpha		=	get_option_double_value(buf_t,l_cnt,"alpha");

	if(check_option_string(buf_t,l_cnt,"dip_vector_prefix") == 1)
		sw->dip_vector_prefix	= 	get_option_single_string(buf_t,l_cnt,"dip_vector_prefix");

	if(check_option_string(buf_t,l_cnt,"frc_vector_prefix") == 1)
		sw->frc_vector_prefix	= 	get_option_single_string(buf_t,l_cnt,"frc_vector_prefix");

	if(check_option_string(buf_t,l_cnt,"bck_vector_prefix") == 1)
		sw->bck_vector_prefix	= 	get_option_single_string(buf_t,l_cnt,"bck_vector_prefix");

	if(check_option_string(buf_t,l_cnt,"fwd_vector_prefix") == 1)
		sw->fwd_vector_prefix	= 	get_option_single_string(buf_t,l_cnt,"fwd_vector_prefix");

	if(check_option_string(buf_t,l_cnt,"charge_prefix") == 1)
		sw->charge_prefix	= 	get_option_single_string(buf_t,l_cnt,"charge_prefix");

	if(check_option_string(buf_t,l_cnt,"file_no") == 1)
		sw->file_no		=	get_option_int_value(buf_t,l_cnt,"file_no");

	if(check_option_string(buf_t,l_cnt,"number_of_vectors") == 1)
		sw->number_of_vectors	=	get_option_int_value(buf_t,l_cnt,"number_of_vectors");

	free_text_area(buf_t,l_cnt);
}

void disp_parameter(SW *sw)
{
	fprintf(stderr,"\n");
	fprintf(stderr,"-------------------------------------------------------------------------\n");
	fprintf(stderr,"Radical distance of the n-th spherical shell: rn = %4d\n",sw->rn);
	fprintf(stderr,"Radical spacing: delta_d = %5.2f\n",sw->delta_d);
	fprintf(stderr,"-------------------------------------------------------------------------\n");
	fprintf(stderr,"Local ensemble radius:    rc = %5.2f\n",sw->rc);
	fprintf(stderr,"Gauss Weight parameter:   alpha = %5.2f\n",sw->alpha);
	fprintf(stderr,"-------------------------------------------------------------------------\n");
	fprintf(stderr,"\n");
}

static char *check_option_file(char **buf,int l_cnt)
{
	char *file_name;
	int i;

	for(i = 0; i < l_cnt; i++){
		if(strcmp(*(buf+i),"-op") == 0){
			i++;
			file_name = (char *)calloc(strlen(*(buf+i))+1,sizeof(char ));
			strcpy(file_name,*(buf+i));
			return file_name;
		}
	}
	return NULL;
}

static int check_option_string(char **buf_t,int l_cnt,char *op)
{
	char *tag,*val;
	int i,flg;

	tag = (char *)calloc(MAXLINE*100,sizeof(char ));
	val = (char *)calloc(MAXLINE*100,sizeof(char ));
	flg = 0;
	for(i = 0; i < l_cnt; i++){
		if(*(*(buf_t+i)) == '#' || *(*(buf_t+i)) == '\n' || *(*(buf_t+i)) == '\0' || strcmp(*(buf_t+i),"//") == 0)
			continue;
		get_tag_val(*(buf_t+i),tag,val);
		if(strcmp(tag,op) == 0){
			flg = 1;
			break;
		}
	}
	free(tag);
	free(val);
	return flg;
}

static char *get_option_single_string(char **buf_t,int l_cnt,char *op)
{
	char *name,*tag,*val;
	int i;

	tag = (char *)calloc(MAXLINE*100,sizeof(char ));
	val = (char *)calloc(MAXLINE*100,sizeof(char ));
	name = NULL;
	for(i = 0; i < l_cnt; i++){
		if(*(*(buf_t+i)) == '#' || *(*(buf_t+i)) == '\n' || *(*(buf_t+i)) == '\0'  || strcmp(*(buf_t+i),"//") == 0)
			continue;
		get_tag_val(*(buf_t+i),tag,val);
		if(strcmp(tag,op) == 0){
			name = (char *)calloc(strlen(val)+1,sizeof(char ));
			strcpy(name,val);
		}
	}
	free(tag);
	free(val);
	return name;
}

static char **get_option_multi_string(char **buf_t,int l_cnt,int *o_cnt,char *op)
{
	char **name,*tag,*val;
	int i,k;

	*o_cnt = get_option_multi_tag_num(buf_t,l_cnt,op);
	name = (char **)calloc(*o_cnt,sizeof(char *));
	tag = (char *)calloc(MAXLINE*100,sizeof(char ));
	val = (char *)calloc(MAXLINE*100,sizeof(char ));
	for(i = k = 0; i < l_cnt; i++){
		if(*(*(buf_t+i)) == '#' || *(*(buf_t+i)) == '\n' || *(*(buf_t+i)) == '\0'  || strcmp(*(buf_t+i),"//") == 0)
			continue;
		get_tag_val(*(buf_t+i),tag,val);
		if(strcmp(tag,op) == 0){
			*(name+k) = (char *)calloc(strlen(val)+1,sizeof(char ));
			strcpy(*(name+k),val);
			k++;
		}
	}
	free(tag);
	free(val);
	return name;
}

static double get_option_double_value(char **buf_t,int l_cnt,char *op)
{
        char *tag,*val;
        int i;
        double x;

        tag = (char *)calloc(MAXLINE*100,sizeof(char ));
        val = (char *)calloc(MAXLINE*100,sizeof(char ));
        x = 0.0;
        for(i = 0; i < l_cnt; i++){
                if(*(*(buf_t+i)) == '#' || *(*(buf_t+i)) == '\n' || *(*(buf_t+i)) == '\0'  || strcmp(*(buf_t+i),"//") == 0)
                        continue;
                get_tag_val(*(buf_t+i),tag,val);
                if(strcmp(tag,op) == 0)
                        x = (double )atof(val);
        }
        free(tag);
        free(val);
        return x;
}

static int get_option_int_value(char **buf_t,int l_cnt,char *op)
{
        char *tag,*val;
        int i,x;

        tag = (char *)calloc(MAXLINE*100,sizeof(char ));
        val = (char *)calloc(MAXLINE*100,sizeof(char ));
        x = 0;
        for(i = 0; i < l_cnt; i++){
                if(*(*(buf_t+i)) == '#' || *(*(buf_t+i)) == '\n'  || *(*(buf_t+i)) == '\0' || strcmp(*(buf_t+i),"//") == 0)
                        continue;
                get_tag_val(*(buf_t+i),tag,val);
                if(strcmp(tag,op) == 0)
                        x = (int )atoi(val);
        }
        free(tag);
        free(val);
        return x;
}

static int get_option_multi_tag_num(char **buf_t,int l_cnt,char *op)
{
	char *tag,*val;
	int i,n;

	tag = (char *)calloc(MAXLINE*100,sizeof(char ));
	val = (char *)calloc(MAXLINE*100,sizeof(char ));
	n = 0;
	for(i = 0; i < l_cnt; i++){
		if(*(*(buf_t+i)) == '#' || *(*(buf_t+i)) == '\n' || *(*(buf_t+i)) == '\0'  || strcmp(*(buf_t+i),"//") == 0)
			continue;
		get_tag_val(*(buf_t+i),tag,val);
		if(strcmp(tag,op) == 0)
			n++;
	}
	free(tag);
	free(val);
	return n;
}

static void get_tag_val(char *buf,char *tag,char *val)
{
	int i,k;

	for(i = 0; *(buf+i) != ' ' && *(buf+i) != '\t' && *(buf+i) != '\0' && *(buf+i) != '\n'; i++)
		*(tag+i) = *(buf+i);
	*(tag+i) = '\0';
	i++;
	for(; *(buf+i) != '\0' && (*(buf+i) == ' ' || *(buf+i) == '\t'); i++)
		;
	for(k = 0; *(buf+i) != '\0' && *(buf+i) != '\n'; k++,i++)
		*(val+k) = *(buf+i);
	*(val+k) = '\0';
}

