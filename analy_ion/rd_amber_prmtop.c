/*                                              */
/*      Copyright : Koji Ogata                  */
/*                                              */
/*      Date      : 2005.08.21                  */
/*                                              */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "cal_vector.h"

#define CHRG_UNIT 18.2223

static int	skip_line(char **,int ,char *);
static void	get_prmtop_pntr(PRMTOP *,char **,int );
static char	**get_label_name_area(int ,int );
static char	**get_label_name(int ,int ,char **,int ,char *);
static void	free_label_name_area(char **,int );
static double	*get_double_value(int ,char **,int ,char *);
static int	*get_int_value(int ,char **,int ,char *);
static void	convert_charge(double *,int );
static void	sum_numex(int *,int );
static void	get_index_area(PRMTOP *);

static void	get_index_values_3(int *,int *,int *,int ,char **,int ,char *);
static void	get_index_values_4(int *,int *,int *,int *,int ,char **,int ,char *);
static void	get_index_values_5(int *,int *,int *,int *,int *,int ,char **,int ,char *);
static void	check_new_line(int ,int ,int *,int *,char **);
static char	*identify_term(PRMTOP *);
static char 	*identify_term_from_bond(PRMTOP *prmtop);
static void	check_connect_residue(PRMTOP *,int *,int ,int );
static int	check_connect_res(PRMTOP *,int ,int ,int ,int );

PRMTOP *rd_amber_prmtop(char *fname)
{
	char **buf;
	int l_cnt,i,n,*w;
	PRMTOP *prmtop;

	prmtop = (PRMTOP *)calloc(1,sizeof(PRMTOP ));
	buf = rd_text_file(fname,&l_cnt);
	get_prmtop_pntr(prmtop,buf,l_cnt);
	prmtop->igraph = get_label_name(prmtop->natom,5,buf,l_cnt,"%FLAG ATOM_NAME");
	prmtop->chrg = get_double_value(prmtop->natom,buf,l_cnt,"%FLAG CHARGE");
	prmtop->amass = get_double_value(prmtop->natom,buf,l_cnt,"%FLAG MASS");

	prmtop->iac = get_int_value(prmtop->natom,buf,l_cnt,"%FLAG ATOM_TYPE_INDEX");
	prmtop->numex = get_int_value(prmtop->natom,buf,l_cnt,"%FLAG NUMBER_EXCLUDED_ATOMS");
	sum_numex(prmtop->numex,prmtop->natom);
	prmtop->ico = get_int_value(prmtop->ntypes*prmtop->ntypes,buf,l_cnt,"%FLAG NONBONDED_PARM_INDEX");
	prmtop->labres = get_label_name(prmtop->nres,5,buf,l_cnt,"%FLAG RESIDUE_LABEL");
	prmtop->ipres = get_int_value(prmtop->nres,buf,l_cnt,"%FLAG RESIDUE_POINTER");
	prmtop->rk = get_double_value(prmtop->numbnd,buf,l_cnt,"%FLAG BOND_FORCE_CONSTANT");
	prmtop->req = get_double_value(prmtop->numbnd,buf,l_cnt,"%FLAG BOND_EQUIL_VALUE");
	prmtop->tk = get_double_value(prmtop->numang,buf,l_cnt,"%FLAG ANGLE_FORCE_CONSTANT");
	prmtop->teq = get_double_value(prmtop->numang,buf,l_cnt,"%FLAG ANGLE_EQUIL_VALUE");

	prmtop->pk = get_double_value(prmtop->nptra,buf,l_cnt,"%FLAG DIHEDRAL_FORCE_CONSTANT");
	prmtop->pn = get_double_value(prmtop->nptra,buf,l_cnt,"%FLAG DIHEDRAL_PERIODICITY");
	prmtop->phase = get_double_value(prmtop->nptra,buf,l_cnt,"%FLAG DIHEDRAL_PHASE");

	prmtop->one_scee = get_double_value(prmtop->nptra,buf,l_cnt,"%FLAG SCEE_SCALE_FACTOR");
	prmtop->one_scnb = get_double_value(prmtop->nptra,buf,l_cnt,"%FLAG SCNB_SCALE_FACTOR");

	prmtop->solty = get_double_value(prmtop->natyp,buf,l_cnt,"%FLAG SOLTY");
	n = prmtop->ntypes*(prmtop->ntypes+1)/2;

	prmtop->cn1 = get_double_value(n,buf,l_cnt,"%FLAG LENNARD_JONES_ACOEF");
	prmtop->cn2 = get_double_value(n,buf,l_cnt,"%FLAG LENNARD_JONES_BCOEF");
	get_index_area(prmtop);

	get_index_values_3(prmtop->ibh,prmtop->jbh,prmtop->icbh,prmtop->nbonh,buf,l_cnt,"%FLAG BONDS_INC_HYDROGEN");
	get_index_values_3(prmtop->ib,prmtop->jb,prmtop->icb,prmtop->nbona,buf,l_cnt,"%FLAG BONDS_WITHOUT_HYDROGEN");

	get_index_values_4(prmtop->ith,prmtop->jth,prmtop->kth,prmtop->icth,prmtop->ntheth,buf,l_cnt,"%FLAG ANGLES_INC_HYDROGEN");
	get_index_values_4(prmtop->it,prmtop->jt,prmtop->kt,prmtop->ict,prmtop->ntheta,buf,l_cnt,"%FLAG ANGLES_WITHOUT_HYDROGEN");

	get_index_values_5(prmtop->iph,prmtop->jph,prmtop->kph,prmtop->lph,prmtop->icph,prmtop->nphih,buf,l_cnt,"%FLAG DIHEDRALS_INC_HYDROGEN");
	get_index_values_5(prmtop->ip,prmtop->jp,prmtop->kp,prmtop->lp,prmtop->icp,prmtop->nphia,buf,l_cnt,"%FLAG DIHEDRALS_WITHOUT_HYDROGEN");

	prmtop->natex = get_int_value(prmtop->nnb,buf,l_cnt,"%FLAG EXCLUDED_ATOMS_LIST");
	prmtop->asol = get_double_value(prmtop->nphb,buf,l_cnt,"%FLAG HBOND_ACOEF");
	prmtop->bsol = get_double_value(prmtop->nphb,buf,l_cnt,"%FLAG HBOND_BCOEF");
	prmtop->hbcut = get_double_value(prmtop->nphb,buf,l_cnt,"%FLAG HBCUT");
	prmtop->isymbol = get_label_name(prmtop->natom,5,buf,l_cnt,"%FLAG AMBER_ATOM_TYPE");
	prmtop->itree = get_label_name(prmtop->natom,5,buf,l_cnt,"%FLAG TREE_CHAIN_CLASSIFICATION");
	prmtop->join = get_int_value(prmtop->natom,buf,l_cnt,"%FLAG JOIN_ARRAY");
	prmtop->irotat = get_int_value(prmtop->natom,buf,l_cnt,"%FLAG IROTAT");

	if(prmtop->ifbox > 0){
		w = get_int_value(3,buf,l_cnt,"%FLAG SOLVENT_POINTERS");
		prmtop->iptres = *w;
		prmtop->nspm = *(w+1);
		prmtop->nspsol = *(w+2);
		free(w);
		prmtop->nsp = get_int_value(prmtop->nspm,buf,l_cnt,"%FLAG ATOMS_PER_MOLECULE");
		w = get_int_value(4,buf,l_cnt,"%FLAG BOX_DIMENSIONS");
		prmtop->beta = *w;
		prmtop->box = (double *)calloc(3,sizeof(double ));
		*prmtop->box = *w;
		*(prmtop->box+1) = *(w+1);
		*(prmtop->box+2) = *(w+2);
		free(w);
		prmtop->term = identify_term(prmtop);
	}
	else
		prmtop->term = identify_term_from_bond(prmtop);

	free_text_area(buf,l_cnt);
	return prmtop;
}

void free_prmtop_area(PRMTOP *prmtop)
{
	free_label_name_area(prmtop->igraph,prmtop->natom);
	free(prmtop->chrg);
	free(prmtop->amass);
	free(prmtop->iac);
	free(prmtop->numex);
	free(prmtop->ico);
	free_label_name_area(prmtop->labres,prmtop->nres);
	free(prmtop->ipres);
	free(prmtop->rk);
	free(prmtop->req);
	free(prmtop->tk);
	free(prmtop->teq);
	free(prmtop->pk);
	free(prmtop->pn);
	free(prmtop->phase);
}

static void get_prmtop_pntr(PRMTOP *prmtop,char **buf,int l_cnt)
{
	int n;

	n = skip_line(buf,l_cnt,"%FLAG POINTERS")+2;
	prmtop->natom  = atoi(*(buf+n));
	prmtop->ntypes = atoi(*(buf+n)+8);
	prmtop->nbonh  = atoi(*(buf+n)+16);
	prmtop->mbona  = atoi(*(buf+n)+24);
	prmtop->ntheth = atoi(*(buf+n)+32);
	prmtop->mtheta = atoi(*(buf+n)+40);
	prmtop->nphih  = atoi(*(buf+n)+48);
	prmtop->mphia  = atoi(*(buf+n)+56);
	prmtop->nhparm = atoi(*(buf+n)+64);
	prmtop->nparm  = atoi(*(buf+n)+72);
	n++;
	prmtop->nnb    = atoi(*(buf+n));
	prmtop->nres   = atoi(*(buf+n)+8);
	prmtop->nbona  = atoi(*(buf+n)+16);
	prmtop->ntheta = atoi(*(buf+n)+24);
	prmtop->nphia  = atoi(*(buf+n)+32);
	prmtop->numbnd = atoi(*(buf+n)+40);
	prmtop->numang = atoi(*(buf+n)+48);
	prmtop->nptra  = atoi(*(buf+n)+56);
	prmtop->natyp  = atoi(*(buf+n)+64);
	prmtop->nphb   = atoi(*(buf+n)+72);
	n++;
	prmtop->ifpert = atoi(*(buf+n));
	prmtop->nbper  = atoi(*(buf+n)+8);
	prmtop->ngper  = atoi(*(buf+n)+16);
	prmtop->ndper  = atoi(*(buf+n)+24);
	prmtop->mbper  = atoi(*(buf+n)+32);
	prmtop->mgper  = atoi(*(buf+n)+40);
	prmtop->mdper  = atoi(*(buf+n)+48);
	prmtop->ifbox  = atoi(*(buf+n)+56);
	prmtop->nmxrs  = atoi(*(buf+n)+64);
	prmtop->ifcap  = atoi(*(buf+n)+72);
	n++;
	prmtop->nextra = atoi(*(buf+n));

/*
	fprintf(stderr,"%8d",prmtop->natom);
	fprintf(stderr,"%8d",prmtop->ntypes);
	fprintf(stderr,"%8d",prmtop->nbonh);
	fprintf(stderr,"%8d",prmtop->mbona);
	fprintf(stderr,"%8d",prmtop->ntheth);
	fprintf(stderr,"%8d",prmtop->mtheta);
	fprintf(stderr,"%8d",prmtop->nphih);
	fprintf(stderr,"%8d",prmtop->mphia);
	fprintf(stderr,"%8d",prmtop->nhparm);
	fprintf(stderr,"%8d\n",prmtop->nparm);

	fprintf(stderr,"%8d",prmtop->nnb);
	fprintf(stderr,"%8d",prmtop->nres);
	fprintf(stderr,"%8d",prmtop->nbona);
	fprintf(stderr,"%8d",prmtop->ntheta);
	fprintf(stderr,"%8d",prmtop->nphia);
	fprintf(stderr,"%8d",prmtop->numbnd);
	fprintf(stderr,"%8d",prmtop->numang);
	fprintf(stderr,"%8d",prmtop->nptra);
	fprintf(stderr,"%8d",prmtop->natyp);
	fprintf(stderr,"%8d\n",prmtop->nphb);

	fprintf(stderr,"%8d",prmtop->ifpert);
	fprintf(stderr,"%8d",prmtop->nbper);
	fprintf(stderr,"%8d",prmtop->ngper);
	fprintf(stderr,"%8d",prmtop->ndper);
	fprintf(stderr,"%8d",prmtop->mbper);
	fprintf(stderr,"%8d",prmtop->mgper);
	fprintf(stderr,"%8d",prmtop->mdper);
	fprintf(stderr,"%8d",prmtop->ifbox);
	fprintf(stderr,"%8d",prmtop->nmxrs);
	fprintf(stderr,"%8d\n",prmtop->ifcap);

	fprintf(stderr,"%8d\n",prmtop->nextra);
*/
	
}

static char **get_label_name(int natom,int len,char **buf,int l_cnt,char *tag)
{
	char **atom_n,wname[MAXLINE];
	int n,i,k;
	
	n = skip_line(buf,l_cnt,tag)+2;
	atom_n = get_label_name_area(natom,len);
	for(i = k = 0; i < natom; i++){

		memset(wname,'\0',sizeof(*(atom_n+i)));
		strncpy(wname,*(buf+n)+k*(len-1),len-1);
		sscanf(wname,"%s",*(atom_n+i));
		k++;
		if(k == 20 && i != natom-1){
			n++;
			if(*(*(buf+n)) == '%'){
				(void)fprintf(stderr,"Invalid prmtop file\n");
				exit(0);
			}
			k = 0;
		}
	}
	return atom_n;
}

static char **get_label_name_area(int natom,int len)
{
	char **atom_n;
	int i;

	atom_n = (char **)calloc(natom,sizeof(char *));
	for(i = 0; i < natom; i++)
		*(atom_n+i) = (char *)calloc(len+1,sizeof(char ));
	return atom_n;
}

static void free_label_name_area(char **atom_n,int natom)
{
	int i;

	for(i = 0; i < natom; i++)
		free(*(atom_n+i));
	free(atom_n);
}

static double *get_double_value(int natom,char **buf,int l_cnt,char *tag)
{
	int i,n,k;
	double *chg;

	chg = (double *)calloc(natom,sizeof(double ));
	n = skip_line(buf,l_cnt,tag)+2;
	for(i = k = 0; i < natom; i++){
		*(chg+i) = atof(*(buf+n)+k*16);
		k++;
		if(k == 5 && i != natom-1){
			n++;
			if(*(*(buf+n)) == '%'){
				(void)fprintf(stderr,"Invalid prmtop file\n");
				exit(0);
			}
			k = 0;
		}
	}
	return chg;
}

static int *get_int_value(int natom,char **buf,int l_cnt,char *tag)
{
	int i,n,k,*a_idx;

	a_idx = (int *)calloc(natom,sizeof(int ));
	n = skip_line(buf,l_cnt,tag)+2;
	for(i = k = 0; i < natom; i++){
		*(a_idx+i) = atoi(*(buf+n)+k*8);
		k++;
		if(k == 10 && i != natom-1){
			n++;
			if(*(*(buf+n)) == '%'){
				(void)fprintf(stderr,"Invalid prmtop file\n");
				exit(0);
			}
			k = 0;
		}
	}
	return a_idx;
}

static int skip_line(char **buf,int l_cnt,char *word)
{
	int i,len;

	len = strlen(word);
	for(i = 0; i < l_cnt; i++){
		if(strncmp(*(buf+i),word,len) == 0)
			return i;
	}
	(void)fprintf(stderr,"Invalid word %s\n",word);
	exit(0);
}

static void convert_charge(double *chg,int natom)
{
	int i;

	for(i = 0; i < natom; i++)
		*(chg+i) /= CHRG_UNIT;
}

static void sum_numex(int *numex,int natom)
{
	int i,n,*w;

	w = (int *)calloc(natom,sizeof(int ));
	for(i = 0; i < natom; i++)
		*(w+i) = *(numex+i);
	n = 1;
	*numex = n;
	for(i = 1; i < natom; i++){
		n += *(w+i-1);
		*(numex+i) = n;
	}
/*
for(i = 0; i < natom; i++){
fprintf(stderr,"%d %d\n",*(w+i),*(numex+i));
getchar();
}
*/
	free(w);
}

static void get_index_area(PRMTOP *prmtop)
{
	prmtop->ibh = (int *)calloc(prmtop->nbonh,sizeof(int ));
	prmtop->jbh = (int *)calloc(prmtop->nbonh,sizeof(int ));
	prmtop->icbh = (int *)calloc(prmtop->nbonh,sizeof(int ));

	prmtop->ib = (int *)calloc(prmtop->nbona,sizeof(int ));
	prmtop->jb = (int *)calloc(prmtop->nbona,sizeof(int ));
	prmtop->icb = (int *)calloc(prmtop->nbona,sizeof(int ));

	prmtop->ith = (int *)calloc(prmtop->ntheth,sizeof(int ));
	prmtop->jth = (int *)calloc(prmtop->ntheth,sizeof(int ));
	prmtop->kth = (int *)calloc(prmtop->ntheth,sizeof(int ));
	prmtop->icth = (int *)calloc(prmtop->ntheth,sizeof(int ));

	prmtop->it = (int *)calloc(prmtop->ntheta,sizeof(int ));
	prmtop->jt = (int *)calloc(prmtop->ntheta,sizeof(int ));
	prmtop->kt = (int *)calloc(prmtop->ntheta,sizeof(int ));
	prmtop->ict = (int *)calloc(prmtop->ntheta,sizeof(int ));

	prmtop->iph = (int *)calloc(prmtop->nphih,sizeof(int ));
	prmtop->jph = (int *)calloc(prmtop->nphih,sizeof(int ));
	prmtop->kph = (int *)calloc(prmtop->nphih,sizeof(int ));
	prmtop->lph = (int *)calloc(prmtop->nphih,sizeof(int ));
	prmtop->icph = (int *)calloc(prmtop->nphih,sizeof(int ));

	prmtop->ip = (int *)calloc(prmtop->nphia,sizeof(int ));
	prmtop->jp = (int *)calloc(prmtop->nphia,sizeof(int ));
	prmtop->kp = (int *)calloc(prmtop->nphia,sizeof(int ));
	prmtop->lp = (int *)calloc(prmtop->nphia,sizeof(int ));
	prmtop->icp = (int *)calloc(prmtop->nphia,sizeof(int ));
}

static void get_index_values_3(int *ib,int *jb,int *icb,int natom,char **buf,int l_cnt,char *tag)
{
	int i,n,k;

	n = skip_line(buf,l_cnt,tag)+2;
	for(i = k = 0; i < natom; i++){
		*(ib+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(jb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(icb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
	}
}

static void get_index_values_4(int *ib,int *jb,int *kb,int *icb,int natom,char **buf,int l_cnt,char *tag)
{
	int i,k,n;

	n = skip_line(buf,l_cnt,tag)+2;
	for(i = k = 0; i < natom; i++){
		*(ib+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(jb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(kb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(icb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
	}
}

static void get_index_values_5(int *ib,int *jb,int *kb,int *lb,int *icb,int natom,char **buf,int l_cnt,char *tag)
{
	int i,k,n;

	n = skip_line(buf,l_cnt,tag)+2;
	for(i = k = 0; i < natom; i++){
		*(ib+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(jb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(kb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(lb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
		*(icb+i) = atoi(*(buf+n)+k*8)/3;
		check_new_line(i,natom,&k,&n,buf);
	}
}

static void check_new_line(int i,int natom,int *k,int *n,char **buf)
{
	(*k)++;
	if(*k == 10 && i != natom-1){
		(*n)++;
		if(*(*(buf+ *n)) == '%'){
			(void)fprintf(stderr,"Invalid prmtop file\n");
			exit(0);
		}
		(*k) = 0;
	}
}

static char *identify_term(PRMTOP *prmtop)
{
	char *term;
	int i,idx;

	term = get_init_char_array(prmtop->natom,'d');
	idx = 0;
	for(i = 0; i < prmtop->nspm; i++){
		idx += *(prmtop->nsp+i);
		*(term+idx-1) = 't';
	}
	return term;
}

static char *identify_term_from_bond(PRMTOP *prmtop)
{
	char *term;
	int *chain,chno,i,j,k,sp,ep,max_i;

	chain = get_init_int_array(prmtop->nres,0);
	chno = 0;
	for(i = 0; i < prmtop->nres; i++){
		if(*(chain+i) != 0)
			continue;
		chno++;
		check_connect_residue(prmtop,chain,i,chno);
	}
	for(k = 1; k <= chno; k++){
		max_i = -1;
		for(i = 0; i < prmtop->nres; i++){
			if(*(chain+i) == k){
				if(max_i < i)
					max_i = i;
			}
		}
		*(chain+max_i) = -1;
	}

	term = get_init_char_array(prmtop->natom+1,'d');
	for(i = 0; i < prmtop->nres; i++){
		if(*(chain+i) == -1){
			get_prmtop_sp_ep_idx(prmtop,i+1,i+1,&sp,&ep);
			*(term+ep) = 't';
		}
	}
	free(chain);
	return term;
}

static void check_connect_residue(PRMTOP *prmtop,int *chain,int no,int chno)
{
	int i,sp0,ep0,sp1,ep1;

	*(chain+no) = chno;
	get_prmtop_sp_ep_idx(prmtop,no+1,no+1,&sp0,&ep0);
	for(i = 0; i < prmtop->nres; i++){
		if(*(chain+i) != 0)
			continue;
		get_prmtop_sp_ep_idx(prmtop,i+1,i+1,&sp1,&ep1);
		if(check_connect_res(prmtop,sp0,ep0,sp1,ep1) == 1)
			check_connect_residue(prmtop,chain,i,chno);
	}
}


static int check_connect_res(PRMTOP *prmtop,int sp0,int ep0,int sp1,int ep1)
{
	int i13,i23,j13,j23,k;

	i13 = sp0;
	i23 = ep0;
	j13 = sp1;
	j23 = ep1;
	for(k = 0; k < prmtop->nbona; k++){
		if(*(prmtop->ib+k) >= i13 && *(prmtop->ib+k) <= i23 && *(prmtop->jb+k) >= j13 && *(prmtop->jb+k) <= j23)
				return 1;
		if(*(prmtop->ib+k) >= j13 && *(prmtop->ib+k) <= j23 && *(prmtop->jb+k) >= i13 && *(prmtop->jb+k) <= i23)
				return 1;
	}
	return 0;
}

void get_prmtop_sp_ep_idx(PRMTOP *prmtop,int res_0,int res_1,int *idx_0,int *idx_1)
{
        *idx_0 = *(prmtop->ipres+res_0-1)-1;
        if(res_1 != prmtop->nres)
                *idx_1 = *(prmtop->ipres+res_1)-2;
        else
                *idx_1 = prmtop->natom-1;
}
