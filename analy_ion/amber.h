#define EU	332.0685563
#define KB	0.001986231278448	/*    1.3806504E-23*1.43862E+20  */

#define SC_VDW	2.0
#define SC_ELE	1.2

#define	PROBE	1.4

typedef struct prmtop_t{
		int natom;
		int ntypes;
		int nbonh;
		int mbona;
		int ntheth;
		int mtheta;
		int nphih;
		int mphia;
		int nhparm;
		int nparm;
		int nnb;
		int nres;
		int nbona;
		int ntheta;
		int nphia;
		int numbnd;
		int numang;
		int nptra;
		int natyp;
		int nphb;
		int ifpert;
		int nbper;
		int ngper;
		int ndper;
		int mbper;
		int mgper;
		int mdper;
		int ifbox;
		int nmxrs;
		int ifcap;
		int nextra;
		char **igraph;
		double *chrg;
		double *amass;
		int *iac;
		int *numex;
		int *ico;
		char **labres;
		int *ipres;
		double *rk;
		double *req;
		double *tk;
		double *teq;
		double *pk;
		double *pn;
		double *phase;
		double *one_scee;
		double *one_scnb;
		double *solty;
		double *cn1;
		double *cn2;

		int *ibh;
		int *jbh;
		int *icbh;

		int *ib;
		int *jb;
		int *icb;

		int *ith;
		int *jth;
		int *kth;
		int *icth;

		int *it;
		int *jt;
		int *kt;
		int *ict;

		int *iph;
		int *jph;
		int *kph;
		int *lph;
		int *icph;

		int *ip;
		int *jp;
		int *kp;
		int *lp;
		int *icp;

		int *natex;
		double *asol;
		double *bsol;
		double *hbcut;
		
		char **isymbol;
		char **itree;
		int *join;
		int *irotat;
		int iptres;
		int nspm;
		int nspsol;
		int *nsp;
		double beta;
		double *box;
		char *term;
	}PRMTOP;

typedef struct xyz_t_t{
		char *aname;
		double x;
		double y;
		double z;
	}XYZ_T;

typedef struct traj_t{
		int n;
		XYZ_T *box;
		XYZ_T *mdcrd;
	}TRAJ;

typedef struct rst_t_t{
		char *title;
		char *natom_s;
		int natom;
		double t;
		ATOM_T *r;
		ATOM_T *v;
		ATOM_T box;
		ATOM_T box_a;
	}RST_T;
