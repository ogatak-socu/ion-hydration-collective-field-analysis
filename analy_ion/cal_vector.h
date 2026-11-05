/*
 *
*/
#define MAXLINE 1024

#define EU	332.0685563
#define KB	0.001986231278448	/*    1.3806504E-23*1.43862E+20  */
#define C2DEBYE 4.803204

#define GRID_SIZE 5.0
#define SC_VDW	2.0
#define SC_ELE	1.2
#define	PROBE	1.4


typedef struct sw_t_t{
		char	*prmtop;
		char	*inpcrd;

		char	*system_size;

		char	*spherical_grid_points_file;

		int	target_atom_index;

		int	tr_cnt;
		int	f_cnt;
		char	**mdcrd;
		char	**mdfrc;

		int	rn;
		double	delta_d;

		double	rc;
		double	alpha;

		double grid_size;

		char    *dip_vector_prefix;
		char	*frc_vector_prefix;
		char	*bck_vector_prefix;
		char	*fwd_vector_prefix;

                char    *charge_prefix;

		int	file_no;

		int	number_of_vectors;

	}SW;

typedef struct xyz_t_t{
		char *aname;
		double x;
		double y;
		double z;
	}XYZ_T;

typedef struct grid_t_t{
                double r;
		int npoint;
                XYZ_T *p;
                XYZ_T *dxyz_t;
		unsigned long int    *n;
		double *weight;
		double *charge;
        }GRID_T;

typedef struct poly_a_t{
                GRID_T *grid_d;		/* Dipole	*/
                GRID_T *grid_f;		/* Force	*/
                GRID_T *grid_m;		/* Forward	*/
                GRID_T *grid_b;		/* Backward	*/
        }POLY_A;

typedef struct wgrid_t_t{
                int layer;
                int idx;
                double x;
                double y;
                double z;
                unsigned long int n;
                double weight;
        }WGRID_T;

typedef struct charge_t_t{
                GRID_T *grid_p;		/* Positive charge	*/
                GRID_T *grid_n;		/* Negative charge	*/
        }CHARGE_T;

typedef struct wcharge_t_t{
                int layer;
                int idx;
                double charge_p;
                double charge_n;
                double weight_p;
                double weight_n;
                unsigned long int n_p;
                unsigned long int n_n;
        }WCHARGE_T;

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


typedef struct traj_t{
		int n;
		XYZ_T *box;
		XYZ_T *mdcrd;
	}TRAJ;

typedef struct rst_t{
		char *title;
		char *natom_s;
		int natom;
		double t;
		XYZ_T *r;
		XYZ_T *v;
		XYZ_T box;
		XYZ_T box_a;
	}RST;

typedef struct xyz_i_t{
		int x;
		int y;
		int z;
		double dx;
		double dy;
		double dz;
        }XYZ_I;

typedef struct poly_p_t{
                int n;
                int *layer;
                int *idx;
        }POLY_P;

typedef struct range_p_t{
                double weight;
                int idx;
                int layer;
        }RANGE_P;


void cal_vector_mdcrd(PRMTOP *,SW *);


SW *rd_options(char **,int );
void disp_parameter(SW *);

PRMTOP *rd_amber_prmtop(char *);
void get_prmtop_sp_ep_idx(PRMTOP *,int ,int ,int *,int *);

RST *rd_amber_restart_file(char *);
void free_amber_rst(RST *);


FILE *rd_file_open(char *);
FILE *rd_binary_file_open(char *);
FILE *wt_file_open(char *);
FILE *wt_binary_file_open(char *);
FILE *append_file_open(char *);
void rd_file_close(char *,FILE *);
void wt_file_close(char *,FILE *);
char **rd_text_file(char *,int *);
char **rd_text_file_without_cr(char *,int *);

int get_text_line_cnt(char *);
char **get_text_area_fix(int );
char **get_text_area_flexible(int ,int );
void free_text_area(char **,int );

char *get_init_char_array(int ,char );
void set_init_char_array(char *,int ,char );
char **get_char_array_2d(int ,int );
char **get_init_char_array_2d(int ,int ,char );
void set_init_char_array_2d(char **,int ,int ,char );
void cp_char_array(char *a,char *,int );
void free_char_array_2d(char **,int );
int *get_init_int_array(int ,int );
void set_init_int_array(int *,int ,int );
int **get_int_array_2d(int ,int );
int **get_init_int_array_2d(int ,int ,int );
void set_init_int_array_2d(int **,int ,int ,int );
void cp_int_array(int *,int *,int );
void free_int_array_2d(int **,int );
unsigned long int *get_init_unsigned_long_int_array(int ,unsigned long int );
void set_init_unsigned_long_int_array(unsigned long int *,int ,unsigned long int );
unsigned long int **get_unsigned_long_int_array_2d(int ,int );
unsigned long int **get_init_unsigned_long_int_array_2d(int ,int ,unsigned long int );
void set_init_unsigned_long_int_array_2d(unsigned long int **,int ,int ,unsigned long int );
void cp_unsigned_long_int_array(unsigned long int *,unsigned long int *,int );
void free_unsigned_long_int_array_2d(unsigned long int **,int );
double *get_init_double_array(int ,double );
void set_init_double_array(double *,int ,double );
double **get_double_array_2d(int ,int );
double **get_init_double_array_2d(int ,int ,double );
void set_init_double_array_2d(double **,int ,int ,double );
void free_double_array_2d(double **,int );
void cp_double_array(double *,double *,int );

XYZ_I *get_box_size(SW *,int *);

POLY_A *get_grid_system(int ,SW *);
GRID_T *get_init_grid_system_file(int ,double ,char *);
GRID_T *cp_init_grid_t(GRID_T *,int ,SW *,int );
void mk_link_grid_coord(GRID_T *,GRID_T *,int );
void free_grid_t_area(GRID_T *,int );

XYZ_T *get_init_xyz_t(int );
void scaling_xyz(XYZ_T *,int ,double );
double distance_xyz(XYZ_T *,XYZ_T *);
double distance_2_xyz(XYZ_T *,XYZ_T *);
double distance0_xyz(XYZ_T *);
double distance0_2_xyz(XYZ_T *);
double vec_angle_xyz(XYZ_T *,XYZ_T *,XYZ_T *);
double vec_angle0_xyz(XYZ_T *,XYZ_T *);
void cal_center_xyz(XYZ_T *,int ,XYZ_T *);
void zero_vector_xyz(XYZ_T *);
void x_axis_vector_xyz(XYZ_T *);
void y_axis_vector_xyz(XYZ_T *);
void z_axis_vector_xyz(XYZ_T *);
void unit_vector_xyz(XYZ_T *,XYZ_T *);
void unit_vector_self_xyz(XYZ_T *);
void pos_vector_xyz(XYZ_T *,XYZ_T *,XYZ_T *);
void pos_unit_vector_xyz(XYZ_T *,XYZ_T *,XYZ_T *);
void cp_vector_xyz(XYZ_T *,XYZ_T *);
void add_vector_xyz(XYZ_T *,XYZ_T *);
void add_vector2_xyz(XYZ_T *,XYZ_T *,XYZ_T *);
void del_vector_xyz(XYZ_T *,XYZ_T *);
void del_vector2_xyz(XYZ_T *,XYZ_T *,XYZ_T *);
void center_vector_xyz(XYZ_T *,XYZ_T *,XYZ_T *);
void scale_vector_xyz(XYZ_T *,double );
void multi_scalar2vector_xyz(XYZ_T *,XYZ_T *,double );
void divid_vector_xyz(XYZ_T *,double );
void divid_scalar2vector_xyz(XYZ_T *,XYZ_T *,double );
double inner_product_xyz(XYZ_T *,XYZ_T *);
void outer_product_xyz(XYZ_T *,XYZ_T *,XYZ_T *);

POLY_P *assign_spherical_system2cube_system(GRID_T *,XYZ_I *,double ,int ,int ,int *);
int cal_total_grid_index(SW *,int );
int g_ijk2idx(XYZ_I *,int ,int ,int );
int g_xyz2idx(XYZ_T *pos,XYZ_I *,XYZ_T *,double );
int g_box_range(int ,int ,int ,XYZ_I *);
void g_xyz2ijk_mass(XYZ_T *,int *,int *,int *,XYZ_T *,double );


TRAJ *rd_amber_traj(char *,int ,int ,int ,int *);
TRAJ *rd_amber_traj_other(char *,int ,int ,int ,int *);
void free_traj_area(TRAJ *,int );
void free_traj_area_other(TRAJ *,int );
TRAJ *get_init_traj_area(int );
void cp_traj_mdcrd(XYZ_T *,XYZ_T *,int );


void shift_traj_center(TRAJ *,XYZ_T *,int );
void shift_traj_center_res(TRAJ *,XYZ_T *,PRMTOP *);

void check_water_hydrogen_all(XYZ_T *,XYZ_T *,PRMTOP *);

void get_delta_xyz_values_array(XYZ_T *,XYZ_T *,XYZ_T *,XYZ_T *,PRMTOP *);

void accumulate_vector_all(SW *);

void cal_charge_mdcrd(PRMTOP *,SW *);
CHARGE_T *get_init_charge_t(int ,SW *);

void accumulate_charge_all(SW *);


RANGE_P *get_init_range_p(int );
void set_init_range_p(RANGE_P *,int );
void get_range_point_polygon_grid_normal(RANGE_P *,GRID_T *,XYZ_T *,POLY_P *,XYZ_I *,double ,int ,double ,double );
