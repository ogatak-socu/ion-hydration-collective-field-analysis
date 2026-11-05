#
# makefile
#
.SUFFIXES : .c .o

DEST    = .
CC      = gcc
CFLAGS   = -I.  -m64 -pg -g -O -c 
LIBS    = -lm
LIB_FLAG =
LMD     =

HEADER	= cal_vector.h

OBJS    = rd_options.o file_open.o rd_file_tools.o array_func.o get_box_size.o get_grid_system.o cal_coord_func.o assign_spherical_system2cube_system.o \
	shift_traj_center.o check_water_hydrogen_all.o  rd_amber_prmtop.o rd_amber_restart_file.o rd_amber_traj.o get_delta_xyz_values_array.o cal_charge_func.o \
	assign_atoms2points.o

CAL_VECTOR	= cal_vector.o cal_vector_mdcrd.o $(OBJS)

ACCUMULATE_VECTOR	= accumulate_vector.o accumulate_vector_all.o $(OBJS)

CAL_CHARGE	= cal_charge.o cal_charge_mdcrd.o $(OBJS)

ACCUMULATE_CHARGE	= accumulate_charge.o accumulate_charge_all.o $(OBJS)

all : cal_vector accumulate_vector cal_charge accumulate_charge

.PHONY : all

cal_vector : $(CAL_VECTOR) $(HEADER)
	$(CC) $(CAL_VECTOR) -o $@ $(LIBS)

accumulate_vector : $(ACCUMULATE_VECTOR) $(HEADER)
	$(CC) $(ACCUMULATE_VECTOR) -o $@ $(LIBS)

cal_charge : $(CAL_CHARGE) $(HEADER)
	$(CC) $(CAL_CHARGE) -o $@ $(LIBS)

accumulate_charge : $(ACCUMULATE_CHARGE) $(HEADER)
	$(CC) $(ACCUMULATE_CHARGE) -o $@ $(LIBS)

.c.o   : $(OBJS) $(HEADER)
	$(CC) $(CFLAGS) $<

clean  :
	rm $(OBJS)

.o   : $(HEADER)

cal_vector.o: $(HEADER)
cal_vector_mdcrd.o: $(HEADER)
accumulate_vector.o: $(HEADER)
accumulate_vecto_all.o: $(HEADER)
cal_charge.o: $(HEADER)
cal_charge_mdcrd.o: $(HEADER)
accumulate_charge.o: $(HEADER)
accumulate_charge_all.o: $(HEADER)
rd_amber_prmtop.o: $(HEADER)
rd_amber_restart_file.o: $(HEADER)
rd_amber_traj.o: $(HEADER)
rd_options.o: $(HEADER)
file_open.o: $(HEADER)
rd_file_tools.o: $(HEADER)
array_func.o: $(HEADER)
get_box_size.o: $(HEADER)
get_grid_system.o: $(HEADER)
cal_coord_func.o: $(HEADER)
assign_spherical_system2cube_system.o: $(HEADER)
check_water_hydrogen_all.o: $(HEADER)
shift_traj_center.o: $(HEADER)
get_delta_xyz_values_array.o: $(HEADER)
cal_charge_func.o: $(HEADER)
assign_atoms2points.o: $(HEADER)
