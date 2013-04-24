#==== These first flags are likely to be quite machine specific. ====#
SRCDIR = /home/stefan/him/source/

#==== Linux/gcc directives ====#
 CC=mpicc
 CFLAGS = -O3
 LDFLAGS = -I`pwd` -I$(SRCDIR) -I/usr/include 
 LIBFLAGS = -lm -lnetcdf -L/usr/lib  

#==== SGI Origin3000 directives ====#
# CFLAGS = -O2
# LDFLAGS = -I`pwd` -I$(SRCDIR) -I/usr/local/include
# LIBFLAGS = -lm -lnetcdf -lsma -lmpi -L/usr/local/lib

#==== SGI Origin 3000 debugging directives ====#
# CFLAGS = -g2
# LDFLAGS = -I`pwd` -I$(SRCDIR) -I/usr/local/include
# TOTALVIEW = -L/opt/totalview/4.1.0-3/irix6-mips/lib -ldbfork_n64
# LIBFLAGS = -lm -lnetcdf -lsma -lmpi -L/usr/local/lib $(TOTALVIEW)


#==== The remainder of the file is probably not machine specific. ====#

OBJS1 = initialize.o initialize_output.o surface_forcing.o \
      barotropic.o continuity.o CoriolisAdv.o diabatic_driver.o \
      diabatic_entrain.o diagnostics.o eqn_of_state.o \
      fit_compressibility.o HIM_driver.o HIM_field_output.o \
      HIM_io.o HIM_parallel.o HIM_parser.o HIM_restart.o \
      HIM_sum_output.o hor_visc.o mixed_layer.o PointAccel.o \
      PressureForce.o set_metrics.o sponge.o thickness_diffuse.o \
      timetype.o tracer.o vertvisc.o

OBJS = HIM.o $(OBJS1)

OBJS_CHECKPAR = HIM_checkparallel.o $(OBJS1)

HEADERS = Makefile init.h $(SRCDIR)/HIM_io.h $(SRCDIR)/HIM.h \
      $(SRCDIR)/hor_visc.h $(SRCDIR)/metrics.h

HIMTOCDFSRC = $(SRCDIR)/HIMtocdf.c $(SRCDIR)/HIMtocdf_read.c \
      $(SRCDIR)/HIMtocdf_poisson.c $(SRCDIR)/set_metrics.c \
      $(SRCDIR)/HIM_io.c $(SRCDIR)/HIM_par_IO.c



HIM: $(HEADERS) $(OBJS) Makefile
	$(CC) -o $@ $(OBJS) $(LIBFLAGS) $(LDFLAGS)

HIM_checkparallel: $(HEADERS) $(OBJS_CHECKPAR) Makefile
	$(CC) -o $@ $(OBJS_CHECKPAR) $(LIBFLAGS) $(LDFLAGS)

HIMtocdf: $(HEADERS) $(HIMTOCDFSRC) Makefile
	$(CC) -DSERIAL_IO_CODE $(HIMTOCDFSRC) -o $@ $(LIBFLAGS) $(LDFLAGS)
	-rm *.u
	-rm *.o

param_suggest: $(SRCDIR)/param_suggest.c init.h
	@$(CC) -o $@ $(SRCDIR)/param_suggest.c $(LIBFLAGS) $(LDFLAGS)
	@param_suggest
	@rm param_suggest


#==== Object files follow. ====#

barotropic.o: $(SRCDIR)/barotropic.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/barotropic.c $(LDFLAGS)

continuity.o: $(SRCDIR)/continuity.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/continuity.c $(LDFLAGS)

continuity_FCT.o: $(SRCDIR)/continuity_FCT.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/continuity_FCT.c $(LDFLAGS)

continuity_mpdata.o: $(SRCDIR)/continuity_mpdata.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/continuity_mpdata.c $(LDFLAGS)

CoriolisAdv.o: $(SRCDIR)/CoriolisAdv.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/CoriolisAdv.c $(LDFLAGS)

diabatic_driver.o: $(SRCDIR)/diabatic_driver.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/diabatic_driver.c $(LDFLAGS)

diabatic_entrain.o: $(SRCDIR)/diabatic_entrain.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/diabatic_entrain.c $(LDFLAGS)

diagnostics.o: $(SRCDIR)/diagnostics.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/diagnostics.c $(LDFLAGS)

eqn_of_state.o: $(SRCDIR)/eqn_of_state.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/eqn_of_state.c $(LDFLAGS)

eqn_of_state_linear.o: $(SRCDIR)/eqn_of_state_linear.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/eqn_of_state_linear.c $(LDFLAGS)

eqn_of_state_UNESCO.o: $(SRCDIR)/eqn_of_state_UNESCO.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/eqn_of_state_UNESCO.c $(LDFLAGS)

fit_compressibility.o: $(SRCDIR)/fit_compressibility.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/fit_compressibility.c $(LDFLAGS)

HIM.o: $(SRCDIR)/HIM.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM.c $(LDFLAGS)

HIM_checkparallel.o: $(SRCDIR)/HIM_checkparallel.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_checkparallel.c $(LDFLAGS)

HIM_driver.o: $(SRCDIR)/HIM_driver.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_driver.c $(LDFLAGS)

HIM_field_output.o: $(SRCDIR)/HIM_field_output.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_field_output.c $(LDFLAGS)

HIM_io.o: $(SRCDIR)/HIM_io.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_io.c $(LDFLAGS)

HIM_parallel.o: $(SRCDIR)/HIM_parallel.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_parallel.c $(LDFLAGS)

HIM_parser.o: $(SRCDIR)/HIM_parser.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_parser.c $(LDFLAGS)

HIM_restart.o: $(SRCDIR)/HIM_restart.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_restart.c $(LDFLAGS)

HIM_sum_output.o: $(SRCDIR)/HIM_sum_output.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/HIM_sum_output.c $(LDFLAGS)

hor_visc.o: $(SRCDIR)/hor_visc.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/hor_visc.c $(LDFLAGS)

initialize.o: initialize.c $(HEADERS)
	$(CC) $(CFLAGS) -c initialize.c $(LDFLAGS)

initialize_output.o: initialize_output.c $(HEADERS)
	$(CC) $(CFLAGS) -c initialize_output.c $(LDFLAGS)

mixed_layer.o: $(SRCDIR)/mixed_layer.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/mixed_layer.c $(LDFLAGS)

PointAccel.o: $(SRCDIR)/PointAccel.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/PointAccel.c $(LDFLAGS)

PressureForce.o: $(SRCDIR)/PressureForce.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/PressureForce.c $(LDFLAGS)

set_metrics.o: $(SRCDIR)/set_metrics.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/set_metrics.c $(LDFLAGS)

sponge.o: $(SRCDIR)/sponge.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/sponge.c $(LDFLAGS)

surface_forcing.o: $(SRCDIR)/surface_forcing.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/surface_forcing.c $(LDFLAGS)

thickness_diffuse.o: $(SRCDIR)/thickness_diffuse.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/thickness_diffuse.c $(LDFLAGS)

timetype.o: $(SRCDIR)/timetype.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/timetype.c $(LDFLAGS)

tracer.o: $(SRCDIR)/tracer.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/tracer.c $(LDFLAGS)

vertvisc.o: vertvisc.c $(HEADERS)
	$(CC) $(CFLAGS) -c vertvisc.c $(LDFLAGS)

USER_tracer_example.o: $(SRCDIR)/USER_tracer_example.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCDIR)/USER_tracer_example.c $(LDFLAGS)
