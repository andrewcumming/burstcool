# Makefile for cool.cc 

# common directories
CDIR = c
ODIR = o

# compiler
CC = g++-9
FORT = gfortran
CFLAGS = -lm -lgsl -lgslcblas -L/usr/local/lib -I/usr/local/include

# main code
OBJS = $(ODIR)/burst.o $(ODIR)/root.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o $(ODIR)/vector.o $(ODIR)/timer.o $(ODIR)/ns.o $(ODIR)/condegin13.o

burstcool : $(OBJS) $(ODIR)/burstcool.o
	$(CC) -o burstcool $(OBJS) $(ODIR)/burstcool.o $(CFLAGS) 

$(ODIR)/burstcool.o : $(CDIR)/burstcool.cc
		$(CC) -c $(CDIR)/burstcool.cc -o $(ODIR)/burstcool.o $(CFLAGS)

$(ODIR)/burst.o : $(CDIR)/burst.cc
		$(CC) -c $(CDIR)/burst.cc -o $(ODIR)/burst.o $(CFLAGS)

$(ODIR)/condegin13.o : $(CDIR)/condegin13.f
	$(FORT) -c $(CDIR)/condegin13.f -o $(ODIR)/condegin13.o $(CFLAGS)

$(ODIR)/makegrid.o : $(CDIR)/makegrid.cc
	$(CC) -c $(CDIR)/makegrid.cc -o $(ODIR)/makegrid.o $(CFLAGS)

# compile routines from the common directory

$(ODIR)/root.o : $(CDIR)/root.c
	$(CC) -c $(CDIR)/root.c -o $(ODIR)/root.o $(CFLAGS)

$(ODIR)/eos.o : $(CDIR)/eos.cc $(ODIR)/condegin13.o
	$(CC) -c $(CDIR)/eos.cc -o $(ODIR)/eos.o $(CFLAGS)

$(ODIR)/odeint.o : $(CDIR)/odeint.cc
	$(CC) -c $(CDIR)/odeint.cc -o $(ODIR)/odeint.o $(CFLAGS)

$(ODIR)/spline.o : $(CDIR)/spline.cc
	$(CC) -c $(CDIR)/spline.cc -o $(ODIR)/spline.o $(CFLAGS)

$(ODIR)/vector.o : $(CDIR)/vector.cc
	$(CC) -c $(CDIR)/vector.cc -o $(ODIR)/vector.o $(CFLAGS)

$(ODIR)/timer.o : $(CDIR)/timer.c
	$(CC) -c $(CDIR)/timer.c -o $(ODIR)/timer.o $(CFLAGS)

$(ODIR)/ns.o : $(CDIR)/ns.c
	$(CC) -c $(CDIR)/ns.c -o $(ODIR)/ns.o $(CFLAGS)


# clean up

clean:
	rm -f $(ODIR)/*.o

grid_sorty: out/grid
	sort out/grid -k1,1 -g >out/grid_sorty

movie:
	ffmpeg -qscale 1 -r 20 -b 9600 -i png/%3d.png movie.mp4

movie2:
	ffmpeg -r 30 -i png/%3d.png  -vcodec libx264 -y -an movie.mp4

cleanpng:
	rm -f png/*.png

