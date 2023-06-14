CC := gcc
LIBS := -lm -lgsl -lgslcblas
CFLAGS := -mcmodel=medium

$(shell mkdir -p bin/)

all : nonlocal_all CPpackage AddCore AddCoreXSF
nonlocal_all: nonlocal nonlocal_pw nonlocal_merge_pw

nonlocal : src/nonlocal.c
	$(CC) $(CFLAGS) $< -o bin/$@ -O3 $(LIBS) -lfftw3 

nonlocal_pw : src/nonlocal_by_pw.c
	$(CC) $(CFLAGS) $< -o bin/$@ -O3 $(LIBS)

nonlocal_merge_pw : src/nonlocal_merge_pw.c
	$(CC) $(CFLAGS) $< -o bin/$@ $(LIBS)

CPpackage : src/CPpackage.c
	$(CC) $(CFLAGS) $< -o bin/$@ -O3 $(LIBS) -lxc

AddCore : src/AddCore.c
	$(CC) $(CFLAGS) $< -o bin/$@ $(LIBS)

AddCoreXSF : src/AddCoreXSF.c
	$(CC) $(CFLAGS) $< -o bin/$@ $(LIBS)

.PHONY : clean
clean :
	rm bin/*
