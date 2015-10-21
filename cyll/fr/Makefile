all : f_par.bin f_ort.bin

f_par.bin : Grad_parallel_SDSS.f
	gfortran -O2 Grad_parallel_SDSS.f -o f_par.bin


f_ort.bin : Grad_ortogonal_SDSS.f
	gfortran -O2 Grad_ortogonal_SDSS.f -o f_ort.bin


