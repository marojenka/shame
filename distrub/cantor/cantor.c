#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

double R = 1.0, R2=0.1;
double d = 2.2, k = 0.01, q=1;


int cantor_square(	double x, double y, double z, // Координаты центра нового кубика
					double p, double h, double lim,  // вероятнось, половина стороны и минимальный размер 
					FILE* out, int ji ) 	{ // выходной файл и номер точки

	double r;
	if( lim >= h ) {
		x += h * (( (double) rand() / RAND_MAX) * 2 - 1) ;
		y += h * (( (double) rand() / RAND_MAX) * 2 - 1) ;
		z += h * (( (double) rand() / RAND_MAX) * 2 - 1) ;		
        if( (sqrt(x*x + y*y + z*z) < R) & 
            ( ((double) rand() / RAND_MAX) < q ) ) {
		    fprintf(out, "%e %e %e\n", x, y, z);
		    return 1;
        } else 
            return 0;
	}
	h /= 2;
	int i=0;
	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x+h,y+h,z+h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x+h,y+h,z-h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x+h,y-h,z-h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x+h,y-h,z+h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x-h,y-h,z-h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x-h,y+h,z-h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x-h,y-h,z+h, p,h,lim, out, i); 

	if ( ((double) rand() / RAND_MAX) <= p ) 
	 	i	+=	cantor_square(x-h,y+h,z+h, p,h,lim, out, i); 
	
	return i;
}

int 
get_cantor( char* output ) { 
	double p, r, x, y, z;
	FILE *out;
	int i;
	
	if( (out = fopen(output, "w")) == NULL ) {
		printf("Can't create file %s \n", output);
		return 1; 
	}

	p = pow(2, d - 3); 
	printf("p = %lf\n", p);
	//i = cantor_square(0, 0, 0, p, R, k, out, 0); 
	double h;
	x = y = z = 0;
	h = R;
    k = R2;
	i	+=	cantor_square(x,y,z, p, h, k, out, 0); 
		
	printf("N = %d\n", i);
	fclose(out);
	return i;
}


int main(int argc, char** argv) {
	char optString[] = {"R:s:r:d:q:"};
    char * data = "cantor.dat";
	FILE *LOG; 

	int opt; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
        switch( opt ) {
			case 's': // Название файла с данным
                data = optarg; 
                break;
                             
            case 'R': // Границы стки 
				sscanf(optarg, "%lf", &R); 
                break;
			
			case 'r':
				sscanf(optarg, "%lf", &R2); 
                break;
			
			case 'd':
				sscanf(optarg, "%lf", &d); 
                break;
			
			case 'q':
				sscanf(optarg, "%lf", &q); 
                break;
			
            default:
				printf("What?");
                /* сюда на самом деле попасть невозможно. */
				break;
        }
        opt = getopt( argc, argv, optString );
    }
	
	printf("generate %s	D=%1.2lf\n", data, d);
	printf("Radius = %lf, R2 = %lf\n", R, R2 );
	srand(time(NULL));

     struct timeval time; 
     gettimeofday(&time,NULL);

     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
     srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

     // The trouble here is that the seed will repeat every
     // 24 days or so.

     // If you use 100 (rather than 1000) the seed repeats every 248 days.

     // Do not make the MISTAKE of using just the tv_usec
     // This will mean your seed repeats every second.
	get_cantor(data);	

	LOG = fopen("distrub.log", "a");
	fprintf(LOG, "generate %s	D=%1.2lf\n", data, d);
	fprintf(LOG, "Radius = %lf, R2 = %lf\n", R, R2 );
	fclose(LOG);
	return 0;
}

