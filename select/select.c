#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#define MAX(a, b)       ((a) > (b) ? (a) : (b))
#define MIN(a, b)       ((a) < (b) ? (a) : (b))

char *f_in_name, *f_out_name; 
int N; // how many points
double q = 1;

void
read_data() {
    FILE *f_in, *f_out; 
    char *line;
    size_t len=0; 
    ssize_t read;
    int i;
    
    f_in = fopen(f_in_name, "r");
    f_out= fopen(f_out_name, "w");

	if( f_in == NULL ) {
		printf("Can't find file %s.\n", f_in_name);
		exit(10);
	}
    
    double toss; 
    char string[1000]; 
    while ((read = getline(&line, &len, f_in)) != -1) {
		toss = ((double) rand() / RAND_MAX);
        if( toss < q ) 
            fprintf(f_out, "%s", line);
    }
    if( line) 
        free(line); 
    fclose(f_in);
    fclose(f_out);
}

int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:q:o:"};
	int opt, got_out = 0; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'd':
				f_in_name = optarg; 
				break;
			case 'q': 
				sscanf(optarg, "%lf", &q); 
				break;
			case 'o':
				got_out = 1;
				f_out_name = optarg; 
				break;

			default:
				printf("What? %d\n", opt);
				break;
    	}
	    opt = getopt( argc, argv, optString );
	}
	if( got_out != 1 ) { 
        char tmp[100];
        f_out_name = tmp;
	    sprintf(f_out_name, "%s_q%d.dat", f_in_name, (int) q*100);
    }
	read_data();
	return(0);	
}

