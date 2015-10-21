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

int SELECT_COORDINATES = 0; // CARTESIAN

char *name; 
int N; // how many points
double **p;  // Array of points. p[i][N]
double **P;  // Array of points. p[i][N]
double *R_par, *R_ort;     // Rmax ofc
int    *R_par_ind, *R_ort_ind; // Rmax compared to grid. grid[R_par_ind] < R_par[i] < grid[R_par_ind + 1] ???
double ***mass; 
int ***par, ***ort; 
//int **par1, **ort1; 
//int **par2, **ort2; 

double h = 1;

int grid_n = 100;
double *grid, grid_h1, grid_R[2] = {0.1, 100}; // grid itseld, h+1 and limits. 

double A[2] = {0, 90}, B[2] = {0, 90}, R[2] = {0, 100};
double q = 1;

double
Angle(double x1, double y1, double z1, double x2, double y2, double z2) {
	double r1, r2;
	r1=sqrt(x1*x1+y1*y1+z1*z1);
	r2=sqrt(x2*x2+y2*y2+z2*z2);
	if( (r1==0) || (r2==0))
		return 0;
	return acos((x1*x2+y1*y2+z1*z2) / r1 / r2 );
}

void
spherical_coordinats(int n, double *x, double *y, double *z, double *a, double *b, double *r) {
	int i; 
	for( i=0; i<n; i++ ) {
		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        a[i] = atan2(y[i], x[i]) ;  
		if( a[i] < 0 ) 
			a[i] += 2*M_PI; 
		if( x[i] == 0 & y[i] == 0 ) {
			if( z[i] > 0 ) b[i] = +M_PI/2.;
			if( z[i] < 0 ) b[i] = -M_PI/2.;
		}
        b[i] = asin(z[i]/r[i]) ; 
	}
}	

void
cartesian_coordinats(int n, double *x, double *y, double *z, double *a, double *b, double *r) {
	int i; 
	for( i=0; i<n; i++ ) {
		x[i] = r[i] * cos(b[i]) * cos(a[i]); 
		y[i] = r[i] * cos(b[i]) * sin(a[i]); 
		z[i] = r[i] * sin(b[i]); 
	}
}	

int
inside_sp(double a,double b,double r) {
	if(r < R[0]) { return 1; }
	if(r > R[1]) { return 2; }
	if(a < A[0]) { return 3; }
	if(a > A[1]) { return 4; }
	if(b < B[0]) { return 5; }
	if(b > B[1]) { return 6; }
	return 0;
}

int
inside_ca(double x,double y, double z) {
	double a,b,r; 
	spherical_coordinats(1,&x,&y,&z,&a,&b,&r);
	return inside_sp(a,b,r);
}

void
cross_product(double *a, double *b, double *res) {
	res[0] = a[1]*b[2]-a[2]*b[1];
	res[1] = a[2]*b[0]-a[0]*b[2];
	res[2] = a[0]*b[1]-a[1]*b[0];
}


rmax(double x, double y, double z) {
	double rm, a, b, r; 
	double da = 2*M_PI, db = 2*M_PI; 
	int i; 
	spherical_coordinats(1,&x,&y,&z,&a,&b,&r);
	if(inside_sp(a,b,r) != 0) return(-1);
	da = 2*M_PI; db = 2*M_PI;
	if( A[0] != A[1] - 2*M_PI) {
		da = MIN(a - A[0], A[1] - a); 
	}
	db = b - B[0];
	db = MIN(db, B[1] - b);
	rm = MIN(r - R[0], R[1] - r); 
	if( da <= M_PI/2 ) {
		rm = MIN(rm, r * cos(b) * sin( da )) ; 
	}
	if( db <= M_PI/2 ) {
		rm = MIN(rm, r * sin(db) ); 
	}
	return( rm );
}

void
count() {
	int i, j, k;
	double rd, rp, x, rh;
	double dx, dy, dz; 
	double vec[3]; 
	int sign; 

	vec[2] = 0;		
	for(i=0; i<N; i++) {
	printf("\ri=%d", i);
		vec[0] = 1. / sqrt(1 + p[i][0]*p[i][0] / (p[i][1]*p[i][1])); 
		vec[1] = - p[i][0]/p[i][1] * vec[0]; 
		for(j=0; j<N; j++) {
			if(i==j) continue;
			dx=p[i][0]-p[j][0];
			dy=p[i][1]-p[j][1];
			dz=p[i][2]-p[j][2];
			rd = sqrt(dx*dx+dy*dy+dz*dz); 

		    rp = (p[i][0]*p[j][0]+p[i][1]*p[j][1]+p[i][2]*p[j][2])/P[i][2]; 
			rh = sqrt(P[j][2]*P[j][2] - rp*rp); 	
			x = sqrt(rd*rd-rh*rh);
			if( ( rh < h ) & (x < R_par[i]) ) { 
				//for(k=MAX(index_grid(x),0); k<MIN(R_par_ind[i], grid_n); k++) {
				//	par[i][k]++; 			
				//}
				if(rp<0) sign=0; else sign=1;
				par[sign][i][index_grid(x)]++; 
			}

			//rd = sqrt(dx*dx+dy*dy+dz*dz);
			rp = dx*vec[0] + dy*vec[1];// + dz*vec[2]; 
			rh = sqrt(rd*rd-rp*rp);
			if( ( rh < h ) & ( rp < R_ort[i]  ) ) { 
				//for(k=MAX(index_grid(rp),1); k<MIN(R_ort_ind[i], grid_n); k++) {
				//	ort[i][k]++; 			
				//}
				if(rp<0) sign=0; else sign=1;
				ort[sign][i][index_grid(fabs(rp))]++; 
			}
			//if(i%10) { printf("%d %lf %lf %d %d\n", i, x, rp, par[i][index_grid(x)], ort[i][index_grid(rp)] ); }
		}
	}
	printf("\n");
}

double 
delta(int k, int *r_i, int ***m) { 
	int i, count;
	double res;
	count = 0;
	res=0;
	for(i=0; i<N; i++) {
		if((r_i[i] > k)) {
			res += m[0][i][k] - m[1][i][k];
			count++;
		}
	}
	if(count == 0) return 0; 
	return res / count;
}


double 
av_n(int k, int *r_i, int ***m) { 
	int i, count;
	double res;
	count = 0;
	res=0;
	for(i=0; i<N; i++) {
		if((r_i[i] > k)) {
			res += m[0][i][k] + m[1][i][k];
			count ++; 
		}
	}
	if(count == 0) return 0; 
	return res / count;
}

double 
av_N(int K, int *r_i, int ***m) { 
	int i, k, count;
	double res;
	count = 0;
	res=0;
	for(i=0; i<N; i++) {
		if((r_i[i] > K)) {
			for(k=0; k<K; k++)
				res += m[0][i][k] + m[1][i][k];
			count ++; 
		}
	}
	if(count == 0) return 0; 
	return res / count;
}

void
write_result() {
	FILE *output; 
	char output_filename[80]; 
	int k, i; 
	
	sprintf(output_filename, "%s_cy.dat", name);
	output = fopen(output_filename, "w");
	for(k=1; k<grid_n; k++) {
		fprintf(output, " %lf", grid[k]);
		fprintf(output, " %lf", av_n(k, R_par_ind, par)/2*(grid[k]-grid[k-1])*M_PI*h*h);
		fprintf(output, " %lf", av_n(k, R_ort_ind, ort)/2*(grid[k]-grid[k-1])*M_PI*h*h);
		fprintf(output, " %lf", av_N(k, R_par_ind, par)/2*grid[k]*M_PI*h*h);
		fprintf(output, " %lf", av_N(k, R_ort_ind, ort)/2*grid[k]*M_PI*h*h);
		fprintf(output, " %lf", delta(k, R_par_ind, par));
		fprintf(output, " %lf", delta(k, R_ort_ind, ort));
		fprintf(output, "\n");
	}
	fclose(output);
	/*
	// printf one fulled sphere 
	double a,b,r;
	for(k=10; k<grid_n; k+=20) 
	{
		sprintf(output_filename, "%s_sp_%1.lf.dat", name, grid[k]);
		output = fopen(output_filename, "w");
			
		for(i=0; i<N; i++) {
			if(R_par_ind[i] > k) {
				fprintf(output, "%lf ", P[i][2] );
				fprintf(output, "%lf ", Angle(p[i][0], p[i][1], p[i][2], mass[i][k][0], mass[i][k][1], mass[i][k][2]));
				fprintf(output, "%lf ", sqrt(pow(mass[i][k][0],2)+pow(mass[i][k][1],2)+pow(mass[i][k][2],2))/(m[i][k]+1));  
				fprintf(output, "%d  ", m[i][k] );
				fprintf(output, "%lf ", (p[i][0]*mass[i][k][0] + p[i][1]*mass[i][k][1] + p[i][2]*mass[i][k][2]) / P[i][2] / (m[i][k]+1) );
				fprintf(output, "%lf %lf %lf\n", mass[i][k][0], mass[i][k][1], mass[i][k][2]);
			}
		}
		fclose(output);
	}
	*/
}

void
boundary() {
	int index; 
	double r1, r2, r3, r4; 
	double da = 2*M_PI, db = 2*M_PI; 
	int i; 
	for(i=0; i<N; i++) { 
		{
			r1 = sqrt(R[1]*R[1] - h*h); 
			da = 2*M_PI, db = 2*M_PI;
			if( A[0] != A[1] - 2*M_PI) {
				da = MIN(P[i][0] - A[0], A[1] - P[i][0]); 
			}
			db = MIN(P[i][1] - B[0], B[1] - P[i][1]); 
			//R_par[i] = MIN(r - R[0], R[1] - r); 
			r2 = R[0]; 
			if( da < M_PI/2 ) 
				r2 = MAX(r2, h/(sin(da)*cos(P[i][1]))) ; 
			if( db < M_PI/2 ) 
				r2 = MAX(r2, h/tan(db) ); 
		
			if( (r2>=P[i][2])||(P[i][2]>r1) ) {
				R_par[i] = 0; 
				R_par_ind[i] = 0; 
				//continue; 
			} else { 
				R_par[i]     = MIN(P[i][2] - r2, r1 - P[i][2]);	
				R_par_ind[i] = index_grid(R_par[i]);
			}
		}
		
		{
			double vec[3], point[3]; 
			double r1, r2;
			double da = 2*M_PI, db = 2*M_PI; 
			vec[0] = 1. / sqrt(1 + p[i][0]*p[i][0] / (p[i][1]*p[i][1])); 
			vec[1] = - p[i][0]/p[i][1] * vec[0]; 
			vec[2] = 0;		
			for(index=0; index<grid_n; index++) {
				point[0] = vec[0] * grid[index]; 
				point[1] = vec[1] * grid[index]; 
				point[2] = vec[2] * grid[index]; 
				
				r1 = rmax(p[i][0]+point[0], p[i][1]+point[1], p[i][2]+point[2]);			
				r2 = rmax(p[i][0]-point[0], p[i][1]-point[1], p[i][2]-point[2]);			
				if(  ( r2 < h ) || ( r1 < h ) ||
					 ( r1 < 0 ) || ( r2 < 0 ) ) break;
			}

			index = MAX(0, index-1);
			if(index == 0) { 
				R_ort[i] = 0; 
				R_ort_ind[i] = 0; 
			} else {
				R_ort[i] = grid[index];
				R_ort_ind[i] = index_grid(grid[index]);
			}
		}
		/*
		{
			//if(((double) rand() / RAND_MAX)>0.0005) continue;
			double circle[400]; 
			double ex[3], ey[3], center[3], point[3], vec[3]; 
			double t; 
			double x,y,z;
			x=p[i][0]; 
			y=p[i][1];
			z=p[i][2];
			//printf("x y z = %lf %lf %lf \n", x, y, z); 	
			vec[0] = 1. / sqrt(1 + p[i][0]*p[i][0] / (p[i][1]*p[i][1])); 
			vec[1] = - p[i][0]/p[i][1] * vec[0]; 
			vec[2] = 0;		
			//printf("vec = %lf %lf %lf \n", vec[0], vec[1], vec[2]);
			
			int l, m;
			FILE *file, *miss;
			double cy[2][1000][3];
			char name[80], bane[80]; 
			sprintf(name, "foo_%d.dat", i);
			sprintf(bane, "boo_%d.dat", i);
			file = fopen(name, "a");
			miss = fopen(bane, "a");
			fprintf(file, "%lf %lf %lf\n", x, y, z);
			for(l=0; l<2; l++) 
			{//i=0;
				//cartesian_coordinats(1,&center[0],&center[1],&center[2],&a,&b,&foo);
				//printf("r_ort : %lf\n", R_ort[i]); 
				center[0] = vec[0] * R_ort[i];
				center[1] = vec[1] * R_ort[i];
				center[2] = vec[2] * R_ort[i];
				if(l==1) {center[0] *= -1; center[1] *= -1;center[2] *= -1;}
				//printf("center : %lf %lf %lf\n", center[0]+x, center[1]+y, center[2]+z);
				// vectors for orthogonal plane
				if((center[1]==0) & (center[2] == 0)) 
					{point[0]=0;point[1]=1;point[2]=0;}
				else
					{point[0]=1;point[1]=0;point[2]=0;}
				cross_product(center, point, ex);
				cross_product(center, ex   , ey);
				// normalization for EX and EY
				t=sqrt(ex[0]*ex[0] + ex[1]*ex[1] + ex[2]*ex[2]); 
					   ex[0]=ex[0]/t;ex[1]=ex[1]/t;ex[2]=ex[2]/t;
				t=sqrt(ey[0]*ey[0] + ey[1]*ey[1] + ey[2]*ey[2]); 
					   ey[0]=ey[0]/t;ey[1]=ey[1]/t;ey[2]=ey[2]/t;
				for(m=0; m<100; m++) {
					//for(t=0; t<2*M_PI; t+=M_PI/100) {
					t = (double) m/100 * 2 * M_PI;
					//center + radius * cos (t) * ex + radius * sin (t) * ey 
					cy[l][m][0] = x+center[0] + h*cos(t)*ex[0] + h*sin(t)*ey[0];
					cy[l][m][1] = y+center[1] + h*cos(t)*ex[1] + h*sin(t)*ey[1];
					cy[l][m][2] = z+center[2] + h*cos(t)*ex[2] + h*sin(t)*ey[2];
				}
			}
			int sorvalos = 0; 
			
			for(m=0; m<100; m++) {
				for(l=0; l<2; l++) {
					fprintf(file, "%lf %lf %lf %lf %lf %lf\n", cy[l][m][0], cy[l][m][1], cy[l][m][2], x, y, z);
					if(inside_ca(cy[l][m][0], cy[l][m][1], cy[l][m][2]) != 0) {
						sorvalos = 10; 
						fprintf(miss, "%le %le %le %d\n", cy[l][m][0], cy[l][m][1], cy[l][m][2], inside_ca(cy[l][m][0], cy[l][m][1], cy[l][m][2]));
					}
				}
			}
			fclose(file);
			fclose(miss);
			if(sorvalos == 10) printf("pew :: %lf %lf %lf\n", x, y, z);
		}*/
		//printf("%lf %lf\n", R_par[i], R_ort[i]);
	}	
}

void
read_data() {
	FILE 	*file; 
	char	input_name[80];
	double 	x, y, z, a, b, r;
	int i, oor = 0, hack; // out of range;	
	int NMAX = N;
	
	sprintf(input_name, "%s.dat", name);
	file = fopen(input_name, "r");
	if( file == NULL ) {
		printf("Can't find file %s.\n", input_name);
		exit(10);
	}
	double foo;
	for(i=0; i<N; i++) {
		foo=((double) rand() / RAND_MAX);
		hack = fscanf(file, "%lf %lf %lf", &(p[i][0]), &(p[i][1]), &(p[i][2]));
		spherical_coordinats(1, &(p[i][0]), &(p[i][1]), &(p[i][2]), &(P[i][0]), &(P[i][1]), &(P[i][2]));
		if( (inside_sp(P[i][0], P[i][1], P[i][2]) != 0) || ( q < foo))  {
			oor++;
			i--; N--;
		} 	//else
			//printf("%lf\n", P[i][1]);
	}
	fclose(file);
	printf("There are %d out of range. %d left.\n", oor, N);
	//
	//file = fopen("/tmp/selected.dat", "w"); 
	//for(i=0; i<N; i++) {
	//	fprintf(file, "%lf %lf %lf\n", a[i], b[i], r[i]);
	//}
	//fclose(file);
}

void
init_begining() {
	A[0]=A[0]*M_PI/180;  
	A[1]=A[1]*M_PI/180;  
	B[0]=B[0]*M_PI/180;  
	B[1]=B[1]*M_PI/180;  
	int i, k;
	// points. N of them. 
	
	p = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		p[i] = (double *) calloc(3, sizeof(double));
	}
	P = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		P[i] = (double *) calloc(3, sizeof(double));
	}
}

void
init_second() {
	int i,k;
	R_ort     = (double *) calloc(N, sizeof(double));
	R_ort_ind = (int    *) calloc(N, sizeof(int   ));
	R_par     = (double *) calloc(N, sizeof(double));
	R_par_ind = (int    *) calloc(N, sizeof(int   ));
	grid = (double *) calloc(grid_n, sizeof(double));
	mass = (double ***) calloc(N, sizeof(double));
	for(i=0; i<N; i++) {
		mass[i] = (double **) calloc(grid_n, sizeof(double));
		for(k=0; k<grid_n; k++) {
			mass[i][k] = (double *) calloc(3, sizeof(double)); 
		}
	}
	par = (int ***) calloc(2, sizeof(int **));
	for(k=0;k<2;k++) {
		par[k] = (int **) calloc(N, sizeof(int *));
		for(i=0; i<N; i++) {
			par[k][i] = (int *) calloc(grid_n, sizeof(int));
		}
	}
	ort = (int ***) calloc(2, sizeof(int **));
	for(k=0;k<2;k++) {
		ort[k] = (int **) calloc(N, sizeof(int *));
		for(i=0; i<N; i++) {
			ort[k][i] = (int *) calloc(grid_n, sizeof(int));
		}
	}	
	/*	
	pos = (int ***) calloc(N, sizeof(int **));
	for(i=0; i<N; i++) {
		pos[i] = (int **) calloc(grid_n, sizeof(int *));
		if(pos[i] == NULL) { printf("I CANT SORRY\n"); exit(42);}
		//for(k=0; k<grid_n; k++) {
		//	pos[i][k] = (int *) calloc(2, sizeof(int));
		//}
	}
	*/
}

int
index_grid(double r) {
	int i;
	i = (int) (log(r / grid[0] ) / log(grid_h1)); 
	if(i<0) i = 0;
	return i;
}

make_grid() {
	grid_h1 = pow(grid_R[1]/grid_R[0], (double) 1./((double) grid_n - 1));
	grid[0] = grid_R[0]; 
	int k;
	for( k=1; k<grid_n; k++)	{
		grid[k] = grid[k-1] * grid_h1;
	}	
}

print_summ() {
	printf("%d points\n.", N);
	printf("A:%lf-%lf B:%lf-%lf R:%lf-%lf\n", A[0]*180/M_PI, A[1]*180/M_PI, B[0]*180/M_PI, B[1]*180/M_PI, R[0], R[1]);
}

void
gogogo() {
	init_begining();
	read_data();
	init_second();
	print_summ();
	make_grid();
	boundary();
	count(); 
	write_result(); 
}


int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:g:n:h:a:b:r:q:"};
	int opt; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'd':
				name = optarg; 
				break;
			case 'N':
				sscanf(optarg, "%d", &N); 
				break;
			case 'g':
				sscanf(optarg, "%lf", &grid_R[0]); 
				sscanf(argv[optind], "%lf", &grid_R[1]);
				optind++;
				break;
			case 'n':
				sscanf(optarg, "%d", &grid_n); 
				break;
			case 'a':
				sscanf(optarg, "%lf", &A[0]); 
				sscanf(argv[optind], "%lf", &A[1]);
				optind++;
				break ;
			case 'b':
				sscanf(optarg, "%lf", &B[0]); 
				sscanf(argv[optind], "%lf", &B[1]);
				optind++;
				break ;
			case 'r':
				sscanf(optarg, "%lf", &(R[0])); 
				sscanf(argv[optind], "%lf", &(R[1]));
				optind++;
				break ;
			case 'q': 
				sscanf(optarg, "%lf", &q); 
				//sscanf(optarg, "%lf", %q);
				break;
			case 'h': 
				sscanf(optarg, "%lf", &h); 
				break;
						
			default:
				printf("What? %d\n", opt);
				break;
        }
        opt = getopt( argc, argv, optString );
    }
	gogogo();
	return(0);	
}

/*
double
calc_nc(int k) {
	int i;
	int nc; 
	for(i=0; i<N; i++) {
		if(R_par_ind[i] > k)
			nc++; 
	}
	return nc;
}

double
calc_shift(int k) {
	// Distance: r(radius of the spere) -- shift of center of mass. 
	// there is no direction btw. Simple like that. 
	double x,y,z;
	double shift; 
	int i, count;
	
	shift = 0;
	count = 0;
	for(i=0; i<N; i++) {
		if((R_par_ind[i] > k)&(m[i][k]!=1)) {
			shift += sqrt(pow(mass[i][k][0],2)+pow(mass[i][k][1],2)+pow(mass[i][k][2],2))/(m[i][k]+1);
			count++; 
		}
	}
	return shift/count;
}
double
vec_shift(int k) {
	double x,y,z;
	int i, count;
	
	x=0;y=0;z=0;
	count = 0;
	for(i=0; i<N; i++) {
		if(R_par_ind[i] > k) {
			x += mass[i][k][0]/(m[i][k]+1);
			y += mass[i][k][1]/(m[i][k]+1);
			z += mass[i][k][2]/(m[i][k]+1);
			count++; 
		}
	}
	return sqrt(x*x+y*y+z*z)/count;
}

double
angle_to_the_Line_Of_Sight(int k) {
	double x,y,z;
	double angle; 
	int i, count;
	
	angle = 0;
	count = 0;
	x=0;y=0;z=0;
	for(i=0; i<N; i++) {
		if((R_par_ind[i] > k)&(m[i][k]!=0)) {
			//x += mass[i][k][0]/2;
			//y += mass[i][k][1]/2;
			//z += mass[i][k][2]/2;
			angle += Angle(p[i][0], p[i][1], p[i][2], mass[i][k][0], mass[i][k][1], mass[i][k][2]);
			count ++; 
		}
	}
	//shift = sqrt(x*x+y*y+z*z); 
	return angle / count;
}

double 
delta( int k) { 
	double d; 
	int i, count;
	double n1,n2,n; 
	double res;
	n1=0;n2=0;n=0;
	for(i=0; i<N; i++) {
		if(R_par_ind[i] > k) {
			n1 += pos[i][k][0];
			n2 += pos[i][k][1]; 
			n  += m[i][k];	
			//printf("%lf %d\n", grid[k], m[i][k]);
		}
	}
	res = (double) (n1-n2)/n;
	if(n==0) res = 0;
	//printf("********************%lf %lf %lf %lf\n", n1, n2, n, res);
	return res;
}
double 
projection( int k) { 
	double d; 
	int i, count;
	double res;
	count = 0;
	res=0;
	for(i=0; i<N; i++) {
		if((R_par_ind[i] > k)&(m[i][k]!=0)) {
			res += (p[i][0]*mass[i][k][0] + p[i][1]*mass[i][k][1] + p[i][2]*mass[i][k][2]) / P[i][2] / (m[i][k]+1);
			count ++; 
		}
	}
	if(count == 0) return 0; 
	return res / count;
}
*/
