/* Program vectrino */

#include <sys/dir.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

#define MAX(A,B) ((A)>(B) ? (A):(B))
#define MIN(A,B) ((A)<(B) ? (A):(B))

// variable definition

double tol = 1E-8; // Tolerance

double Karman = 0.41;   // Von Karman constant
double C = 5.1;         // Constant log-law
double mu = 1E-3;       // Dynamic velocity
double rho = 1E+3;      // Density
double nu = 1E-6;       // Kinematic velocity

float b;

int i, j, k, n;

FILE *fpin;
FILE *fpout;

char outfile[512];

int d, sd, f;

struct direct **dirs;
int ndirs;

struct direct **subdirs;
int nsubdirs;

struct direct **files;
int nfiles;

int dir_filter (const struct direct *entry);
char file[512];

int dnull;
float fnull;

int direction;

float Re; 				// Reynolds number
float W, H;  			// flume width; water flow depth
float h;  				// step height
float dh;  				// hydraulic diameter

// Vectrino

float *Vx_vec, *Vy_vec, *Vz_vec, *Vz1_vec, *Vz2_vec;  	
float *snrx_vec, *snry_vec, *snrz_vec, *snrz1_vec, *snrz2_vec;
float *corx_vec, *cory_vec, *corz_vec, *corz1_vec, *corz2_vec;  
float snrx_vec_av, snry_vec_av, snrz_vec_av;
float snr_vec_av;
float corx_vec_av, cory_vec_av, corz_vec_av;
float cor_vec_av;

// Filter data

int nav;
float Vx_lmin;
float Vx_lmax;
  
float Vx_vec_av, Vy_vec_av, Vz_vec_av;

// Backstep data  

float x, y, z;  			// x, y, z: backstep
float *xh, yh, zh;  			// positions divided by the step height
float *Vx, *Vy, *Vz;			// velocity  	
float *Vx_f, *Vy_f, *Vz_f;  		// fluctuations of the velocities	
float Vx_av;
float Vy_av;
float Vz_av;
float Vx_std;
float Vy_std;
float Vz_std;

float Vx_min_old;
float Vx_max_old;
float Vy_min_old;
float Vy_max_old;
float Vz_min_old;
float Vz_max_old;

float Vx_min;
float Vx_max;
float Vy_min;
float Vy_max;
float Vz_min;
float Vz_max;

float sumx;
float sumy;
float sumz;
float sum;
float Sx;  		// skewness
float Sy;
float Sz;
float Fx;  		// flatness
float Fy;
float Fz;
float ReS;  		// Reynolds stresses (Reynolds shear stresses)
float K;  		// turbulent kinetic energy
float *utau;  		// shear or friction velocity
float Ue;  		// maximum free-stream velocity
float *cf;  		// skin friction
float yplus, Uplus;	// dimensionless distance from the wall / velocity
float *delta_sl;	// sublayer thikness
float *tauw;    	// wall shear stress
float tau11;		// normal Reynolds stresses in x direction
float tau22;		// normal Reynolds stresses in y direction
float TIx;		// turbulent intensity in x direction
float TIy;		// turbulent intensity in y direction

float *xh_aux;

float offset;  	// offset y from wall
int type;		// backstep type
int phase_space_filter; // phase space filter
float yhmin;		// yhmin for utau

float ymin, ymax;
int first;
int last;

float m;
float rlength;

float db = 3.0;        // snr 
float cq = 0.85;       // cor

int nt;

void
sort(float *A, float *B, int n)
{

	int i, j, increment;
	float tmp1, tmp2;

	for (increment = n / 2; increment > 0; increment /= 2)  
		for (i = increment; i < n; i++)
		{
			tmp1 = A[i];
			tmp2 = B[i];
			for(j = i; j >= increment; j -= increment)
				if(tmp1 < A[j - increment])
				{
					A[j] = A[j - increment];
					B[j] = B[j - increment];
				}
				else
					break;
			A[j] = tmp1;
			B[j] = tmp2;
		}
}

int
compare (const void *a, const void *b)
{     
	const float *da = (const float *) a;
	const float *db = (const float *) b;

	return (*da > *db) - (*da < *db);
}

int
file_filter (const struct direct *entry)
{
	if (strcmp (entry->d_name + strlen (entry->d_name) - 3, "dat") != 0)
		return 0;

	return 1;
}

int
dir_filter (const struct direct *entry)
{
	int isdir;

	if (strcmp (entry->d_name, ".") == 0 || strcmp (entry->d_name, "..") == 0)
		return 0;

	isdir = chdir (entry->d_name);

	if (isdir == -1)
		return 0;

	chdir ("..");

	return 1;
}

float average(float *u, int n)
{
	int i;
	float sum = 0.0f;

	for( i = 0; i < n; i++)
		sum += u[i];

	if (n > 0)
		return (sum / n);
	else 
		return 0.0f;

}

float min(float *u, int n)
{
	
	float u_min = +1E+30;

	for (i = 0; i < n; i++)
	{
		u_min = MIN (u_min, u[i]);
	}
	
	return u_min;
	
}

float max(float *u, int n)
{
	
	float u_max = -1E+30;

	for (i = 0; i < n; i++)
	{
		u_max = MAX (u_max, u[i]);
	}
	
	return u_max;
	
}

float standard_dev(float *u, int n)
{

	float sum = 0.0f;
	float mean = average(u, n);
	int i;

	for(i = 0; i < n; i++)
		sum += pow(u[i] - mean, 2);

	if (n > 0)
		return sqrt(sum / n);
	else
		return 0.0f;

}

float calc_teta(float *u, float *d2u, int n)
{

	int i;
	float sum1 = 0.0f;
	float sum2 = 0.0f;

	for (i = 0; i < n; i++)
	{
		sum1 += u[i] * d2u[i];
		sum2 += u[i] * u[i];
	}

	return atan(sum1 / sum2);

}

int point_inside_ellipse(float x, float y, float a, float b)
{

	float f = (x * x) / (a * a) + (y * y) / (b * b); 

	if (f <= 1) 
		return 1;
	else 
		return 0;

}

void smooth_spike(float *u, float u_av, int i, int n)
{

	if (i > 1 && i < n - 1)
		u[i] = (u[i+1] + u[i-1]) * 0.5;
	else
		u[i] = (u[n - 1] + u[0]) * 0.5;
	
	//u[i] = u_av;

	//u[i] = u_max;
	
	/*
	if (i > 0)
		u[i] = u[i-1];
	else
		u[i] = u[n - 1];
	*/
	
}

int remove_spikes(float *u, float *v, float *w, int n, float u_min, float u_max)
{
	
	int m;
	float *u_aux, *v_aux, *w_aux;
	int *delete;

	u_aux = malloc (n * sizeof (float));
	v_aux = malloc (n * sizeof (float));
	w_aux = malloc (n * sizeof (float));

	delete = malloc (n * sizeof (int));

	for (i = 0; i < n; i++)
		delete[i] = 0;
	
	printf("u_min: %f\n", u_min);
	printf("u_max: %f\n", u_max);
	
	for (i = 0; i < n; i++)
	{
		if(u[i] < u_min || u[i] > u_max)
			delete[i] = 1;
	}
	
	m = 0;
	
	for (i = 0; i < n; i++)
	{
		if(delete[i] == 0)
		{
			u_aux[m] = u[i];
			v_aux[m] = v[i];
			w_aux[m] = w[i];
			m++;
		}
	}

	printf ("\nNumber of non-physical spikes removed: %d\n", n - m);
	
	n = m;
	
	for (i = 0; i < n; i++)
	{
		u[i] = u_aux[i];
		v[i] = v_aux[i];
		w[i] = w_aux[i];
	}
	
	free(u_aux);
	free(v_aux);
	free(w_aux);
	
	free(delete);
	
	return n;
	
}

int read_raw_data(char *dir_name)
{

	chdir(dir_name);

	nfiles = scandir (".", &files, file_filter, alphasort);

	i = 0;

	for (f = 0; f < nfiles; f++)
	{
		
		strcpy (file, files[f]->d_name);

		if (strncmp (file, "exp", 3) != 0)
			continue;
	
		fpin = fopen (file, "r");

		//printf ("\n************ Reading file: %s\n", file);

		do
		{

			fnull = -1;

			fscanf (fpin, "%f", &fnull);

			if (fnull == -1)
				break;

			//printf("%f\n", fnull); 

			fscanf (fpin,
				"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
				&fnull, 
				&Vx_vec[i], &Vy_vec[i], &Vz1_vec[i], &Vz2_vec[i], 
				&fnull, &fnull, &fnull, &fnull, 
				&snrx_vec[i], &snry_vec[i], &snrz1_vec[i], &snrz2_vec[i], 
				&corx_vec[i], &cory_vec[i], &corz1_vec[i], &corz2_vec[i]);  

			Vz_vec[i] = (Vz1_vec[i] + Vz2_vec[i]) * 0.5;

			//if (Vx_vec[i] >= Vx_lmin && Vx_vec[i] <= Vx_lmax)	
			
			i++;

		}
		while (!feof (fpin));

		printf ("\n************ Done reading file: %s\n", file);

	}

	chdir("..");
	
	printf ("\nTotal number of samples in raw data: %d\n", i);
	
	return i;

}

int save_vector(char *file_name, float *u, int n)
{
	
	FILE *fp = fopen (file_name, "w");
	
	for (i = 0; i < n; i++)
	{
		fprintf(fp, "%f\n", u[i]);
	}

	fclose(fp);
	
}

/*
int calc_power_spectrum(float *u, int n, float u_av, float sample_rate)
{

	float *ul;	
	float *frequencies;
	
	int halfpoint;

	frequencies = malloc (n * sizeof (float));
	ul = malloc (n * sizeof (float));

	float nyquistfreq = sample_rate * 0.5;

	// Fluctuating vel u'
	for (i = 0; i < n; i++)
		ul[i] = u[i] - u_av;
	
	halfpoint = n / 2;

	for (i = 0; i < halfpoint; i++)
		frequencies[i] = i / halfpoint * nyquistfreq;
			
	thefft=fft(wvel);
	Power=thefft.*conj(thefft);
	
	free(frequencies);
	free(ul);
	
}*/

int apply_phase_space_filter(float *u, int n)
{

	int ip, im;

	int n_spikes;
	float perc_spikes;

	int iter = 0;

	int *spike;
	
	float u_min_old = +1E+30;
	float u_max_old = -1E+30;
		
	float u_min, u_max;

	float *du, *d2u;

	float u_av;
	
	float sigma_u;
	float sigma_du;
	float sigma_d2u;

	float Emax_u;	
	float Emax_du;
	float Emax_d2u;
	
	float f1, f2;
	float fa, fb;
	float a, b;

	float teta;

	int continue_loop;
	
	du = malloc (n * sizeof (float));
	d2u = malloc (n * sizeof (float));

	spike = malloc (n * sizeof (int));
	
	// Apply Phase-space thresholding method filter	(Goring and Nikora, 2002)

	for (i = 0; i < n; i++)
	{
		spike[i] = 0;
	}

	printf ("\nBegin phase-space thresholding method\n");

	u_av = average(u, n);

	// Fluctuating vel u'
	for (i = 0; i < n; i++)
		u[i] -= u_av;
	
	/*
	for (i = 0; i < 30; i++)
		printf("u = %f\n", u[i]);
	*/
	
	do
	{

		continue_loop = 0;

		u_max = max(u, n);
		u_min = min(u, n);
		
		if (fabs(u_min - u_min_old) < 1E-8 &&
			fabs(u_max - u_max_old) < 1E-8)
			break;

		u_min_old = u_min;
		u_max_old = u_max;

		printf ("iter=%d\n", iter);
		
		printf ("u_min=%f\n", u_min);
		printf ("u_max=%f\n", u_max);

		// Step 1
		for (i = 0; i < n; i++)
		{

			ip = i + 1;
			im = i - 1;

			if (im < 0) im = n;
			if (ip > n - 1) ip = 0;

			du[i] = (u[ip] - u[im]) * 0.5;
		}

		for (i = 0; i < n; i++)
		{

			ip = i + 1;
			im = i - 1;

			if (im < 0) im = n;
			if (ip > n - 1) ip = 0;

			d2u[i] = (du[ip] - du[im]) * 0.5;
		}

		// Step 2
		printf("u_av = %f\n", u_av);
		
		sigma_u = standard_dev(u, n);
		sigma_du = standard_dev(du, n);
		sigma_d2u = standard_dev(d2u, n);

		printf("sigma_u = %f\n", sigma_u);
		
		Emax_u = sqrt(2 * log(n)) * sigma_u;
		Emax_du = sqrt(2 * log(n)) * sigma_du;
		Emax_d2u = sqrt(2 * log(n)) * sigma_d2u;

		printf("Emax_u = %f\n", Emax_u);
		
		// Step 3
		teta = calc_teta(u, d2u, n);

		// Step 4
		//Ellipse x1
		f1 = Emax_u;
		f2 = Emax_d2u;

		fa = ((f1*cos(teta))*(f1*cos(teta)) - (f2*sin(teta))*(f2*sin(teta))) / (cos(teta)*cos(teta) - sin(teta)*sin(teta));
		fb = ((f2*cos(teta))*(f2*cos(teta)) - (f1*sin(teta))*(f1*sin(teta))) / (cos(teta)*cos(teta) - sin(teta)*sin(teta));

		a = sqrt(fa);
		b = sqrt(fb);
				
		/*
		printf ("iter=%d\n", iter);
		printf ("fa=%f\n", fa);
		printf ("fb=%f\n", fb);
		printf ("teta=%f\n", teta);
		printf ("x -> a=%f, b=%f\n", a, b);
		*/

		for (i = 0; i < n; i++)
		{

			if (spike[i] == 1) continue;
			
			if (!point_inside_ellipse(u[i], d2u[i], a, b))
			{
				//printf("e1 -> u[i]=%f, d2u[i]=%f, a=%f, b=%f\n", u[i], d2u[i], a, b);
				
				// Smooth spike
				spike[i] = 1;
				smooth_spike(u, u_av, i, n);
				continue_loop = 1;
			}
		}
	
		//Ellipse x2
		a = Emax_u;
		b = Emax_du;

		for (i = 0; i < n; i++)
		{

			if (spike[i] == 1) continue;

			if (!point_inside_ellipse(u[i], du[i], a, b))
			{
				//printf("e2 -> u[i]=%f, du[i]=%f, a=%f, b=%f\n", u[i], du[i], a, b);

				// Smooth spike
				spike[i] = 1;
				smooth_spike(u, u_av, i, n);
				continue_loop = 1;					
			}
		}

		//Ellipse x3
		a = Emax_du;
		b = Emax_d2u;

		for (i = 0; i < n; i++)
		{

			if (spike[i] == 1) continue;

			if (!point_inside_ellipse(du[i], d2u[i], a, b))
			{
				//printf("e3 -> du[i]=%f, d2u[i]=%f, a=%f, b=%f\n", du[i], d2u[i], a, b);
			
				// Smooth spike
				spike[i] = 1;
				smooth_spike(u, u_av, i, n);
				continue_loop = 1;					
			}
		}

		iter++;

	}while (continue_loop);

	// Statistics 
	// Number of spikes smoothed
	n_spikes = 0;

	for (i = 0; i < n; i++)
	{
		if (spike[i] == 1) 
			n_spikes++;
	}

	free(du);
	free(d2u);

	free(spike);
	
	// Percentage filtered
	perc_spikes = ((float) n_spikes) / ((float) n) * 100;

	printf ("\nNumber of spikes smoothed: %d\n", n_spikes);
	printf ("\nPercentage filtered: %.2f\n", perc_spikes);

	// Fluctuating vel u'
	for (i = 0; i < n; i++)
		u[i] += u_av;
	
	/*
	if (n_spikes == n)
	{

		sigma_u = standard_dev(u, n);
		printf("sigma_u = %f\n", sigma_u);
	
		for (i = 0; i < 30; i++)
			printf("u = %f\n", u[i]);
			
		exit(0);
	}
	*/
	
}

int transform_coord()
{

	// Transform vectrino coordinate system to backstep coordinate system
	for (i = 0; i < n; i++)
	{
		Vx[i] = Vx_vec[i];
		Vy[i] = Vz_vec[i];
		Vz[i] = Vy_vec[i];
	}

}

int
main (int argc, char *argv[])
{

	// Commands
	printf ("\n");
	printf ("******************* Program: vectrino *******************\n");
	printf ("\n");

	if (argc == 8)
	{
		Re = atof (argv[1]);
	    Vx_lmin = atof (argv[2]);
		Vx_lmax = atof (argv[3]);
		type = atoi (argv[4]);
		phase_space_filter = atoi (argv[5]);
		offset = atof (argv[6]);
		yhmin = atof (argv[7]);
	}
	else
	{
		// Re: estimated Reynolds
		// Vx_lmin: minimum x velocity threshold
		// Vx_lmax: maximum x velocity threshold
		// type: 0- without step, 1- with step, 2- with sand, 3- with step and tap
		// phase_space_filter: 1 - enable, 0 - disable
		// offset: offset in y from wall
		// yhmin: ymin for utau
		printf("./vectrino Re Vx_lmin Vx_lmax type phase_space_filter offset yhmin\n\n");
		return 0;
	} 

	printf ("Re: %g\n", Re);

	/***************************** START DEFINE GEOMETRY *****************************/

	// Flume geometry
	W = 0.304;

	// Calculate free stream U   

	// Type: 0- without step; 1- with step; 2- with sand; 3- with step and tap

	if (type == 0)     
	{
		H = 0.24;
		dh = 2 * W * H / (2*W + 2*H);    
		Ue = Re * nu / dh;
		h = 0.0;
	}

	if (type == 1)     
	{
		H = 0.24;
		h = 0.04;
		Ue = Re * nu / h;
	}

	if (type == 2)     
	{
		H = 0.186;
		dh = 2 * W * H / (2*W + 2*H);    
		Ue = Re * nu / dh;
		h = 0.0;
	}

	if (type == 3)     
	{
		H = 0.289;
		h = 0.04;
		Ue = Re * nu / h;
	}

	/***************************** END DEFINE GEOMETRY *****************************/

	printf ("Minimum limit Vx: %g\n", Vx_lmin);
	printf ("Maximum limit Vx: %g\n", Vx_lmax);
	printf ("Type: %d\n", type);
	printf ("Phase Space Filter: %d\n", phase_space_filter);	
	printf ("Offset: %f\n", offset);
	printf ("Value yhmin: %f\n", yhmin);

	/***************************** START MEMORY ALLOCATION *****************************/

	n = 100000;

	Vx_vec    = malloc (n * sizeof (float));
	Vy_vec    = malloc (n * sizeof (float));
	Vz_vec    = malloc (n * sizeof (float));
	Vz1_vec   = malloc (n * sizeof (float));
	Vz2_vec   = malloc (n * sizeof (float));
	snrx_vec  = malloc (n * sizeof (float));
	snry_vec  = malloc (n * sizeof (float));
	snrz_vec  = malloc (n * sizeof (float));
	snrz1_vec = malloc (n * sizeof (float));
	snrz2_vec = malloc (n * sizeof (float));
	corx_vec  = malloc (n * sizeof (float));
	cory_vec  = malloc (n * sizeof (float));
	corz_vec  = malloc (n * sizeof (float));
	corz1_vec = malloc (n * sizeof (float));
	corz2_vec = malloc (n * sizeof (float));

	Vx   = malloc (n * sizeof (float));
	Vy   = malloc (n * sizeof (float));
	Vz   = malloc (n * sizeof (float));
	Vx_f = malloc (n * sizeof (float));
	Vy_f = malloc (n * sizeof (float));
	Vz_f = malloc (n * sizeof (float));

	xh = malloc (30 * sizeof (float));
	xh_aux = malloc (30 * sizeof (float));
	cf = malloc (30 * sizeof (float));
	utau = malloc (30 * sizeof (float));
	delta_sl = malloc (30 * sizeof (float));
	tauw = malloc (30 * sizeof (float));

	/***************************** END MEMORY ALLOCATION *****************************/

	ndirs = scandir (".", &dirs, dir_filter, alphasort);

	for (d = 0; d < ndirs; d++)
	{

		if (strcmp(dirs[d]->d_name, "end") == 0) continue;
		if (strcmp(dirs[d]->d_name, "start") == 0) continue;

		printf ("\n* Directory: %s\n", dirs[d]->d_name);
		xh[d] = atof(dirs[d]->d_name + 3);

		chdir (dirs[d]->d_name);

		nsubdirs = scandir (".", &subdirs, dir_filter, alphasort);

		sprintf (outfile, "results_%.3f.dat", xh[d]);

		fpout = fopen (outfile, "w");

		fprintf (fpout, "#         1              2");
		fprintf (fpout, "               3              4");
		fprintf (fpout, "               5               6                7               8             9               10");
		fprintf (fpout, "               11              12               13              14            15              16");
		fprintf (fpout, "               17              18               19              20            21              22");
		fprintf (fpout, "               23           24              25  ");
		fprintf (fpout, "            26             27               28");
		fprintf (fpout, "            29             30          ");
		fprintf (fpout, "       31");
		fprintf (fpout, "              32              33 ");
		fprintf (fpout, "              34   ");
		fprintf (fpout, "           35 ");
		fprintf (fpout, "           36  ");
		fprintf (fpout, "           37   ");
		fprintf (fpout, "            38                    39    ");
		fprintf (fpout, "         40"); 
		fprintf (fpout, "         41 ");  
		fprintf (fpout, "         42             ");
		fprintf (fpout, "         43             ");
		fprintf (fpout, "         44             ");
		fprintf (fpout, "\n");  

		fprintf (fpout, "#       x/h              y/h");
		fprintf (fpout, "               x              y");
		fprintf (fpout, "              Vx_min         Vx_max           Vx_av           Vx_std          Sx              Fx");
		fprintf (fpout, "              Vy_min         Vy_max           Vy_av           Vy_std          Sy              Fy");
		fprintf (fpout, "              Vz_min         Vz_max           Vz_av           Vz_std          Sz              Fz");
		fprintf (fpout, "           Vx_av/Ue      Vy_av/Ue       Vz_av/Ue");
		fprintf (fpout, "           ReS              K              utau");
		fprintf (fpout, "        ReS/(Ue*Ue)   K/(Ue*Ue)");
		fprintf (fpout, "       cf");
		fprintf (fpout, "             yplus          Uplus");
		fprintf (fpout, "           delta_sl");
		fprintf (fpout, "          tauw");
		fprintf (fpout, "          tau11");
		fprintf (fpout, "           tau22");
		fprintf (fpout, "     tau11/(Ue*Ue) tau22/(Ue*Ue)");
		fprintf (fpout, "      TIx"); 
		fprintf (fpout, "         TIy");  
		fprintf (fpout, "         ReS/(utau*utau)");
		fprintf (fpout, "   tau11/(utau*utau)");
		fprintf (fpout, "    tau22/(utau*utau)");
		fprintf (fpout, "\n");  

		/***************************** START READ RESULTS *****************************/

		printf ("\nSTART READ RESULTS\n");
	
		ymin = +1E+30;      
		ymax = -1E+30;      
		first = 0;
		last = 0;

		for (sd = 0; sd < nsubdirs; sd++)
		{

			printf ("\n*** Sub-Directory: %s\n", subdirs[sd]->d_name);
			yh = atof(subdirs[sd]->d_name + 3);

			if (h != 0.0)	  
			{
				x = xh[d] * h;	  
				y = yh * h + offset;
			}
			else
			{
				x = xh[d] * 0.04;	  
				y = yh * 0.04 + offset;	  
			}

			if (x > 0 && y < 0.004) continue;
			if (x < 0 && y < h + 0.004) continue;

			if (h != 0.0)	  
				yh = y / h;
			else
				yh = y;	  

			// calculate closest y to the bottom 
			//(first= folder having closest y to the bottom; last = folder having further away y from the bottom)

			if (y < ymin)
			{
				ymin = y;
				first = sd;
			}

			if (y > ymax)
			{
				ymax = y;
				last = sd;
			}

		}

		/***************************** START CALCULATE DIRECTION *****************************/

		// First calculate direction
		printf ("\nCALCULATE DIRECTION\n");
		printf ("\n*** Sub-Directory: %s\n", subdirs[first]->d_name);
		yh = atof(subdirs[first]->d_name + 3);

		x = xh[first] * h;
		y = yh * h + offset;

		n = read_raw_data(subdirs[first]->d_name);
		
		n = remove_spikes(Vx_vec, Vy_vec, Vz_vec, n, Vx_lmin, Vx_lmax);
		
		if (n == 0)
		{		
			printf ("Error, n = 0\n");
			exit (0);
		}
		
		if (phase_space_filter == 1)
		{
			apply_phase_space_filter(Vx_vec, n);
			apply_phase_space_filter(Vy_vec, n);
			apply_phase_space_filter(Vz_vec, n);
		}

		transform_coord();

		for (i = 0; i < n; i++)
		{	      
			snrz_vec[i] = (snrz1_vec[i] + snrz2_vec[i]) * 0.5;
			corz_vec[i] = (corz1_vec[i] + corz2_vec[i]) * 0.5;
		}

		// Determination of average - snr
		snrx_vec_av = average(snrx_vec, n);
		snry_vec_av = average(snry_vec, n);
		snrz_vec_av = average(snrz_vec, n);

		snr_vec_av = (snrx_vec_av + snry_vec_av + snrz_vec_av) / 3.0;

		// Determination of average - cor
		corx_vec_av = average(corx_vec, n);
		cory_vec_av = average(cory_vec, n);
		corz_vec_av = average(corz_vec, n);

		cor_vec_av = (corx_vec_av + cory_vec_av + corz_vec_av) / 3.0;

		Vx_av = average(Vx, n);

		if (Vx_av > 0)
			direction = +1;
		else
			direction = -1;      

		/***************************** END CALCULATE DIRECTION *****************************/

		/***************************** START CALCULATE FRICTION VELOCITY USING REYNOLDS SHEAR STRESS *****************************/

		printf ("\nSTART CALCULATE FRICTION VELOCITY USING REYNOLDS SHEAR STRESS\n");
	
		utau[d] = 0.0;
		nt = 0;

		for (sd = 0; sd < nsubdirs; sd++)
		{

			printf ("\n*** Sub-Directory: %s\n", subdirs[sd]->d_name);
			yh = atof(subdirs[sd]->d_name + 3);

			if (h != 0.0)	  
			{
				x = xh[d] * h;	  
				y = yh * h + offset;
			}
			else
			{
				x = xh[d] * 0.04;	  
				y = yh * 0.04 + offset;	  
			}

			//if (x > 0 && y < 0.004) continue;
			//if (x < 0 && y < h + 0.004) continue;

			if (h != 0.0)	  
				yh = y / h;
			else
				yh = y;	  

			//printf ("\n x = %f\n", x);
			//printf ("\n y = %f\n", y);

			n = read_raw_data(subdirs[sd]->d_name);

			n = remove_spikes(Vx_vec, Vy_vec, Vz_vec, n, Vx_lmin, Vx_lmax);

			if (phase_space_filter == 1)
			{
				apply_phase_space_filter(Vx_vec, n);
				apply_phase_space_filter(Vy_vec, n);
				apply_phase_space_filter(Vz_vec, n);
			}

			transform_coord();

			for (i = 0; i < n; i++)
			{	      
				snrz_vec[i] = (snrz1_vec[i] + snrz2_vec[i]) * 0.5;
				corz_vec[i] = (corz1_vec[i] + corz2_vec[i]) * 0.5;
			}

			// Determination of average - snr
			snrx_vec_av = average(snrx_vec, n);
			snry_vec_av = average(snry_vec, n);
			snrz_vec_av = average(snrz_vec, n);

			snr_vec_av = (snrx_vec_av + snry_vec_av + snrz_vec_av) / 3.0;

			// Determination of average - cor
			corx_vec_av = average(corx_vec, n);
			cory_vec_av = average(cory_vec, n);
			corz_vec_av = average(corz_vec, n);

			cor_vec_av = (corx_vec_av + cory_vec_av + corz_vec_av) / 3.0;

			nav = 0;	  
			sumx = 0.0;
			sumy = 0.0;
			sumz = 0.0;

			// Determination of average - velocity
			for (i = 0; i < n; i++)
			{
			  if (snrx_vec[i] > snr_vec_av - db && snry_vec[i] > snr_vec_av - db && snrz_vec[i] > snr_vec_av - db &&
				  corx_vec[i] > cor_vec_av * cq && cory_vec[i] > cor_vec_av * cq && corz_vec[i] > cor_vec_av * cq)
				{
					sumx += Vx[i];
					sumy += Vy[i];
					sumz += Vz[i];
					nav++;
				}
			}

			if (nav != 0)
			{
				Vx_av = sumx / nav;
				Vy_av = sumy / nav;
				Vz_av = sumz / nav;
			  
				//printf("n - nav: %d\n", n - nav);
			}
			else
			{
				printf ("Error, nav = 0\n");
				exit (0);
			}

			// Determination of Reynolds shear stresses
			sum = 0.0;

			for (i = 0; i < n; i++)
			{
				Vx_f[i] = (Vx[i] - Vx_av);
				Vy_f[i] = (Vy[i] - Vy_av);
				Vz_f[i] = (Vz[i] - Vz_av);
				sum += Vx_f[i] * Vy_f[i];
			}

			if (n > 0)
				ReS = -sum / n;
			else
			{
				printf ("Error, n = 0\n");
				exit (0);
			}	    

			if (yh > yhmin && yh < 2.0)   
			{
				
				utau[d] += sqrt(fabs(1.0/(1.0 - y / (2.0 * H)) * (nu * Vx_av / y - ReS)));				
				nt++;
			}

		}

		// Average friction velocity

		if (nt != 0)
			utau[d] /= nt;
		else
		{
			printf ("Error, nt = 0\n");
			exit(0);
		}

		printf("\nFriction velocity using expression for open-channel flows\n");  
		printf("--> utau = %f\n", utau[d]);

		/***************************** END CALCULATE FRICTION VELOCITY USING REYNOLDS SHEAR STRESS *****************************/

		/***************************** START CALCULATE ALL OTHER VALUES *****************************/
		
		printf ("\nSTART CALCULATE ALL OTHER VALUES\n");
	
		for (sd = 0; sd < nsubdirs; sd++)
		{

			printf ("\n*** Sub-Directory: %s\n", subdirs[sd]->d_name);
			yh = atof(subdirs[sd]->d_name + 3);

			if (h != 0.0)	  
			{
				x = xh[d] * h;	  
				y = yh * h + offset;
			}
			else
			{
				x = xh[d] * 0.04;	  
				y = yh * 0.04 + offset;	  
			}

			//if (x > 0 && y < 0.004) continue;
			//if (x < 0 && y < h + 0.004) continue;

			if (h != 0.0)	  
				yh = y / h;
			else
				yh = y;	  

			//printf ("\n x = %f\n", x);
			//printf ("\n y = %f\n", y);

			n = read_raw_data(subdirs[sd]->d_name);

			chdir(subdirs[sd]->d_name);
	
			sprintf (outfile, "results_Vx_raw.dat");
			save_vector(outfile, Vx_vec, n);
			
			sprintf (outfile, "results_Vy_raw.dat");
			save_vector(outfile, Vy_vec, n);

			sprintf (outfile, "results_Vz_raw.dat");
			save_vector(outfile, Vz_vec, n);
			
			n = remove_spikes(Vx_vec, Vy_vec, Vz_vec, n, Vx_lmin, Vx_lmax);
						
			if (phase_space_filter == 1)
			{
				apply_phase_space_filter(Vx_vec, n);
				sprintf (outfile, "results_Vx_filtered.dat");
				save_vector(outfile, Vx_vec, n);
			
				apply_phase_space_filter(Vy_vec, n);
				sprintf (outfile, "results_Vy_filtered.dat");
				save_vector(outfile, Vy_vec, n);
				
				apply_phase_space_filter(Vz_vec, n);
				sprintf (outfile, "results_Vz_filtered.dat");
				save_vector(outfile, Vz_vec, n);
			}

			chdir("..");
			
			transform_coord();

			for (i = 0; i < n; i++)
			{	      
				snrz_vec[i] = (snrz1_vec[i] + snrz2_vec[i]) * 0.5;
				corz_vec[i] = (corz1_vec[i] + corz2_vec[i]) * 0.5;
			}

			// Determination of average - snr
			snrx_vec_av = average(snrx_vec, n);
			snry_vec_av = average(snry_vec, n);
			snrz_vec_av = average(snrz_vec, n);

			snr_vec_av = (snrx_vec_av + snry_vec_av + snrz_vec_av) / 3.0;

			// Determination of average - cor
			corx_vec_av = average(corx_vec, n);
			cory_vec_av = average(cory_vec, n);
			corz_vec_av = average(corz_vec, n);

			cor_vec_av = (corx_vec_av + cory_vec_av + corz_vec_av) / 3.0;

			// Determination of average - velocity
			Vx_av = average(Vx, n);
			Vy_av = average(Vy, n);
			Vz_av = average(Vz, n);

			// Calculate power spectrum
			// Kolmogorov 5/3 law
			//calc_power_spectrum(Vx, n, Vx_av, 200);
			//calc_power_spectrum(Vy, n, Vy_av, 200);
			//calc_power_spectrum(Vz, n, Vz_av, 200);
			
			// Determination of limits - velocity
			Vx_min = min(Vx, n);
			Vx_max = max(Vx, n);

			Vy_min = min(Vy, n);
			Vy_max = max(Vy, n);

			Vz_min = min(Vz, n);
			Vz_max = max(Vz, n);
						
			// Determination of standard deviation
			Vx_std = standard_dev(Vx, n);
			Vy_std = standard_dev(Vy, n);
			Vz_std = standard_dev(Vz, n);

			// Determination of turbulence intensity

			TIx = Vx_std / Ue * 100;
			TIy = Vy_std / Ue * 100; 

			// Determination of Reynolds shear stresses
			sum = 0.0;

			for (i = 0; i < n; i++)
			{
				Vx_f[i] = (Vx[i] - Vx_av);
				Vy_f[i] = (Vy[i] - Vy_av);
				Vz_f[i] = (Vz[i] - Vz_av);
				sum += Vx_f[i] * Vy_f[i];
			}

			if (n > 0)
				ReS = -sum / n;
			else
			{
				printf ("Error, n = 0\n");
				exit (0);
			}	    

			// Determination of Skewness
			Sx = Vx_av / (Vx_std * Vx_std * Vx_std);
			Sy = Vy_av / (Vy_std * Vy_std * Vy_std);
			Sz = Vz_av / (Vz_std * Vz_std * Vz_std);

			// Determination of Flatness
			Fx = Vx_av / (Vx_std * Vx_std * Vx_std * Vx_std);
			Fy = Vy_av / (Vy_std * Vy_std * Vy_std * Vy_std);
			Fz = Vz_av / (Vz_std * Vz_std * Vz_std * Vz_std);

			// Determination of Reynolds normal stresses and turbulent kinetic energy
			sumx = 0.0;
			sumy = 0.0;

			for (i = 0; i < n; i++)
			{
				sumx += Vx_f[i] * Vx_f[i];
				sumy += Vy_f[i] * Vy_f[i];
			}

			if (n > 0)
			{

				tau11 = sumx / n;
				tau22 = sumy / n;

				K = 0.5 * (tau11 + tau22);
			}
			else
			{
				printf ("Error, n = 0\n");
				exit (0);
			}

			//printf ("Ue: %f\n", Ue);

			// Determination of skin friction
			cf[d] = +2 * direction * utau[d] * utau[d] / (Ue * Ue);

			// Determination of yplus and Uplus
			yplus = y * utau[d] / nu;
			Uplus = Vx_av / utau[d];

			// Sub-layer thickness
			delta_sl[d] = 30 * nu / utau[d];

			// Wall shear stress
			tauw[d] = rho * utau[d] * utau[d];

			/***************************** END CALCULATE ALL OTHER VALUES *****************************/

			/***************************** START PRINT VALUES *****************************/

			printf ("\n");
			printf ("Number of entries:   %d\n", n);
			printf ("\n");
			printf ("Minimum Vx [m/s]: %f\n", Vx_min);
			printf ("Maximum Vx [m/s]: %f\n", Vx_max);
			printf ("Average Vx [m/s]: %f\n", Vx_av);
			printf ("St Dev  Vx [m/s]: %f\n", Vx_std);
			printf ("Skewness x   [-]: %f\n", Sx);
			printf ("Flatness x   [-]: %f\n", Fx);
			printf ("\n");
			printf ("Minimum Vy [m/s]: %f\n", Vy_min);
			printf ("Maximum Vy [m/s]: %f\n", Vy_max);
			printf ("Average Vy [m/s]: %f\n", Vy_av);
			printf ("St Dev  Vy [m/s]: %f\n", Vy_std);
			printf ("Skewness y   [-]: %f\n", Sy);
			printf ("Flatness y   [-]: %f\n", Fy);
			printf ("\n");
			printf ("Minimum Vz [m/s]: %f\n", Vz_min);
			printf ("Maximum Vz [m/s]: %f\n", Vz_max);
			printf ("Average Vz [m/s]: %f\n", Vz_av);
			printf ("St Dev  Vz [m/s]: %f\n", Vz_std);
			printf ("Skewness z   [-]: %f\n", Sz);
			printf ("Flatness z   [-]: %f\n", Fz);
			printf ("\n");
			printf ("Friction velocity [m/s]: %f\n", utau[d]);
			printf ("Skin friction       [-]: %f\n", cf[d]);
			printf ("Sub-layer thickness [m]: %f\n", delta_sl[d]);
			printf ("Wall shear stress  [Pa]:  %f\n", tauw[d]);
			printf ("\n");	  
			printf ("Reynolds shear stress x,y [m^2/s^2]: %f\n", ReS);
			printf ("Turb. kin. ener. x,y      [m^2/s^2]: %f\n", K);	  
			printf ("Normal Reynolds stress x,x [m^2/s^2]: %f\n", tau11);
			printf ("Normal Reynolds stress y,y [m^2/s^2]: %f\n", tau22);
			printf ("Turbulence intensity x          [%%]: %f\n", TIx);
			printf ("Turbulence intensity y          [%%]: %f\n", TIy);

			/***************************** END PRINT VALUES *****************************/

			/***************************** START SAVE VALUES TO FILE *****************************/

			fprintf (fpout, " %+.8E %+.8E", xh[d], yh);
			fprintf (fpout, " %+.8E %+.8E", x, y);
			fprintf (fpout, " %+.8E %+.8E %+.8E %+.8E %+.8E %+.8E", Vx_min, Vx_max, Vx_av, Vx_std, Sx, Fx);
			fprintf (fpout, " %+.8E %+.8E %+.8E %+.8E %+.8E %+.8E", Vy_min, Vy_max, Vy_av, Vy_std, Sy, Fy);
			fprintf (fpout, " %+.8E %+.8E %+.8E %+.8E %+.8E %+.8E", Vz_min, Vz_max, Vz_av, Vz_std, Sz, Fz);
			fprintf (fpout, " %+.8E %+.8E %+.8E", Vx_av / Ue, Vy_av / Ue, Vz_av / Ue);
			fprintf (fpout, " %+.8E %+.8E %+.8E", ReS, K, utau[d]);
			fprintf (fpout, " %+.8E %+.8E", ReS / (Ue * Ue), K / (Ue * Ue));
			fprintf (fpout, " %+.8E", cf[d]);
			fprintf (fpout, " %+.8E %+.8E", yplus, Uplus);
			fprintf (fpout, " %+.8E", delta_sl[d]);
			fprintf (fpout, " %+.8E", tauw[d]);
			fprintf (fpout, " %+.8E", tau11);
			fprintf (fpout, " %+.8E", tau22);
			fprintf (fpout, " %+.8E %+.8E", tau11 / (Ue * Ue), tau22 / (Ue * Ue));  
			fprintf (fpout, " %+.8E", TIx);
			fprintf (fpout, " %+.8E", TIy);
			fprintf (fpout, " %+.8E", ReS / (utau[d] * utau[d]));
			fprintf (fpout, " %+.8E", tau11 / (utau[d] * utau[d]));
			fprintf (fpout, " %+.8E", tau22 / (utau[d] * utau[d]));
			fprintf (fpout, "\n");

			/***************************** END SAVE VALUES TO FILE *****************************/

		}

		fclose(fpout);

		chdir ("..");

		/***************************** END READ RESULTS *****************************/

	}
	
	/*amelia_inicio*/

	sprintf (outfile, "results_utau.dat");

	fpout = fopen (outfile, "w");

	fprintf (fpout, "# x/h      utau\n");

	for (d = 0; d < ndirs; d++)
	{
		xh_aux[d] = xh[d];
	}

	sort(xh_aux, utau, ndirs);

	for (d = 0; d < ndirs; d++)
	{
		fprintf (fpout, " %+.8E %+.8E\n", xh_aux[d], utau[d]);    
	}    

	fclose(fpout);
		
	/*amelia_fim*/

	sprintf (outfile, "results_tauw.dat");

	fpout = fopen (outfile, "w");

	fprintf (fpout, "# x/h      tauw\n");

	for (d = 0; d < ndirs; d++)
	{
		xh_aux[d] = xh[d];
	}

	sort(xh_aux, tauw, ndirs);

	for (d = 0; d < ndirs; d++)
	{
		fprintf (fpout, " %+.8E %+.8E\n", xh_aux[d], tauw[d]);    
	}    

	fclose(fpout);

	sprintf (outfile, "results_cf.dat");

	fpout = fopen (outfile, "w");

	fprintf (fpout, "# x/h      cf\n");

	for (d = 0; d < ndirs; d++)
	{
		xh_aux[d] = xh[d];
	}

	sort(xh_aux, cf, ndirs);

	for (d = 0; d < ndirs; d++)
	{
		fprintf (fpout, " %+.8E %+.8E\n", xh_aux[d], cf[d]);    
	}    

	fclose(fpout);
	
	for (d = 0; d < ndirs - 1; d++)
	{
		if(cf[d+1] > 0 && cf[d] < 0)
		{
			// line passing through two points
			m = (cf[d+1] - cf[d]) / (xh_aux[d+1] - xh_aux[d]);
			b = cf[d] - m * xh_aux[d];

			if (m != 0)
				rlength = -b / m;
			else
				rlength = 0.0;

			break;
		}
	}    

	printf ("\nReattachment length = %f\n", rlength);
	
	printf ("\n");

	// Exit program
	return (0);

}
