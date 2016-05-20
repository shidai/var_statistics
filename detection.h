#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
////#include <gsl/gsl_rng.h>
////#include <gsl/gsl_randist.h>
//#include <fftw3.h>
#include "T2toolkit.h"
//#include "cpgplot.h"
//#include <omp.h>

typedef struct acfStruct {
	int n; // number of dynamic spectrum

	double cFlux;  // pulsar flux density
	double whiteLevel;  // white noise level

	int nchn;
	int nsubint;

	float **dynPlot; // dynamic spectrum for pgplot
	float *psrMean;
	float *psrVar; 

	float snr_mean;
	float snr_var;
	float ratio;
} acfStruct;

typedef struct noiseStruct {
	int npixel;
	int n; // number of realization
	int nchn;
	int nsubint;

	float **noisePlot; // npixel*(nchn*nsubint)
	float *noiseMean; // mean of noise averaged over npixel for each realization
	float *noiseVar; // variance of noise averaged over npixel for each realization
	double whiteLevel;  // white noise level
} noiseStruct;

typedef struct controlStruct {
	int npixel; 
	int n; // number of dynamic spectrum
	//char oname[1024]; // output file name
	char dname[1024]; // graphics device

	int nsub;  // number of subintegrations
	int nchan; // number of subchannels

	double whiteLevel;   // white noise level, mJy
	double cFlux;        // flux density of pulsars, mJy

	int noplot;
}controlStruct;

void deallocateMemory (acfStruct *acfStructure, noiseStruct *noiseStructure);
float find_peak_value (int n, float *s);
int readParams(char *fname, char *dname, int n, controlStruct *control);
void initialiseControl(controlStruct *control);

//void heatMap (acfStruct *acfStructure, char *dname);
//void palett(int TYPE, float CONTRA, float BRIGHT);
//int plotDynSpec (char *pname, char *dname);

int qualifyVar (acfStruct *acfStructure, noiseStruct *noiseStructure, controlStruct *control);
float chiSquare (float *data, int n, float noise);
float moduIndex (float *data, int n);
int histogram (float *data, int n, float *x, float *val, float low, float up, int step);
float find_max_value (int n, float *s);
float find_min_value (int n, float *s);

int calNoise (noiseStruct *noiseStructure, controlStruct *control);
int simNoise (noiseStruct *noiseStructure, long seed);
int calNoiseMean (noiseStruct *noiseStructure, int n);
float variance (float *data, int n);
float mean (float *data, int n);
int calPsr (acfStruct *acfStructure, controlStruct *control);
int simPsr (acfStruct *acfStructure, long seed, int k);

void readNoise (char *Fname, double *fdiss);
int readNum (char *Tname);

int calNoise (noiseStruct *noiseStructure, controlStruct *control)
{
	int nchn, nsubint;
	long seed;
	int i;
	int n, npixel;

	noiseStructure->npixel = control->npixel; 
	noiseStructure->n = control->n; 
	noiseStructure->nchn = control->nchan; 
	noiseStructure->nsubint = control->nsub; 

	nchn = noiseStructure->nchn;
	nsubint = noiseStructure->nsubint;

	noiseStructure->whiteLevel = sqrt(nchn*nsubint)*control->whiteLevel; // mJy
	n = control->n;  // number of realization
	npixel = control->npixel; // number of pixel to average
	// allocate memory
	noiseStructure->noisePlot = (float **)malloc(sizeof(float *)*npixel);
	noiseStructure->noiseMean = (float *)malloc(sizeof(float)*n);
	noiseStructure->noiseVar = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < npixel; i++)
	{
		noiseStructure->noisePlot[i] = (float *)malloc(sizeof(float)*nsubint*nchn);
	}

	// simulate noise
	for (i=0; i<noiseStructure->n; i++)
	{
		seed = TKsetSeed();
		simNoise (noiseStructure, seed);
		calNoiseMean (noiseStructure, i);
	}

	return 0;
}

int simNoise (noiseStruct *noiseStructure, long seed)
{
	int nchn = noiseStructure->nchn;
	int nsubint = noiseStructure->nsubint;
	int npixel = noiseStructure->npixel;

	int i, j, k;

	for (k = 0; k < npixel; k++)
	{
		for (i = 0; i < nchn; i++)
		{
			for (j = 0; j < nsubint; j++)
			{
				//noiseStructure->noisePlot[k][i*nsubint+j] = 16*(float)(TKgaussDev(&seed));   // create noise image pixels
				noiseStructure->noisePlot[k][i*nsubint+j] = (float)(noiseStructure->whiteLevel*TKgaussDev(&seed));   // create noise image pixels
			}
		}
	}

	return 0;
}

int calNoiseMean (noiseStruct *noiseStructure, int n)
{
	int nchn = noiseStructure->nchn;
	int nsubint = noiseStructure->nsubint;
	int npixel = noiseStructure->npixel;

	int i;
	float meanVal, varVal;

	meanVal = 0.0;
	varVal = 0.0;
	for (i = 0; i < npixel; i++)
	{
		meanVal += mean (noiseStructure->noisePlot[i], nchn*nsubint);   // create noise image pixels
		varVal += variance (noiseStructure->noisePlot[i], nchn*nsubint);   // create noise image pixels
	}
	meanVal = meanVal/npixel;
	varVal = varVal/npixel;

	noiseStructure->noiseMean[n] = meanVal;
	noiseStructure->noiseVar[n] = varVal;

	return 0;
}

int calPsr (acfStruct *acfStructure, controlStruct *control)
{
	//FILE *fin;
	int nsubint, nchn;
	int n; // number of realization
	long seed;
	int i;

	//printf ("Starting simulating dynamic spectrum\n");
	acfStructure->n = control->n; 
	acfStructure->cFlux = control->cFlux; // mJy
	acfStructure->nchn = control->nchan;
	acfStructure->nsubint = control->nsub;

	nchn = acfStructure->nchn;
	nsubint = acfStructure->nsubint;
	n = acfStructure->n;

	acfStructure->whiteLevel = sqrt(nchn*nsubint)*control->whiteLevel; // mJy
	// allocate memory
	acfStructure->dynPlot = (float **)malloc(sizeof(float *)*n);
	acfStructure->psrMean = (float *)malloc(sizeof(float)*n);
	acfStructure->psrVar = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		acfStructure->dynPlot[i] = (float *)malloc(sizeof(float)*nsubint*nchn);
	}


	for (i=0; i<acfStructure->n; i++)
	{
		seed = TKsetSeed();
		simPsr (acfStructure, seed, i);
		acfStructure->psrMean[i] = mean (acfStructure->dynPlot[i], nsubint*nchn);
		acfStructure->psrVar[i] = variance (acfStructure->dynPlot[i], nsubint*nchn);
	}

	return 0;
}

int simPsr (acfStruct *acfStructure, long seed, int k)
{
	int nchn = acfStructure->nchn;
	int nsubint = acfStructure->nsubint;
	float re, im, noise;

	int i, j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsubint; j++)
		{
			re = TKgaussDev(&seed);
			im = TKgaussDev(&seed);
			noise = TKgaussDev(&seed);
			//printf ("%f %f %f\n", re, im, noise);
			//acfStructure->dynPlot[k][i*nsubint+j] = (float)(2.0*(pow(re,2.0)+pow(im,2.0))+16*noise);  
			//acfStructure->dynPlot[k][i*nsubint+j] = (float)(2.0*(pow(re,2.0)+pow(im,2.0))+acfStructure->whiteLevel*noise);  
			acfStructure->dynPlot[k][i*nsubint+j] = (float)(0.5*(pow(re,2.0)+pow(im,2.0))+acfStructure->whiteLevel*noise);  
		}
	}

	return 0;
}

void deallocateMemory (acfStruct *acfStructure, noiseStruct *noiseStructure)
{
	int npixel = noiseStructure->npixel; // number of dynamic spectrum
	int n = acfStructure->n; // number of dynamic spectrum
	
	free(acfStructure->psrMean);
	free(acfStructure->psrVar);

	int i;

	for (i = 0; i < n; i++)
	{
		free(acfStructure->dynPlot[i]);
	}
	//free(acfStructure->dynPlot);
	
	// free noise memory
	free(noiseStructure->noiseMean);
	free(noiseStructure->noiseVar);
	for (i = 0; i < npixel; i++)
	{
		free(noiseStructure->noisePlot[i]);
	}
}

float find_peak_value (int n, float *s)
{
	int i;
	float temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	return temp[n-1];
}
						
int readParams(char *fname, char *dname, int n, controlStruct *control)
//int readParams(char *fname, char *oname, char *dname, int n, controlStruct *control)
{
	FILE *fin;
	char param[1024];
	int endit=-1;
	int finished=0;

	// define the output file name
	//strcpy(control->oname,oname);
	strcpy(control->dname,dname);

	control->n = n;

	///////////////////////////////////////
	if ((fin=fopen(fname,"r"))==NULL)
	{
		printf ("Can't open file!\n");
		exit(1);
	}

	//printf("Reading parameters...\n");

      	// Find the start observation
	while (!feof(fin))
	{
		if (fscanf(fin,"%s",param)==1)
		{
			if (strcasecmp(param,"START_OBS")==0)
			{
				endit=0;
				break;
			}
		}
		else 
			return 1;
	}

	if (endit==-1)
		return 1;

	do
	{
		fscanf(fin,"%s",param);
		if (strcasecmp(param,"END_OBS")==0)
			endit=1;
		else
		{
			if (strcasecmp(param,"NCHAN")==0)
			      	fscanf(fin,"%d",&(control->nchan));
			else if (strcasecmp(param,"NSUB")==0)
      				fscanf(fin,"%d",&(control->nsub));
			//else if (strcasecmp(param,"WHITE_LEVEL")==0)
      			//	fscanf(fin,"%lf",&(control->whiteLevel));
			else if (strcasecmp(param,"CFLUX")==0)
			      	fscanf(fin,"%lf",&(control->cFlux));
			else if (strcasecmp(param,"NPIXEL")==0)
			      	fscanf(fin,"%d",&(control->npixel));
		}
	} while (endit==0);

	if (fclose(fin))
	{
		printf ("Can't close file.\n");
		exit(1);
	}

	return finished;
}

void initialiseControl(controlStruct *control)
{
	strcpy(control->dname,"1/xs");
	
	control->n = 1; // simulate 1 dynamic spectrum by default

	// Standard defaults
	control->nchan = 32;
	control->nsub = 8;
	control->whiteLevel = 0;   // mJy
	control->cFlux = 0.0;   // mJy
}

//void heatMap (acfStruct *acfStructure, char *dname)
//{
//	//int i,j;                     
//	//int dimx = acfStructure.ns;
//	//int dimy = acfStructure.nf; // dimensions 
//	//float tab[dimx*dimy];       // value
//	char caption[1024];
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", acfStructure->cFreq, "MHz", "BW:", acfStructure->bw, "MHz", "Length:", acfStructure->tint, "s");
//  
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	int dimx = acfStructure->nsubint;
//	int dimy = acfStructure->nchn;
//	double bw = acfStructure->bw;
//  
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//
//	double f1 = acfStructure->cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = acfStructure->cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//printf ("f1 f2: %lf %lf\n", f1, f2);
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, acfStructure->dynPlot[0]);
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
//
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,dname,1,1);
//	//cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	//cpgswin(0,dimx,164,132); // set window
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//	//cpgbox("BCTSIN",10,5,"BCTSIN",50,5);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//      	//cpglab("Subintegration","Frequency (MHz)","Freq: 150.0 MHz BW: -32.000 MHz Length: 960.0 s");
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(acfStructure->dynPlot[0],dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(acfStructure->dynPlot,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//} 
//
//void palett(int TYPE, float CONTRA, float BRIGHT)
//{
////-----------------------------------------------------------------------
//// Set a "palette" of colors in the range of color indices used by
//// PGIMAG.
////-----------------------------------------------------------------------
//	float GL[] = {0.0, 1.0};
//	float GR[] = {0.0, 1.0};
//	float GG[] = {0.0, 1.0};
//	float GB[] = {0.0, 1.0};
//	float RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
//	float RR[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
//	float RG[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
//	float RB[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
//	float HL[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float HR[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float HG[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float HB[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	float WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
//	float WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
//	float WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
//	float WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
//	float AL[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
//	float AR[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//	float AG[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
//	float AB[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//      
//	if (TYPE == 1)
//	{   
//		//-- gray scale
//		cpgctab(GL, GR, GG, GB, 2, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 2) 
//	{
//		//-- rainbow
//		cpgctab(RL, RR, RG, RB, 9, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 3) 
//	{
//		//-- heat
//		cpgctab(HL, HR, HG, HB, 5, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 4) 
//	{
//		//-- weird IRAF
//		cpgctab(WL, WR, WG, WB, 10, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 5) 
//	{
//		//-- AIPS
//		cpgctab(AL, AR, AG, AB, 20, CONTRA, BRIGHT);
//	}
//}

//int plotDynSpec (char *pname, char *dname)
//{
//	FILE *fin;
//	char start[128];
//	double tint,bw,cFreq;
//	int nsub,nchn;
//	float val;
//	int n1, n2;
//	float *dynSpec;
//
//	int i;
//
//	int dimx;
//	int dimy;
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	char caption[1024];
//
//	if ((fin=fopen(pname,"r"))==NULL)
//	{
//		printf ("Can't open dynamic spectrum!\n");
//		exit(1);
//	}
//
//      	// Find the start of dynamic spectrum, which contains basic info
//	while (!feof(fin))
//	{
//		if (fscanf(fin,"%s %d %d %lf %lf %lf",start,&nsub,&nchn,&bw,&tint,&cFreq)==6)
//		{
//			if (strcasecmp(start,"START")==0)
//			{
//				break;
//			}
//		}
//		//else 
//		//	return 1;
//	}
//
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", cFreq, "MHz", "BW:", bw, "MHz", "Length:", tint, "s");
//
//	dynSpec = (float*)malloc(sizeof(float)*nsub*nchn);
//
//	i = 0;
//	while (fscanf(fin,"%d %d %f", &n1, &n2, &val)==3)
//	{
//		dynSpec[i] = val;
//		i++;
//	}
//
//	if (fclose(fin))
//	{
//		printf ("Can't close dynamic spectrum!\n");
//		exit(1);
//	}
//
//	dimx = nsub;
//	dimy = nchn;
//  
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, dynSpec);
//
//	double f1 = cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
// 
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,dname,1,1);
//	//cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//
//	free(dynSpec);
//
//	return 0;
//} 

int qualifyVar (acfStruct *acfStructure, noiseStruct *noiseStructure, controlStruct *control)
{
	int i;
	int n = acfStructure->n;

	float *source_amp, *source_var;

	float meanM, varM;
	float meanV, varV;

	source_amp = (float *)malloc(sizeof(float)*n);
	source_var = (float *)malloc(sizeof(float)*n);

	for (i=0; i<n; i++)
	{
		source_amp[i] = acfStructure->psrMean[i] - noiseStructure->noiseMean[i];
		source_var[i] = acfStructure->psrVar[i] - noiseStructure->noiseVar[i];
	}

	meanM = mean (source_amp, n);
	varM = variance (source_amp, n);
	meanV = mean (source_var, n);
	varV = variance (source_var, n);

	acfStructure->snr_mean = meanM/sqrt(varM);
	acfStructure->snr_var = meanV/sqrt(varV);
	acfStructure->ratio = (meanM/sqrt(varM))/(meanV/sqrt(varV));
	/*
	maxV = find_max_value(n,var);
	minV = find_min_value(n,var);
	maxVN = find_max_value(n_n,var_n);
	minVN = find_min_value(n_n,var_n);
	////////////////////////////
	// make histogram
	if (control->noplot != 1)
	{
		cpgbeg(0,control->dname,1,1);

		xHis = (float*)malloc(sizeof(float)*step);
		val = (float*)malloc(sizeof(float)*step);
		xHisN = (float*)malloc(sizeof(float)*step);
		valN = (float*)malloc(sizeof(float)*step);
		histogram (var, n, xHis, val, minV-2*(maxV-minV)/step, maxV+2*(maxV-minV)/step, step);
		histogram (var_n, n_n, xHisN, valN, minVN-2*(maxVN-minVN)/step, maxVN+2*(maxVN-minVN)/step, step);

      		cpgsch(1); // set character height
      		cpgscf(1); // set character font

		// find the max and min
		max1 = find_max_value(step,val);
		max2 = find_max_value(step,valN);
		max = (max1 >= max2 ? max1 : max2);
      		//cpgenv(-5,5,0,4500,0,1); // set window and viewport and draw labeled frame
      		cpgenv(minVN-1, maxV+1, 0, max+0.1*max, 0, 1); // set window and viewport and draw labeled frame

		sprintf(caption, "%s", "Variance histogram");
      		cpglab("Variance","Number",caption);
		cpgbin(step,xHis,val,0);
		cpgsci(2);
		cpgbin(step,xHisN,valN,0);
		cpgsci(3);
		cpgline(100, xThreshold, yThreshold);
		///////////////////////////////////////////////////////
		cpgend();

		free(xHis);
		free(val);
		free(xHisN);
		free(valN);
	}
	*/

	free(source_amp);
	free(source_var);

	return 0;
}

float chiSquare (float *data, int n, float noise)
{
	int i;

	float ave;
	float chiS;
	
	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	chiS = 0.0;
	for (i=0; i<n; i++)
	{
		chiS += pow(data[i]-ave,2)/pow(noise,2);
	}

	return chiS/(n-1);
}

float moduIndex (float *data, int n)
{
	int i;

	float ave, devi;
	float m;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	devi = 0.0;
	for (i=0; i<n; i++)
	{
		devi += pow(data[i]-ave,2);
	}
	devi = sqrt(devi/n);

	m = devi/ave;

	return m;
}

float variance (float *data, int n)
{
	int i;

	float ave, devi;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	devi = 0.0;
	for (i=0; i<n; i++)
	{
		devi += pow(data[i]-ave,2);
	}
	devi = devi/n;

	return devi;
}

float mean (float *data, int n)
{
	int i;

	float ave;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	return ave;
}

int histogram (float *data, int n, float *x, float *val, float low, float up, int step)
{
	int i,j,count;
	float width;
	float *temp;

	temp = (float*)malloc(sizeof(float)*(step+1));

	width = (up-low)/step;
	for (i=0; i<step; i++)
	{
		x[i] = low + i*width + width/2.0;
	}

	for (i=0; i<=step; i++)
	{
		temp[i] = low + i*width;
	}

	for (i=0; i<step; i++)
	{
		count = 0;
		for (j=0; j<n; j++)
		{
			if (data[j]>=temp[i] && data[j]<temp[i+1])
			{
				count += 1;
			}
		}
		//val [i] = count;
		val[i] = (float)(count)/n;
	}

	free(temp);
	return 0;
}

float find_max_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

float find_min_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a <= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

int readNum (char *Tname)
{
	int nt;
	FILE *fint;
	double temp;

	if ((fint=fopen(Tname, "r"))==NULL)
	{
		printf ("Can't open scintillation time-scale...\n");
		exit(1);
	}

	nt = 0;
	while (fscanf(fint, "%lf", &temp) == 1)
	{
		nt++;
	}

	if (fclose(fint))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}

	return nt;
}



void readNoise (char *Fname, double *fdiss)
{
	int i;
	FILE *finf;
	double temp;

	///////////////////////////////////////////////////////////////////////
	if ((finf=fopen(Fname, "r"))==NULL)
	{
		printf ("Can't open flux...\n");
		exit(1);
	}

	i = 0;
	while (fscanf(finf, "%lf", &temp) == 1)
	{
		fdiss[i] = temp;
		i++;
	}

	if (fclose(finf))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}
}
