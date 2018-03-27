/*

SMFORMULA.CPP, based on

 HR2.C
 V1.02

 A program to calculate elemental compositions for a given mass.
 See the file README for details.

--------------------------------------------------------------------
 Copyright (c) 2001...2005 Joerg Hau <joerg.hau(at)dplanet.ch>.

 mail: joerg.hau@dplanet.ch
 www:  http://www.mysunrise.ch/users/joerg.hau/

 *changed version by Tobias Kind (TK), 2006 , Fiehnlab,
 *added extended valencies, added implementation of
  seven golden rules of molecular formula filtering
 
 *changed version by Robert Winkler (RW), CINVESTAV, 2013-2017,
 robert.winkler@cinvestav.mx, robert.winkler@bioprocess.org 
 *added csv export of results.
 *integration into SpiderMass for generation of molecular formula.
 *Heavy isotope names changed for compatibility with smisotope.c

 This program is free software; you can redistribute it and/or
 modify it under the terms of version 2 of the GNU General Public
 License as published by the Free Software Foundation. See the
 file LICENSE for details.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
--------------------------------------------------------------------

 Creation:	somewhere in 1992 by JHa.
 Revision:  2001-04-18, GPL'd, first public release (JHa)
            2001-04-21, improved help text (JHa)
            2002-06-27, added sodium (JHa)
            2002-10-09, added 15N (JHa)
            2005-02-25, added -v option; license now GPL *v2* (JHa)
            2005-02-27, optimised code in calc loop (JHa)
            2005-02-28, verified and updated atomic masses (JHa)
            2005-06-17, added GPL text when "-h" is used (JHa)
			2006-01-01, extended version for BMC Bioinformatics publication - HR2 (TK)
			2006-03-03, added element ratio checks, extended valencies, only even electrons - HR2 (TK)
			2006-09-09,	1000x-10000x speedup hand optimized hehe. - HR2 (TK)
						-->special version for CHNSOP-F-Cl-Br-Si 
			2013-10-16, added csv export for integration in SpiderMass (RW)
			2013-10-17, accurate masses for elements revised (NIST,2013), as well as CNOPS/ element ratios (RW)
			2014-02-21, heavy isotope names changed to 1-letter code: 13C->X, 15N->M
			2017-09-04, revision of formula generation with 2H
 This is ANSI C and should compile with any C compiler; use
 something along the lines of "gcc -Wall -O3 -o hr hr.c".
 Optimize for speed, you may gain factor 3!
 NOW compiled under Visual C++ Express (faster than GCC) in C++ mode for boolean type.


 ---------------------------------------------------------------------
 Example arguments:
 1) -m 1 -t 100000 -C 1-100 -H 1-220 -N 0-10 -O 0-10 -P 0-10 -S 0-10 -L 0-10 -B 0-10
 2) -m 500 -t 1 -C 50-100 -H 10-220 -N 0-10 -O 0-10 -P 0-10 -S 0-10 -L 0-10 -B 0-10
 3) hr2-all-res -c "Hexaflumuron_459.982882Da_3ppm" -m 459.982882 -t 1.37995 -C 0-39 -H 0-98 -N 0-34 -O 0-30 -P 0-12 -S 0-12 -F 0-12 -L 0-14 -B 0-7 -I 0-0 
	945 formulas found in    253 seconds. (now 4 seconds, before eternal)
 4) hr2 -m 459.982882 -t 1.37995 -C 10-39 -H 28-98 -N 4-34 -O 0-30 -P 1-12 -S 1-12 -F 0-12 -L 1-14 -B 2-6 -I 0-0
	1 formula in 0 seconds (former eternal time)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sstream>
#include <math.h>
#include <fstream>
#include <unistd.h>
#include <iostream>
using namespace std; //RW

#define VERSION "20170904"	/* String ! */
#define TRUE 	1
#define FALSE 	0
#define MAXLEN  181          /* max. length of input string */

#define _CRT_SECURE_NO_DEPRECATE 1

typedef struct 	{
		const char *sym;	/* symbol */
		const double mass;	/* accurate mass */
		const float val;	/* to calculate unsaturations */
		const int key;		/* used for decoding cmd line */
		int min,		/* atom count min */
		    max,		/* atom count max */
		    cnt,		/* atom count actual */
		    save;		/* atom count old  - for loop exiting*/
		} Element;


/* --- atomic masses as published by IUPAC, 2002-10-02 ---------- */

Element el[]=
/* array of the elements used here:
   Symbol, exact mass, dbe, keycode, number, default-min, default-max
   RW: Actualization of exact masses, values from NIST, 2013
*/
// ele |    mass   |  dbe | key| min | max | cnt | save
{
{ "C",  12.000000000,   +2.0, 'C', 0, 41, 0,0 },
{ "X", 13.0033548378, +2.0, '1', 0, 0, 0 ,0 }, //13C
{ "H",   1.0078250321,  -1.0, 'H', 0, 72, 0 ,0},
{ "D",   2.0141017778,  -1.0, 'D', 0, 0, 0 ,0}, //2H
{ "N",  14.0030740048,  +1.0, 'N', 0, 34, 0,0 },		//org +1 = valence = 3: now +3 for valence = 5
{ "M", 15.0001088982,   +1.0, 'M', 0, 0, 0,0 }, //15N
{ "O",  15.9949146196,   0.0, 'O', 0, 30, 0 ,0},
{ "F",  18.99840322,    -1.0, 'F', 0, 0, 0 ,0},
{ "Na", 22.9897692809,    -1.0, 'A', 0, 0, 0 ,0},
{ "Si", 27.9769265325,  +2.0, 'I', 0, 0, 0 ,0},	
{ "P",  30.97376163,    +3.0, 'P', 0, 0, 0 ,0},		//org +1 valence = 3: now +3 for valence = 5
{ "S",  31.972071,    +4.0, 'S', 0, 0, 0 ,0},		//org 0 = valence = 2; now +4 for valence = 6
{ "Cl", 34.96885268,    -1.0, 'L', 0, 0, 0 ,0},
{ "Br", 78.9183371,     -1.0, 'B', 0, 0, 0 ,0},
};

const double electron = 0.000549;	/* mass of the electron in amu */


/* --- global variables --- */

double  charge,		/* charge on the molecule */
        tol;		/* mass tolerance in mmu */
char    comment[MAXLEN]="";	/* some text ;-) */
int     single;		/* flag to indicate if we calculate only once and exit */
int     nr_el;		/* number of elements in array (above) */


int     input(char *text, double *zahl);
int     readfile(char *whatfile);
double  calc_mass(void);
float   calc_rdb(void);
long     do_calculations(double mass, double tolerance);
int     clean (char *buf);
//you have to compile with C++ or define yourself this bool type (C99 compiler definition)
bool calc_element_ratios(bool element_probability);

/* --- threading ------------------- */
/* mass and RDB calculation could be in several other threads
however the loop is very fast, context switching takes longer and
slows down the whole process. Single-Thread multiprocessor (cluster) implementation
seems superior here.
*/
int hThread2,hThread3; //RW: HANDLE changed to int for cross-compilation
unsigned threadID2,threadID3;


/* --- main --- */

int main (int argc, char *argv[])  //RW: optionsS type changed for cross-compilation

{
double mz;	/* mass */
char buf[MAXLEN];
int i, tmp;

static const char *id =
"hr version %s. Copyright (C) by Joerg Hau 2001...2005, Tobias Kind 2006 :-) & Robert Winkler 2013...2017 ;-).\n";

static const char *disclaimer =
"\nD=2H, X=13C, M=15N \n"
"\nThis program is free software; you can redistribute it and/or modify it under\n"
"the terms of version 2 of the GNU General Public License as published by the\n"
"Free Software Foundation.\n\n"
"This program is distributed in the hope that it will be useful, but WITHOUT ANY\n"
"WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n"
"PARTICULAR PURPOSE. See the GNU General Public License for details.\n";

static const char *msg =
"Calculates possible elemental compositions for a given mass.\n\n"
"usage: hr [options] file\n\nValid command line options are:\n"
"-h      This Help screen.\n"
"-v      Display version information.\n"
"-t tol  Set tolerance to 'tol' mmu (default 5).\n"
"-m mz   Set mass to 'mz'.\n"
"-c txt  Set comment to 'txt' (only useful together with '-m').\n"
"-p      Positive ions; electron mass is removed from the formula.\n"
"-n      Negative ions; electron mass is added to the formula.\n"
"-X a-b  For element X, use atom range a to b. List of valid atoms:\n\n"
"           X    key   mass (6 decimals shown)\n"
"        -------------------------------------\n";

/* initialise variables */

single = FALSE;			/* run continuously */
charge = 0.0;	       	 	/* default charge is neutral */
tol = 5.0;			/* default tolerance in mmu */
nr_el = sizeof(el)/sizeof(el[0]);	/* calculate array size */


/* decode and read the command line */

while ((tmp = getopt(argc, argv, "hvpnt:m:c:C:H:N:M:O:D:1:S:F:L:B:P:I:A:")) != EOF)
	switch (tmp)
		{
		case 'h':     	  		/* help me */
			printf (id, VERSION);
            printf ("%s",msg);
			for (i=0; i < nr_el; i++)
				printf ("        %4s     -%c %15.6lf\n",
				el[i].sym, el[i].key, el[i].mass);
			printf(disclaimer, "\n");
			return 0;
		case 'v':     	  		/* version */
			printf (id, VERSION);
			return 0;
		case 'p':    			/* positive charge */
			charge = +1.0;
			continue;
		case 'n':			/* negative carge */
			charge = -1.0;
			continue;
		case 't':       		/* tolerance */
			strcpy(buf, optarg);
			sscanf(buf, "%lf", &tol);
			continue;
		case 'm':			/* single mass */
			strcpy(buf, optarg);
		        sscanf(buf, "%lf", &mz);
		        single = TRUE;
	        	continue;
		case 'c':			/* comment for single mass */
	   		strcpy(comment, optarg);
		        continue;
		case 'C':      		/* C12 */
 		case 'H':      		/* 1H */
		case 'N':      		/* 14N */
		case 'M':      		/* 15N */
		case 'O':      		/* 16O */
 		case 'D':      		/* 2H */
 		case '1':      		/* 13C */
		case 'A':      		/* Na ('N' is taken for Nitrogen!) */
		case 'S':      		/* 32S */
		case 'F':      		/* 19F */
		case 'L':      		/* 35Cl ('C' is taken!) */
		case 'B':      		/* 79Br */
		case 'P':      		/* 31P */
		case 'I':      		/* 28Si ('S' is taken!) */
			i = 0;
			/* compare keys until found */
			while ((i < nr_el) && (el[i].key != tmp))
				i++;
			strcpy(buf, optarg);
			sscanf(buf, "%d-%d", &el[i].min, &el[i].max);	/* copy over */
			if (el[i].min > el[i].max)			/* swap them */
				{
				tmp = el[i].min;
				el[i].min = el[i].max;
				el[i].max = tmp;
				}

			// printf ("\n %c = %c ... %s (%d-%d)", tmp, el[i].key, el[i].sym, el[i].min, el[i].max);

			continue;
		case '~':    	  	/* invalid arg */
		default:
			printf ("'%s -h' for help.\n", argv[0]);
			return 1;
		}

if (argv[optind] != NULL)	 /* remaining parameter on cmd line? */
	/* must be a file -- treat it line by line */
	return (readfile (argv[optind]));

if (single == TRUE)  	   	 	/* only one calculation requested? */
	do_calculations(mz, tol);       /* do it, then exit ... */
else
	{				/* otherwise run a loop */
	while (input(comment, &mz))
		{
		tmp = do_calculations(mz, tol);
		printf("\n");
		}
	}

return 0;
}


/***************************************************************************
* INPUT:	reads a dataset in "dialog mode".			   *
* Input: 	Pointer to comment text, pointer to mass.		   *
* Returns:	Result of sscanf(), 0 if prog is to be finished. 	   *
* Note:		This is a fairly primitive way to do it, but it works ;-)  *
****************************************************************************/
int input(char *txt, double *zahl)
{
char buf[MAXLEN];				/* input line */

*zahl = 0.0;				    /* reset */

printf("\n\nComment             : ");	/* display prompt */

char* checkreadline = fgets(buf, MAXLEN-1, stdin);	/* read line */
if(checkreadline == ""){
  // Reading the line failed
}

buf[MAXLEN] = 0x0;				/* terminate string */
clean (buf);				    /* remove linefeed */
strcpy(txt, buf);			    /* copy text over */

printf("Mass (ENTER to quit): ");	/* display prompt */

char* checkreadline2 = fgets(buf, MAXLEN-1, stdin);	/* read line */
if(checkreadline2 == ""){
  // Reading the line failed
}

buf[MAXLEN] = 0x0;			    /* terminate string */
if (!clean (buf))			    /* only a CR ? --> quit */
	return 0;
sscanf(buf,"%lf", zahl);		/* scan string */

return 1;
}

/***************************************************************************
* READFILE:	reads dataset from file.				   *
* Input: 	Pointer to comment text, pointer to mass.		   *
* Returns:	0 if OK, 1 if error.					   *
****************************************************************************/
int readfile(char *whatfile)
{
double mz;		/* measured mass */
char buf[MAXLEN];		/* input line */
FILE *infile;

infile = fopen(whatfile, "r");
if (NULL == infile)
	{
	fprintf (stderr, "Error: Cannot open %s.", whatfile);
	return 1;
	}

while (fgets(buf, MAXLEN-1, infile))
	{
	buf[MAXLEN] = 0x0;			/* terminate string */
	if (*buf == ';')		/* comment line */
		continue;
	if (!clean (buf))		/* only a CR ? --> quit */
		return 0;
	sscanf(buf,"%s %lf", comment, &mz);	/* scan string */
	do_calculations(mz, tol);
	mz = 0.0;				/* reset */
	}
return 0;
}


/************************************************************************
* CALC_MASS:	Calculates mass of an ion from its composition.	  	*
* Input: 	nothing (uses global variables) 	      		*
* Returns. 	mass of the ion.	  				*
* Note:		Takes care of charge and electron mass!   		*
* 		(Positive charge means removal of electrons).	 	*
*************************************************************************/
double calc_mass(void)
{
int i;
double sum = 0.0;

for (i=0; i < nr_el; i++)
	sum += el[i].mass * el[i].cnt;

return (sum - (charge * electron));
}


/************************************************************************
* CALC_RDB:	Calculates rings & double bond equivalents.    		*
* Input: 	nothing (uses global variables)			   	*
* Returns. 	RDB.				       			*
*************************************************************************/
float calc_rdb(void)
{
int i;
float sum = 2.0;

for (i=0; i < nr_el; i++)
	sum += el[i].val * el[i].cnt;

return (sum/2.0);
}
/************************************************************************
* Calculates element ratios , CH2 (more than 8 electrons needed is not handled)  		
* Calculations element probabilities if element_probability = true 
* Input: 	nothing (uses global variables)			   	
* Returns. true/false.				       			
*************************************************************************/
bool calc_element_ratios(bool element_probability)
{
bool CHNOPS_ok;	
float HC_ratio;
float NC_ratio;
float OC_ratio;
float PC_ratio;
float SC_ratio;

float C_count = (float)el[0].cnt+(float)el[1].cnt; //RW added isotopes
float H_count = (float)el[2].cnt+(float)el[3].cnt; //RW added isotopes
float N_count = (float)el[4].cnt+(float)el[5].cnt; //RW added isotopes
float O_count = (float)el[6].cnt;
float P_count = (float)el[10].cnt;
float S_count = (float)el[11].cnt;


//RW ELEMENT RATIOS and CNOPS adjusted, according to Kind & Fiehn, 2007

		/* ELEMENT RATIOS allowed
			MIN		MAX (99.99%)
		H/C	0.1		6.00
		N/C	0.00	4.00
		O/C	0.00	3.00
		P/C	0.00	2.00
		S/C	0.00	3.00
		*/	

//RW Probability check for common range (covering 99.7%)

	// set CHNOPS_ok = true and assume all ratios are ok
	CHNOPS_ok = true;	
	
	
	if (C_count && H_count >0)					// C and H  must have one count anyway (remove for non-organics//
	{	
		HC_ratio = H_count/C_count;
		if (element_probability)
		{
			if ((HC_ratio <  0.2) || (HC_ratio >  3.1)) // this is the H/C probability check ;
			CHNOPS_ok = false;
		}
		else if (HC_ratio >  6.0) // this is the normal H/C ratio check - type cast from int to float is important
			CHNOPS_ok = false;
	}

	if (N_count >0)	// if positive number of nitrogens then thes N/C ratio else just calc normal
	{
		NC_ratio = N_count/C_count;
		if (element_probability)
		{
			if (NC_ratio >  1.3) // this is the N/C probability check ;
			CHNOPS_ok = false;
		}
		else if (NC_ratio >  4.0)
			CHNOPS_ok = false;
	}	
	
	if (O_count >0)	// if positive number of O then thes O/C ratio else just calc normal
	{	
		OC_ratio = O_count/C_count;
		if (element_probability)
		{
			if (OC_ratio >  1.2) // this is the O/C  probability check ;
			CHNOPS_ok = false;		
		}
		else if (OC_ratio >  3.0)
				CHNOPS_ok = false;
	}	


	if (P_count >0)	// if positive number of P then thes P/C ratio else just calc normal
	{	
		PC_ratio = 	P_count/C_count;
		if (element_probability)
		{
			if (PC_ratio >  0.3) // this is the P/C  probability check ;
			CHNOPS_ok = false;	
		
		}
		else if (PC_ratio >  2.0)
			CHNOPS_ok = false;
	}	

	if (S_count >0)	// if positive number of S then thes S/C ratio else just calc normal
	{	
		SC_ratio = 	S_count/C_count;
		if (element_probability)
		{
			if (SC_ratio >  0.8) // this is the S/C  probability check ;
			CHNOPS_ok = false;	
		}
		else if (SC_ratio >  3.0)
			CHNOPS_ok = false;
	}	

//-----------------------------------------------------------------------------	
		
	// check for multiple element ratios together with probability check 
	//if N<10, O<20, P<4, S<3 then true
	if (element_probability && (N_count > 10) && (O_count > 20) && (P_count > 4) && (S_count > 3))
		CHNOPS_ok = false;	
	
	// NOP check for multiple element ratios together with probability check
	// NOP all > 3 and (N<11, O <22, P<6 then true)
	if (element_probability && (N_count > 3) && (O_count > 3) && (P_count > 3))
		{
		if (element_probability && (N_count > 11) && (O_count > 22) && (P_count > 6))
			CHNOPS_ok = false;	
		}
	
	// OPS check for multiple element ratios together with probability check
	// O<14, P<3, S<3 then true
	if (element_probability && (O_count > 14) && (P_count > 3) && (S_count > 3))
		CHNOPS_ok = false;	

	// PSN check for multiple element ratios together with probability check
	// P<3, S<3, N<4 then true
	if (element_probability && (P_count > 3) && (S_count > 3) && (N_count >4))
		CHNOPS_ok = false;	

	
	// NOS check for multiple element ratios together with probability check
	// NOS all > 6 and (N<19 O<14 S<8 then true)
	if (element_probability && (N_count >6) && (O_count >6) && (S_count >6))
	{
		if (element_probability && (N_count >19) && (O_count >14) && (S_count >8))
			CHNOPS_ok = false;	
	}	


	// function return value;
	if (CHNOPS_ok == true)
		return true;
	else 
		return false;
}

/************************************************************************
* DO_CALCULATIONS: Does the actual calculation loop.			*
* Input: 	   measured mass (in amu), tolerance (in mmu)	    	*
* Returns. 	   number of hits.	       				*
*************************************************************************/
long do_calculations (double measured_mass, double tolerance)
{
time_t start, finish;
double elapsed_time;

double mass;			/* calc'd mass */
double limit_lo, limit_hi;	/* mass limits */
float rdb, lewis;			/* Rings & double bonds */
long i;
long long hit;		/* counts the hits, with long declaration, overflow after 25h with all formulas < 2000 Da
							    long = FFFFFFFFh = 4,294,967,295d*/
long long counter;
bool elementcheck;
bool set_break;


time( &start );		// start time
printf("\n");		/* linefeed */

/* calculate limits */

limit_lo = measured_mass - (tolerance / 1000.0);
limit_hi = measured_mass + (tolerance / 1000.0);

// if (strlen(comment))	/* print only if there is some text to print */
// 	printf ("Text      \t%s\n", comment);

// printf ("Composition\t");
// for (i=0; i < nr_el; i++)
// 	if (el[i].max > 0)
// 		printf("%s:%d-%d ", el[i].sym, el[i].min, el[i].max);
// printf ("\n");

// printf ("Tol (mmu)\t%.1f\n",tolerance);
// printf ("Measured\t%.4lf\n", measured_mass);
// printf ("Charge  \t%+.1lf\n", charge);



/*
RW defining the csv file name and writing the header
*/


// ofstream denovofile; //RW define output file variable
// denovofile.open ("HR3.csv"); //RW define output file name
stringstream hroutstream;   //RW string stream used for the conversion to string and file output
hroutstream << "Formula" << ";" << "RDB" << ";" << "LEWIS"  << ";"  << "Mass_Da" << ";" << "Mass_Error_mDa" << " \n"; //RW
string stringResult;          //RW resulting string variable
stringResult = hroutstream.str(); //RW conversion of the stream to a string
cout << stringResult; //RW writing the string to the file

hit = 0;			/* Reset counter */
counter = 0;
set_break = false;	/* set breaker for element counts to false */

/* Now let's run the big big loop ... I'd like to do that
   recursively but did not yet figure out how ;-) 
   TK Adds: the loop is just fine.
*/

/* now comes the "COOL trick" for calculating all formulae:
sorting the high mass elements to the outer loops, the small weights (H)
to the inner loops;

This will reduce the computational time by factor ~10-60-1000
OLD HR: Cangrelor at 1ppm  4465 formulas found in   5866 seconds.
NEW HR2: Cangrelor at 1ppm 4465 formulas found in     96 seconds.
NEW2 HR2: Cangrelor at 1ppm 4465 formulas found in     60 seconds.
NEW3 HR2: Cangrelor at 1ppm 4465 formulas found in     59 seconds.
HR2 Fast: Cangrelor at 1ppm 4465 formulas found in     41 seconds by evaluating 2,003,436,894 formulae.
hr2 -c "Cangrelor" -m  774.948 -t 0.77 -C 1-64 -H 1-112 -N 0-30 -O 0-80 -P 0-12 -S 0-9 -F 0-10 -L 0-10

Another additional trick is to end the 2nd.. 3rd.. 4th.. xth innermost loop
to prevent loops which are just higher and higher in mass.
*/

el[13].cnt = el[13].min - 1;  el[13].save = el[13].cnt; 
while (el[13].cnt++ < el[13].max) /* "Br"*/ { 

el[12].cnt = el[12].min - 1;  el[12].save = el[12].cnt; 
while (el[12].cnt++ < el[12].max) /*"Cl"*/ { 
	 
el[11].cnt = el[11].min - 1;  el[11].save = el[11].cnt; 
while (el[11].cnt++ < el[11].max) /*"S"*/ { 
	 
el[10].cnt = el[10].min - 1;  el[10].save = el[10].cnt; 
while (el[10].cnt++ < el[10].max) /*"P"*/ { 
	 
el[9].cnt = el[9].min - 1;  el[9].save = el[9].cnt; 
while (el[9].cnt++ < el[9].max) /*"Si"*/ { 

el[8].cnt = el[8].min - 1;  el[8].save = el[8].cnt; 
while (el[8].cnt++ < el[8].max) /*"Na"*/{ 

el[7].cnt = el[7].min - 1;  el[7].save = el[7].cnt; 
while (el[7].cnt++ < el[7].max) /*"F"*/ { 
 
el[6].cnt = el[6].min - 1;  el[6].save = el[6].cnt; 
while (el[6].cnt++ < el[6].max) /*"O"*/ { 
	 
el[5].cnt = el[5].min - 1;  el[5].save = el[5].cnt; 
while (el[5].cnt++ < el[5].max) /*"15N"*/{ 

el[4].cnt = el[4].min - 1; el[4].save = el[4].cnt; 
while (el[4].cnt++ < el[4].max) /*"N"*/{ 
	 
el[1].cnt = el[1].min - 1; el[1].save = el[1].cnt; 
while (el[1].cnt++ < el[1].max) /*"13C"*/ { 

el[0].cnt = el[0].min - 1; el[0].save = el[0].cnt; 
while (el[0].cnt++ < el[0].max) /* "C"*/ { 

el[3].cnt = el[3].min - 1; 	el[3].save = el[3].cnt; 
while (el[3].cnt++ < el[3].max) /*"D"*/{ 

el[2].cnt = el[2].min - 1; el[2].save = el[2].cnt; 
while (el[2].cnt++ < el[2].max) /*"H"*/{ 

	mass = calc_mass();
	counter++;

	//just for debug purposes
	//if (mass > limit_hi)  
	//printf("mass: %f\tC: %d  H: %d  N: %d O: %d P: %d S: %d Cl: %d Br: %d\n",mass,el[0].cnt,el[2].cnt,el[4].cnt,el[6].cnt,el[10].cnt,el[11].cnt,el[12].cnt,el[13].cnt);
    	
	/* if we exceed the upper limit, we can stop the calculation
       for this particular element (JHa 20050227). <-- comment TK that will only bust the innermost while loop, which is "H"*/

	// break H loop 	if (mass > limit_hi)  break;

    //************************************************************************************************************/	
	//Calculus loop with print out
	//************************************************************************************************************/	
	


	if ((mass >= limit_lo) && (mass <= limit_hi)) /* within limits? */
	{	
		// element check will be performed always, if variable bool element_probability is true also probabilities will be calculated
		// not an elegant implementation, but fast.
		 elementcheck = calc_element_ratios(true);
		 if (elementcheck)
	{ 
	rdb = calc_rdb();	/* get RDB */
	lewis = (float)(fmod(rdb, 1)); /*calc remainder*/
	if ((rdb >= 0) && (lewis != 0.5) && (lewis !=-0.5))/* less than -0.5 RDB does not make sense */

		{													/* NO(!) CH3F10NS2 exists , RDB =  -4.0   M= 282.9547*/

		
		hit ++;
		for (i = 0; i < nr_el; i++)	 /* print composition */
		    if (el[i].cnt > 0)	/* but only if useful */
		 
			  {
			  // printf("%s%d.", el[i].sym, el[i].cnt);	//print formula to screen
			  ostringstream hroutstream;   //RW string stream used for the conversion to string and file output
			  hroutstream << el[i].sym << el[i].cnt; //RW generation of formula string
			  string stringResult;          //RW resulting string variable
			  stringResult = hroutstream.str(); //RW conversion of the stream to a string
			  cout << stringResult; //RW writing the string to the stdout
		      }
			  // printf("\t\t%.1f\t%.4lf\t%+.1lf mmu \n", rdb, mass, 1000.0 * (measured_mass - mass));
			  
			  double mass_out = mass;
  		
	  		  float rdb_out = rdb;
	  		  
	  		  float lewis_out = lewis;
	  		  
	  		  ostringstream hroutstream;   //RW stream used for the conversion to string and file output
	  		  hroutstream << ";" << rdb_out << ";" << lewis_out  << ";"  << mass_out << ";" << 1000.0 * (measured_mass - mass) << " \n"; //RW
			  
			  string stringResult;          //RW resulting string variable
			  stringResult = hroutstream.str(); //RW conversion of the stream to a string
			  cout << stringResult; //RW writing the string to the file
			  
  
		}	/* end of 'rdb' loop */

	}	// end of elementcheck loop
	
	}	/* end of 'limit' loop */
	//************************************************************************************************************/
		

	/*
	TK: if the current mass is larger than the limit the loop can be exited.
	Each element must point to the element which is in use and before.
	This is a static implementation which can be enhanced with a pointer chain to the lower element.
	Actually now its only allowed for CHNSOP-Fl-Cl-Br-Si !!! Brute-force <> elegance :-)
	*/
		} /*"H"*/
		
		} /*"D"*/
		
		if ((mass >= limit_lo) && (el[2].save == el[2].cnt-1)) break;
		} /* "C"*/
		
		} /*"13C"*/

		if ((mass >= limit_lo) && (el[0].save == el[0].cnt-1)) break;
		} /*"N"*/
		
		} /*"15N"*/

        if ((mass >= limit_lo) && (el[4].save == el[4].cnt-1)) break;
		} /*"O"*/
		
	    if ((mass >= limit_lo) && (el[6].save == el[6].cnt-1)) break;
		} /*"F"*/
		
		} /*"Na"*/
		
	    if ((mass >= limit_lo) && (el[7].save == el[7].cnt-1)) break;
		}  /*"Si"*/
		
		if ((mass >= limit_lo) && (el[9].save == el[9].cnt-1)) break;
		} /*"P"*/
		
		if ((mass >= limit_lo) && (el[10].save == el[10].cnt-1)) break;
		} /*"S"*/
		
		if ((mass >= limit_lo) && (el[11].save == el[11].cnt-1)) break;
		} /*"Cl"*/
		
		if ((mass >= limit_lo) && (el[12].save == el[12].cnt-1)) break;
		} /*"Br" ends*/
/* close that giant loop thing started above */

// denovofile.close(); //RW
//return 0; //RW
	

time(&finish);		// stop timer
elapsed_time = difftime(finish , start);	// calulate time difference



return hit;
}


/************************************************************************
* CLEAN:	"cleans" a buffer obtained by fgets() 			*
* Input: 	Pointer to text buffer					*
* Returns:	strlen of buffer.   					*
*************************************************************************/
int clean (char *buf)
{
int i;

for(i = 0;  i < strlen(buf);  i++)		/* search for CR/LF */
	{
	if(buf[i] == '\n' || buf[i] == '\r')
		{
		buf[i] = 0;			/* stop at CR or LF */
		break;
		}
	}
return (strlen(buf));
}


