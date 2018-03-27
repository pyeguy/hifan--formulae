/*

SMISOTOPE.C, based on ISOTOPE.C

modified by Dr. Robert Winkler, robert.winkler@ira.cinvestav.mx, robert.winkler@bioprocess.org, 2014
for integration in SpiderMass

main changes:

2013-10-21 Isotope tables revised, according to NIST
2014-02-21 Heavy isotopes added as D (2H), X (13C) and N (15N) 
--------------------------------------------------------------------


 I S O T O P E . C

 A program to calculate isotope patterns from a given formula.

--------------------------------------------------------------------
 Copyright (c) 1996...2005 Joerg Hau <joerg.hau(at)dplanet.ch>.

 mail: joerg.hau at dplanet.ch
 www:  http://www.mysunrise.ch/users/joerg.hau/

 This program is free software; you can redistribute it and/or
 modify it under the terms of version 2 of the GNU General Public
 License as published by the Free Software Foundation. See the
 file LICENSE for details.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
--------------------------------------------------------------------

 Acknowledgement: This program is based on a file named MASS.C,
 which one of my (former) colleagues found around 1994 "somewhere
 on Internet". He could not recall the source, and the file did
 neither carry any copyright nor was the author identified
 somehow.
 As the original file was already publicly available, I put
 this modified version under the GNU Public License.

 Note: Version code is reflected by date in ISO writing (20020125)

 Modification history:

 1996-xx-xx, put element and isotope tables IN the code (JHa).
 1998-03-23, extended comments and output (JHa).
 2002-01-25, some minor (stylistic) fixes (JHa)
 2002-01-26, fix error message when entering new element
 2005-06-17, added command line calculation (JHa)

To run the program interactively, enter formula when asked. Element symbols
must be correctly typed, with upper and lower case as usual. Elements can be
repeated in a formula, but brackets are not understood: CH3OH is ok, but
Ni(CO)4 is not; this has to be typed as NiC4O4.

If you run the program, the output is by defaut normalized to 100.00
percent for the bigest peak. If you launch it with option '-f' instead, the
output peaks will be printed as fractions of the total intensity.

You can redirect the output to a file or printer:

    isotope > filename
or
    isotope > printername.

You can specify a formula on the command line, but only the first formula
on the line is taken into account. Example:

    isotope C12H22O11

Calculation of a series of spectra can be done in batch mode. Create an
input file of all the formula, one per line, with a 'q' (quit) as the first
character of the last line. Then type

    isotope < infile > outfile.

(... but don't make any mistakes in the formulas, if you do it this way!).

 This is ANSI C and should compile with any C compiler; use
 something along the lines of "gcc -Wall -O3 -o isotope isotope.c".
 Optimize for speed!

*/

#define VERSION "20140221"	/* string! */
#define CUTOFF 1e-7

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

typedef struct {int m; float fr; }                  isotope;    /* mass, abundance */
typedef struct {char *sym; int niso;  isotope *p; } element;    /* symbol, no. of isotopes, ptr */
typedef struct {int atno; int count; }              atom;
typedef struct {int mass; float intens; }           peak;

element el[] =
 {
 { "H",  2, NULL },
 { "D", 1, NULL}, //RW Deuterium, 2H
 { "He", 2, NULL },
 { "Li", 2, NULL },
 { "Be", 1, NULL },
 { "B",  2, NULL },
 { "C",  2, NULL },
 { "X", 1, NULL }, //RW 13C
 { "N",  2, NULL },
 { "M", 1, NULL }, //RW 15N
 { "O",  3, NULL },
 { "F",  1, NULL },
 { "Ne", 3, NULL },
 { "Na", 1, NULL },
 { "Mg", 3, NULL },
 { "Al", 1, NULL },
 { "Si", 3, NULL },
 { "P",  1, NULL },
 { "S",  4, NULL },
 { "Cl", 2, NULL },
 { "Ar", 3, NULL },
 { "K",  3, NULL },
 { "Ca", 6, NULL },
 { "Sc", 1, NULL },
 { "Ti", 5, NULL },
 { "V",  2, NULL },
 { "Cr", 4, NULL },
 { "Mn", 1, NULL },
 { "Fe", 4, NULL },
 { "Co", 1, NULL },
 { "Ni", 5, NULL },
 { "Cu", 2, NULL },
 { "Zn", 5, NULL },
 { "Ga", 2, NULL },
 { "Ge", 5, NULL },
 { "As", 1, NULL },
 { "Se", 6, NULL },
 { "Br", 2, NULL },
 { "Kr", 6, NULL },
 { "Rb", 2, NULL },
 { "Sr", 4, NULL },
 { "Y",  1, NULL },
 { "Zr", 5, NULL },
 { "Nb", 1, NULL,  },
 { "Mo", 7, NULL },
 { "Tc", 1, NULL },
 { "Ru", 7, NULL },
 { "Rh", 1, NULL },
 { "Pd", 6, NULL },
 { "Ag", 2, NULL },
 { "Cd", 8, NULL },
 { "In", 2, NULL },
 { "Sn", 10, NULL },
 { "Sb", 2, NULL },
 { "Te", 8, NULL },
 { "I" , 1, NULL },
 { "Xe", 9, NULL },
 { "Cs", 1, NULL },
 { "Ba", 7, NULL },
 { "La", 2, NULL },
 { "Ce", 4, NULL },
 { "Pr", 1, NULL },
 { "Nd", 7, NULL },
 { "Pm", 1, NULL },
 { "Sm", 7, NULL },
 { "Eu", 2, NULL },
 { "Gd", 7, NULL },
 { "Tb", 1, NULL },
 { "Dy", 7, NULL },
 { "Ho", 1, NULL },
 { "Er", 6, NULL },
 { "Tm", 1, NULL },
 { "Yb", 7, NULL },
 { "Lu", 2, NULL },
 { "Hf", 6, NULL },
 { "Ta", 2, NULL },
 { "W" , 5, NULL },
 { "Re", 2, NULL },
 { "Os", 7, NULL },
 { "Ir", 2, NULL },
 { "Pt", 5, NULL },
 { "Au", 1, NULL },
 { "Hg", 7, NULL },
 { "Tl", 2, NULL },
 { "Pb", 4, NULL },
 { "Bi", 1, NULL },
 { "Po", 1, NULL },
 { "At", 1, NULL },
 { "Rn", 1, NULL },
 { "Fr", 1, NULL },
 { "Ra", 1, NULL },
 { "Ac", 1, NULL },
 { "Th", 1, NULL },
 { "Pa", 1, NULL },
 { "U" , 3, NULL },
 { "Np", 1, NULL },
 { "Pu", 1, NULL },

	};

isotope iso[] =
{
{  1,  .999885   /* H */ },
{  2,  .000115 },
{  2, 1.0 /* D */ },
{  3,  .00000134  /* He */ },
{  4,  .99999866   },
{  6,  .0759    /* Li */ },
{  7,  .9241 },
{  9,  1.00     /* Be */ },
{ 10,  .199     /* B */ },
{ 11,  .801 },
{ 12,  .9893   /* C */ },
{ 13,  .0107 },
{ 13, 1.0 /* 13C */ },
{ 14,  .99636   /* N */ },
{ 15,  .00364 },
{ 15, 1.0 /* 15N */},
{ 16,  .99757   /* O */ },
{ 17,  .00038 },
{ 18,  .00205 },
{ 19,  1.00     /* F */ },
{ 20,  .9048    /* Ne */ },
{ 21,  .0027 },
{ 22,  .0925 },
{ 23,  1.00     /* Na */ },
{ 24,  .7899    /* Mg */ },
{ 25,  .1000 },
{ 26,  .1101 },
{ 27,  1.00     /* Al */ },
{ 28,  .92223    /* Si */ },
{ 29,  .04685 },
{ 30,  .03092 },
{ 31,  1.00     /* P */ },
{ 32,  .9499     /* S */ },
{ 33,  .0075 },
{ 34,  .0425 },
{ 36,  .0001 },
{ 35,  .7576   /* Cl */ },
{ 37,  .2424 },
{ 36,  .003365   /* Ar */ },
{ 38,  .000632 },
{ 40,  .996003 },
{ 39,  .932581    /* K */ },
{ 40,  .000117 },
{ 41,  .067302 },
{ 40,  .96941    /* Ca */ },
{ 42,  .00647 },
{ 43,  .00135 },
{ 44,  .02086 },
{ 46,  .00004 },
{ 48,  .00187 },
{ 45,  1.00     /* Sc */ },
{ 46,  .0825    /* Ti */ },
{ 47,  .0744 },
{ 48,  .7372 },
{ 49,  .0541 },
{ 50,  .0518 },
{ 50,  .00250    /* V */ },
{ 51,  .99750 },
{ 50,  .04345    /* Cr */ },
{ 52,  .83789 },
{ 53,  .09501 },
{ 54,  .02365 },
{ 55,  1.00     /* Mn */ },
{ 54,  .05845    /* Fe */ },
{ 56,  .91754 },
{ 57,  .02119 },
{ 58 , .00282 },
{ 59,  1.00     /* Co */ },
{ 58,  .680769    /* Ni */ },
{ 60,  .262231 },
{ 61,  .011399 },
{ 62,  .036345 },
{ 64,  .009256 },
{ 63,  .6915,    /* Cu */ },
{ 65,  .3085 },
{ 64,  .48268    /* Zn */ },
{ 66,  .27975 },
{ 67,  .04102 },
{ 68,  .19024 },
{ 70,  .00631 },
{ 69,  .60108     /* Ga */ },
{ 71,  .39892 },
{ 70,  .2038    /*Ge */ },
{ 72,  .2731 },
{ 73,  .0776 },
{ 74,  .3672 },
{ 76,  .0783 },
{ 75,  1.00     /* As */ },
{ 74,  .0089    /* Se */ },
{ 76,  .0937 },
{ 77,  .0763 },
{ 78,  .2377 },
{ 80,  .4961 },
{ 82,  .0873 },
{ 79,  .5069    /* Br */ },
{ 81,  .4931 },
{ 78,  .00355    /* Kr */ },
{ 80,  .02286 },
{ 82,  .11593 },
{ 83,  .11500 },
{ 84,  .56987 },
{ 86,  .17279 },
{ 85,  .7217    /* Rb */ },
{ 87,  .2783  },
{ 84,  .0056    /* Sr */ },
{ 86,  .0986 },
{ 87,  .0700 },
{ 88,  .8258   },
{ 89,  1.00     /* Y */ },
{ 90,  .5145    /* Zr */ },
{ 91,  .1122 },
{ 92,  .1715 },
{ 94,  .1738 },
{ 96,  .0280 },
{ 93,  1.00     /* Nb */ },
{ 92,  .1477    /* Mo */ },
{ 94,  .0923 },
{ 95,  .1590 },
{ 96,  .1668 },
{ 97,  .0956 },
{ 98,  .2419 },
{100,  .0967 },
{ 98,  1.00     /* Tc */ },
{ 96,  .0554    /* Ru */ },
{ 98,  .0187 },
{ 99,  .1276 },
{100,  .1260 },
{101,  .1706 },
{102,  .3155 },
{104,  .1862 },
{103,  1.00     /* Rh */ },
{102,  .0102    /* Pd */ },
{104,  .1114 },
{105,  .2233 },
{106,  .2733 },
{108,  .2646 },
{109,  .1172 },
{107,  .51839    /* Ag */ },
{109,  .48161 },
{106,  .0125    /* Cd */ },
{108,  .0089 },
{110,  .1249 },
{111,  .1280 },
{112,  .2413 },
{113,  .1222 },
{114,  .2873 },
{116,  .0749 },
{113,  .0429    /* In */ },
{115,  .9571 },
{112,  .0097    /* Sn */ },
{114,  .0066 },
{115,  .0034 },
{116,  .1454 },
{117,  .0768 },
{118,  .2422 },
{119,  .0859 },
{120,  .3258 },
{122,  .0463               },
{124,  .0579 },
{121,  .5721    /* Sb */ },
{123,  .4279 },
{120,  .0009   /* Te */ },
{122,  .0255 },
{123,  .0089 },
{124,  .0474 },
{125,  .0707 },
{126,  .1884 },
{128,  .3174 },
{130,  .3408 },
{127,  1.00     /* I */ },
{124,  .000952    /* Xe */ },
{126,  .000890 },
{128,  .019102 },
{129,  .264006 },
{130,  .040710 },
{131,  .212324 },
{132,  .269086 },
{134,  .104357 },
{136,  .088573 },
{133,  1.00     /* Cs */ },
{130,  .00106   /* Ba */ },
{132,  .00101 },
{134,  .02417 },
{135,  .06592 },
{136,  .07854 },
{137,  .11232 },
{138,  .71698 },
{138,  .00090   /* La */ },
{139,  .99910 },
{136,  .00185   /* Ce */ },
{138,  .00251 },
{140,  .88450 },
{142,  .11114 },
{141,  1.00     /* Pr */ },
{142,  .272    /* Nd */ },
{143,  .122 },
{144,  .238 },
{145,  .083 },
{146,  .172 },
{148,  .057 },
{150,  .056 },
{145,  1.00     /* Pm */ },
{144,  .0307    /* Sm */ },
{147,  .1499 },
{148,  .1124 },
{149,  .1382 },
{150,  .0738 },
{152,  .2675 },
{154,  .2275 },
{151,  .4781    /* Eu */ },
{153,  .5219 },
{152,  .0020    /* Gd */ },
{154,  .0218 },
{155,  .1480 },
{156,  .2047 },
{157,  .1565 },
{158,  .2484 },
{160,  .2186 },
{159,  1.00     /* Tb */ },
{156,  .00056    /* Dy */ },
{158,  .00095 },
{160,  .02329 },
{161,  .18889 },
{162,  .25475 },
{163,  .24896 },
{164,  .28260 },
{165,  1.00     /* Ho */ },
{162,  .00139     /* Er */ },
{164,  .01601 },
{166,  .33503 },
{167,  .22869 },
{168,  .26978 },
{170,  .14910 },
{169,  1.00     /* Tm */ },
{168,  .0013    /* Yb */ },
{170,  .0304 },
{171,  .1428 },
{172,  .2183 },
{173,  .1613 },
{174,  .3183 },
{176,  .1276 },
{175,  .9741     /* Lu */ },
{176,  .0259 },
{174,  .0016    /* Hf */ },
{176,  .0526 },
{177,  .1860 },
{178,  .2728 },
{179,  .1362 },
{180,  .3508 },
{180,  .00012  /* Ta */ },
{181,  .99988 },
{180,  .0012    /* W */ },
{182,  .2650 },
{183,  .1431 },
{184,  .3064 },
{186,  .2843 },
{185,  .3740    /* Re */ },
{187,  .6260 },
{184,  .0002   /* Os */ },
{186,  .0159 },
{187,  .0196 },
{188,  .1324 },
{189,  .1615 },
{190,  .2626 },
{192,  .4078 },
{191,  .373     /* Ir */ },
{193,  .627 },
{190,  .00014,  /* Pt */ },
{192,  .00782 },
{194,  .32967 },
{195,  .33832 },
{196,  .25242 },
{198,  .07163 },
{197,  1.00     /* Au */ },
{196,  .0015   /* Hg */ },
{198,  .0997 },
{199,  .1687 },
{200,  .2310 },
{201,  .1318 },
{202,  .2986 },
{204,  .0687 },
{203,  .2952     /* Tl */ },
{205,  .7048 },
{204,  .014,    /* Pb */ },
{206,  .241 },
{207,  .221 },
{208,  .524 },
{209,  1.00     /* Bi */ },
{209,  1.00     /* Po */ },
{210,  1.00     /* At */ },
{222,  1.00     /* Rn */ },
{223,  1.00     /* Fr */ },
{226,  1.00     /* Ra */ },
{227,  1.00     /* Ac */ },
{231,  1.00     /* Pa */ },
{234,  .000054    /* U */ },
{235,  .007204 },
{238,  .992742 },
{237,  1.00      /* Np */ },
{244,  1.00      /* Pu */ },
};

#define ADDBASE 100        /* above the natural elements */
#define MAXADD 10          /* space for user-defined 'elements' */
#define MAXAT  50

element	addel[MAXADD];      	/* e.g. isotope enriched ones, etc. */
isotope	addiso[5 * MAXADD];
atom	atoms[MAXAT];    	/* up to MAXAT different atoms in one formula */
int 	natoms;                 /* number of atoms */
int 	eadd=0, iadd=0;        	/* added elements & isotopes */

int nel = sizeof(el) / sizeof(element);
int niso = sizeof(iso) / sizeof(isotope);

/* --- some variables needed for reading the cmd line --- */

char   *optarg;		/* global: pointer to argument of current option */
static int optind = 1;	/* global: index of which argument is next. Is used
                as a global variable for collection of further
                arguments (= not options) via argv pointers. */


/* --- FILE OUTPUT --- (rw) */

FILE *isotopefile;


/* --- Function prototypes --- */

void 	setpointers (void);
void 	addelement(void);
int 	squob(char *s);
void 	foutput(int a, int n);
isotope	imin(int atno);
isotope	imax(int atno);
int 	atno(char *s);
int 	formula(char *in);
int     getopt(int argc, char *argv[], char *optionS);

void setpointers (void)      /* set pointers in el entries to start of isotopes */
{                            /* for that element in the iso table */
int i;
isotope *p;

p = iso;
for (i = 0;  i < nel;  i++)
	{
	el[i].p = p;
	p += el[i].niso;
	}
}


int squob(char *buf)      /* squeeze out the blanks in s and return the length */
			  /* of the resultant string not counting terminating null*/
{
int i,j;
char c;

i=j=0;
while ( (c = buf[i++]) != '\0')
	if (c != ' ') buf[j++]=c;
buf[j]='\0';
return(j);
}


void foutput(int oz, int nr)
{
atoms[natoms].atno = oz;		/* Ordnungszahl */
atoms[natoms].count = nr;		/* No. of atoms of this element */
natoms++;
if (natoms >= MAXAT)
	{
	printf("Formula too long\n");
	exit(1);
	}
}


isotope imin(int atno)        /* isotope of lowest mass of element atno */
{
if (atno < ADDBASE)
      return (*(el[atno-1].p));  /* element 1 is in table entry 0, etc. */
else return(*(addel[atno-ADDBASE].p));
}


isotope imax(int atno)        /* isotope of highest mass of element atno */
{
if (atno < ADDBASE)
      return( *(el[atno - 1].p + el[atno -1].niso - 1));
else return(*(addel[atno-ADDBASE].p + addel[atno -ADDBASE].niso - 1));
}


int atno(char *str)  			/* return atom # or 0 if str not valid element symbol */
{
int i;

for (i = 0;  i < nel;  i++)             /* try natural elements first */
    if (0 == strcmp(el[i].sym, str))	/* if symbol is found in elem. table */
	return(i + 1);			/* return correct atom number */

for (i = 0;  i < eadd;  i++)		/* try 'user-def.' elements */
    if (0 == strcmp(addel[i].sym, str))	/* if symbol is found in 'user-def.'elem. table */
	return(i + ADDBASE);            /* return a high 'atomic #' on user-defined elements */

return 0;				/* this is 'else' */
}


int formula(char *in)
/* Determine if input is valid formula. Set composition in table *atoms.
   return 1 if OK, 0 if not call foutput() */
{
int 	a, n,      /* to handle (element, count) pairs. a is El.No., n is the number */
	state;     /* number of chars read */
char 	ch, buf[3];

state = a = n= 0;				/* here, no char was read */

while(0 != (ch = *(in++)))		/* while not end of string */
	{
	if (ch == '\n' || ch == '\r')	/* junk at end from fgets ? */
		continue;
	switch(state)
		{
		case 0:
		    if(isupper(ch))		/* uppercase letter ? */
			    {
			    buf[0] = ch;		/* take first char of element symbol */
			    buf[1] = 0;			/* ASCII zero to terminate string */
			    state = 1;			/* i.e. 1st char was read */
			    }
		    else
			goto error;
		    break;

		case 1:
		    if (isdigit(ch))			/* is it a number ? */
			    {
			    if(0 == (a = atno(buf)))	/* if NOT an element */
				goto error;
			    n = ch - '0';		/* just like 'atoi()' */
			    state = 2;			/* 2nd char was read */
			    }
		    else if(islower(ch))
			    {
			    buf[1] = ch;       	 	/* take char as 2nd letter */
			    buf[2] = 0;			/* terminate string */
			    if(0 == (a = atno(buf)))	/* see above */
				goto error;
			    state = 3;
			    }
		    else if(isupper(ch))
			    {
			    if(0 == (a = atno(buf)))	/* see above */
				goto error;
			    n = 1;                      /* 1st char read */
			    foutput(a, n);
			    buf[0] = ch;		/* take first char of element symbol */
			    buf[1] = 0;			/* ASCII zero to terminate string */
			    state = 1;			/* i.e. 1st char was read */
			    }
		    else
			goto error;
		    break;

		case 2:
		    if (isdigit(ch))                    /* is it a number */
			n = 10 * n + ch - '0';		/* YES -> get value (tens) */
		    else if(isupper(ch))
			    {
			    foutput(a, n);
			    buf[0] = ch;		/* take first char of element symbol */
			    buf[1] = 0;			/* ASCII zero to terminate string */
			    state = 1;			/* i.e. 1st char was read */
			    }
		    else
			goto error;
		    break;

		case 3:
		    if (isdigit(ch))      		/* is it a number */
			{				/* YES */
			n = ch - '0';			/* get value */
			state = 2;
			}
		    else if (isupper(ch))
			{
			if(0 == (a = atno(buf)))
				goto error;
			n = 1;
			foutput(a, n);
			buf[0] = ch;			/* take first char of element symbol */
			buf[1] = 0;			/* ASCII zero to terminate string */
			state = 1;			/* i.e. 1st char was read */
			}
		    else
			goto error;
		    break;

		}  /* end of case */
	}          /* end of while */

if (state == 1 || state == 3) n = 1;
if (state != 0)
	{
	if(0 == (a = atno(buf))) goto error;
	foutput(a, n);
	return(1);
	}

error:
	printf("Bad formula\n\n");
	return(0);
}



void addelement(void)            /* user-defined element */
{
int i, j, k, first;
float f, sumpc;
char buf1[81], buf2[81], *p;

while(1)
	{
	printf("Enter symbol for element: ");
	fgets(buf1, 80, stdin);
	for(i = 0;  i < strlen(buf1);  i++)
		if (buf1[i] == '\n' || buf1[i] == '\r')
			buf1[i] = ' ';
	squob(buf1);
	buf1[2] = 0;
	if (!isupper(buf1[0]) || (buf1[1] != 0 && (!islower(buf1[1]))))
		{
		printf("Bad symbol\n\n");
		continue;
		}
	if(0 == atno(buf1))     /* not already known symbol */
		{
		addel[eadd].sym = (char *)malloc(3);
		strcpy(addel[eadd].sym, buf1);
		break;
		}
	printf("This symbol is already in use.\n\n");
	}
printf("\nEnter mass - percent abundance pairs, one pair per line, with space between.\n");
printf("Enter an empty line to finish.\n");

startiso:
    addel[eadd].niso = 0;
    first = 1;
    sumpc = 0.0;
    while(1)
	{
	printf(">");
	fgets(buf1, 80, stdin);
	if (strlen(buf1) <= 1)
		{
		if(sumpc > 99.5 && sumpc < 100.5) break;
		else if(sumpc <= 99.5)
			{
			printf("Sum of isotope percentages is %.1f, press Y if ok, N to enter more:", sumpc);
			fgets(buf2, 80, stdin);
			if (buf2[0] == 'y' || buf2[0] == 'Y') break;
			else continue;
			}
		else
			{
			printf("Sum of isotope percentages is %.1f, press Y if ok, N to start over:", sumpc);
			fgets(buf2, 80, stdin);
			if (buf2[0] == 'y' || buf2[0] == 'Y') break;
			else goto startiso;
			}
		}
	j = strtol(buf1, &p, 10);

	if (j > 0 && j < 300)
		addiso[iadd].m = j;
	else if(j <= 0)
		{
		printf("Negative or zero mass. Re-enter line.");
		continue;
		}
	else
		{
		printf("Mass > 300. Type Y to confirm, or N to re-enter line. ");
		fgets(buf2, 80, stdin);
		k = buf2[0];
		if (k != 'y' && k != 'Y')
			continue;
		addiso[iadd].m = j;
		}

	f = atof(p);
	if (f > 0 && f <= 100)
		{
		addiso[iadd].fr = f / 100;
		sumpc += f;
		}
	else
		{
		printf("Impossible percentage. Re-enter line\n");
		continue;
		}

	addel[eadd].niso++;
	if(first)
		{
		addel[eadd].p = addiso + iadd;
		first = 0;
		}
	iadd++;
   }
eadd++;
}



int main (int argc, char *argv[])
{
int i, j, m, q, nold, nnew, oldmin, oldmax, newmin, newmax, ii, csvi, ix, fraction, ns, read_cmd, tmp;
register int k;
char buf[81], stars[71];
float fr, maxintens, sumintens;
isotope s;
element e;
peak *old, *new;
register peak *pp;

static char *id =
"isotope version %s. Copyright (C) by Joerg Hau 1996...2005, modified by Robert Winkler (RW) 2014 \n";

static char *msg =
"\nusage: isotope [-h] [-v] [-f] [formula]\n\nValid command line options are:\n"
"    -h       This Help screen.\n"
"    -v       Display version information.\n"
"    -f       Print fractional intensities (default: scaled to 100%).\n"
"    formula  Chemical formula, e.g. 'C12H11O11'. 2H=D, 13C=C, 15N=M.\n";

static char *disclaimer =
"\nThis program is free software; you can redistribute it and/or modify it under\n"
"the terms of version 2 of the GNU General Public License as published by the\n"
"Free Software Foundation.\n\n"
"This program is distributed in the hope that it will be useful, but WITHOUT ANY\n"
"WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n"
"PARTICULAR PURPOSE. See the GNU General Public License for details.\n";

fraction = 0;   /* normalize to 100 max, or print fractions on -f cmd switch */
read_cmd = 0;   /* != 0 if formula is read via cmd line */
setpointers();
nnew = 0;
new = NULL;

/* decode and read the command line */

while ((tmp = getopt(argc, argv, "hvf")) != EOF)
	switch (tmp)
		{
		case 'h':     	  		/* help me */
			printf (id, VERSION);
            printf (msg);
			printf(disclaimer, "\n");
			return 0;
		case 'v':     	  		/* version */
			printf (id, VERSION);
			return 0;
		case 'f':    			/* print fractional intensities */
			fraction = 1;
			continue;
		case '~':    	  	/* invalid arg */
		default:
			printf ("'%s -h' for help.\n", argv[0]);
			return 1;
		}

if (argv[optind] != NULL)	 /* remaining parameter on cmd line? */
    {
	strcpy(buf, argv[optind]);     /* read it */
    read_cmd = 1;           /* set flag */
    }

while(1)
	{
    if (!read_cmd)          /* if NOT read via cmd line, use interactive mode */
        {
        printf("Enter a formula (no brackets). Q to quit, E to define an extra element.\n:");
	    fgets(buf, 80, stdin);				/* read line from stdin */
	    buf[80] = 0;					/* terminate string */
        }

	squob(buf);					/* clean up */
	for(k = 0;  k < strlen(buf);  k++)		/* search for CR/LF */
		{
		if(buf[k] == '\n' || buf[k] == '\r')
			{
			buf[k] = 0;			/* stop at CR or LF */
			break;
			}
		}

	if (buf[0] == 'q' || buf[0] == 'Q')
		exit(0);

	if ((buf[0] == 'e' || buf[0] == 'E') && strlen(buf) == 1)	/* add -e-lement */
		{
		addelement();
		continue;
		}

	natoms = 0;					/* init. */
	k = formula(buf);
	if (k == 0)                 /* problem ? */
        {
        if (read_cmd)           /* if formula was read via cmd line, */
            exit (1);           /* quit here */
        continue;
        }

	old = (peak *)malloc(sizeof(peak));		/* init. */
	old->mass = 0;
	old->intens = 1;
	nold = 1;
	oldmin = oldmax = 0;

	for (i = 0;  i < natoms;  i++)				/* for all elements */
		{
		for(j = 0;  j < atoms[i].count;  j++)		/* for all atoms of an element */
			{
			s = imin(atoms[i].atno);
			newmin = oldmin + s.m;			/* min. mass */
			s = imax(atoms[i].atno);
			newmax = oldmax + s.m;			/* max. mass */
			nnew = newmax - newmin + 1;		/* number */
			new = (peak *)malloc(nnew * sizeof(peak));
			if (new == NULL)
				{
				printf("\nOut of Memory!\n");
				exit(1);
				}
			for (k = 0;  k < nnew;  k++)
				new[k].intens = 0;		/* init. */

			if (atoms[i].atno < ADDBASE)
				e = el[atoms[i].atno - 1];
			else
				e = addel[atoms[i].atno - ADDBASE];
				
				
			for (k = 0;  k < e.niso;  k++)          /* for all isotopes */
				{
				m = (k + e.p)->m;		/* mass */
				fr = (k + e.p)->fr;		/* inty */
				for (q = 0;  q < nold;  q++)
					{
					ix = m + old[q].mass - newmin;
					new[ix].mass = m + old[q].mass;		/* shift mass */
					new[ix].intens += fr * old[q].intens; 	/* add inty */
					}
				}	/* end of 'k' loop (isotopes) */

		       /* normalize to maximum intensity of 1.0 */
			maxintens = 0;
			for (ii = 0;  ii < nnew;  ii++)
				if (new[ii].intens > maxintens)		/* find max. value */
					maxintens = new[ii].intens;
			for(ii = 0;  ii < nnew;  ii++)
					new[ii].intens /= maxintens;

			/* throw away very small peaks */
			for (ii = 0;  ii < nnew;  ii++)
				{
				if(new[ii].intens < CUTOFF)
					{
					for (k = ii, pp = new + ii;  k < nnew-1;  k++, pp++)
						{
						pp->mass = (pp+1)->mass;
						pp->intens = (pp+1)->intens;   /* inner loop...*/
						}                            /* avoid structure copy */
					nnew--;
					ii--;
					}
				}
			free(old);
			old = new;
			nold = nnew;
			oldmin = newmin;
			oldmax = newmax;
			}	/* end of 'j' loop (atoms) */
	       }		/* end of 'i' loop (elements) */

	maxintens = sumintens = 0;
	for (ii=0;  ii<nnew;  ii++)		/* find max. */
		{
		sumintens += new[ii].intens;
		if (new[ii].intens > maxintens)
			maxintens = new[ii].intens;
		}
	for (ii = 0;  ii < nnew;  ii++)		/* calculate fraction/percent */
		{
		if (fraction)
			new[ii].intens /= sumintens;
		else
			new[ii].intens *= 100;  /* they are already normalized to max=1 */
		}
	if (fraction)
		maxintens /= sumintens;
	else
		maxintens *= 100;

	
	for(ii = 0;  ii < nnew;  ii++)
		{
		ns = .5 + 60.0 * new[ii].intens / maxintens;	/* no. of stars */
		for (k = 0;  k < ns;  k++)
			stars[k] = '*';
		stars[ns] = 0;
		printf(fraction? "%5d%8.4f  |%s\n" : "%5d%8.2f  |%s\n",
				 new[ii].mass, new[ii].intens, stars);
		}
	printf("\n");
	
	isotopefile = fopen("isotopes.csv", "w"); //RW open outputfile for writing
	fprintf(isotopefile,"TM0;TM1;TM2;TM3; \n");
	
	for(ii = 0;  ii < 4;  ii++) //RW
		{
		fprintf(isotopefile,"%.10f;",new[ii].intens); //RW write intensities to output file
		}
	fclose(isotopefile); //RW close the output file
	
	free(new);
    if (read_cmd)          /* if formula was read via cmd line, quit here */
        exit (0);
	}		/* end of 'while (1)...' */
}


/***************************************************************************
* GETOPT: Command line parser, system V style.
*
*  This routine is widely (and wildly) adapted from code that was
*  made available by Borland International Inc.
*
*  Standard option syntax is:
*
*    option ::= SW [optLetter]* [argLetter space* argument]
*
*  where
*    - SW is '-'
*    - there is no space before any optLetter or argLetter.
*    - opt/arg letters are alphabetic, not punctuation characters.
*    - optLetters, if present, must be matched in optionS.
*    - argLetters, if present, are found in optionS followed by ':'.
*    - argument is any white-space delimited string.  Note that it
*      can include the SW character.
*    - upper and lower case letters are distinct.
*
*  There may be multiple option clusters on a command line, each
*  beginning with a SW, but all must appear before any non-option
*  arguments (arguments not introduced by SW).  Opt/arg letters may
*  be repeated: it is up to the caller to decide if that is an error.
*
*  The character SW appearing alone as the last argument is an error.
*  The lead-in sequence SWSW ("--") causes itself and all the rest
*  of the line to be ignored (allowing non-options which begin
*  with the switch char).
*
*  The string *optionS allows valid opt/arg letters to be recognized.
*  argLetters are followed with ':'.  Getopt () returns the value of
*  the option character found, or EOF if no more options are in the
*  command line. If option is an argLetter then the global optarg is
*  set to point to the argument string (having skipped any white-space).
*
*  The global optind is initially 1 and is always left as the index
*  of the next argument of argv[] which getopt has not taken.  Note
*  that if "--" or "//" are used then optind is stepped to the next
*  argument before getopt() returns EOF.
*
*  If an error occurs, that is an SW char precedes an unknown letter,
*  then getopt() will return a '~' character and normally prints an
*  error message via perror().  If the global variable opterr is set
*  to false (zero) before calling getopt() then the error message is
*  not printed.
*
*  For example, if
*
*    *optionS == "A:F:PuU:wXZ:"
*
*  then 'P', 'u', 'w', and 'X' are option letters and 'A', 'F',
*  'U', 'Z' are followed by arguments. A valid command line may be:
*
*    aCommand  -uPFPi -X -A L someFile
*
*  where:
*    - 'u' and 'P' will be returned as isolated option letters.
*    - 'F' will return with "Pi" as its argument string.
*    - 'X' is an isolated option.
*    - 'A' will return with "L" as its argument.
*    - "someFile" is not an option, and terminates getOpt.  The
*      caller may collect remaining arguments using argv pointers.
***************************************************************************/
int getopt(int argc, char *argv[], char *optionS)
{
static char *letP	= NULL;		/* remember next option char's location */
static char SW		= '-';		/* switch character */

int opterr = 1;				/* allow error message	*/
unsigned char ch;
char *optP;

if (argc > optind)
	{
	if (letP == NULL)
		{
		if ((letP = argv[optind]) == NULL || *(letP++) != SW)
			goto gopEOF;

		if (*letP == SW)
			{
			optind++;
			goto gopEOF;
			}
		}
	if (0 == (ch = *(letP++)))
		{
		optind++;
		goto gopEOF;
		}
	if (':' == ch  ||  (optP = strchr(optionS, ch)) == NULL)
		goto gopError;
	if (':' == *(++optP))
		{
		optind++;
		if (0 == *letP)
			{
			if (argc <= optind)
				goto  gopError;
			letP = argv[optind++];
			}
		optarg = letP;
		letP = NULL;
	}
	else
	{
	if (0 == *letP)
		{
		optind++;
		letP = NULL;
		}
	optarg = NULL;
	}
	return ch;
}

gopEOF:
	optarg = letP = NULL;
	return EOF;

gopError:
	optarg = NULL;
	errno  = EINVAL;
	if (opterr)
		perror ("Command line option");
	return ('~');
}
