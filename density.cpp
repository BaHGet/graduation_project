#include "bin/files.h"
#include "bin/plot.h"
#include "bin/algebra.h"
int flag_opposite = 1;
void convert(char *namei, int ls,
             int colx, int coly, int colI, int *colII, int I_flag,
             numero x0, numero y0, numero x1, numero y1,
             numero Imin, numero Imax, int logf,
             int dimx, int dimy, int limits,
             numero *xx, numero *yy, int *cc, int ncurve, int color,
             int four, char *dtext)
{
	char nameo[200];
	strcpy(nameo, namei);
	int l = strlen(namei);
	if (l < 5) strcat(nameo, ".gif");
	else                       \
		if (nameo[l-4] == '.')
		{
			nameo[l-3] = 'g';
			nameo[l-2] = 'i';
			nameo[l-1] = 'f';
		}
	int nc = number_of_columns(ls, namei);
	int nr = number_of_rows(ls, namei);
	numero *intensity, *intensityF;
	intensity = new numero [nr];
	if (nc < 3) on_error("density", "input file must contain at least 3 columns");
	if (colx < 0 || colx >= nc || coly < 0 || coly >= nc)
		on_error("density", "colx/coly must be within the range 1 -", nc, "");
	FILE *fin;
	fin = fopen(namei, "r");
	skip(fin, ls);
	numero x, y, xp, val, vv, xx0, xx1, yy0, yy1, IImin = infinite, IImax = -infinite;
	numero valr, vali, q, II;
	int i, j, m, nx = -1, ny;
	for (i = 0; i < nr; i++)
	{
		val = 0;
		for (j = m = 0; j < nc; j++)
		{
			vv = alread_numero(fin);
			if (j == coly && colI > nc) val = vv; // this is to plot the color scale
			if (j == colx) x = vv;
			else
				if (j == coly) y = vv;
				else
					if (I_flag)  val += vv * vv;
					else
						if (j == colI && I_flag == 0) val = flag_opposite * vv;
						else
							if (colI < 0) if (j == colII[m])
								{
									val += vv * vv;
									m++;
								}
		}
		if (i == 0)
		{
			xx0 = xp = x;
			yy0 = y;
		}
		if (nx < 0 && xp != x)
		{
			ny = i;
			nx = nr / ny;
		}
		if (val < IImin) IImin = val;
		if (val > IImax) IImax = val;
		if (logf)
		{
			if (val < 0) val = -infinite;
			else val = log(val);
		}
		intensity[i] = val;
	}
	xx1 = x;
	yy1 = y;
	fclose(fin);
	if (four)
	{
		intensityF = new numero [nr];
		IImin = infinite;
		IImax = -infinite;
		if (four == 2)
		{
			for (i = 0; i < nx; i++)
			{
				for (m = 0; m < ny; m++)
				{
					q = y0 + m * (y1 - y0) / (ny - 1);
					valr = vali = 0;
					for (j = 0; j < ny; j++)
					{
						y = yy0 + j * (yy1 - yy0) / (ny - 1);
						II = intensity[i*ny+j];
						valr += cos(q * y) * II;
						vali += sin(q * y) * II;
					}
					intensityF[i*ny+m] = val = sqrt(valr * valr + vali * vali);
					if (val < IImin) IImin = val;
					if (val > IImax) IImax = val;
				}
			}
		}
		else
		{
			for (j = 0; j < ny; j++)
			{
				for (m = 0; m < ny; m++)
				{
					q = x0 + m * (x1 - x0) / (nx - 1);
					valr = vali = 0;
					for (i = 0; i < nx; i++)
					{
						x = xx0 + i * (xx1 - xx0) / (nx - 1);
						II = intensity[i*ny+j];
						valr += cos(q * x) * II;
						vali += sin(q * x) * II;
					}
					intensityF[m*ny+j] = val = sqrt(valr * valr + vali * vali);
					if (val < IImin) IImin = val;
					if (val > IImax) IImax = val;
				}
			}
		}
		delete [] intensity;
		intensity = intensityF;
	}
	if (four != 1)
	{
		if (x0 < xx0) x0 = xx0;
		if (x1 > xx1) x1 = xx1;
	}
	else
	{
		xx0 = x0;
		xx1 = x1;
	}
	if (four != 2)
	{
		if (y0 < yy0) y0 = yy0;
		if (y1 > yy1) y1 = yy1;
	}
	else
	{
		yy0 = y0;
		yy1 = y1;
	}
	if (Imax > IImax) Imax = IImax;
	if (Imin < IImin) Imin = IImin;
	printf("%s -> %s, skip=%d, colx=%d coly=%d colI=%d log=%d\n",
	       namei, nameo, ls, colx + 1, coly + 1, colI + 1, logf);
	printf("   nx=%d  ny=%d\n", nx, ny);
	printf("   physical data: %g <x< %g,", xx0, xx1);
	printf("  %g <y< %g,", yy0, yy1);
	printf("  %g <I< %g\n", IImin, IImax);
	printf("   plotted data:  %g <x< %g,", x0, x1);
	printf("  %g <y< %g,", y0, y1);
	printf("  %g <I< %g\n", Imin, Imax);
	printf("   dimx=%d  dimy=%d\n\n", dimx, dimy);
	if (logf)
	{
		if (Imin > 0) Imin = log(Imin);
		Imax = log(Imax);
	}
	plot_density(nameo, intensity,
	             x0, x1, y0, y1, xx0, xx1, yy0, yy1, nx, ny, dimx, dimy,
	             Imin, Imax, limits, xx, yy, cc, ncurve, 1, color, dtext);
	delete [] intensity;
}
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("This program yields a contour plot of a function I=f(x,y).\n");
		printf("The format is\n\n");
		printf("   density input [- skip=0 colx=1 coly=2 colI=3 ...\n");
		printf("   ... [x0=a y0=a x1=a y1=a ...\n");
		printf("   ... [Imin=a Imax=a log=0 [dimx=500 dimy=500 ...\n");
		printf("   ... density-color=0\n");
		printf("   ... [limits=0 [curve [cls=0 ix=1 iy=2 color=0 [text]]]]]]]]\n\n");
		printf("[Be aware that the '-' sign is separated by spaces.]\n");
		printf("Type 'density h' for further explanation of input parameters.\n");
		printf("(c) 2007 F. Javier Garcia de Abajo, Instituto de Optica, CSIC, Spain.\n");
		return 0;
	}
	if (!strcmp(argv[1], "h") || !strcmp(argv[1], "-"))
	{
		printf("input: input file name in ascii format;\n");
		printf("       wild chars are allowed in this parameter, in which\n");
		printf("       case the output-file parameter is ignored and the\n");
		printf("       output-file names are constructed from the input\n");
		printf("       files by adding a .gif extension.\n");
		printf("       The output file in gif format, and its name\n");
		printf("       will be taken as the input-file parameter with\n");
		printf("       a .gif extension added.\n");
		printf("skip: number or rows skiped from input file (e.g.,comments).\n");
		printf("colx,coly,colI: file columns taken as x, y, and intensity.\n");
		printf("    if colI<0 then -colI is the number of columns to be used\n");
		printf("    to plot col1^2+col2^2+..., and col1, col2, ... are given\n");
		printf("    immediatly after colI;\n");
		printf("    if colI>1000, then -I is plotted for column colI-1000\n");
		printf("(x0,y0), (x1,y1): corners of plotted rectangular area in x-y.\n");
		printf("Imin, Imax: minimum and maximum values of plotted intensity.\n");
		printf("log: log (linear) scale if set to 1 (0).\n");
		printf("dimx, dimy: number of pixels in x and y directions.\n");
		printf("density-color:  0 for color plot\n");
		printf("                1 for black&white plot (min&max)\n");
		printf("                2 for black&red plot\n");
		printf("                3 for black&blue plot\n");
		printf("                4 for white&black plot (min&max)\n");
		printf("                5 for white&red plot\n");
		printf("                6 for white&blue plot\n");
		printf("                7 for blue-black-red-white\n");
		printf("                8 for blue-black-red\n");
		printf("                9 for blue-white-red\n");
		printf("               10 for yellow-black-red\n");
		printf("               11 for yellow-white-red\n");
		printf("               12 for green-white-red\n");
		printf("               13 for green-black-red\n");
		printf("               14 for black&yellow\n");
		printf("               15 for black-red-yellow\n");
		printf("               16 for black-yellow-white\n");
		printf("               20 for black-red-white\n");
		printf("               21 for black-green-white\n");
		printf("               22 for black-blue-white\n");
		printf("               23 for black-pink-white\n");
		printf("               24 for black-cyan-white\n");
		printf("               25 for white-pink\n");
		printf("               26 for white-gray\n");
		printf("               50 for black-blue-red-white\n");
		printf("               60 for red-black in (-1,0); white-blue in(0,1); yellow outside\n");
		printf("limits: set to 1 to plot in white (black) values above\n");
		printf("       (below) the range specified by (Imin,Imax);\n");
		printf("       set to 2 to plot in black (white) values above\n");
		printf("       (below) the range specified by (Imin,Imax);\n");
		printf("       set to 3 to plot in white values above range;\n");
		printf("       set to 4 to plot in black values below range;\n");
		printf("       set to 5 to plot in black values above range;\n");
		printf("       set to 6 to plot in white values below range;\n");
		printf("       otherwise, these values ar plotted in white (black).\n");
		printf("curve: input file from where the first cls rows are ignored\n");
		printf("       and the columns ix and iy of which define a curve to\n");
		printf("       be superimposed on the density plot with color\n");
		printf("       0 (black), 1 (red), 2 (blue), or 3 (white).\n");
		printf("       Use curve='none' or '0' to ignore the curve.\n");
		printf("text:  text to be shown with the plot\n");
		printf("       (set to 000 to show the name of the input data file.\n\n");
		printf("Default values are shown above. 'a' means that the maximum\n");
		printf("possible range is selected. 'a' can be used as parameter.\n");
		printf("The number of columns in the file is detected automatically.\n");
		printf("The input file must have been generated by a\n");
		printf("'loop-x of loop-y' structure.\n");
		printf("\n");
		printf("Miscellaneous:\n");
		printf("\n");
		printf("Use 'density input -Fx ...' rather than 'density input - ...'\n");
		printf("to produce the image of the Fourier transform of the input\n");
		printf("along the x direction, in which case x0 and x1 are\n");
		printf("interpreted as the range of momenta qx, with the same.\n");
		printf("number of qx's as the number of x'x\n");
		printf("(x0 and x1 are mandatory arguments in that case).\n");
		printf("\n");
		printf("Use 'density input -Fy ...' rather than 'density input - ...'\n");
		printf("to produce the image of the Fourier transform of the input\n");
		printf("along the y direction, the the same number of qy's as y'x\n");
		printf("(y0 and y1 are mandatory args., otherwise expect meaningless results.\n");
		return 0;
	}
	int ls = 0, colx = 0, coly = 1, colI = 2, *colII, I_flag = 0, logf = 0, dimx = 500, dimy = 500;
	numero x0 = -infinite, y0 = -infinite, x1 = infinite, y1 = infinite;
	numero Imin = -2 * infinite, Imax = 2 * infinite, val;
	int limits = 0, i, j, l = 1, l0 = 1, cls = 0, ix = 0, iy = 1, cvcol = 0, color = 0, four = 0;
	char *dtext = NULL, *ddtext;
	while (argc > l + 1 && strcmp(argv[l], "-")
	        && strcmp(argv[l], "-Fx") && strcmp(argv[l], "-Fy")) l++;
	if (!strcmp(argv[l], "-"))   four = 0;
	else
		if (!strcmp(argv[l], "-Fx")) four = 1;
		else
			if (!strcmp(argv[l], "-Fy")) four = 2;
	if (!strcmp(argv[l], "-") || !strcmp(argv[l], "-Fx")
	        || !strcmp(argv[l], "-Fy"))
	{
		l--;
		l0 = l;
		if (argc > l + 2) ls =  atoi(argv[ l+2]);
		if (argc > l + 3) colx = atoi(argv[ l+3]) - 1;
		if (argc > l + 4) coly = atoi(argv[ l+4]) - 1;
		if (argc > l + 5) if (strcmp(argv[ l+5], "a")) colI = atoi(argv[l+5]) - 1;
			else                       I_flag = 1;
		if (colI < 0)
		{
			colI++;
			colII = new int [-colI+1];
			for (i = 0; i < -colI; i++)
			{
				l++;
				colII[i] = atoi(argv[l+5]) - 1;
			}
		}
		else
			if (colI > 1000)
			{
				flag_opposite = -1;
				colI = colI - 1000;
			}
		if (argc > l + 6) if (strcmp(argv[ l+6], "a")) x0 =  alread_numero(argv[ l+6]);
		if (argc > l + 7) if (strcmp(argv[ l+7], "a")) y0 =  alread_numero(argv[ l+7]);
		if (argc > l + 8) if (strcmp(argv[ l+8], "a")) x1 =  alread_numero(argv[ l+8]);
		if (argc > l + 9) if (strcmp(argv[ l+9], "a")) y1 =  alread_numero(argv[ l+9]);
		if (argc > l + 10) if (strcmp(argv[l+10], "a")) Imin = alread_numero(argv[l+10]);
		if (argc > l + 11) if (strcmp(argv[l+11], "a")) Imax = alread_numero(argv[l+11]);
		if (argc > l + 12) logf = atoi(argv[l+12]);
		if (argc > l + 13) dimx = atoi(argv[l+13]);
		if (argc > l + 14) dimy = atoi(argv[l+14]);
		if (argc > l + 15) color = atoi(argv[l+15]);
		if (argc > l + 16) limits = atoi(argv[l+16]);
		if (argc > l + 18) cls =   atoi(argv[l+18]);
		if (argc > l + 19) ix =    atoi(argv[l+19]) - 1;
		if (argc > l + 20) iy =    atoi(argv[l+20]) - 1;
		if (argc > l + 21) cvcol = atoi(argv[l+21]);
		if (argc > l + 22) dtext = argv[l+22];
	}
//  if(cvcol==0) cvcol=0; else if(cvcol==1) cvcol=254; else
//  if(cvcol==2) cvcol=1; else if(cvcol==3) cvcol=255;
	int ncur = 0, nco, *ccur = NULL;
	numero *xcur = NULL, *ycur = NULL;
	FILE *fin;
	if (argc > l + 17 && strcmp(argv[l+17], "none") && strcmp(argv[l+17], "0"))
	{
		nco = number_of_columns(cls, argv[l+17]);
		if (ix < 0 || ix >= nco || iy < 0 || iy >= nco)
			on_error("density", "ix/iy must be within the range 1 -", nco, "");
		ncur = number_of_rows(cls, argv[l+17]);
		fin = fopen(argv[l+17], "r");
		skip(fin, cls);
		ccur = new int [ncur];
		xcur = new numero [ncur];
		ycur = new numero [ncur];
		for (i = 0; i < ncur; i++)
		{
			ccur[i] = cvcol;
			for (j = 0; j < nco; j++)
			{
				val = alread_numero(fin);
				if (j == ix) xcur[i] = val;
				if (j == iy) ycur[i] = val;
			}
		}
		fclose(fin);
	}
	else ncur = 0;
	for (i = 1; i <= l; i++)
	{
		ddtext = dtext;
		if (dtext != NULL) if (!strcmp(dtext, "000")) ddtext = argv[i];
		convert(argv[i], ls, colx, coly, colI, colII, I_flag,
		        x0, y0, x1, y1, Imin, Imax, logf, dimx, dimy, limits,
		        xcur, ycur, ccur, ncur, color, four, ddtext);
	}
	return 0;
}

