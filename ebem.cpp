#define BEM_HELP                                                             "\
**************************************************************************  \n\
ebem.cpp                                             19-xi-04 --> 25-iv-06  \n\
**************************************************************************  \n\
                                                                            \n\
SOLUTION OF SCHRODINGER\'S EQUATION IN 2D FOR TWO DIFFERENT METHODS OF      \n\
SOLUTION:                                                                   \n\
                                                                            \n\
      - ELECTRON BOUNDARY ELEMENT METHOD (EBEM, WHICH GIVES THE NAME        \n\
        TO THIS CODE)                                                       \n\
      - ELECTRON PLANE WAVE APPROCH (EPW FOR PERIODIC SYSTEMS v             \n\
                                                                            \n\
(c)  Dr. F. Javier Garcia de Abajo,                   2004-2006             \n\
     Instituto de Optica, CSIC, Madrid, Spain, and                          \n\
     Donostia International Physics Center (DIPC), San Sebastian, Spain.    \n\
     E-mail: jga@sw.ehu.es                                                  \n\
                                                                            \n\
I - GENERAL DESCRIPTION                                                     \n\
                                                                            \n\
Solution of Schrodinger\'s equation in 2D including customizable barriers   \n\
and using the electron boundary element method (EBEM) for finite systems    \n\
or electron plane wave (EPW) expansions for periodic systems.               \n\
                                                                            \n\
This program is prepared to accept commands from an input file to define    \n\
the boundary geometry and instructs it to calculate the local density of    \n\
states (LDOS), photoemission (PE) intensities, and wave functions. This     \n\
program is well suited to deal with surface states in surfaces like         \n\
Cu(111), Au(111), etc. decorated with steps, molecules, adatoms, etc.       \n\
                                                                            \n\
The program uses atomic units (a.u., e=m=hbar=1) both internally and        \n\
on input/output (distances, energies, etc.). Angles are in degrees          \n\
on input and in radians inside the program. Conversion factors:             \n\
                                                                            \n\
      Energy   -> 1 a.u. = 27.2113834 eV = 4.35974381e-18 J.                \n\
      Length   -> 1 a.u. = 0.529177208 Angstroms.                           \n\
      Time     -> 1 a.u. = 2.4188843e-17 s.                                 \n\
      Mass     -> 1 a.u. = 9.1093819e-31 kg = electron mass.                \n\
      Velocity -> 1 a.u. = 2.18769125e6 m/s.                                \n\
                                                                            \n\
Additionally, the program allows one to calculate band structures for       \n\
periodic systems. This is done by using the commands described in Sec.      \n\
II.4a below, or also via an inline command as described in Sec. III.        \n\
                                                                            \n\
The EBEM code can be invoked in two different modes:                        \n\
                                                                            \n\
      1) 'ebem input-file': execution based upon an input file.             \n\
      2) 'ebem epw <parameters>': inline calculation of band structures.    \n\
                                                                            \n\
Mode 1 is the most general one and requires an input file. Mode 2 can       \n\
yield for band structures in simple periodic systems using plane waves,     \n\
as described in Sec. III below. Mode 1 includes all features of mode 2.     \n\
                                                                            \n\
II - INPUT FILE                                                             \n\
                                                                            \n\
The input file consists of a list of commands that are executed in          \n\
consecutive order. It has to have the following structure:                  \n\
                                                                            \n\
      1.- Define boundary geometry.                                         \n\
      2.- Define regions properties (potential), boundary properties,       \n\
          effective electron mass, etc.                                     \n\
      3.- Run calculations of LDOS, PE, band structures, or wave functions, \n\
          as many as needed, in consecutive order.                          \n\
      4.- Exit the program.                                                 \n\
                                                                            \n\
The 2D x-y plane is divided into regions by the boundaries defined in part  \n\
1 of the input file. Each region has a numerical label (0, 1, ...) and a    \n\
value of the potential as defined in part 2 of the input file (0 by         \n\
default). Each boundary must separe different regions, and even if they     \n\
have identical potential, their mumerical labels must be different.         \n\
                                                                            \n\
The maximum number of regions, segments (i.e., straight-line boundaries     \n\
or arcs, as defined below), and parametrization points is controlled        \n\
by variables defined inside the code by the following directives:           \n\
                                                                            \n\
      #define n_regions n1   // (regions from 0 to n1-1)                    \n\
      #define nl_max    n2   // (up to n2 boundary segments)                \n\
                                                                            \n\
and these values (n1,n2) can be changed by the user (just look              \n\
for #define ... inside the source code and change the values as needed,     \n\
but keep in mind that higher values will involve more RAM memory, and       \n\
remember to re-compile the code before using it again).                     \n\
                                                                            \n\
II.1 DESCRIPION OF COMMANDS USED TO DEFINE THE BOUNDARY GEOMETRY            \n\
                                                                            \n\
 -->  add-polygon-boundary                                                  \n\
             N  mu1 mu2  x0 y0  x1 y1 n1  x2 y2 n2 ... xN yN nN             \n\
                                                                            \n\
      This adds N straight-line segments representing potential barriers    \n\
      that separate different regions. The first segment is described by n1 \n\
      discretization points and goes from point (x0,y0) to point (x1,y1).   \n\
      The second segment, with n2 points, goes from (x1,y1) to (x2,y2).     \n\
      And so on up to the Nth segment, from (xN-1,yN-1) to (xN,yN) with nN  \n\
      points. The relative position of regions mu1 and mu2 is shown in the  \n\
      following diagrams                                                    \n\
                                                                            \n\
                    j+1             j                                       \n\
                    /              /                                        \n\
               mu1 /          mu2 /                                         \n\
                  / mu2          / mu1                                      \n\
                 /              /                                           \n\
                j             j+1                                           \n\
                                                                            \n\
      for the segment going from point (xj,yj) to point (xj+1,yj+1).        \n\
                                                                            \n\
 -->  add-arc-boundary  mu1 mu2  x0 y0  a b  alpha  th0 th1  n              \n\
                                                                            \n\
      This adds an arc of an ellipse described by n discretization points.  \n\
      The ellipse is centered at (x0,y0) and has a principal axis of radius \n\
      a forming an angle alpha degrees with respect to the x direction, and \n\
      the other principal axis of radius b forms an angle alpha degrees     \n\
      with respect to the y direction. The region inside/outside the        \n\
      ellipse is mu1/mu2:                                                   \n\
                                                                            \n\
                     |                                                      \n\
                   b |   mu1     th1                                        \n\
                     |            |   mu2                                   \n\
                     *------------|                                         \n\
              (x0,y0)     a       |                                         \n\
                                 /                                          \n\
                            th0 /                                           \n\
                                                                            \n\
      The ellipse spans the angular region (th0,th1) in degrees, with the   \n\
      angle measured counter-clockwise with respect to principal axis a.    \n\
                                                                            \n\
 -->  duplicate-geometry  dx dy                                             \n\
                                                                            \n\
      Create a copy of the boundaries that have been defined in the         \n\
      input file up to this command, but displaced by a vector (dx,dy)      \n\
      with respect to the original one.                                     \n\
                                                                            \n\
 -->  rotate-geometry xr yr phi                                             \n\
                                                                            \n\
      Rotate the geometry so far defined by an angle phi degrees (>0 for    \n\
      counter-clockwise rotation) with respect to the point (xr,yr).        \n\
      (This does not rotate the unit cell vectors in periodic systems.)     \n\
                                                                            \n\
 -->  translate-geometry Dx Dy                                              \n\
                                                                            \n\
      Translate geometry by a vector (Dx,Dy).                               \n\
                                                                            \n\
 -->  repeat-periodic-geometry a1 a2 alpha beta Rmax                        \n\
                                                                            \n\
      Use the geometry so far defined as a basis to construct a periodic 2D \n\
      lattice consisting of copies of that geometry. The 2D lattice has     \n\
      vectors of lengths a1 and a2. The first vector forms an angle alpha   \n\
      degrees with respect to the x axis. The second vector forms an angle  \n\
      beta degrees with respect to the first vector. All lattice sites up   \n\
      to a distance Rmax from (0,0) will be included.                       \n\
                                                                            \n\
 -->  repeat-periodic-geometry-1D a phi n m                                 \n\
                                                                            \n\
      Create n+m additional copies of the geometry defined so far along the \n\
      direction given by the polar angle phi in degrees with respect to the \n\
      x axis. The spacing between copies is a. m copies go in the positive  \n\
      direction of phi and n copies go in the negative direction.           \n\
                                                                            \n\
 -->  periodicity  x0 y0  a1 a2  alpha beta  n1 n2                          \n\
                                                                            \n\
      Set up lattice parameters to be used with EPW calculations assuming   \n\
      that the system is periodic, with the origin of the first unit        \n\
      cell at (x0,y0) and with lattice parameters a1, a2, alpha1, alpha2    \n\
      as defined in the \"repeat-periodic-geometry ...\" command.           \n\
      n1 and n2 control the number of points used to obtain                 \n\
      potential matrices for EPW calculations (see below).                  \n\
      These two parameters have to be large enough to guarantee             \n\
      convergence whenever a periodic system is considered.                 \n\
      Important note: the user has to make sure that the boundaries are all \n\
      defined inside the first unit cell.                                   \n\
                                                                            \n\
 -->  begin-virtual-structure                                               \n\
                                                                            \n\
      After this command is invoked, further boundary elements are only     \n\
      used to determine internally the region in which (x,y) points are     \n\
      located, but they are not employed in the final parametrization for   \n\
      solving Schrodinger equations. This command is useful when open       \n\
      structures are needed (e.g., an infinite cone), since the internal    \n\
      determination of the regions require that the boundaries are          \n\
      closed (no open ends).                                                \n\
                                                                            \n\
 -->  write-geometry output-file                                            \n\
                                                                            \n\
      Print out information on discretization points as defined so far      \n\
      using the above commands. The output of this command is written       \n\
      in file output-file. Each line of the output contains the following   \n\
      information: i, xs, ys, ths, mu1, mu2, nxs, nys, coef, and ns, where  \n\
      (xs,ys) are the coordinates of point i, (nxs,nys) is the boundary     \n\
      normal, which must point towards medium mu2 (Important warning: the   \n\
      geometry, and in particular the orientation of the normal, must be    \n\
      checked by the user with this command before calculating LDOS, PE     \n\
      intensities, and bands), medium mu1 is on the other side of the       \n\
      boundary, ths is the boundary parameter [x=xs(ths), y=ys(ths)],       \n\
      ns is the Jacobian of this parametrization, and coef is the           \n\
      parameter increment used to integrate along the boundaries.           \n\
                                                                            \n\
 -->  write-region mu n output-file                                         \n\
                                                                            \n\
      Print out file output-file with those points (x,y) of a grid of       \n\
      n by n points that are contained inside region mu. Each output line   \n\
      contains \"x y\".                                                     \n\
                                                                            \n\
 -->  write-regions n output-file                                           \n\
                                                                            \n\
      Print out file output-file with points covering the entire            \n\
      geometry using a grid of n by n points. Each output line has the      \n\
      format \"x y mu\", where mu is the region where (x,y) is located.     \n\
                                                                            \n\
 -->  clear-all                                                             \n\
                                                                            \n\
      Clear the geometry and disregard all previous fragments. Also,        \n\
      set all parameters to their default values (non-periodic system       \n\
      with effective mass =1, potentials and barriers =0, etc.).            \n\
                                                                            \n\
II.2 DESCRIPION OF COMMANDS USED TO SET Im{E}, meff, V, AND BARRIERS        \n\
                                                                            \n\
 -->  Ei value                                                              \n\
                                                                            \n\
      Assign Im{E}=value (common to all regions), where E is the            \n\
      energy of the surface electron (Im{E}=0 by default). Notice that      \n\
      one must have Ei>=0 in order to represent physical systems where      \n\
      the number of electrons either stays constant (Ei=0) or decreases     \n\
      (Ei>0) due to coupling to non-elastic channels, etc.                  \n\
                                                                            \n\
 -->  meff value                                                            \n\
                                                                            \n\
      Assign meff=value, the effective electron mass (meff=1 by default).   \n\
                                                                            \n\
 -->  potential mu Vr Vi                                                    \n\
                                                                            \n\
      Assign V_mu=Vr+i*Vi, the complex potential experienced by the         \n\
      electron in region mu (V_mu=0 by default). Vi<=0 !!!                  \n\
                                                                            \n\
 -->  barrier mu mup ar ai                                                  \n\
                                                                            \n\
      Introduce a delta-function barrier (ar+i*ai)*delta_s                  \n\
      at the boundary between regions mu and mup, where delta_s is          \n\
      a delta function in the direction normal to the boundary.             \n\
      ar=ai=0 for all boundaries by default.                                \n\
                                                                            \n\
 -->  infinite-barriers                                                     \n\
                                                                            \n\
      Set infinite barriers at all boundaries. Do not use the commands      \n\
      \"potential ...\" and \"barrier ...\" in the same input file where    \n\
      \"infinite-barriers\" are used. The boundaries have to be defined     \n\
      so that mu=0 is the region where the wave function is calculated.     \n\
      Other regions will be excluded from the calculation of LDOS and PE.   \n\
      This command cannot be used in combination with periodic systems.     \n\
                                                                            \n\
II.3 CALCULATIONS OF LDOS, PE, AND WAVE FUNCTIONS BASED UPON THE EBEM       \n\
                                                                            \n\
The EBEM has been designed to calculate LDOS, PE, and wave functions in     \n\
finite structures. Important: the calculated LDOS does not include spin     \n\
degeneracy, so it must be multiplied by a factor of 2 for electrons.        \n\
                                                                            \n\
II.3a Calculation of the local density of states (LDOS)                     \n\
                                                                            \n\
For finite systems, the LDOS is calculated directly for a probe             \n\
at a given position in 2D real space, in a way similar to the mode of       \n\
operation of scanning tunneling microscopes. Important: LDOS without spin.  \n\
                                                                            \n\
 -->  LDOS-rE-BEM x0 y0 x1 y1 nr E0 E1 nE output-file                       \n\
                                                                            \n\
      Calculate LDOS along the segment going from (x0,y0) to (x1,y1)        \n\
      (nr points along this segment) for electron energies (real part)      \n\
      from E0 to E1 (a total of nE energies), and write the result in       \n\
      file output-file with the format \"r E LDOS\" in each row (nr*nE      \n\
      rows in total), where r is the distance along the segment.            \n\
                                                                            \n\
 -->  LDOS-xy-BEM x0 y0 dx dy Er nx ny output-file                          \n\
                                                                            \n\
      Calculate LDOS in a rectangle of lower-left corner                    \n\
      (x0,y0) and of sides dx and dy along the x and y directions,          \n\
      respectively. The results are stored in file output-file              \n\
      with the format \"x y LDOS\" in each line. The points                 \n\
      (x,y) form a grid of nx points along x and ny points along y.         \n\
      The surface electron energy is Er (real part).                        \n\
                                                                            \n\
 -->  LDOS-xyE-BEM x0 y0 dx dy nx ny E0 E1 nE output-file                   \n\
                                                                            \n\
      Calculate LDOS as a function of (x,y,E). Each output line contains    \n\
      E x y LDOS.                                                           \n\
                                                                            \n\
II.3b Calculation of photoemission intensities (PE)                         \n\
                                                                            \n\
PE intensities are calculated in the EBEM by integrating over a real-space  \n\
area determined via the command \"PE-area ...\".                            \n\
                                                                            \n\
The emission occurs at discrete energies for confined systems.              \n\
Important note: It it is convenient to                                      \n\
introduce either a positive imaginary part in the final PE energy (Ei) or   \n\
a negative imaginary part in the electron potential (Vi) at some or all     \n\
regions in order to capture the bulk of the PE intensity at the points      \n\
of the energy grids used in the following commands.                         \n\
                                                                            \n\
 -->  PE-area-BEM  x0 y0  a1 a2  alpha beta  n1 n2                          \n\
                                                                            \n\
      Define an integration cell with the same parameters as in the         \n\
      command \"periodicity ...\", but without assuming that the            \n\
      system is periodic. This is used to calculate PE from finite          \n\
      structures with localized states using EBEM. n1 and n2 are            \n\
      integrationg grid parameters.                                         \n\
                                                                            \n\
 -->  PE-qE-BEM qx0 qy0 qx1 qy1 nQ E0 E1 nE output-file                     \n\
                                                                            \n\
      Calculate PE intensities along the segment going from                 \n\
      (qx0,qy0) to (qx1,qy1) in photoelectron 2D parallel-momentum          \n\
      space q (nQ points along this segment) for electron                   \n\
      energies (real part) from E0 to E1 (a total of nE energies),          \n\
      and write the result in file output-file with the                     \n\
      format \"q E PE-intensity\" in each row (nQ*nE rows in total).        \n\
                                                                            \n\
 -->  PE-qxqy-BEM qx0 qy0 dqx dqy Er nx ny output-file                      \n\
                                                                            \n\
      Calculate PE intensities in a rectangle of lower-left corner          \n\
      (qx0,qy0) and of sides dqx and dqy along the x and y directions,      \n\
      respectively, in photoelectron 2D parallel-momentum space q.          \n\
      The results are saved in file output-file with the                    \n\
      format \"qx qy PE-intensity\" in each line. The points                \n\
      (qx,qy) form a grid with nx points along x and ny points along y.     \n\
      The surface electron energy is Er (real part).                        \n\
                                                                            \n\
II.3c Calculation of wave functions                                         \n\
                                                                            \n\
The following command calculate wave functions                              \n\
corresponding to specific excitations of a finite system                    \n\
(e.g., a localized point source, like what one would have when an           \n\
electron is launched from a STM tip onto the surface).                      \n\
                                                                            \n\
 -->  WF-xy-1 xa ya x0 y0 dx dy Er nx ny output-file                        \n\
                                                                            \n\
      Calculate wave function in a rectangle of lower-left corner           \n\
      (x0,y0) and of sides dx and dy along the x and y directions,          \n\
      respectively, corresponding to an electron point source at (xa,ya).   \n\
      The results are stored in file output-file with the format            \n\
      \"x y |wf|^2 argument{wf}\" in each line. The points                  \n\
      (x,y) form a grid of nx points along x and ny points along y.         \n\
      The surface electron energy is Er (real part).                        \n\
      This command works only for finite systems.                           \n\
                                                                            \n\
 -->  WF-xy-2 xa ya phi D n mu0 x0 y0 dx dy Er nx ny output-file            \n\
                                                                            \n\
      Same as WF-xy-1, but for an incoming Gaussian beam moving             \n\
      towards the azumithal direction phi, with approximate focal amplitude \n\
      width D, with focus at (xa,ya), and constructed from a superposition  \n\
      of 2*n+1 plane waves (n=50 is typically enough). Use n=0 to have      \n\
      just a plane wave. The incoming wave function takes a value of 1      \n\
      at (xa,ya). Notice that the actual focal width cannot be smaller      \n\
      than approximately half a wavelength (diffraction limit).             \n\
      The beam is assumed to come from medium mu0.                          \n\
                                                                            \n\
 -->  WF-bound xa ya Er output-file                                         \n\
                                                                            \n\
      Same as WF-xy-1, with the wave function given at the boundary points. \n\
      Output: \"xi yi dli nxi nyi |wf|^2 arg{wf}\", where                   \n\
      dli is the length of interval i centered at (xi,yi) and with normal   \n\
      vector (nxi,nyi) in the parametrization.                              \n\
                                                                            \n\
 -->  WF-rbound x0 y0 x1 y1 nr Er output-file                               \n\
                                                                            \n\
      Same as WF-bound-1, with the excitation point along a segment, like   \n\
      in WF-rE...  Output: \"r xi yi dli nxi nyi |wf|^2 arg{wf}\".          \n\
                                                                            \n\
II.4 CALCULATIONS OF 2D BAND STRUCTURES, LDOS, AND PE BASED UPON EPW        \n\
                                                                            \n\
The commands described in this section can be used to calculate band        \n\
structures, LDOS, and photoemission intensitites using a plane-wave         \n\
expansion of the electronic wave functions (EPW).                           \n\
They can be used only for periodic systems and finite barriers.             \n\
                                                                            \n\
II.4a Calculation of band strcuture for periodic systems (BANDS)            \n\
                                                                            \n\
Band structures can be calculated with this code using a plane-wave         \n\
representation of the 2D potential and wave functions. Bands for            \n\
hexagonal or square lattices of simple objects can be calculated using      \n\
an inline command, as shown in Sec. III. Other lattices or more complex     \n\
unit cells can be calculated as well, also using plane waves, with the      \n\
geometry and lattice description introduced via an input file using the     \n\
above commands.                                                             \n\
                                                                            \n\
 -->  BANDS-hexagonal gmax nQ Emax output-file                              \n\
                                                                            \n\
      Calculate band structure for a periodic system with hexagonal         \n\
      lattice. The electron eigen-energies are given for nQ 2D momentum     \n\
      points along an excursion of representative points within the first   \n\
      Brillouin zone: Gamma-M-K-Gamma. gmax controls the number of plane    \n\
      waves used in the calculation (gmax lies typically in the 4-12        \n\
      range). The output file consists of lines with the format             \n\
      \"Q Re{E} Im{E}\", where E is the eigen-energy corresponding          \n\
      to point Q, and Q varies from 0 to 3^0.5 to 1+3^0.5 to 3+3^0.5        \n\
      in the excursion from Gamma to M to K to Gamma.                       \n\
                                                                            \n\
 -->  BANDS-square gmax nQ Emax output-file                                 \n\
                                                                            \n\
      Similar to the previous command, but for square lattices, with nQ     \n\
      points along the excursion Gamma-X-M-Gamma, and with Q running from   \n\
      0 to 1 to 2 to 2+2^0.5 when going from Gamma to X to M to Gamma.      \n\
                                                                            \n\
 -->  BANDS gmax qx0 qy0 qx1 qy1 nQ Emax output-file                        \n\
                                                                            \n\
      Similar to the previous commands, but for arbitrary lattices, with    \n\
      nQ points sampled from (qx0,qy0) to (qx1,qy1) in momentum space.      \n\
                                                                            \n\
II.4b Calculation of the local density of states (LDOS) using EPW           \n\
                                                                            \n\
LDOS can be calculated for periodic systems using EPW. This method is       \n\
fast, but a large number of cell points (n1 and n2 parameters in command    \n\
\"periodicity ...\" above) have to be given, and either the potential       \n\
or the electron energy (Ei command) have to have a non-zero imaginary part  \n\
in order to obtain non-zero results for LDOS. Those imaginary parts have    \n\
to be sufficiently large as to define bands with a width in momentum larger \n\
than the grid spacing defined by n1 and n2. Just run calculations with      \n\
increasing values of n1 and n2 until a converged result is obtained.        \n\
                                                                            \n\
The execution of these commands can take considerable time. A line is       \n\
printed with the number of Q points that have been calculated already as    \n\
the program proceeds (the total number of Q points is n1*n2).               \n\
                                                                            \n\
 -->  LDOS-rE-EPW gmax x0 y0 x1 y1 nr E0 E1 nE output-file                  \n\
                                                                            \n\
      Calculate LDOS in the same way as with the command \"LDOS-rE-BEM      \n\
      ...\", but using EPW and assuming that the system is periodic. See    \n\
      commands \"BANDS ...\" for the meaning of parameter gmax.             \n\
                                                                            \n\
 -->  LDOS-xy-EPW gmax x0 y0 dx dy Er nx ny output-file                     \n\
                                                                            \n\
      Calculate LDOS in the same way as with the command \"LDOS-xy-BEM      \n\
      ...\", but using EPW and assuming that the system is periodic. See    \n\
      commands \"BANDS ...\" for the meaning of parameter gmax.             \n\
                                                                            \n\
 -->  DOS-EPW gmax E0 E1 nE output-file                                     \n\
                                                                            \n\
      Calculate total DOS in the same way as with the command \"LDOS-rE-BEM \n\
      ...\" (i.e., the integral of the LDOS over the first unit cell).      \n\
                                                                            \n\
II.4c Calculation of photoemission intensitites using EPW                   \n\
                                                                            \n\
 -->  PE-qE-hexagonal-EPW gmax nQ E0 E1 nE output-file                      \n\
                                                                            \n\
      Calculate PE intensities along the Gamma-M-K-Gamma excursion within   \n\
      the first Brillouin zone in photoelectron 2D parallel-momentum space  \n\
      (the system must have the periodicity of a hexabonal lattice and      \n\
      nQ points along the mentioned excursion are calculated) for electron  \n\
      energies (real part) ranging from E0 to E1 (a total of nE energies),  \n\
      and write the result in file output-file with the                     \n\
      format \"q E PE-intensity\" in each row (a total of nQ*nE rows),      \n\
      where q is 0, 3^0.5, 1+3^0.5, and 3+3^0.5 in Gamma, M, K, and         \n\
      Gamma, respectively.                                                  \n\
      The parameter gmax controls the number of plane waves. Check out      \n\
      the description of \"BANDS-hexagonal ...\" above (Sec. II.4a) for     \n\
      further details on the use of of plane waves and the meainig of gmax. \n\
                                                                            \n\
 -->  PE-qE-square-EPW gmax nQ E0 E1 nE output-file                         \n\
                                                                            \n\
      Calculate PE intensities along the Gamma-X-M-Gamma excursion within   \n\
      the first Brillouin zone in photoelectron 2D parallel-momentum space  \n\
      (the system must have the periodicity of a square lattice and         \n\
      nQ points along the mentioned excursion are calculated) for electron  \n\
      energies (real part) ranging from E0 to E1 (a total of nE energies),  \n\
      and write the result in file output-file with the                     \n\
      format \"q E PE-intensity\" in each row (a total of nQ*nE rows),      \n\
      where q is 0, 1, 2, and 2+2^0.5 in Gamma, X, M, and Gamma,            \n\
      respectively. See above for a description of gmax.                    \n\
                                                                            \n\
 -->  PE-qE-EPW gmax qx0 qy0 qx1 qy1 nQ E0 E1 nE output-file                \n\
                                                                            \n\
      Calculate PE intensities along the segment going from                 \n\
      (qx0,qy0) to (qx1,qy1) in photoelectron 2D parallel-momentum          \n\
      space q (nQ points along this segment) for electron                   \n\
      energies (real part) from E0 to E1 (a total of nE energies),          \n\
      and write the result in file output-file with the                     \n\
      format \"q E PE-intensity\" in each row (nQ*nE rows in total).        \n\
                                                                            \n\
 -->  PE-qxqy-EPW gmax qx0 qy0 dqx dqy Er nx ny output-file                 \n\
                                                                            \n\
      Calculate PE intensities in a rectangle of lower-left corner          \n\
      (qx0,qy0) and of sides dqx and dqy along the x and y directions,      \n\
      respectively, in photoelectron 2D parallel-momentum space q.          \n\
      The results are saved in file output-file with the                    \n\
      format \"qx qy PE-intensity\" in each line. The points                \n\
      (qx,qy) form a grid with nx points along x and ny points along y.     \n\
      The surface electron energy is Er (real part).                        \n\
                                                                            \n\
II.5 EXIT THE EXECUTION OF THE PROGRAM AND INSERT INPUT FILE                \n\
                                                                            \n\
 -->  end                                                                   \n\
 -->  exit                                                                  \n\
 -->  quit                                                                  \n\
                                                                            \n\
      Any of these commands will finish the execution of the program.       \n\
                                                                            \n\
 -->  include file-name                                                     \n\
                                                                            \n\
      This interprets the contents of input file file-name. Notice          \n\
      that file-name has to end with an end, exit, or quit command. This    \n\
      is useful to store pre-defined shapes as separate input files.        \n\
                                                                            \n\
III - CALCULATION OF BAND STRUCTURES USING AN INLINE COMMAND                \n\
                                                                            \n\
The EBEM can also calculate one-electron bands                              \n\
of electrons subject to a 2D periodic potential within a plane wave         \n\
approach using an inline command with the following format:                 \n\
                                                                            \n\
      'ebem epw output-file a R gmax nQ Re{U1} Re{U0}                       \n\
               [object=0 lattice=0 Im{U1}=0 Im{U0}=0                        \n\
                meff=1 Emax=infinity]',                                     \n\
                                                                            \n\
where parameters inside square brackets are optional and their default      \n\
values are shown above. This mode of operation allows one to calculate 2D   \n\
electron band structures for one electron subject to hexagonal or square    \n\
lattices (lattice equals 0 or 1, respectively) with one circle or triangle  \n\
(object equals 0 or 1, respectively) inside each unit cell. Here, a is the  \n\
lattice constant, R is the circle radius or triangle side in each case,     \n\
U1 and U0 are the values of the complex electron potential inside and       \n\
outside the circle/trianle, respectively. meff is the effective electron    \n\
mass. The output-file format, the role of gmax, and the momenta for which   \n\
eigen-energies are calculated are the same as with commands                 \n\
BANDS-hexagonal and BANDS-square (see above) for lattice=0 and 1,           \n\
respectively.                                                               \n\
                                                                            \n\
**************************************************************************  \n"

#include "bin/jga.h"        // Basic tools (constants, numero, etc.).
#include "bin/complex.h"    // Complex numbers.
#include "bin/matrix.h"     // Matrices (handling, inversion, etc.).
#include "bin/bessel.h"     // Some Bessel functions (Jm, Ym, Km, etc.).
#include "bin/lattice.h"    // 2D lattice (reciprocal lattice, etc.).
#include "bin/eispack.h"    // Matrix diagonalization.
#include "bin/algebra.h"    // Input algebra.

// **************************************************************************
// *** Parameters controling allocation of memory for EBEM.               ***
// **************************************************************************

#define n_regions   20  // Maximum number of regions (from 0 to n_regions-1).
#define nl_max   100000  // Maximum number of boundary segments.

// **************************************************************************
// *** Boundary geometry.                                                 ***
// **************************************************************************

// --- Homogeneous regions limited by boundaries.

numero    meff=1;                    // Effective electron mass.
complex   Vmu[n_regions];            // Potential of each region.
complex   kmu[n_regions];            // Electron momentum in each region.
complex   aa[n_regions][n_regions];  // Delta-function barrier.
int       finite_barriers=1;         // 0/1 for infinite/finite barriers.

// --- Boundary fragments (segments and arcs) limiting the regions.

numero xx1[nl_max],yy1[nl_max],xx2[nl_max],yy2[nl_max]; // Fragment
numero th0[nl_max],th1[nl_max],galph[nl_max];           // parametrization.
int    gnn[nl_max];                  // Number of points in each fragment.
int    gmu1[nl_max],gmu2[nl_max];    // Regions on either side of boundary.
int    gtype[nl_max];                // Type of fragment: 0=segment, 1=arc.
int    nl=0;                         // Current number of fragments.
int    mumin,mumax;                  // Range of mu values actually used.
int    mu_env;                       // Region outside geometry.

numero xth(numero th)
{
  int j=0;  while(th>1) {th--; j++;}
  if(gtype[j]==0)  return xx1[j]+(xx2[j]-xx1[j])*th;
  if(gtype[j]==1)
    return  xx1[j] + xx2[j]*cos(th0[j]+th*(th1[j]-th0[j]))*cos(galph[j])
                   - yy2[j]*sin(th0[j]+th*(th1[j]-th0[j]))*sin(galph[j]);
  return 0;
}

numero yth(numero th)
{
  int j=0;  while(th>1) {th--; j++;}
  if(gtype[j]==0)  return yy1[j]+(yy2[j]-yy1[j])*th;
  if(gtype[j]==1)
    return  yy1[j] + xx2[j]*cos(th0[j]+th*(th1[j]-th0[j]))*sin(galph[j])
                   + yy2[j]*sin(th0[j]+th*(th1[j]-th0[j]))*cos(galph[j]);
  return 0;
}

numero xthp(numero th)
{
  int j=0;  while(th>1) {th--; j++;}
  if(gtype[j]==0)  return xx2[j]-xx1[j];
  if(gtype[j]==1)
  return -xx2[j]*(th1[j]-th0[j])*sin(th0[j]+th*(th1[j]-th0[j]))*cos(galph[j])
         -yy2[j]*(th1[j]-th0[j])*cos(th0[j]+th*(th1[j]-th0[j]))*sin(galph[j]);
  return 0;
}

numero ythp(numero th)
{
  int j=0;  while(th>1) {th--; j++;}
  if(gtype[j]==0)  return yy2[j]-yy1[j];
  if(gtype[j]==1)
  return -xx2[j]*(th1[j]-th0[j])*sin(th0[j]+th*(th1[j]-th0[j]))*sin(galph[j])
         +yy2[j]*(th1[j]-th0[j])*cos(th0[j]+th*(th1[j]-th0[j]))*cos(galph[j]);
  return 0;
}

void add_segment(int n, int mu1, int mu2,
                 numero x1, numero y1, numero x2, numero y2)
{
  if(nl>=nl_max) on_error("EBEM: init geometry","too many segments =",nl+1,"");
  gnn[nl]=n;  xx1[nl]=x1;  yy1[nl]=y1;  xx2[nl]=x2;  yy2[nl]=y2;
  th0[nl]=th1[nl]=galph[nl]=0;
  gmu1[nl]=mu1;  gmu2[nl]=mu2;  gtype[nl]=0;
  if(mu1<0||n_regions<=mu1) on_error("EBEM: regions:","out of range:",mu1,"");
  if(mu2<0||n_regions<=mu2) on_error("EBEM: regions:","out of range:",mu2,"");
  nl++;
}

void add_arc(int n, int mu1, int mu2, numero x0, numero y0,
             numero a, numero b, numero th0_, numero th1_, numero galph_)
{
  if(nl>=nl_max) on_error("EBEM: init geometry","too many segments =",nl+1,"");
  gnn[nl]=n;  xx1[nl]=x0;  yy1[nl]=y0;  xx2[nl]=a;  yy2[nl]=b;
  th0[nl]=th0_;  th1[nl]=th1_;  galph[nl]=galph_;
  gmu1[nl]=mu1;  gmu2[nl]=mu2;  gtype[nl]=1;
  if(mu1<0 || n_regions<=mu1) on_error("EBEM: regions","out of range:",mu1,"");
  if(mu2<0 || n_regions<=mu2) on_error("EBEM: regions","out of range:",mu2,"");
  nl++;
}

// --------------------------------------------------------------------------

int      nth=0;                 // Number of parametrization points.
numero   *xs,*ys;               // (x,y) for the different points.
numero   *nxs,*nys;             // Boundary normal (it points towards mu2).
numero   *ths,*ns;              // Parameter th, Jacobian (R')^2+(z')^2)^0.5.
numero   *coef;                 // Increment in th.
int      *mu1,*mu2;             // Regions types.
complex  *aamm;                 // 2*aa*meff.
int nth_real_structure=-1,nthh; // real/virtual structures

numero   period_x0,period_y0;   // Origin of first unit cell.
numero   period_a1x,period_a1y,period_a2x,period_a2y;  // Unit cell vectors.
int      period_n1=0,period_n2=0;              // Param. pts. of cell bound.
int      period_flag=0;         // 0/1 for non-periodic/periodic systems.

void init_periodicity(FILE *fin)
{
  period_x0=alread_numero(fin);     period_y0=alread_numero(fin);
  a0=alread_numero(fin);            b0=alread_numero(fin);
  alpha=alread_numero(fin)*pi/180;  beta=alread_numero(fin)*pi/180;
  period_n1=alread_int(fin);        period_n2=alread_int(fin);
  period_a1x=a0*cos(alpha);       period_a2x=b0*cos(alpha+beta);
  period_a1y=a0*sin(alpha);       period_a2y=b0*sin(alpha+beta);
  lattice(1);  lattice_free();  // Reciprocal vectors (see bin/lattice.h).
}

void init_geometry(int &i, int j)   // Initialize fragment j, starting with
{                                   // point i (on input) and returning
  int ii;                           // starting point i of next fragment
                                    // (on output).
  if(gmu1[j]==gmu2[j])
    on_error("EBEM: init geometry","barriers must separate different regions");

  for(ii=0; ii<gnn[j]; ii++, i++) {
    ths[i]=j+(ii+0.5)/gnn[j];  coef[i]=1.0/gnn[j];
    xs[i]=xth(ths[i]);    ys[i]=yth(ths[i]);
    nxs[i]=ythp(ths[i]);  nys[i]=-xthp(ths[i]);
    ns[i]=sqrt(nxs[i]*nxs[i]+nys[i]*nys[i]);
    nxs[i]=nxs[i]/ns[i];  nys[i]=nys[i]/ns[i];
    mu1[i]=gmu1[j];       mu2[i]=gmu2[j];
    if(mu1[i]>=0 && mu1[i]<mumin) mumin=mu1[i];  if(mu1[i]>mumax) mumax=mu1[i];
    if(mu2[i]>=0 && mu2[i]<mumin) mumin=mu2[i];  if(mu2[i]>mumax) mumax=mu2[i];
    if(mu1[i]<0 || mu2[i]<0)  aamm[i]=0;  else
    aamm[i]=2*meff*aa[mu1[i]][mu2[i]];
} }

void free_geometry(void)
{
  mumin=n_regions;  mumax=-1;  nth_real_structure=-1;
  if(nth>0) {
    delete [] xs;  delete [] ys;  delete [] nxs;  delete [] nys;
    delete [] ths;  delete [] ns;  delete [] coef;
    delete [] mu1;  delete [] mu2;  delete [] aamm;
    nth=0;
} }

void init_geometry(void)
{
  numero x1,y1,x2,y2,x3,y3,x4,y4;
  int i,j;

  nthh=nth_real_structure;  free_geometry();  nth_real_structure=nthh;
  for(j=0; j<nl; j++)  nth+=gnn[j];
  if(nth<=0) on_error("EBEM: init_geometry", "a structure must be defined");
  xs=new numero [nth];  ys=new numero [nth];
  nxs=new numero [nth];  nys=new numero [nth];
  ths=new numero [nth];  ns=new numero [nth];  coef=new numero [nth];
  mu1=new int [nth];  mu2=new int [nth];  aamm=new complex [nth];

  for(i=j=0; j<nl; j++)  init_geometry(i,j);         // Actual geometry.

  numero xx=-infinite;
  for(i=0; i<nth; i++) if(xs[i]>xx) {xx=xs[i]; j=i;}
  if(nxs[j]>0) mu_env=mu2[j]; else mu_env=mu1[j];    // Environment region.

  if(nth_real_structure>0) {nthh=nth;  nth=nth_real_structure;}
  else nthh=nth;
}

void end_real_structure(void) {init_geometry();  nth_real_structure=nth;}

int where_is_the_point(int mu0, numero x, numero y)  // If mu0>=0, this routine
{                                                    // returns 1 (0) if (x,y)
  complex point_a[n_regions], val0;                  // is (is not) inside
  numero v0,vv;                                      // region mu0. If mu0<0,
  int i,j,mu,mmu;                                    // this routine returns
                                                     // the value of the region
  for(mu=mumin; mu<=mumax; mu++) point_a[mu]=0;      // where (x,y) is located.

  for(i=0; i<nthh; i++) {
    val0=ns[i]*coef[i]*complex(-nys[i],nxs[i])/complex(xs[i]-x,ys[i]-y);
    if(mu1[i]>=0) point_a[mu1[i]]+= val0;
    if(mu2[i]>=0) point_a[mu2[i]]+=-val0;
  }

  for(mu=mumin, v0=1, mmu=mu_env; mu<=mumax; mu++) {
    vv=real(point_a[mu]/(i_c*pi));  if(vv>v0) {v0=vv; mmu=mu;}
  }

  if(mu0>=0)  return (mu0==mmu);

  return mmu;
}

int where_is_the_point(numero x, numero y) {return where_is_the_point(-1,x,y);}

// --------------------------------------------------------------------------

int n_fragments=0;

void copy_geometry(numero dx, numero dy)
{
  int i;

  for(i=0; i<n_fragments; i++)
    if(gtype[i]==1) add_arc(gnn[i],gmu1[i],gmu2[i],xx1[i]+dx,yy1[i]+dy,
                            xx2[i],yy2[i],th0[i],th1[i],galph[i]);
    else            add_segment(gnn[i],gmu1[i],gmu2[i],xx1[i]+dx,yy1[i]+dy,
                                xx2[i]+dx,yy2[i]+dy);
}

void rotate_point(numero xr, numero yr, numero phi, numero &x, numero &y)
{
  numero r=sqrt(sqr(x-xr)+sqr(y-yr));
  numero ff=atan2(y-yr,x-xr);
  x=xr+r*cos(ff+phi);
  y=yr+r*sin(ff+phi);
}

void rotate_geometry(numero xr, numero yr, numero phi)
{
  int i;

  for(i=0; i<n_fragments; i++)
    if(gtype[i]==1) {
      rotate_point(xr,yr,phi,xx1[i],yy1[i]);
      galph[i]+=phi;
    } else {
      rotate_point(xr,yr,phi,xx1[i],yy1[i]);
      rotate_point(xr,yr,phi,xx2[i],yy2[i]);
}   }

void translate_geometry(numero Dx, numero Dy)
{
  int i;

  for(i=0; i<n_fragments; i++)
    if(gtype[i]==1) {
      xx1[i]=xx1[i]+Dx;  yy1[i]=yy1[i]+Dy;
    } else {
      xx1[i]=xx1[i]+Dx;  yy1[i]=yy1[i]+Dy;
      xx2[i]=xx2[i]+Dx;  yy2[i]=yy2[i]+Dy;
}   }

void write_geometry(char *name)
{
  FILE *fout; fout=fopen(name,"w");
  int i;
  fprintf(fout, "# i       xs       ys     ths  m1  m2      nxs      nys");
  fprintf(fout, "     coef       ns\n");
  for(i=0; i<nth; i++)
    fprintf(fout,"%3d %8.3f %8.3f %7.2f  %2d  %2d  %7.2f  %7.2f %8.3f %8.2f\n",
                 i, xs[i], ys[i], ths[i], mu1[i], mu2[i],
           nxs[i], nys[i], coef[i], ns[i]);
  fclose(fout);
}

void write_region(int mu, int ll, char *name)
{
  numero x0=infinite,y0=infinite,x1=-infinite,y1=-infinite, x,y;
  int i,j;

  for(i=0; i<nth; i++) if(mu1[i]==mu || mu2[i]==mu || mu<0) {
    if(xs[i]<x0) x0=xs[i];  if(ys[i]<y0) y0=ys[i];
    if(xs[i]>x1) x1=xs[i];  if(ys[i]>y1) y1=ys[i];
  }

  FILE *fout; fout=fopen(name,"w");
  for(i=0; i<ll; i++)
  for(j=0; j<ll; j++) {
    if(ll<=1) x=x0; else x=x0+i*(x1-x0)/(ll-1);
    if(ll<=1) y=y0; else y=y0+j*(y1-y0)/(ll-1);
    if(mu>=0) if(where_is_the_point(mu,x,y)) fprintf(fout, "%g %g\n", x, y);
    if(mu<0) fprintf(fout, "%g %g %d\n", x, y, where_is_the_point(x,y));
  }
  fclose(fout);
}

void clear_all(void)
{
  meff=1;
  finite_barriers=1;  nl=0;  n_fragments=0;  free_geometry();
  period_n1=period_n2=0;  period_flag=0;
  int i,j;
  for(i=0; i<n_regions; i++) {
    Vmu[i]=0;
    for(j=0; j<n_regions; j++) aa[i][j]=0;
} }

// **************************************************************************
// *** Calculation of Green functions.                                    ***
// **************************************************************************

complex besselHm(int m, complex z)
{
  numero J,Y,Jp,Yp, mz=mod(z);
  complex Hm;

  if(mz<1e-30) on_error("EBEM (internal)",
    "either E-V takes a vanishing value or the electron wave function is",
    "evaluated at a discretization point");

  if(ABS(imag(z))<1e-5*mz) {
    besselJY(ABS(m),real(z),J,Y,Jp,Yp);  Hm=complex(J,Y);
    if(m<0 && (-m)%2!=0) Hm=-Hm;
  } else Hm=2/pi/i_l(m+1)*besselK(m,-i_c*z);

  return Hm;
}

complex GG(int mu, numero xi, numero yi, numero nxi, numero nyi,
           numero x0, numero y0, int index)
{
  if(mu<mumin || mumax<mu)
    on_error("EBEM: Green function", "region", mu, "is out of range");

  numero  x,y,R;
  complex kk=kmu[mu];

  x=xi-x0;  y=yi-y0;

  R=sqrt(x*x+y*y);   // Set index=1 or 2 for G or dG/dn, respectively.

  if(index==1)  return  besselHm(0,R*kk) * (-i_c*meff/2);
  else          return -besselHm(1,R*kk) * (-i_c*meff/2)
                        * kk * (nxi*x+nyi*y)/R;
}

complex GG(int mu, int i, numero x0, numero y0, int index)
  {return GG(mu,xs[i],ys[i],nxs[i],nys[i],x0,y0,index);}

int muG,iG,jG,indexG;

complex GG_thp(numero thp)
{
  if(thp==ths[jG]) return ns[jG]*GG(muG,iG,xs[jG],  ys[jG],  indexG);
                   return ns[jG]*GG(muG,iG,xth(thp),yth(thp),indexG);
}

// --------------------------------------------------------------------------

int nG=50, ijmax=4;     // Convergence parameters.
complex eqa1,eqa2;      // Phase factors for periodic structures.

complex GGj(int mu, int i, int j, int index)
{
  int nnG;  if(i==j)  nnG=5*nG;  else  nnG=nG;

  muG=mu;  iG=i;  jG=j;  indexG=index;

  if(ABS(i-j)>ijmax)  return coef[j]*GG_thp(ths[j]);

  int l;                                          // Integration for
  complex sum=0;                                  // neighboring intervals
  numero  a, b, h=coef[j]/2/nnG;                  // with open integration.

  a=ths[j]-coef[j]/2;  b=ths[j];
  sum+=(109*(GG_thp(a+h)+GG_thp(b-h))-5*(GG_thp(a+2*h)+GG_thp(b-2*h))
       +63*(GG_thp(a+3*h)+GG_thp(b-3*h))+49*(GG_thp(a+4*h)+GG_thp(b-4*h)))/48;
  for(l=5; l<=nnG-5; l++)  sum+=GG_thp(a+l*h);

  a=ths[j];  b=ths[j]+coef[j]/2;
  sum+=(109*(GG_thp(a+h)+GG_thp(b-h))-5*(GG_thp(a+2*h)+GG_thp(b-2*h))
       +63*(GG_thp(a+3*h)+GG_thp(b-3*h))+49*(GG_thp(a+4*h)+GG_thp(b-4*h)))/48;
  for(l=5; l<=nnG-5; l++)  sum+=GG_thp(a+l*h);

  return sum*h;
}

// **************************************************************************
// *** Eigensystem.                                                       ***
// **************************************************************************

matrix G1_1,G2_1,S1,S2,D_1;  // Full matrices of the EBEM.
matrix sigma1,sigma2;        // Boundary charges.
matrix psi,psip;             // Boundary sources.

// The photoelectron parallel momentum for PE is defined in
// bin/lattice.h as variables 'numero Qx,Qy;'.

#define mumu(m1,m2,i) ((m1==m2 && m1>=0) || (m1<0 && m2==mu_env))

void init(complex E)
{
  int i,j;  init_geometry();

  // regions initialization

  for(i=0; i<n_regions; i++)  kmu[i]=sqrt(2*meff*(E-Vmu[i])+i_c*1e-20);

  // initialize eigensystem

  matrix G1,G2,H1,H2;

  G1.alloc(nth,nth);
  if(finite_barriers) {
    G2.alloc(nth,nth);  H1.alloc(nth,nth);  H2.alloc(nth,nth);
  }

  for(j=0; j<nth; j++)
  for(i=0; i<nth; i++) if(finite_barriers) {

    if(mumu(mu1[i],mu1[j],i)) {G1.a(i,j, GGj(mu1[j],i,j,1));
                               H1.a(i,j, GGj(mu1[j],i,j,2)
                                           +G1.a(i,j)*aamm[i]);}  else
    if(mumu(mu2[i],mu1[j],i)) {G1.a(i,j,-GGj(mu1[j],i,j,1));
                               H1.a(i,j,-GGj(mu1[j],i,j,2));}

    if(mumu(mu2[i],mu2[j],i)) {G2.a(i,j, GGj(mu2[j],i,j,1));
                               H2.a(i,j, GGj(mu2[j],i,j,2));}   else
    if(mumu(mu1[i],mu2[j],i)) {G2.a(i,j,-GGj(mu2[j],i,j,1));
                               H2.a(i,j,-GGj(mu2[j],i,j,2)
                                           +G2.a(i,j)*aamm[i]);}
                                           // ^^^^... delta-function barriers.
    if(i==j) {H1.add(i,j,-meff);  H2.add(i,j,meff);}  // H1 and H2 for i=j.

  } else G1.a(i,j,GGj(0,i,j,1));

    G1_1=1/G1;  G1.free();  if(finite_barriers) {
                            S1=H1*G1_1;  H1.free();
    G2_1=1/G2;  G2.free();  S2=H2*G2_1;  H2.free();
    D_1=1/(S1-S2);
} }

void solve(void)
{
  matrix eta1;
  if(finite_barriers) {
    eta1=D_1*(psip-S2*psi);  sigma1=G1_1*eta1;  sigma2=G2_1*(eta1-psi);
  } else sigma1=-(G1_1*psi);
}

complex reflected_wave(numero x, numero y)
{
  int i, mu=where_is_the_point(x,y);
  complex val=0, ss;

  if(finite_barriers || mu==0)         // medium mu=0 is the wavefunction
  for(i=0; i<nth; i++) {               // domain for infinite barriers
    if(finite_barriers==0)  ss=sigma1.a(i);  else
    if(mu1[i]==mu)          ss=sigma1.a(i);  else
    if(mu2[i]==mu)          ss=sigma2.a(i);
    if(mu1[i]==mu || mu2[i]==mu)  val+=ns[i]*coef[i]*GG(mu,i,x,y,1)*ss;
//    printf("%d %g %g %g %g %d %d %d\n", i, x, y, mod(sigma1.a(i)), mod(sigma2.a(i)), mu, mu1[i], mu2[i]);
  }

  return val;
}

// **************************************************************************
// *** sources                                                            ***
// **************************************************************************

int outgoing_source(numero x0, numero y0)
{
  int i,s, mu=where_is_the_point(x0,y0);

  psi.alloc(nth,1);  if(finite_barriers) psip.alloc(nth,1);

  if(finite_barriers || mu==0)
  for(i=0; i<nth; i++)  if(mu==mu1[i] || mu==mu2[i]) {
    if(mu==mu1[i] && finite_barriers) s=-1; else s=1;
    psi.a(i, s*GG(mu,i,x0,y0,1));
    if(finite_barriers) {
      psip.a(i, s*GG(mu,i,x0,y0,2));
      if(mu==mu1[i]) psip.add(i,aamm[i]*psi.a(i));  // s=-1 in this case
  } }

  solve();  return mu;
}

complex Gaussian_beam(int mu, numero x, numero y,
                      numero x0, numero y0, numero phi, numero D, int n,
                      numero nx, numero ny, int index)
{
  if(mu<mumin || mumax<mu)
    on_error("EBEM: Gaussian beam", "region", mu, "is out of range");

  numero kk=real(kmu[mu]);  // Set index=1 or 2 for G or dG/dn, respectively.
  if(ABS(imag(kmu[mu]))>1e-3*kk)
    on_error("EBEM: Gaussian beam", "cannot be defined in lossy region");

  numero cfi=cos(phi),sfi=sin(phi);
  x=x-x0;  y=y-y0;                     // (x,y) is now relative to focus.
  numero xp=x*cfi+y*sfi, yp=-x*sfi+y*cfi, val;          // Beam || x now.

  numero alpha=sqr(D/4);              // D is the focal width in amplitude.
  numero fct,Qx,Qy,Qy2;  complex psi;  int i;

  for(i=-n, psi=0, val=0; i<=n; i++) {
    if(n==0) Qy=0; else Qy=i*kk/n;  Qy2=Qy*Qy;
    Qx=kk*kk-Qy2;  if(Qx<0) Qx=0; else Qx=sqrt(Qx);
    fct=exp(-alpha*Qy2);
    if(index==1) psi+=fct*euler(1,Qx*xp+Qy*yp);
    else         psi+=fct*euler(1,Qx*xp+Qy*yp)
                      * i_c*((Qx*cfi-Qy*sfi)*nx+(Qx*sfi+Qy*cfi)*ny);
    val+=fct;
  }

  return psi/val;
}

int Gaussian_beam_source(numero x0, numero y0, numero phi, numero D,
                         int n, int mu)
{
  int i,s;

  psi.alloc(nth,1);  if(finite_barriers) psip.alloc(nth,1);

  if(finite_barriers || mu==0)
  for(i=0; i<nth; i++)  if(mu==mu1[i] || mu==mu2[i]) {
    if(mu==mu1[i] && finite_barriers) s=-1; else s=1;
    psi.a(i, s*Gaussian_beam(mu,xs[i],ys[i],x0,y0,phi,D,n,0,0,1));
    if(finite_barriers) {
      psip.a(i, s*Gaussian_beam(mu,xs[i],ys[i],x0,y0,phi,D,n,nxs[i],nys[i],2));
      if(mu==mu1[i]) psip.add(i,aamm[i]*psi.a(i));  // s=-1 in this case
  } }

  solve();
}

// **************************************************************************
// *** LDOS                                                               ***
// **************************************************************************

numero LDOS(numero x0, numero y0)
{
  int     i,j,mu;
  complex val=0,E;

  mu=outgoing_source(x0,y0);  E=sqr(kmu[mu])/(2*meff);
  val=reflected_wave(x0,y0);

  return  meff/(2*pi*pi)*(pi/2+atan(real(E)/ABS(imag(E)))) - (1/pi)*imag(val);
}

// **************************************************************************
// *** Photoemission (PE)                                                 ***
// **************************************************************************

complex initial_wave(int mu, numero x, numero y, int i)
{
  if(mu<0) on_error("EBEM: initial_wave", "internal error");
  complex PE_coef=4*i_c/(kmu[mu]*kmu[mu]-Qx*Qx-Qy*Qy+i_c*1e-20)*sqr(meff/2);

  if(i>=0) PE_coef=PE_coef*i_c*(Qx*nxs[i]+Qy*nys[i]);

  return  PE_coef * euler(1,Qx*x+Qy*y);
}

complex initial_wave(numero x, numero y)
  {return initial_wave(where_is_the_point(x,y),x,y,-1);}

void init_initial_wave(numero Qx_, numero Qy_)
{
  complex psi1;
  int i;  Qx=Qx_;  Qy=Qy_;

  psi.alloc(nth,1);  if(finite_barriers) psip.alloc(nth,1);

  for(i=0; i<nth; i++) if(finite_barriers) {
    psi1=initial_wave(mu1[i],xs[i],ys[i],-1);
    psi.a(i, initial_wave(mu2[i],xs[i],ys[i],-1) - psi1);
    psip.a(i, initial_wave(mu2[i],xs[i],ys[i], i)
             -initial_wave(mu1[i],xs[i],ys[i], i) - aamm[i]*psi1);
  } else if(mu1[i]==0 || mu2[i]==0) psi.a(i, initial_wave(0,xs[i],ys[i],-1));

  solve();
}

numero PE(numero Qx, numero Qy, complex E)
{
  init_initial_wave(Qx,Qy);

  int i,j;
  numero x,y;
  complex val=0;

  for(i=0; i<period_n1; i++)
  for(j=0; j<period_n2; j++) {
    x=period_x0+(i+0.5)*period_a1x/period_n1+(j+0.5)*period_a2x/period_n2;
    y=period_y0+(i+0.5)*period_a1y/period_n1+(j+0.5)*period_a2y/period_n2;
    val+=(initial_wave(x,y)+reflected_wave(x,y)) * euler(1,-Qx*x-Qy*y)
         *area/(period_n1*period_n2);
  }

  return mod2(val);
}

// **************************************************************************
// *** Electron plane wave (EPW) calculations in 2D                       ***
// **************************************************************************

class epw_def {  // This object is self-contained, except for matrix.h,
public:          // lattice.h, bessel.h, complex.h, jga.h, and eispack.h

  int     obj;         // Object type: -1/0/1 for arbitrary/circle/triangle.
  int     lat;         // Lattice type: -2/-1/0/1, with lat=0/1 for
                       // hexagonal/square, and lat<0 for arb. lattice.
  // Common parameters.

  numero  meff;        // Effective electron mass.
  numero  gmax;        // Radius of accepted 2D G vectors in reciprocal space.

  // Only for obj=0 or 1.

  numero  a,R;         // Lattice constant and object size (radius, side).
  complex U1,U0;       // Potential inside (U1) and outside (U0) object.

  // Only for ojb=-1, 0, or 1.
  int     nQ;          // Number of sampling Q points.

  // Only for obj=-1.
  numero  Q0x,Q0y,Q1x,Q1y;   // Segment in Q-space for lat=-1.

  // Only for obj=-2.
  numero  dQx,dQy;     // Sides of rectangular grid in Q for lat=-2.
  int     nx,ny;       // Number of subdivisions in rectangular grid in Q.

  // Only for bands (opt=0 in get_band).
  numero  Emax;        // Maximum energy under consideration for bands.

  // Only for photoemission (opt=1 in get_band).
  numero  Er,Ei;       // Photoemission: Re{E} and Im{E}.
  numero  E0;          // Photoemission: initial energy.
  numero  E1;          // Photoemission: final energy.
  int     nE;          // Photoemission: number of energies.

  // Only for DOS and LDOS (opt=3-5 in get_band).
  numero x0,y0,x1,y1,dx,dy;   // nx,ny will define a grid in x-y in this case.
  int    nr;                  // Number of points in real-space segment.
  int    n1,n2;               // Grid in momentum space.
  int    nDOS;                // Number of DOS values.
  numero *DOS;                // Actual values of DOS or LDOS.

  // Internal variables.
  matrix  mU,mE,mfi;   // Potential matrix, energy eigenvalues, and eigenfunct.
  int     set_lattice; // Set to 1 to define triang. or square lattice in EPW.

  // Internal routines.
  void solve(numero Qx_, numero Qy_);

  complex eexp(numero G, numero d) {
    numero Gd=G*d;  if(ABS(Gd)>=1e-8) return (1-euler(1,Gd))/G;
    return -i_c*d;
  }

  complex wf(int i, numero x, numero y)
  {
    complex val=0;  int j;
    for(j=0; j<n_g; j++)  val+=mfi.a(j,i)*euler(1,Qgx[j]*x+Qgy[j]*y);
    return val;    
  }

  numero PE(void) {                               // Photoemission probability.
    int i;  complex val=0;
    for(i=0; i<n_g; i++) val+=mod2(mfi.a(0,i))/(complex(Er,Ei)-mE.a(i));
    return (-1/pi)*imag(val);
  }

  numero TDOS(void) {                             // Total DOS.
    int i;  complex val=0;
    for(i=0; i<n_g; i++) val+=1/(complex(Er,Ei)-mE.a(i));
    return (-1/pi)*imag(val);
  }

  numero LDOS(numero x, numero y) {               // LDOS.
    int i;  complex val=0;
    for(i=0; i<n_g; i++) val+=mod2(wf(i,x,y))/(complex(Er,Ei)-mE.a(i));
    return (-1/pi)*imag(val);
  }

  // Main external routine.

  void get_band(char *name, int opt);

  // Default values.

  epw_def(void) {
    obj=0;    // Circles.
    lat=0;    // Hexagonal lattice.
    meff=1;
    Emax=infinite;  set_lattice=0;
    Q0x=Q0y=0;  nQ=1;
} };

void epw_def::solve(numero Qx_, numero Qy_)
{
  matrix ma;  int i;        // Secular matrix and solution of the eigensystem.

  Qx=Qx_;  Qy=Qy_;  lattice(gmax);

  ma=mU;
  for(i=0; i<n_g; i++) ma.add(i,i,(sqr(Qgx[i])+sqr(Qgy[i]))/(2*meff));
  eispack(ma,mE,mfi);       // mE is an array with eigen-energies.
  eispack_order(mE,mfi,8);  // The columns of mfi are the corresponding
}              // normalized eigenfunctions: <Qj|Q'j'>=(2*pi)^2 delta(Q-Q').

void epw_def::get_band(char *name, int opt)
{                               // 'name' is the output file.
  complex val;                  // opt=0 -> band structure
  numero  x,y,Gx,Gy,G,QQ,JJ,YY; // opt=1 -> photoemission (PE) q-E map.
  int     l,i,j,ii,jj,mu;       // opt=2 -> photoemission qx-qy map.
  FILE    *fout;                // opt=3 -> LDOS rE map
                                // opt=4 -> LDOS xy map
  if(set_lattice)               // opt=5 -> DOS
  if(lat==0) {a0=b0=a;  alpha=-pi/6;  beta=pi/3;}  // Hexagonal lattice.
  if(lat==1) {a0=b0=a;  alpha=0;      beta=pi/2;}  // Square lattice.

  if(lat==0) printf("Hexagonal lattice");  else
  if(lat==1) printf("Square lattice");      else
  if(lat<0)  printf("Manually-defined lattice");
  if(obj==0) printf(" of circles");         else
  if(obj==1) printf(" of triangles");       else
  if(obj<0)  printf(" of manually-defined cells");
  lattice(gmax);  printf(" described by %d plane waves.\n", n_g);

  // --- potential matrix

  mU.alloc(n_g,n_g);

  if(obj<0) {                                      // Manually-defined cell.
    #ifdef BEM_HELP
    if(finite_barriers==0)  on_error("EBEM: EPW",
      "periodic systems cannot be used in combination with infinite barriers");
    init_geometry();
    for(i=0; i<period_n1; i++)
    for(j=0; j<period_n2; j++) {
      x=period_x0+(i+0.5)*period_a1x/period_n1+(j+0.5)*period_a2x/period_n2;
      y=period_y0+(i+0.5)*period_a1y/period_n1+(j+0.5)*period_a2y/period_n2;
      mu=where_is_the_point(x,y);
      if(mu>=0)  val=Vmu[mu]/(period_n1*period_n2);  else val=0;
      for(ii=0; ii<n_g; ii++)
      for(jj=0; jj<n_g; jj++)
        mU.add(ii,jj, val*euler(1,(Qgx[ii]-Qgx[jj])*x+(Qgy[ii]-Qgy[jj])*y));
    }
    #endif

  } else                                           // Predefined shapes.
  for(i=0; i<n_g; i++)
  for(j=0; j<n_g; j++) {
    if(i==j) {
      if(obj==0)  mU.a(i,i,U0+(U1-U0)/area*pi*R*R);            // Circle.
      else        mU.a(i,i,U0+(U1-U0)/area*sqrt(3)/4*R*R);     // Triangle.
    } else {
      Gx=Qgx[i]-Qgx[j];  Gy=Qgy[i]-Qgy[j];  G=sqrt(Gx*Gx+Gy*Gy);
      if(obj==0) {
         besselJY(1,G*R,JJ,YY,x,y);  // old -> val=2*pi/(G*G)*JJ;
         val=2*pi/(G*G)*G*R*JJ;                    // Circle.
      } else {                                     // Triangle with one side
        if(ABS(Gy*a0)<1e-8) Gy=1e-8/a0;            // parallel to the y axis.
        val=( eexp(Gx+Gy/sqrt(3),sqrt(3)*R/2)
             -eexp(Gx-Gy/sqrt(3),sqrt(3)*R/2))/Gy;
      }
      val=(U1-U0)/area*val;  mU.a(i,j,val);
  } }

  // --- 1st Brillouin zone sampling, etc.

  if(opt==3) nDOS=nE*nr; else
  if(opt==4) nDOS=nx*ny; else
  if(opt==5) nDOS=nE;    else  nDOS=0;   if(nDOS>0) DOS=new numero [nDOS];

  for(i=0; i<nDOS; i++) DOS[i]=0;

  fout=fopen(name,"w");  fclose(fout);             // Clear output file.
  if(lat==-2) nQ=nx*ny;
  if(lat==-3) nQ=n1*n2;

  for(l=0; l<nQ; l++) {

    if(l%10==0 && l>0)
      printf("Momentum sampling: first %d Q points completed.\n", l);

    if(lat==-1) {                                  // Manually-defined lattice:
        if(nQ<=1) Qx=Q0x; else Qx=Q0x+l*(Q1x-Q0x)/(nQ-1);  // segment in
        if(nQ<=1) Qy=Q0y; else Qy=Q0y+l*(Q1y-Q0y)/(nQ-1);  // Q space.
        QQ=sqrt(sqr(Qx-Q0x)+sqr(Qy-Q0y));
    } else
    if(lat==-2) {                                  // Rectangular grid in Q.
        if(nx==1) Qx=Q0x; else Qx=Q0x+(l/ny)*dQx/(nx-1);
        if(ny==1) Qy=Q0y; else Qy=Q0y+(l%ny)*dQy/(ny-1);
    } else
    if(lat==-3) {                                  // Grid over first Q cell.
        i=l/n2;  j=l%n2;
        Qx=(i+0.5)*Gax/n1+(j+0.5)*Gbx/n2;
        Qy=(i+0.5)*Gay/n1+(j+0.5)*Gby/n2;
    } else
    if(excursion_1BZ(l,nQ,QQ)==0)                  // Assign (Qx,Qy) in 1BZ
      on_error("EBEM: EPW-hexagonal/square",       // if lattice is triang. or
               "the lattice is neither hexagonal nor square");  // square

    // Find eigenenergies and eigenstates of periodic system (with EPW)
    // for wave vector (Qx,Qy). The sign of Qx and/or Qy can be flipped
    // by hand in this call to explore other orientations of the 1BZ
    // excursion, which can be useful for non-symmetric unit cells.
    #define exc_case 1
    if(exc_case==1) solve( Qx, Qy);
    if(exc_case==2) solve( Qx,-Qy);
    if(exc_case==3) solve(-Qx, Qy);
    if(exc_case==4) solve(-Qx,-Qy);

    if(opt<3 || 5<opt) {

      fout=fopen(name,"a");
      if(opt==0) {                                 // Band structure.
        for(i=0; i<n_g; i++) if(mod(mE.a(i))<Emax)
          fprintf(fout, "%g %g %g\n", QQ, real(mE.a(i)), imag(mE.a(i)));
      }
      if(opt==1 && lat!=-2) {                      // Photoemission q-E map.
        for(j=0; j<nE; j++) {
          if(nE<=1) Er=E0; else Er=E0+j*(E1-E0)/(nE-1);
          fprintf(fout, "%g %g %g\n", QQ, Er, PE());
      } }
      if(opt==2 && lat==-2)                        // Photoemission qx-qx map.
        fprintf(fout, "%g %g %g\n", Qx, Qy, PE());
      fclose(fout);

    } else

    if(opt==3) {                                   // Calculation for LDOS.
      for(i=ii=0; i<nE; i++)
      for(j=0; j<nr; j++,ii++) {
        if(nE==1) Er=E0; else Er=E0+i*(E1-E0)/(nE-1);
        if(nr==1) {x=x0; y=y0;} else {x=x0+j*(x1-x0)/(nr-1);
                                      y=y0+j*(y1-y0)/(nr-1);}
        DOS[ii]+=LDOS(x,y);
    } } else

    if(opt==4) {                                   // Calculation for LDOS.
      for(i=ii=0; i<nx; i++)
      for(j=0; j<ny; j++,ii++) {
        if(nx==1) x=x0; else x=x0+i*dx/(nx-1);
        if(ny==1) y=y0; else y=y0+j*dy/(ny-1);
        DOS[ii]+=LDOS(x,y);
    } } else

    if(opt==5) {                                   // Calculation for TDOS.
      for(i=0; i<nE; i++) {
        if(nE==1) Er=E0; else Er=E0+i*(E1-E0)/(nE-1);
        DOS[i]+=TDOS();
  } } }

  if(nDOS>0) {                                     // Print out of DOS.
    for(i=0; i<nDOS; i++) DOS[i]=DOS[i]/(area*n1*n2);   // A^*/(2*pi)^2 = 1/A

    fout=fopen(name,"a");

    if(opt==3) {
      for(i=ii=0; i<nE; i++)
      for(j=0; j<nr; j++,ii++) {
        if(nE==1) Er=E0; else Er=E0+i*(E1-E0)/(nE-1);
        if(nr==1) {x=x0; y=y0;} else {x=x0+j*(x1-x0)/(nr-1);
                                      y=y0+j*(y1-y0)/(nr-1);}
        fprintf(fout, "%g %g %g\n", sqrt(sqr(x-x0)+sqr(y-y0)), Er,
                DOS[ii]);
    } } else

    if(opt==4) {
      for(i=ii=0; i<nx; i++)
      for(j=0; j<ny; j++,ii++) {
        if(nx==1) x=x0; else x=x0+i*dx/(nx-1);
        if(ny==1) y=y0; else y=y0+j*dy/(ny-1);
        fprintf(fout, "%g %g %g\n", x, y, DOS[ii]);
    } } else

    if(opt==5) {
      for(i=0; i<nE; i++) {
        if(nE==1) Er=E0; else Er=E0+i*(E1-E0)/(nE-1);
        fprintf(fout, "%g %g\n", Er, DOS[i]);
    } }

    fclose(fout);
    delete [] DOS;
  }

  mU.free();  mE.free();  mfi.free();
}

// **************************************************************************
// *** Main routine: input-file interpreter                               ***
// **************************************************************************

void program(char *fin_name)
{
  char name[200], command[100];  // File names and commands.
  FILE *fin, *fout;              // Input, output files.
  epw_def epw;                   // EPW calculations.

// --- Try to open input file.

  strcpy(name,fin_name);
  if((fin=fopen(name,"r"))==NULL) {          // Try to open file name, and
    strcat(name,".bem");                     // if this is not available,
    if((fin=fopen(name,"r"))==NULL) {        // open name.bem, ... or exit.
      printf("*** error in 'EBEM': files %s and %s.bem cannot be opened\n",
             fin_name, fin_name);
      exit(1);
  } }

// --- A few internal variables.

  int i,j,n,m,nr,nE,nx,ny,ll,mu1,mu2,mu;
  numero x,y,x0,y0,x1,y1,dx,dy, Q0x,Q0y,Q1x,Q1y,dQx,dQy, E0,E1,Er,Ei=0;
  numero a,b,al1,be1,Rmax, th0,th1,galph,qq, xr,yr,phi, xa,ya, D;
  complex wf;

// --- Input-file interpreter.

  do {
    read_name(fin,command);

// --- definitions based upon = sign, with no blank spaces in between

    if((i=algebra_locate(command,'='))>0) {   // = sign at position i.
      command[i]=0;
      algebra_define(command,command+i+1);
    } else

/// --- add-polygon-boundary: add a polygon to the boundary geometry

    if(!strcmpC(command,"add-polygon-boundary")) {
      ll=alread_int(fin);  n_fragments+=ll;
      mu1=alread_int(fin);    mu2=alread_int(fin);
      x1=alread_numero(fin);  y1=alread_numero(fin);
      for(i=0; i<ll; i++) {
        x0=x1;  y0=y1;  x1=alread_numero(fin);  y1=alread_numero(fin);        
        n=alread_int(fin);
        add_segment(n,mu1,mu2,x0,y0,x1,y1);
    } } else

// --- add-arc-boundary: add an arc of an ellipse to the boundary geometry

    if(!strcmpC(command,"add-arc-boundary")) {
      n_fragments++;
      mu1=alread_int(fin);     mu2=alread_int(fin);
      x0=alread_numero(fin);   y0=alread_numero(fin);
      a=alread_numero(fin);    b=alread_numero(fin);
      galph=alread_numero(fin)*pi/180;
      th0=alread_numero(fin)*pi/180;  th1=alread_numero(fin)*pi/180;
      n=alread_int(fin);
      add_arc(n,mu1,mu2,x0,y0,a,b,th0,th1,galph);
    } else

// --- duplicate-geometry: repeat geometry displaced by a vector

    if(!strcmpC(command,"duplicate-geometry"))  {
      dx=alread_numero(fin);  dy=alread_numero(fin);
      copy_geometry(dx,dy);  n_fragments=2*n_fragments;
    } else

// --- rotate-geometry: rotate the geometry by an angle

    if(!strcmpC(command,"rotate-geometry"))  {
      xr=alread_numero(fin);  yr=alread_numero(fin);
      phi=alread_numero(fin)*pi/180;
      rotate_geometry(xr,yr,phi);
    } else

// --- translate-geometry: translate the geometry by (xr,yr)

    if(!strcmpC(command,"translate-geometry"))  {
      xr=alread_numero(fin);  yr=alread_numero(fin);
      translate_geometry(xr,yr);
    } else

// --- repeat-periodic-geometry: repeat geometry in a lattice

    if(!strcmpC(command,"repeat-periodic-geometry"))  {
      a=alread_numero(fin);  b=alread_numero(fin);
      al1=alread_numero(fin)*pi/180;
      be1=alread_numero(fin)*pi/180;
      Rmax=alread_numero(fin);
      n=int((Rmax/a+Rmax/b)*(1+1e-8));  ll=0;
      for(i=-n; i<n; i++)
      for(j=-n; j<n; j++) if(i!=0 || j!=0) {
        dx=i*a*cos(al1)+j*b*cos(al1+be1);
        dy=i*a*sin(al1)+j*b*sin(al1+be1);
        if(dx*dx+dy*dy<=Rmax*Rmax*(1+1e-8)) {
          ll++;  copy_geometry(dx,dy);
        }
      }
      n_fragments=(ll+1)*n_fragments;
    } else

// --- repeat-periodic-geometry-1D: repeat geometry in a 1D lattice

    if(!strcmpC(command,"repeat-periodic-geometry-1D"))  {
      a=alread_numero(fin);  phi=alread_numero(fin)*pi/180;
      n=alread_int(fin);  m=alread_int(fin);  dx=a*cos(phi);  dy=a*sin(phi);
      for(i=-n; i<=m; i++) if(i!=0)  copy_geometry(i*dx,i*dy);
      n_fragments=(n+m+1)*n_fragments;
    } else

// --- periodicity: make the system periodic

    if(!strcmpC(command,"periodicity"))  {
      init_periodicity(fin);  period_flag=1;
    } else

// --- begin-virtual-structure: start virtual structure definition

    if(!strcmpC(command,"begin-virtual-structure"))
      end_real_structure();  else

/// --- write-geometry: print out of the boundary geometry

    if(!strcmpC(command,"write-geometry")) {
      read_name(fin,name);  init_geometry();  write_geometry(name);
    } else

// --- write-region: print out points inside a given region

    if(!strcmpC(command,"write-region")) {
      init_geometry();  mu=alread_int(fin);  ll=alread_int(fin);
      read_name(fin,name);  write_region(mu,ll,name);
    } else

// --- write-regions: print out grid with the whole region structure

    if(!strcmpC(command,"write-regions")) {
      init_geometry();  ll=alread_int(fin);
      read_name(fin,name);  write_region(-1,ll,name);
    } else

// --- clear-all: recover default values and discard geometry

    if(!strcmpC(command,"clear-all"))  clear_all();  else

// --- Ei: set the value of Im{E} (0 by default)

    if(!strcmpC(command,"Ei")) {
      Ei=alread_numero(fin);
      if(Ei<0)  on_error("EBEM: Ei", "Ei cannot be negative");
    } else

// --- meff: set the value of the effective electron mass (1 by default)

    if(!strcmpC(command,"meff"))  meff=alread_numero(fin);  else

// --- potential: assign potential of a region (0 by default)

    if(!strcmpC(command,"potential")) {
      i=alread_int(fin);
      a=alread_numero(fin);  b=alread_numero(fin);
      if(b>0)  on_error("EBEM: potential", "Im{V} cannot be positive");
      Vmu[i]=complex(a,b);
    } else

// --- barrier: assign a complex coefficient to a delta-function barrier

    if(!strcmpC(command,"barrier")) {
      i=alread_int(fin);     j=alread_int(fin);
      a=alread_numero(fin);  b=alread_numero(fin);
      aa[i][j]=aa[j][i]=complex(a,b);
    } else

// --- infinite-barriers: infinite-barriers in all boundaries

    if(!strcmpC(command,"infinite-barriers"))  finite_barriers=0;
    else

// --- LDOS-rE-BEM: spectrum for a range of E's along a segment using EBEM

    if(!strcmpC(command,"LDOS-rE-BEM")) {
      x0=alread_numero(fin);  y0=alread_numero(fin);
      x1=alread_numero(fin);  y1=alread_numero(fin);
      nr=alread_int(fin);
      E0=alread_numero(fin);  E1=alread_numero(fin);
      nE=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      for(i=0; i<nE; i++) {
        if(nE==1) Er=E0; else Er=E0+i*(E1-E0)/(nE-1);
        init(complex(Er,Ei));
        for(j=0; j<nr; j++) {
          if(nr==1) {x=x0; y=y0;} else {x=x0+j*(x1-x0)/(nr-1);
                                        y=y0+j*(y1-y0)/(nr-1);}
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g\n", sqrt(sqr(x-x0)+sqr(y-y0)), Er,
                  LDOS(x,y));
          fclose(fout);
	}
    } } else

// --- LDOS-xy-BEM: LDOS map for a single E value using EBEM

    if(!strcmpC(command,"LDOS-xy-BEM")) {
      x0=alread_numero(fin);  y0=alread_numero(fin);
      dx=alread_numero(fin);  dy=alread_numero(fin);
      Er=alread_numero(fin);
      nx=alread_int(fin);
      ny=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      init(complex(Er,Ei));
      for(i=0; i<nx; i++) {
        if(nx==1) x=x0; else x=x0+i*dx/(nx-1);
        for(j=0; j<ny; j++) {
          if(ny==1) y=y0; else y=y0+j*dy/(ny-1);
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g\n", x, y, LDOS(x,y));
          fclose(fout);
	}
    } } else

// --- LDOS-xyE-BEM: LDOS as function of (E,x,y) using EBEM

    if(!strcmpC(command,"LDOS-xyE-BEM")) {
      x0=alread_numero(fin);  y0=alread_numero(fin);
      dx=alread_numero(fin);  dy=alread_numero(fin);
      nx=alread_int(fin);
      ny=alread_int(fin);
      E0=alread_numero(fin);  E1=alread_numero(fin);
      nE=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      for(ll=0; ll<nE; ll++) {
        if(nE==1) Er=E0; else Er=E0+ll*(E1-E0)/(nE-1);
        init(complex(Er,Ei));
        for(i=0; i<nx; i++) {
          if(nx==1) x=x0; else x=x0+i*dx/(nx-1);
          for(j=0; j<ny; j++) {
            if(ny==1) y=y0; else y=y0+j*dy/(ny-1);
            fout=fopen(name,"a");
            fprintf(fout, "%g %g %g %g\n", Er, x, y, LDOS(x,y));
            fclose(fout);
	  }
    } } } else

// --- LDOS-rE-EPW: LDOS for a range of E's along a segment using EPW

    if(!strcmpC(command,"LDOS-rE-EPW")) {
      if(period_flag==0)  on_error("EBEM: LDOS-rE-EPW",
                                   "this applies only to periodic systems");
      epw.gmax=alread_numero(fin);
      epw.x0=alread_numero(fin);  epw.y0=alread_numero(fin);
      epw.x1=alread_numero(fin);  epw.y1=alread_numero(fin);
      epw.nr=alread_int(fin);
      epw.E0=alread_numero(fin);  epw.E1=alread_numero(fin);
      epw.nE=alread_int(fin);     epw.Ei=Ei;
      epw.meff=meff;            epw.n1=period_n1;  epw.n2=period_n2;
      epw.obj=-1;             // Manually-defined cell.
      epw.lat=-3;             // Manually-defined lattice with full 1BZ cover.
      read_name(fin,name);  epw.get_band(name,3);   // Calculate LDOS.
    } else

// --- LDOS-xy-EPW: LDOS map for a single E value using EPW

    if(!strcmpC(command,"LDOS-xy-EPW")) {
      if(period_flag==0)  on_error("EBEM: LDOS-xy-EPW",
                                   "this applies only to periodic systems");
      epw.gmax=alread_numero(fin);
      epw.x0=alread_numero(fin);  epw.y0=alread_numero(fin);
      epw.dx=alread_numero(fin);  epw.dy=alread_numero(fin);
      epw.Er=alread_numero(fin);  epw.Ei=Ei;
      epw.nx=alread_int(fin);
      epw.ny=alread_int(fin);
      epw.meff=meff;            epw.n1=period_n1;  epw.n2=period_n2;
      epw.obj=-1;             // Manually-defined cell.
      epw.lat=-3;             // Manually-defined lattice with full 1BZ cover.
      read_name(fin,name);  epw.get_band(name,4);   // Calculate LDOS.
    } else

// --- DOS-EPW: DOS for a range of E's using EPW

    if(!strcmpC(command,"DOS-EPW")) {
      if(period_flag==0)  on_error("EBEM: LDOS-rE-EPW",
                                   "this applies only to periodic systems");
      epw.gmax=alread_numero(fin);
      epw.E0=alread_numero(fin);  epw.E1=alread_numero(fin);
      epw.nE=alread_int(fin);     epw.Ei=Ei;
      epw.meff=meff;            epw.n1=period_n1;  epw.n2=period_n2;
      epw.obj=-1;             // Manually-defined cell.
      epw.lat=-3;             // Manually-defined lattice with full 1BZ cover.
      read_name(fin,name);  epw.get_band(name,5);   // Calculate LDOS.
    } else

// --- PE-area-BEM: define a cell for PE integration in a non-periodic system

    if(!strcmpC(command,"PE-area-BEM"))  init_periodicity(fin);
    else

// --- PE-qE-BEM: PE spectrum for a range of q's and E's

    if(!strcmpC(command,"PE-qE-BEM")) {
      Q0x=alread_numero(fin);  Q0y=alread_numero(fin);
      Q1x=alread_numero(fin);  Q1y=alread_numero(fin);
      nr=alread_int(fin);
      E0=alread_numero(fin);  E1=alread_numero(fin);
      nE=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      for(i=0; i<nE; i++) {
        if(nE==1) Er=E0; else Er=E0+i*(E1-E0)/(nE-1);
        init(complex(Er,Ei));
        for(j=0; j<nr; j++) {
          if(nr==1) {Qx=Q0x; Qy=Q0y;} else {Qx=Q0x+j*(Q1x-Q0x)/(nr-1);
                                            Qy=Q0y+j*(Q1y-Q0y)/(nr-1);}
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g\n", sqrt(sqr(Qx-Q0x)+sqr(Qy-Q0y)), Er,
                  PE(Qx,Qy,complex(Er,Ei)));
          fclose(fout);
	}
    } } else

// --- PE-qxqy-BEM: PE map as a function of (Qx,Qy) for fixed E

    if(!strcmpC(command,"PE-qxqy-BEM")) {
      Q0x=alread_numero(fin);  Q0y=alread_numero(fin);
      dQx=alread_numero(fin);  dQy=alread_numero(fin);
      Er=alread_numero(fin);
      nx=alread_int(fin);
      ny=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      init(complex(Er,Ei));
      for(i=0; i<nx; i++) {
        if(nx==1) Qx=Q0x; else Qx=Q0x+i*dQx/(nx-1);
        for(j=0; j<ny; j++) {
          if(ny==1) Qy=Q0y; else Qy=Q0y+j*dQy/(ny-1);
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g\n", Qx, Qy, PE(Qx,Qy,complex(Er,Ei)));
          fclose(fout);
	}
    } } else

// --- WF-xy-1: wave function for a point source

    if(!strcmpC(command,"WF-xy-1")) {
      xa=alread_numero(fin);  ya=alread_numero(fin);
      x0=alread_numero(fin);  y0=alread_numero(fin);
      dx=alread_numero(fin);  dy=alread_numero(fin);
      Er=alread_numero(fin);
      nx=alread_int(fin);
      ny=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      init(complex(Er,Ei));
      mu=outgoing_source(xa,ya);
      for(i=0; i<nx; i++) {
        if(nx==1) x=x0; else x=x0+i*dx/(nx-1);
        for(j=0; j<ny; j++) {
          if(ny==1) y=y0; else y=y0+j*dy/(ny-1);
          wf=reflected_wave(x,y);
          if(mu==where_is_the_point(x,y))  wf+=GG(mu,x,y,0,0,xa,ya,1);
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g %g\n", x, y, mod2(wf), arg(wf));
          fclose(fout);
	}
    } } else

// --- WF-xy-2: wave function for a Gaussian beam

    if(!strcmpC(command,"WF-xy-2")) {
      xa=alread_numero(fin);  ya=alread_numero(fin);
      phi=alread_numero(fin)*pi/180;  D=alread_numero(fin);
      n=alread_int(fin);  mu=alread_int(fin);
      x0=alread_numero(fin);  y0=alread_numero(fin);
      dx=alread_numero(fin);  dy=alread_numero(fin);
      Er=alread_numero(fin);
      nx=alread_int(fin);
      ny=alread_int(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      init(complex(Er,Ei));
      Gaussian_beam_source(xa,ya,phi,D,n,mu);
      for(i=0; i<nx; i++) {
        if(nx==1) x=x0; else x=x0+i*dx/(nx-1);
        for(j=0; j<ny; j++) {
          if(ny==1) y=y0; else y=y0+j*dy/(ny-1);
          wf=reflected_wave(x,y);
          if(mu==where_is_the_point(x,y))
            wf+=Gaussian_beam(mu,x,y,xa,ya,phi,D,n,0,0,1);
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g %g\n", x, y, mod2(wf), arg(wf));
          fclose(fout);
    } } } else

// --- WF-bound: wave function at a boundary for a point source

    if(!strcmpC(command,"WF-bound")) {
      xa=alread_numero(fin);  ya=alread_numero(fin);
      Er=alread_numero(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      init(complex(Er,Ei));
      mu=outgoing_source(xa,ya);
      for(i=0; i<nth; i++) {
        wf=psi.a(i,0);
        fout=fopen(name,"a");
        fprintf(fout, "%g %g %g %g %g %g %g\n",
                xs[i], ys[i], coef[i]*ns[i],
                nxs[i], nys[i], mod2(wf), arg(wf));
        fclose(fout);
    } } else

// --- WF-rbound: wave function at a boundary for a point source

    if(!strcmpC(command,"WF-rbound")) {
      x0=alread_numero(fin);  y0=alread_numero(fin);
      x1=alread_numero(fin);  y1=alread_numero(fin);
      nr=alread_int(fin);
      Er=alread_numero(fin);
      read_name(fin,name);  fout=fopen(name,"w"); fclose(fout);
      init(complex(Er,Ei));
      for(j=0; j<nr; j++) {
        if(nr==1) {x=x0; y=y0;} else {x=x0+j*(x1-x0)/(nr-1);
                                      y=y0+j*(y1-y0)/(nr-1);}
        mu=outgoing_source(x,y);
        for(i=0; i<nth; i++) {
          wf=psi.a(i,0);
          fout=fopen(name,"a");
          fprintf(fout, "%g %g %g %g %g %g %g %g\n", sqrt(sqr(x-x0)+sqr(y-y0)),
                  xs[i], ys[i], coef[i]*ns[i],
                  nxs[i], nys[i], mod2(wf), arg(wf));
          fclose(fout);
    } } } else

// --- PE-qE-hexagonal/square-EPW: PE for triang./square lattice using EPW.

    if(!strcmpC(command,"PE-qE-hexagonal-EPW") ||
       !strcmpC(command,"PE-qE-square-EPW")) {
      if(period_flag==0)  on_error("EBEM: PE-qE-hexagonal/square-EPW",
                                   "this applies only to periodic systems");
      if(ABS(a0/b0-1)<=1e-5 && ABS(ABS(beta)-pi/3)<=1e-5)  epw.lat=0;  else
      if(ABS(a0/b0-1)<=1e-5 && ABS(ABS(beta)-pi/2)<=1e-5)  epw.lat=1;  else
        on_error("EBEM: PE-qE-hexagonal/square-EPW",
                 "the lattice must be hexagonal/square");
      epw.gmax=alread_numero(fin);
      epw.nQ=alread_int(fin);
      epw.E0=alread_numero(fin);    epw.E1=alread_numero(fin);
      epw.nE=alread_int(fin);       epw.Ei=Ei;
      epw.a=a0;
      epw.meff=meff;
      epw.obj=-1;             // Manually-defined cell.
      read_name(fin,name);  epw.get_band(name,1);   // Calculate PE.
    } else

// --- PE-qE-EPW: PE for periodic system using EPW with arbitrary lattice.

    if(!strcmpC(command,"PE-qE-EPW")) {
      if(period_flag==0)  on_error("EBEM: PE-qE-EPW",
                                   "this applies only to periodic systems");
      epw.gmax=alread_numero(fin);
      epw.Q0x=alread_numero(fin);   epw.Q0y=alread_numero(fin);
      epw.Q1x=alread_numero(fin);   epw.Q1y=alread_numero(fin);
      epw.nQ=alread_int(fin);
      epw.E0=alread_numero(fin);    epw.E1=alread_numero(fin);
      epw.nE=alread_int(fin);       epw.Ei=Ei;
      epw.meff=meff;
      epw.obj=-1;             // Manually-defined cell.
      epw.lat=-1;             // Manually-defined lattice.
      read_name(fin,name);  epw.get_band(name,1);   // Calculate PE.
    } else

// --- PE-qxqy-EPW: PE for periodic system using EPW (qx-qy map).

    if(!strcmpC(command,"PE-qxqy-EPW")) {
      if(period_flag==0)  on_error("EBEM: PE-qxqy-EPW",
                                   "this applies only to periodic systems");
      epw.gmax=alread_numero(fin);
      epw.Q0x=alread_numero(fin);   epw.Q0y=alread_numero(fin);
      epw.dQx=alread_numero(fin);   epw.dQy=alread_numero(fin);
      epw.Er=alread_numero(fin);    epw.Ei=Ei;
      epw.nx=alread_int(fin);       epw.ny=alread_int(fin);
      epw.meff=meff;
      epw.obj=-1;             // Manually-defined cell.
      epw.lat=-2;             // Manually-def. latt. with rectangular Q grid.
      read_name(fin,name);  epw.get_band(name,2);   // Calculate PE.
    } else

// --- BANDS-hexagonal/square: band structure for a hexagonal/square lattice

    if(!strcmpC(command,"BANDS-hexagonal") ||
       !strcmpC(command,"BANDS-square")) {
      if(period_flag==0)  on_error("EBEM: BANDS-hexagonal/square",
                                   "this applies only to periodic systems");
      if(ABS(a0/b0-1)<=1e-5 && ABS(ABS(beta)-pi/3)<=1e-5)  epw.lat=0;  else
      if(ABS(a0/b0-1)<=1e-5 && ABS(ABS(beta)-pi/2)<=1e-5)  epw.lat=1;  else
        on_error("EBEM: BANDS-hexagonal/square",
                 "the lattice must be hexagonal/square");
      epw.gmax=alread_numero(fin);
      epw.nQ=alread_int(fin);
      epw.Emax=alread_numero(fin);
      epw.a=a0;
      epw.meff=meff;
      epw.obj=-1;             // Manually-defined cell.
      read_name(fin,name);  epw.get_band(name,0);   // Calculate bands.
    } else

// --- BANDS: band structure for manually-defined lattice

    if(!strcmpC(command,"BANDS")) {
      if(period_flag==0)  on_error("EBEM: BANDS",
                                   "this applies only to periodic systems");
      epw.gmax=alread_numero(fin);
      epw.Q0x=alread_numero(fin);   epw.Q0y=alread_numero(fin);
      epw.Q1x=alread_numero(fin);   epw.Q1y=alread_numero(fin);
      epw.nQ=alread_int(fin);
      epw.Emax=alread_numero(fin);
      epw.meff=meff;
      epw.obj=-1;             // Manually-defined cell.
      epw.lat=-1;             // Manually-defined lattice.
      read_name(fin,name);  epw.get_band(name,0);   // Calculate bands.
    } else

// --- include: insert input file and interpret it at this point

    if(!strcmpC(command,"include") || !strcmpC(command,"insert")) {
      read_name(fin,name);
      program(name);
    } else

// --- end, exit, quit: terminate the execution of the program

    if(strcmpC(command,"end") && strcmpC(command,"exit")
       && strcmpC(command,"quit"))
      on_error("EBEM", "command not found:", command);

  } while(strcmpC(command,"end") && strcmpC(command,"exit")
          && strcmpC(command,"quit"));

  fclose(fin);
}

// --------------------------------------------------------------------------
//   Main file.
// --------------------------------------------------------------------------

int main(int argc, char **argv)
{

// --- Help on the use of the program if no parameters are given.

  if(argc<2) {
    printf("The EBEM code requires either an input file or some arguments.\n");
    printf("The correct call syntax to use an input file is as follows:\n\n");
    printf("        ebem input\n");
    printf("or      ebem input.bem\n\n");
    printf("The first one of these commands tries to open file \"input\",\n");
    printf("and if it cannot be opened, it tries to open \"input.bem\".\n");
    printf("The second one tries to open file \"input.bem\" directly.\n");
    printf("\n");
    printf("Alternatively, the code can be used to calculate one-electron\n");
    printf("2D band structures using the inline command\n\n");
    printf("        ebem epw output-file a R gmax nQ Re{U1} Re{U0}\n");
    printf("            [obj=0 lat=0 Im{U1}=0 Im{U0}=0 meff=1 Emax=inf]\n\n");
    printf("where a self-contained parameter list must be supplied.\n\n");
    printf("Use the commands\n\n");
    printf("        ebem help\n");
    printf("or      ebem h\n\n");
    printf("for a detailed description on the input file format or\n");
    printf("'ebem epw ...' parameters.\n");
    return 0;
  }

// --- Read first argument and print out help if this is what it wanted.

  if(!strcmp(argv[1],"h") || !strcmp(argv[1],"help"))
    {printf(BEM_HELP);  return 0;}

  clear_all();  // Global initialization.

// --- Inline calculation of band structures.

  if(!strcmp(argv[1],"epw") || !strcmp(argv[1],"EPW")) {
    if(argc<9) on_error("*** syntax error in EBEM: epw",
      "more arguments are needed.\n    Type 'ebem h' for more details.");
    char name[100];
    epw_def epw;
    strcpy(name, argv[2]);                   // Output file.
    epw.a=atof(argv[3]);                     // Lattice constant.
    epw.R=atof(argv[4]);                     // Object size.
    epw.gmax=atof(argv[5]);                  // Number of plane waves.
    epw.nQ=atoi(argv[6]);                    // Number of sampling points.
    epw.U1=atof(argv[7]);                    // Pot. in. obj. (real part).
    epw.U0=atof(argv[8]);                    // Pot. out. obj. (real part).
    if(argc>9)  epw.obj=atoi(argv[9]);       // Object type.
    if(argc>10) epw.lat=atoi(argv[10]);      // Lattice type.
    if(argc>11) epw.U1+=i_c*atof(argv[11]);  // Pot. in. obj. (imag part).
    if(argc>12) epw.U0+=i_c*atof(argv[12]);  // Pot. out. obj. (imag part).
    if(argc>13) epw.meff=atof(argv[13]);     // Eff. electon mass.
    if(argc>14) epw.Emax=atof(argv[14]);     // Max. energy to be considered.
    epw.set_lattice=1;
    epw.get_band(name,0);
    return 0;
  }

  program(argv[1]);
  return 0;
}

// **************************************************************************
