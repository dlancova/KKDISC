# KKDISC
Kluzniak-Kita disk in KORAL code

---- WORK IN PROGRESS -- 


Short manual for KORAL from Miki:

MC040422--additions for KK00 disk setup in Koral
-in problems.h are definitions of problems and, under each problem
loop, the relevant choices of used modules.
-in mnemonics.h are additional defs, shortcuts.
-in choices.h are some physics choices
-in finite.c are routines related to grid & finite difference
-in physics.c is problem independent physics
-in bc.c are boundary conditions
-no dynamo in Newtonian means no MRI

Units? Everything is done in code units, which are geometrized units,
defined in 

MC200422, 210322--basic manual for Koral use
-in Koral/ directory, where you have a running test version, create one
for your problem e.g. KK00disk, by copying the running one as a template
into new name.
-copy the most similar of available problems from PROBLEMS directory into
a new problem, I made a KKdisk, so I give this as an example here.
-modify PROBLEM/KKdisk/define.h and PROBLEM/KKdisk/init.c
-you can add some part of the setup into tools.c, here I left it empty.
-clean KK00disk/dumps and KK00disk/analysis
-add the paths to your gsl and silo libraries in the file Makefile
-be careful with dumps directory, to have it there (or create it if not
there), as Koral will not create it by itself, will crash if it does not
find it.
-for serial compilation, use ./mser.sh, then run with ./ko (you need silo
for serial runs)
-for serial run, you can go to /KK00disk/dumps and plot using silo output
with visit
-for parallel compilation (you do not need silo), use ./mpar.sh, then run
with mpirun -np 4 ./ko for e.g. 4 processors. Take care that NT(XYZ) has
to be divisible with the number of processors (there will be a message if not).
-for plotting after serial run, just go to "dumps" folder and plot with
visit the *.silo files.
-for plotting after parallel run, you need to compile in serial /mser.sh and
then under /KK00disk execute ./ana 0 12 1 where numbers are n_0 n_max n_step
(12 for 12 timesteps), and you get inside /KK00disk/analysis (create the
analysis directory if there is no, it will not be created automatycally!)
the silo output which you can plot with visit. You do not need to recompile
ana with change of number, just run it again for e.g. 150 files of output...
-for restart, in define.h change to have #define restart and put after the
RESTARTNUM X instead X the number of file from which you want to restart.
Interesting that I changed the lines after run, not before, and they were
executed, so this change does not need recompilation. Other changes need.
//restart
/************************************/ 
#define RESTART
#define RESTARTNUM 52
======================== 


