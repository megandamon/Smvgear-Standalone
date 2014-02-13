Smvgear-Standalone
==================
This driver reads the input state from the following entry* files:
 - smv2chem2_entry.proc0133
 - smv2chem1_entry.proc0133
 - physproc_entry.proc0133
*proc0133 is for process 133


To set up the input and exit data for specific process block
==================
For example, obtain the data files for process 133 from the NCCS in this directory:
scp discover:/archive/anon/pub/gmidata2/users/mrdamon/Smvgear-Standalone-Data/proc_archives/0133.tar .

Extract files:
tar xvf 0133.tar

Link the input data:
./LinkData.bash 0133

Building and running the component with specific processor data
==================
make 
./Do_Smv2_Solver.exe -r 133

Directions for validating output data
==================
 ./Diff.bash 0133

Cleaning the installation
==================
make clean
or
make distclean (to remove exit files)



Developer notes - 

Feb 12 2014
An optimization made changed the species concentrations by: 
   Max  =   8.819316628546657623811789175258837E-0010
   Mean =   3.695324632725830753419521841913095E-0013
   Std  =   5.240100518479359891714028047851416E-0012
