Smvgear-Standalone
==================
This driver reads the input state from the following entry files:
smv2chem2_entry.proc0017
smv2chem1_entry.proc0017
physproc_entry.proc0017



To set up the input and exit data
==================
Obtain these files from the NCCS in this directory:
scp discover:/archive/anon/pub/gmidata2/users/mrdamon/Smvgear-Standalone-Data.tar.gz .

Extract files:
gunzip Smvgear-Standalone-Data.tar.gz
tar xvf Smvgear-Standalone-Data.tar

Link the input data:
./LinkData.bash


Building and running the component
==================
make 
./Do_Smv2_Solver.exe 


Directions for validating output data
==================
./Diff.bash


Cleaning the installation
==================
make clean
or
make distclean (to remove exit files)
 


