#!/bin/bash

procString=$1

diff Smvgear-Standalone-Data/exitFiles_O0/smv2chem2_exit.proc${procString} . 
diff Smvgear-Standalone-Data/exitFiles_O0/smv2chem1_exit.proc${procString} . 
diff Smvgear-Standalone-Data/exitFiles_O0/physproc_exit.proc${procString} . 
