#!/bin/bash

procString=$1

diff Smvgear-Standalone-Data/exitFiles/smv2chem2_exit.proc${procString} . 
diff Smvgear-Standalone-Data/exitFiles/smv2chem1_exit.proc${procString} . 
diff Smvgear-Standalone-Data/exitFiles/physproc_exit.proc${procString} . 
