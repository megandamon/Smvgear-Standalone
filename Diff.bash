#!/bin/bash

procString=$1

diff ${procString}/smv2chem2_exit.proc${procString} . 
diff ${procString}/smv2chem1_exit.proc${procString} . 
diff ${procString}/physproc_exit.proc${procString} . 
