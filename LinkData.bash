#!/bin/bash

export procNumber=$1
echo $procNumber

ln -s Smvgear-Standalone-Data/smv2chem2_entry.proc${procNumber} 
ln -s Smvgear-Standalone-Data/smv2chem1_entry.proc${procNumber} 
ln -s Smvgear-Standalone-Data/physproc_entry.proc${procNumber} 
