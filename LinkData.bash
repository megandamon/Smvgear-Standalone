#!/bin/bash

export procNumber=$1
echo $procNumber

ln -s ${procNumber}/smv2chem2_entry.proc${procNumber} 
ln -s ${procNumber}/smv2chem1_entry.proc${procNumber} 
ln -s ${procNumber}/physproc_entry.proc${procNumber} 
