#!/bin/bash

export outputDirectory=$1
echo $outputDirectory

diff smv2chem2_exit.proc0017 $outputDirectory/
diff smv2chem1_exit.proc0017 $outputDirectory/
diff physproc_exit.proc0017 $outputDirectory/
