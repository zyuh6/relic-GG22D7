#!/bin/sh
cmake -DCHECK=off -DARITH=gmp -DFP_PRIME=381 -DFP_QNRES=on -DFP_METHD="BASIC;COMBA;COMBA;MONTY;JMPDS;JMPDS;SLIDE" -DFPX_METHD="INTEG;INTEG;LAZYR" -DPP_METHD="LAZYR;OATEP" -DCFLAGS="-O2 -funroll-loops -fomit-frame-pointer" $1