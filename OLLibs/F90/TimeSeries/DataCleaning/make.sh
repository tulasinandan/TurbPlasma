#!/bin/bash
f2py3 --opt=-O3 -c qsort.f90 nan_algebra.f90 fill_gaps.f90 autobad.f90 hampel.f90 -m dc
