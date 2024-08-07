# CdftT

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

CdftT is set of a library and a program to compute CDFT descriptors. It is written in C++.\
It is still under developpement. However,several tools are already implemented. New tools will be implemented soon.

## Requirements

A C++ compiler.\
Tested with g++. However you can use any recent version of any C++ compiler.

## Installation

- Download
	- Using git, under a terminal, type : 
		```console
		git clone https://github.com/allouchear/CdftT
		```
	- You can also download the .zip file of CdftT :
		**Click on Code and Download ZIP**

- Compilation

		Edit <path_to_CdftT_dir>/CONFIG
		Set LIBCDFTTDIR corresponding to your machine.
	- For **Linux or MacOS**:\
		Under a terminal, type :
		```console
		cd <path_to_CdftT_dir>/src/applications/cdftt
		./cleancdftt.sh
		./compcdftt.sh
		```
	- For **Windows**:\
		Edit <path_to_CdftT_dir>/src/applications/cdftt/Makefile\
		Set LIBCDFTTDIR corresponding to your machine.\
		Under a terminal, type :
		```console
		cd <path_to_CdftT_dir>/src/applications/cdftt
		make
		```

## How to use it 

After compilation, type :
```console
 ./cdftt
```
You obtain a list of all implemented methods.\
A input file is requirement to run it.\
See examples folder. 

## knowns bugs
*Nothing yet.*

## Contributors
 - [Abdul-Rahman ALLOUCHE](https://sites.google.com/site/allouchear/Home)
 - [Dimitri BUFFAT](https://github.com/dbuffat) (Master 1/Physics/Univ Lyon 1, Intern supervised by A.R. Allouche)
 - [Tetautahi MAAMAATUAIAHUTAPU](https://github.com/tmaamaatua) (Master 1/Physics/Univ Lyon 1, Intern supervised by A.R. Allouche)
