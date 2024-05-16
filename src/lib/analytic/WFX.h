#ifndef CDFTT_WFX_H_INCLUDED
#define CDFTT_WFX_H_INCLUDED

#include<iostream>
#include<string>
#include<fstream>

using namespace std;

string get_one_block_from_wfx_file(istream&, string, int);
int get_one_block_int_from_wfx_file(istream&, string,  int);
double get_one_block_real_from_wfx_file(istream&, string, int);
bool get_one_int_from_wfx_file(istream&, string, int);
double get_one_orbital_from_wfx_file(istream&, int, int);
string readFile(string);

#endif