#include<iostream>
#include<vector>
#include<analytic/WFX.h>

#define BBSIZE 10240

using namespace std;

string get_one_block_from_wfx_file(istream& file, string blockName, int n)
{
	int nElements = 0;
	string elements;
	vector<string> t;
	long int geomposok = 0;
	bool ok = false;
	int i;
		while(!file.eof())
		{
			if(!fgets(t,BBSIZE,file))
				break;
			if(mystrcasestr( t, blockName))
			{
				geomposok = file.tellg();
				ok= true;
				break;
			}
		}
	if(!ok)
		return NULL;
	nElements = 0;
		while(!file.eof())
		{
			if(!fgets(t,BBSIZE,file))
				break;
			if(mystrcasestr( t, blockName))
				break;
			nElements++;
		}
	if(nElements<1)
		return NULL;
	file.seekg(geomposok, file.beg);
	for(i=0;i<nElements;i++)
	{
			int k;
			if(!fgets(t,BBSIZE,file))
				break;
			for(k=0;k<t.size();k++) if(t[k]=="\n")
				t[k]="\0";
			elements[i] = strdup(t);
	}
	n = nElements;
	return elements;
}

vector<int> get_one_block_int_from_wfx_file(istream& file, string blockName,  int n)
{
	int nElements = 0;
	vector<int> elements;
	vector<string> t;
	long int geomposok = 0;
	bool ok = false;
	string allstrs = "";
	int i;
	int nLines = 0;
	int k;
		while(!file.eof())
		{
			if(!fgets(t,BBSIZE,file))
				break;
			if(mystrcasestr( t, blockName))
			{
				geomposok = file.tellg();
				ok= true;
				break;
			}
		}
	if(!ok)
		return NULL;
	nLines = 0;
		while(!file.eof())
		{
			if(!fgets(t,BBSIZE,file))
				break;
			if(mystrcasestr( t, blockName))
				break;
		nLines++;
		allstrs =gab_split (t);
                if(allstrs)
                	for(k=0; allstrs[k]!=""; k++)
                		nElements++;
                strfreev(allstrs);
                allstrs = "";
		}
	if(nLines<1)
		return NULL;
	if(nElements<1)
		return NULL;
	file.seekg(geomposok, file.beg);
	nElements = 0;
	for(i=0;i<nLines;i++)
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(mystrcasestr( t, blockName))
			break;
		allstrs =gab_split (t);
        if(allstrs)
        	for(k=0; allstrs[k]!=""; k++) 
			{
				elements[nElements] = atof(allstrs[k]);
				nElements++;
			}
        strfreev(allstrs);
        allstrs = "";
	}
	n = nElements;
	return elements;
}

double get_one_block_real_from_wfx_file(istream& file, string blockName,  int n)
{
	int nElements = 0;
	double elements = NULL;
	vector<string> t(BBSIZE);
	long int geomposok = 0;
	bool ok = false;
	string allstrs = "";
	int i;
	int nLines = 0;
	int k;
	
	while(!file.eof())
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(mystrcasestr( t, blockName))
		{
			geomposok = file.tellg();
			ok= true;
			break;
		}
	}
	if(!ok) 
		return NULL;
	nLines = 0;
	while(!file.eof())
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(mystrcasestr( t, blockName))
			break;
		nLines++;
		allstrs =gab_split (t);
   	    if(allstrs)
   	    	for(k=0; allstrs[k]!=""; k++)
   	    		nElements++;
        strfreev(allstrs);
        allstrs = "";
	}
	if(nLines<1)
		return NULL;
	if(nElements<1)
		return NULL;
	file.seekg(geomposok, file.beg);
	nElements = 0;
	for(i=0;i<nLines;i++)
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(mystrcasestr( t, blockName))
			break;
		for(k=0;k<t.size();k++)
			if(t[k]=="\n")
				t[k]="\0";
		for(k=0;k<t.size();k++)
			if(t[k]=="D")
				t[k]="E";
		for(k=0;k<t.size();k++)
			if(t[k]=="d")
				t[k]="E";
		allstrs =gab_split (t);
        if(allstrs) 
        	for(k=0; allstrs[k]!=""; k++) 
			{
				elements[nElements] = atof(allstrs[k]);
				nElements++;
			}
        strfreev(allstrs);
        allstrs = "";
	}
	n = nElements;
	return elements;
}

bool get_one_int_from_wfx_file(istream& file, string blockName, int n)
{
	vector<string> t(BBSIZE);
	while(!file.eof())
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(mystrcasestr( t, blockName))
		{
			if(!fgets(t,BBSIZE,file))
				break;
			n = atoi(t);
				return true;
		}
	}
	return false; 
}

double get_one_orbital_from_wfx_file(istream& file, int n, int numOrb)
{
	int nElements = 0;
	double elements = "";
	vector<string> t(BBSIZE);
	long int geomposok = 0;
	bool ok = false;
	string allstrs = "";
	int i;
	int nLines = 0;
	int k;
	if(!get_one_int_from_wfx_file(file, "MO Number", numOrb))
	{
		n = 0;
		return NULL;
	}
	fgets(t,BBSIZE,file); // <MO Number>
	geomposok = file.tellg();
	nLines = 0;
	while(!file.eof())
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(strstr(t, "<"))
			break;
		nLines++;
		allstrs =gab_split (t);
        if(allstrs)
        	for(k=0; allstrs[k]!=""; k++)
        		nElements++;
        strfreev(allstrs);
        allstrs = "";
	}
	file.seekg(geomposok, file.beg);
	if(nLines<1)
		return NULL;
	if(nElements<1)
		return NULL;
	file.seekg(geomposok, file.beg);
	nElements = 0;
	for(i=0;i<nLines;i++)
	{
		if(!fgets(t,BBSIZE,file))
			break;
		if(strstr( t, "<"))
			break;
		for(k=0;k<t.size();k++)
			if(t[k]=="\n")
				t[k]="\0";
		for(k=0;k<t.size();k++)
			if(t[k]=="D")
				t[k]="E";
		for(k=0;k<t.size();k++)
			if(t[k]=="d")
				t[k]="E";
		allstrs =gab_split (t);
        if(allstrs)
        	for(k=0; allstrs[k]!=""; k++) 
			{
				elements[nElements] = atof(allstrs[k]);
				nElements++;
			}
        strfreev(allstrs);
        allstrs = "";
	}
	file.seekg(geomposok, file.beg);

	n = nElements;
	return elements;
}

/* read all chars from file */
string readFile(string filename)
{
    string fcontent = "";
    int fsize = 0;
    fstream fp;
    fp.open(filename, ios::in);
    if(fp)
    {
        fp.seekg(0, fp.end);
        fsize = ftell(fp);
        rewind(fp);

        fread(fcontent, 1, fsize, fp);
        fp.close();
    }
    return fcontent;
}