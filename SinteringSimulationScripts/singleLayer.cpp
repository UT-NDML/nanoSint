#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <new>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <mpi.h>
#include <ctime>

using namespace std;
//using namespace MPI; 

//grouping relevant data together
struct variables
{
int xlength, ylength, zlength;
int noOfparticles;
double rhoVap, etaVap,rhoSolid, etaSolid;
} dim;

struct CircleInfo
{
int midX,midY, midZ, particleNo;
double value, radius;
} circle;

struct constants
{
double surfdif, gbdif, voldif,vapdif;
double surfenergy, gbenergy, gbmobility;
double randR, randE, deltaT;
double A,B,Cgbe;
}con;

//memory allocation for 1D
double *store1D(int maxArray)
{
double *k;
k = new double [maxArray];

if (k==NULL)
{
cout<<"Failed to allocate 2D memory";
exit(EXIT_FAILURE);
}
return k;
}

int *store1Dint(int maxArray)
{
int *k;
k = new int [maxArray];

if (k==NULL)
{
cout<<"Failed to allocate 2D memory";
exit(EXIT_FAILURE);
}
return k;
}

//memory allocation for 2D
double **store2D(int maxArray1,int maxArray2)
{
double **m;
m = new double *[maxArray1];

for (int a = 0; a<maxArray1; a++)
{
    m[a]= new double [maxArray2];
}

if (m==NULL)
{
cout<<"Failed to allocate 2D memory";
exit(EXIT_FAILURE);
}
 return m;
}

int **store2Dint(int maxArray1,int maxArray2)
{
int **m;
m = new int *[maxArray1];

for (int a = 0; a<maxArray1; a++)
{
    m[a]= new int [maxArray2];
}

if (m==NULL)
{
cout<<"Failed to allocate 2D memory";
exit(EXIT_FAILURE);
}
 return m;
}


//clear 1D array
void clear1Dint(int* j1)
{
delete [] j1;
return;
}

void clear1D(double* j1)
{
delete [] j1;
return;
}

//clear 2D array
void clear2D(double **j2, int a)
{
for (int i=0;i<a;i++)
{
delete [] j2 [i];
}
delete [] j2;
return;
}

void clear2Dint(int **j2, int a)
{
for (int i=0;i<a;i++)
{
delete [] j2 [i];
}
delete [] j2;
return;
}

/*Draw circle to assign the field variables in the circle (rho and eta for solid) and the variables outside the circle in the vapor stage */
void CreateParticles (double *g, struct variables dim, struct CircleInfo circle, int prank, int zheight, int nopy, int ysize)
{
int x,y,z, i, zh, zpos, ypos, yh;

ypos = prank%nopy;
zpos = floor(prank/nopy);

for (z = 0; z< zheight; z++)
{
    for (y=0;y<ysize;y++)
    {
        for (x=0;x<dim.xlength;x++)
        {
        zh = z + (zheight*zpos);
        yh = y + (ysize*ypos);

            if(((x-circle.midX)*(x-circle.midX))+((yh-circle.midY)*(yh-circle.midY))+((zh-circle.midZ)*(zh-circle.midZ))<=(circle.radius*circle.radius))
            {
                i = (dim.xlength*ysize*z)+(dim.xlength*y) +x;

                g[i]= circle.value;
            }
        }
    }
}

return;
}

/*print results to a file */


void fileres (double t/*like the actual time t*/, double timet/*time you want*/, double deltatime, double **E, string namestart, string nameend /*include .txt*/, int ySize, int xSize, int zSize, int NoOfParticles, int prank)

{

double f1 = timet - deltatime;
double f2 = timet + deltatime;

if ((t>f1)&&(t<f2))
{
ostringstream oss;
string var;
char *name;
oss<<namestart<<prank<<nameend;

var = oss.str();

name = new char [var.length()];

strcpy(name, var.c_str());
ofstream myfile;

myfile.open(name);

myfile<<"x\ty\tz\trho\t";

for (int i = 0; i<NoOfParticles; i++)
{
    myfile<<"eta"<<i<<"\t";
}

myfile<<"\n";

for (int z = 0; z<zSize; z++)
{
int zh = z + (zSize*prank);
for (int y = 0; y<ySize; y++)
{
    for (int x = 0; x<xSize; x++)
    {
        int l = (xSize*ySize*z)+(y*xSize) + x;

        myfile<<x<<"\t"<<y<<"\t"<<zh<<"\t"<<E[0][l]<<"\t";

        for (int i = 1; i<=NoOfParticles; i++)
        {
            myfile<<E[i][l]<<"\t";
        }

        myfile<<"\n";
    }
}
}

myfile.close();

}

}

void fileprintsum/*matlab form*/ (double t/*like the actual time t*/, double timet/*time you want*/, double deltatime, double **E, string namestart, string nameend /*include .txt*/, int ySize, int xSize, int zSize, int NoOfParticles, int prank)

{

double f1 = timet - deltatime;
double f2 = timet + deltatime;

if ((t>f1)&&(t<f2))
{
double sum;
double sump;

ostringstream oss;
string var;
char *name;
oss<<namestart<<prank<<nameend;

var = oss.str();

name = new char [var.length()];

strcpy(name, var.c_str());
ofstream myfile;

myfile.open(name);

for (int z = 0; z<zSize; z++)
{
int zh = z + (zSize*prank);
for (int y = 0; y<ySize; y++)
{
    for (int x = 0; x<xSize; x++)
    {
        int l = (xSize*ySize*z)+(y*xSize) + x;

        sump = 0.0;

        for(int i = 1; i<=NoOfParticles; i++)
        {
            sump = sump + E[i][l];
        }

        sum = E[0][l]*sump;

        myfile<<sum<<"\n";

    }
}
}

myfile.close();

}

}

void fileprintsumpr/*matlab form*/ (double t/*like the actual time t*/, double timet/*time you want*/, double deltatime, double **E, string namestart, string nameend /*include .txt*/, int procSize, int NoOfParticles, int prank)

{

double f1 = timet - deltatime;
double f2 = timet + deltatime;

if ((t>f1)&&(t<f2))
{
double sum;
double sump;

ostringstream oss;
string var;
char *name;
oss<<namestart<<prank<<nameend;

var = oss.str();

name = new char [var.length()];

strcpy(name, var.c_str());
ofstream myfile;

myfile.open(name);

for (int l = 0; l<procSize; l++)
{
        sump = 0.0;

        for(int i = 1; i<=NoOfParticles; i++)
        {
            sump = sump + E[i][l];
        }

        sum = E[0][l]*sump;

        myfile<<sum<<"\n";
}

myfile.close();
}

}

void MPIfileprintsum (double t, double timet, double deltatime, double **E, /*const char *filename,*/ string namestart, string nameend,double eval_time, int SNo /*serial number of run*/, int procSize, int prank,struct variables dim, int zpos, int ypos, int zheight, int ysize)
{

double f1 = timet - deltatime;
double f2 = timet + deltatime;

if ((t>f1)&&(t<f2))
{
ostringstream oss;
string var;
char *filename;
oss<<namestart<<eval_time<<"SN"<<SNo<<nameend;

var = oss.str();

filename = new char [var.length()+1];

strcpy(filename, var.c_str());

double sump, *sum;
double size_dob;
int amode,t,area;

MPI_Status status;
MPI_File fh;
MPI_Offset disp;

amode = (MPI_MODE_CREATE|MPI_MODE_WRONLY);
sum = store1D(procSize);

MPI_File_open(MPI_COMM_WORLD,filename,amode,MPI_INFO_NULL,&fh);

for (int l = 0; l<procSize; l++)
{
    sump = 0.0;

    for(int i = 1; i<=dim.noOfparticles; i++)
    {
        sump = sump + E[i][l];
    }

    sum[l] = E[0][l]*sump;

}

area = dim.xlength*ysize;

for (int z = 0; z<zheight;z++)
{
t = dim.xlength*ysize*z;
disp = ((z*dim.ylength*dim.xlength)+(ypos*ysize*dim.xlength)+(zpos*dim.xlength*dim.ylength*zheight))*sizeof(size_dob);
MPI_File_write_at(fh,disp,&sum[t],area,MPI_DOUBLE,&status);
}

MPI_File_close(&fh);

delete [] filename;
clear1D(sum);

}

}

void MPIprintvalue (double t, string namestart, string nameend, double value, int prank)
{
ostringstream oss;
string var;
char *filename;
oss<<namestart<<t<<nameend;

var = oss.str();

filename = new char [var.length()+1];

strcpy(filename, var.c_str());

double size_dob;
int amode;

MPI_Status status;
MPI_File fh;
MPI_Offset disp;

amode = (MPI_MODE_CREATE|MPI_MODE_WRONLY);

MPI_File_open(MPI_COMM_WORLD,filename,amode,MPI_INFO_NULL,&fh);

disp = prank*sizeof(size_dob);
MPI_File_write_at(fh,disp,&value,1,MPI_DOUBLE,&status);

MPI_File_close(&fh);

delete [] filename;
}

void MPIfileprintrho(double t, double timet, double deltatime, double **E, /*const char *filename,*/ string namestart, string nameend,double eval_time, int SNo, int procSize, int prank,struct variables dim, int zpos, int ypos, int zheight, int ysize)
{

double f1 = timet - deltatime;
double f2 = timet + deltatime;

if ((t>f1)&&(t<f2))
{
ostringstream oss;
string var;
char *filename;
oss<<namestart<<eval_time<<"SN"<<SNo<<nameend;

var = oss.str();

filename = new char [var.length()+1];

strcpy(filename, var.c_str());

double *rhoP;
double size_dob;
int amode,t,area;

MPI_Status status;
MPI_File fh;
MPI_Offset disp;

amode = (MPI_MODE_CREATE|MPI_MODE_WRONLY);
rhoP = store1D(procSize);

MPI_File_open(MPI_COMM_WORLD,filename,amode,MPI_INFO_NULL,&fh);

for (int l = 0; l<procSize; l++)
{

    rhoP[l] = E[0][l];

}

area = dim.xlength*ysize;

for (int z = 0; z<zheight;z++)
{
t = dim.xlength*ysize*z;
disp = ((z*dim.ylength*dim.xlength)+(ypos*ysize*dim.xlength)+(zpos*dim.xlength*dim.ylength*zheight))*sizeof(size_dob);
MPI_File_write_at(fh,disp,&rhoP[t],area,MPI_DOUBLE,&status);
}

MPI_File_close(&fh);

delete [] filename;
clear1D(rhoP);

}

}

void MPIprintfullE (double t, double timet, double deltatime, double **E, /*const char *filename,*/ string namestart, string nameend,double eval_time, int SNo, int procSize, int prank,struct variables dim, int zpos, int ypos, int zheight, int ysize)
{

double f1 = timet - deltatime;
double f2 = timet + deltatime;

if ((t>f1)&&(t<f2))
{
ostringstream oss;
string var;
char *filename;
oss<<namestart<<eval_time<<"SN"<<SNo<<nameend;

var = oss.str();

filename = new char [var.length()+1];

strcpy(filename, var.c_str());

double *Eprt;
double size_dob;
int amode,t,area;

MPI_Status status;
MPI_File fh;
MPI_Offset disp;

amode = (MPI_MODE_CREATE|MPI_MODE_WRONLY);
Eprt = store1D(procSize*(dim.noOfparticles+1));

MPI_File_open(MPI_COMM_WORLD,filename,amode,MPI_INFO_NULL,&fh);

int a = 0;
for (int i = 0; i<procSize; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    Eprt[a] = E[g][i];
    a = a+1;
    }

}

area = dim.xlength*ysize;

for (int z = 0; z<zheight;z++)
{
t = dim.xlength*ysize*z*(dim.noOfparticles+1);
disp = ((z*dim.ylength*dim.xlength)+(ypos*ysize*dim.xlength)+(zpos*dim.xlength*dim.ylength*zheight))*sizeof(size_dob)*(dim.noOfparticles+1);
MPI_File_write_at(fh,disp,&Eprt[t],area*(dim.noOfparticles+1),MPI_DOUBLE,&status);
}

MPI_File_close(&fh);

delete [] filename;
clear1D(Eprt);

}

}

void MPIprint1D (double *E, const char *filename, int arraySize, int prank)
{
double size_dob;
int amode;

MPI_Status status;
MPI_File fh;
MPI_Offset disp;

amode = (MPI_MODE_CREATE|MPI_MODE_WRONLY);

MPI_File_open(MPI_COMM_WORLD,filename,amode,MPI_INFO_NULL,&fh);

disp = prank*sizeof(size_dob)*arraySize;

MPI_File_set_view(fh,disp,MPI_INT,MPI_INT,"native",MPI_INFO_NULL);

MPI_File_write(fh,E,arraySize,MPI_DOUBLE,&status);

MPI_File_close(&fh);
}

void print2D (double **E/*2D quantity to print to file*/, string namestart, string nameend /*include .txt*/, int rowSize, int colSize, int prank)
{
ostringstream oss;
string var;
char *name;
oss<<namestart<<prank<<nameend;

var = oss.str();

name = new char [var.length()];

strcpy(name, var.c_str());
ofstream myfile;

myfile.open(name);

if (rowSize > colSize)
{
myfile<<"row\t";

for (int i = 0; i<colSize; i++)
{
    myfile<<"of["<<i<<"]\t";
}

myfile<<"\n";

for (int i = 0; i<rowSize; i++)
{

        myfile<<i<<"\t";

        for (int g = 0; g<colSize; g++)
        {
            myfile<<E[i][g]<<"\t";
        }

        myfile<<"\n";
}

}

else
{
myfile<<"row\t";

for (int i = 0; i<rowSize; i++)
{
    myfile<<"of["<<i<<"]\t";
}

myfile<<"\n";

for (int i = 0; i<colSize; i++)
{

        myfile<<i<<"\t";

        for (int g = 0; g<rowSize; g++)
        {
            myfile<<E[g][i]<<"\t";
        }

        myfile<<"\n";
}

}

myfile.close();

}

void print1D (double *E/*1D quantity to print to file*/, string namestart, string nameend /*include .txt*/, int Size, int prank)

{
ostringstream oss;
string var;
char *name;
oss<<namestart<<prank<<nameend;

var = oss.str();

name = new char [var.length()];

strcpy(name, var.c_str());
ofstream myfile;

myfile.open(name);

myfile<<"row\n";

for (int i = 0; i<Size; i++)
{
        myfile<<i<<"\t"<<E[i]<<"\n";
}

myfile.close();

}

void timeME(int prank, int chosenrank, clock_t start)
{

if (prank == chosenrank)
{
	cout<<"Time elapsed: "<<(clock()-start)/((double)CLOCKS_PER_SEC)<<"secs\n";
}

}

int main( int argc, char* argv[])
{
clock_t start;
start = clock();

double **circleProp, *BFE_partialEta, *Etabrcvup, *Etabrcvdown, *Etabrcvleft, *Etabrcvright;
double *rhobounds2down, *rhobounds2up, *rhobounds2right, *rhobounds2left;
double **NewEta, **Eta, **Etaboundsup, **Etaboundsdown, **Etaboundsright, **Etaboundsleft;
double *Eta1, *Etatr, **Etas, *EtaR, *Eta1y;
double BFE_partialRho, intveta, intvrho,*procInit, *procrho, sumrhototInit, sumrhotot,intrvl,temprho;
double *NewRho, *Dif, *h, *rhot, *rhob, *rhor, *rhol;
double *ht, *hb, *Dift, *Difb, *hr, *hl, *Difr, *Difl;
double *rhorup,*rholup,*rhordn,*rholdn;
double RaND, etaSqSum, etaCubeSum, rhoSq, rhoQuad, A2, B12;
double phi,  sumrho,sumeta, sumrhoInit, sumetaInit, BFE_partialRhot, limR, cutoffR;
double t, finalt, rlfinalt, temp, tempn, tempE, tempEN, evat, evatf, countt, counttf;
int Size, N_size, SNo;
int pixelRadius, rad, ffile, sizebv,*procbv,sizebvtot, pid, tempbv, tempbvtot ;
int **N, **N_proc, **Ny_proc, simcheckV;
int runno, noruns, maxruns,startruns,markct,fileno,simno,progno,tempI;

double intvs, intvf; //intvs for printsum interval intvf printfull interval

double rho0,rho1,rho2,rho3,rho4,rho5,eta1,eta0;
double Dif4,Dif5,Dif2,Dif3,Dif1,Dif0,h4,h5,h2,h3,h0,h1;

//defining MPI related terms
int prank, psize, sendcount, psizey, psizez, sendcounty;
int procSize;
int zheight, area,ysize,zpos,ypos, areay;

//Modification eta_product instead of etasum
double etaij = 0.0;

//ffile 0-no input from file 1-input from file
ffile = 0;

//assigning values
t = 0.0;
finalt = 330.0001;
N_size = 6;
dim.etaSolid = 1.0;
dim.etaVap = 0.0;
//dim.noOfparticles = 31;
dim.rhoSolid = 0.9998;
dim.rhoVap = 0.000000089;
//dim.xlength = 51;
dim.ylength = 192;
dim.zlength = 192;

//command line arguments
stringstream strg(argv[1]), strg2(argv[2]);
strg >> dim.noOfparticles;
strg2 >> dim.xlength;

con.randE = 0.000000;
con.randR = 0.0000001;
con.deltaT = 0.000095;
con.gbmobility = 10.0;
con.vapdif = 0.0;

/*
con.A =12.03;
con.B =4.1;
con.Cgbe =7;
con.gbenergy =6;
con.surfenergy =26.12;
con.gbdif =9.53;
con.surfdif =105.3;
con.voldif =0.0753;
*/

con.A =8.03;
con.B =4.1;
con.Cgbe =7;
con.gbenergy =6;
con.surfenergy =23.12;
con.gbdif =7.13;
con.surfdif =71.3;
con.voldif =0.0713;

intvs = 0.5;
intvf = 100;
int printEND = 1;
int noisy = 0; //variable to control console output  

SNo = 504; //sequence number -------> change with cutOffR <---------

pid = 0; //processor to do sumation analysis
limR = 0.0041; 
cutoffR = 0.00002; //0.0001;

pixelRadius = 10;
Size = dim.xlength*dim.ylength*dim.zlength;

vector <vector<double> > eta;

MPI_Init (&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD,&prank);
MPI_Comm_size(MPI_COMM_WORLD,&psize);

MPI_Status Status;

psizey = 32; //set size/number of processors per z slice
psizez = 32; //set number of z slices

if (psize!=(psizey*psizez))
{
cout<<"wrong number of proessors check the number of processors being used.\nExiting.";
MPI_Finalize ();
return 0;
}

A2 = 2*con.A;
B12 = 12*con.B;

procSize = (int) Size/psize;
zheight = (int) dim.zlength/psizez;
ysize = (int) dim.ylength/psizey;

area = dim.xlength*ysize;
areay = dim.xlength*zheight;

//positions in y and z
zpos = floor(prank/psizey);
ypos = prank%psizey;

N = store2Dint(procSize, N_size);
N_proc = store2Dint(area, N_size);
Ny_proc = store2Dint(areay, N_size);

//Creating nearest neighbors for numerical differentiation using partial derivatives
for (int z=0; z<zheight; z++)
{
for (int y=0;y<ysize;y++)
{
    for (int x=0;x<dim.xlength;x++)
    {
    int k;

    k = (area*z)+(y*dim.xlength) + x;

    N[k][0] = ((x+dim.xlength+1)%dim.xlength) + (y*dim.xlength) + (ysize*dim.xlength*z);//f(x+dx,y,z) dx=dy=dz=di=1
    N[k][1] = ((x+dim.xlength-1)%dim.xlength) + (y*dim.xlength) + (ysize*dim.xlength*z);//f(x-dx,y,z)
    N[k][2] = x + ((y+ysize+1)%ysize)*dim.xlength + (ysize*dim.xlength*z);//f(x,y+dy,z)
    N[k][3] = x + ((y+ysize-1)%ysize)*dim.xlength + (ysize*dim.xlength*z);//f(x,y-dy,z)
    N[k][4] = x + (y*dim.xlength)+ (ysize*dim.xlength*((z+zheight+1)%zheight));//f(x,y,z+dz)
    N[k][5] = x + (y*dim.xlength) + (ysize*dim.xlength*((z+zheight-1)%zheight));//f(x,y,z-dz)

    }
}
}

for (int y=0;y<ysize;y++)
{
    for (int x=0;x<dim.xlength;x++)
    {
    int k,z;

    z = 1;

    k = (y*dim.xlength) + x;

    N_proc[k][0] = ((x+dim.xlength+1)%dim.xlength) + (y*dim.xlength) + (ysize*dim.xlength*z);//f(x+dx,y,z) dx=dy=dz=di=1 at z = 1 for all N_proc just for rhot and rhob creating ht and hb
    N_proc[k][1] = ((x+dim.xlength-1)%dim.xlength) + (y*dim.xlength) + (ysize*dim.xlength*z);//f(x-dx,y,z)
    N_proc[k][2] = x + ((y+ysize+1)%ysize)*dim.xlength + (ysize*dim.xlength*z);//f(x,y+dy,z)
    N_proc[k][3] = x + ((y+ysize-1)%ysize)*dim.xlength + (ysize*dim.xlength*z);//f(x,y-dy,z)
    N_proc[k][4] = x + (y*dim.xlength) + (ysize*dim.xlength*((z+3+1)%3));//f(x,y,z+dz)
    N_proc[k][5] = x + (y*dim.xlength) + (ysize*dim.xlength*((z+3-1)%3));//f(x,y,z-dz)

    }
}

for (int z=0;z<zheight;z++)
{
    for (int x=0;x<dim.xlength;x++)
    {
    int y = 1;
    int k = (z*dim.xlength) + x;

    Ny_proc[k][0] = ((x+dim.xlength+1)%dim.xlength) + (y*dim.xlength) + (3*dim.xlength*z);//f(x+dx,y,z) dx=dy=dz=di=1
    Ny_proc[k][1] = ((x+dim.xlength-1)%dim.xlength) + (y*dim.xlength) + (3*dim.xlength*z);//f(x-dx,y,z)
    Ny_proc[k][2] = x + ((y+3+1)%3)*dim.xlength + (3*dim.xlength*z);//f(x,y+dy,z)
    Ny_proc[k][3] = x + ((y+3-1)%3)*dim.xlength + (3*dim.xlength*z);//f(x,y-dy,z)
    Ny_proc[k][4] = x + (y*dim.xlength)+ (3*dim.xlength*((z+zheight+1)%zheight));//f(x,y,z+dz)
    Ny_proc[k][5] = x + (y*dim.xlength) + (3*dim.xlength*((z+zheight-1)%zheight));//f(x,y,z-dz)

    }
}

sendcount = area*(dim.noOfparticles+1);
sendcounty = areay*(dim.noOfparticles+1);

eta.resize(procSize);

for(int i = 0; i<procSize; i++)
{
eta[i].resize(dim.noOfparticles+1);
}

//allocating memory
Etas = store2D(dim.noOfparticles+1,procSize);

if(ffile == 0)
{
circleProp = store2D(dim.noOfparticles,N_size);

//start by setting the entire sample area to the vapor phase
for (int i = 0; i<procSize; i++)
{
RaND = ((double)rand()/RAND_MAX);
Etas[0][i]=dim.rhoVap +(0.5-RaND)*con.randR;
 //density field variable
    for (int g = 1; g<=dim.noOfparticles; g++)
    {
    Etas[g][i] = dim.etaVap+(0.5-RaND)*con.randE; //orientation parameter
    }
}

/* help set circle properties, "in pixels"
LEGEND:- (:,0)- cirlce radius; (:,1)-x middle of circle; (:,2)-ymiddle of circle; (:,3)-zmiddle of circle; (0 or 1,:) particle number */


ifstream particleInfofile ("bed.txt");

for (int i = 1; i <= dim.noOfparticles; i++)
{
int j = i-1;

particleInfofile>>circleProp[j][0]>>circleProp[j][3]>>circleProp[j][2]>>circleProp[j][1];
}


//filling the circle with the right field variable value for rho and eta
for (int g =1; g<=dim.noOfparticles;g++)
{
    circle.radius = circleProp[g-1][0];
    circle.midX = circleProp[g-1][1];
    circle.midY = circleProp[g-1][2];
    circle.midZ = circleProp[g-1][3];
    circle.value = dim.rhoSolid;
    circle.particleNo = g;
    CreateParticles(Etas[0],dim,circle,prank,zheight,psizey,ysize);
    circle.value = dim.etaSolid;
    CreateParticles(Etas[circle.particleNo],dim,circle,prank,zheight,psizey,ysize);
}

clear2D(circleProp,dim.noOfparticles);
}

else
{
double size_dob;
int amode,ti;

MPI_Status status;
MPI_File fh;
MPI_Offset disp;

amode = (MPI_MODE_RDONLY);
EtaR = store1D(procSize*(dim.noOfparticles+1));

MPI_File_open(MPI_COMM_WORLD,"fullT2600SN504000002.dat",amode,MPI_INFO_NULL,&fh);

for (int z = 0; z<zheight;z++)
{
ti = dim.xlength*ysize*z*(dim.noOfparticles+1);
disp = ((z*dim.ylength*dim.xlength)+(ypos*ysize*dim.xlength)+(zpos*dim.xlength*dim.ylength*zheight))*sizeof(size_dob)*(dim.noOfparticles+1);
MPI_File_read_at(fh,disp,&EtaR[ti],area*(dim.noOfparticles+1),MPI_DOUBLE,&status);
}

MPI_File_close(&fh);

int counta = 0;

for (int i = 0; i< procSize; i++)
{
    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    tempE = EtaR[counta];
    Etas[g][i] = tempE;

    counta = counta + 1;

    }

}

clear1D(EtaR);

//where simulation starts
t = 260.0;
}

//check initial circle microstructure
//MPIfileprintsum(t,t,con.deltaT,Etas,"3p","0.0001.dat",t,procSize,prank,dim,zpos,ypos,zheight,ysize);
MPIfileprintrho(t,t,con.deltaT,Etas,"rho",".dat",t,SNo,procSize,prank,dim,zpos,ypos,zheight,ysize);
//MPIprintfullE(t,t,con.deltaT,Etas,"fullTb",".dat",10*t,procSize,prank,dim,zpos,ypos,zheight,ysize);


sumrhoInit = 0.0;

for (int i = 0; i<procSize; i++)
{
    sumrhoInit = sumrhoInit +Etas[0][i];
}

//MPIprintvalue (t, "Rhosum", ".dat", sumrhoInit, prank);

for (int i = 0; i<procSize; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {
   	eta[i][g] = 0.0;
    	eta[i][g] = Etas[g][i];
    }

}

clear2D(Etas,dim.noOfparticles+1);

if (ffile==0)
{
 countt = 0;
 counttf = 0;
}

else
{
 countt = t/intvs;
 counttf = t/intvf;
}

MPI_Reduce(&sumrhoInit,&sumrhototInit,1,MPI_DOUBLE,MPI_SUM,pid,MPI_COMM_WORLD);
MPI_Bcast(&sumrhototInit,1,MPI_DOUBLE,pid,MPI_COMM_WORLD);

int chosenrank = 20;

//if (prank==chosenrank){cout<<"1. Entering While\n";}
//timeME(prank, chosenrank, start);

//===================================//
//start sintering simulation analysis//
//===================================//

while (t < finalt)
{

//>>allocating memory for received data from other processors<<//

rhobounds2down = store1D(area);
rhobounds2up = store1D(area);
rhobounds2right = store1D(areay);
rhobounds2left = store1D(areay);

rhorup = store1D(dim.xlength);
rholup = store1D(dim.xlength);
rhordn = store1D(dim.xlength);
rholdn = store1D(dim.xlength);

if (zpos!=0)
{
Etabrcvdown = store1D(sendcount+area);
}
if (zpos!=(psizez-1))
{
Etabrcvup = store1D(sendcount+area);
}

if (ypos!=0)
{
Etabrcvleft = store1D(sendcounty+areay);
}
if (ypos!=(psizey-1))
{
Etabrcvright = store1D(sendcounty+areay);
}

//----------------------------------------------------//
//>>sending and receiving data from other processors<<//
//----------------------------------------------------//

Etatr = store1D(4*dim.xlength);
Eta1 = store1D(2*(sendcount+area));

int a = 0;
for (int i = 0; i<area; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    Eta1[a] = eta[i][g];
    a = a+1;
    }

}

for (int i = area; i<2*area; i++)
{

    Eta1[a] = eta[i][0];
    a++;

}

for (int i = procSize-area; i<procSize; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    Eta1[a] = eta[i][g];
    a = a+1;
    }

}

for (int i = procSize-2*area; i<procSize-area; i++)
{

    Eta1[a] = eta[i][0];
    a++;

}

a = 0;
for (int i = 0; i<dim.xlength; i++)
{
    Etatr[a] = eta[i][0];
    a++;
}

for (int i = dim.xlength*(ysize-1); i<dim.xlength*(ysize-1)+dim.xlength; i++)
{
    Etatr[a] = eta[i][0];
    a++;
}

for (int i = (zheight-1)*area; i<(zheight-1)*area+dim.xlength; i++)
{
    Etatr[a] = eta[i][0];
    a++;
}

for (int i = (zheight-1)*area + dim.xlength*(ysize-1); i<procSize; i++)
{
    Etatr[a] = eta[i][0];
    a++;
}


//even number of z processors
if(psizez%2==0)
{

tempI = sendcount+area;
if (zpos%2!=0)
{
MPI_Send(&Eta1[0],tempI,MPI_DOUBLE,(prank-psizey), 4 + (prank-psizey),MPI_COMM_WORLD);
}

else
{
MPI_Recv(Etabrcvup,tempI,MPI_DOUBLE,(prank+psizey),4+prank,MPI_COMM_WORLD,&Status);
}

if(zpos%2==0)
{

if(zpos!=0)
{
MPI_Send(&Eta1[0],tempI,MPI_DOUBLE,(prank-psizey),4 + (prank-psizey),MPI_COMM_WORLD);
}

}

else
{
if(zpos !=(psizez-1))
{
MPI_Recv(Etabrcvup,tempI,MPI_DOUBLE,(prank+psizey), 4+prank,MPI_COMM_WORLD,&Status);
}
}

if(zpos%2==0)
{
a = sendcount+area;
MPI_Send(&Eta1[a],tempI,MPI_DOUBLE,(prank+psizey),204+(prank+psizey),MPI_COMM_WORLD);
}
else
{
MPI_Recv(Etabrcvdown,tempI,MPI_DOUBLE,(prank-psizey),204+prank,MPI_COMM_WORLD,&Status);
}

if(zpos%2!=0)
{
a = sendcount+area;
if(zpos !=(psizez-1))
{
MPI_Send(&Eta1[a],tempI,MPI_DOUBLE,(prank+psizey),204+(prank+psizey),MPI_COMM_WORLD);
}
}

else
{
if(zpos!=0)
{
MPI_Recv(Etabrcvdown,tempI,MPI_DOUBLE,(prank-psizey),204+prank,MPI_COMM_WORLD,&Status);
}
}

//sending corners
//receiving top right corner
if (zpos%2!=0)
{ if (ypos!=0){
MPI_Send(&Etatr[0],dim.xlength,MPI_DOUBLE,(prank-1-psizey), 2004 + (prank-1-psizey),MPI_COMM_WORLD);
}}

else
{
if (ypos!=(psizey-1))
{
MPI_Recv(rhorup,dim.xlength,MPI_DOUBLE,(prank+psizey+1),2004+prank,MPI_COMM_WORLD,&Status);
}
}

if(zpos%2==0)
{
if((zpos!=0)&&(ypos!=0))
{
MPI_Send(&Etatr[0],dim.xlength,MPI_DOUBLE,(prank-psizey-1),2004 + (prank-psizey-1),MPI_COMM_WORLD);
}
}

else
{
if((zpos !=(psizez-1))&&(ypos!=psizey-1))
{
MPI_Recv(rhorup,dim.xlength,MPI_DOUBLE,(prank+psizey+1), 2004+prank,MPI_COMM_WORLD,&Status);
}
}

//receiving top left corner
if (zpos%2!=0)
{if (ypos!=(psizey-1)){
a = dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank-psizey+1), 2204 + (prank-psizey+1),MPI_COMM_WORLD);
}}

else
{
if (ypos!=0)
{
MPI_Recv(rholup,dim.xlength,MPI_DOUBLE,(prank+psizey-1),2204+prank,MPI_COMM_WORLD,&Status);
}
}

if(zpos%2==0)
{
if ((zpos!=0)&&(ypos!=(psizey-1)))
{
a = dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank-psizey+1),2204 + (prank-psizey+1),MPI_COMM_WORLD);
}
}

else
{
if ((zpos !=(psizez-1))&&(ypos!=0))
{
MPI_Recv(rholup,dim.xlength,MPI_DOUBLE,(prank+psizey-1), 2204+prank,MPI_COMM_WORLD,&Status);
}
}

//receiving bottom right corner
if(zpos%2==0)
{
if (ypos!=0)
{
a = 2*dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey-1), 2404+(prank+psizey-1),MPI_COMM_WORLD);
}}

else
{if (ypos!=(psizey-1)){
MPI_Recv(rhordn,dim.xlength,MPI_DOUBLE,(prank-psizey+1),2404+prank,MPI_COMM_WORLD,&Status);
}}

if(zpos%2!=0)
{
a = 2*dim.xlength;
if ((zpos != (psizez-1))&&(ypos!=0))
{
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey-1),2404+(prank+psizey-1),MPI_COMM_WORLD);
}
}

else
{
if((zpos!=0)&&(ypos!=(psizey-1)))
{
MPI_Recv(rhordn,dim.xlength,MPI_DOUBLE,(prank-psizey+1),2404+prank,MPI_COMM_WORLD,&Status);
}
}

//receiving bottom left corner
if(zpos%2==0)
{
if(ypos!=(psizey-1)){
a = 3*dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey+1), 2604+(prank+psizey+1),MPI_COMM_WORLD);
}}

else
{
if (ypos!=0){
MPI_Recv(rholdn,dim.xlength,MPI_DOUBLE,(prank-psizey-1),2604+prank,MPI_COMM_WORLD,&Status);
}}

if(zpos%2!=0)
{
a = 3*dim.xlength;
if((zpos != (psizez-1))&&(ypos!=(psizey-1)))
{
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey+1),2604+(prank+psizey+1),MPI_COMM_WORLD);
}
}

else
{
if((zpos!=0)&&(ypos!=0))
{
MPI_Recv(rholdn,dim.xlength,MPI_DOUBLE,(prank-psizey-1),2604+prank,MPI_COMM_WORLD,&Status);
}
}

}

//for odd number of z processors
else
{
tempI = sendcount+area;
//sending etabounds
if (zpos%2!=0)
{
MPI_Send(&Eta1[0],tempI,MPI_DOUBLE,(prank-psizey), 4 + (prank-psizey),MPI_COMM_WORLD);
}

else
{
if(zpos!=(psizez-1))
{
MPI_Recv(Etabrcvup,tempI,MPI_DOUBLE,(prank+psizey),4+prank,MPI_COMM_WORLD,&Status);
}
}

if(zpos%2==0)
{
if(zpos!=0)
{
MPI_Send(&Eta1[0],tempI,MPI_DOUBLE,(prank-psizey),4 + (prank-psizey),MPI_COMM_WORLD);
}
}

else
{
MPI_Recv(Etabrcvup,tempI,MPI_DOUBLE,(prank+psizey), 4+prank,MPI_COMM_WORLD,&Status);
}

if(zpos%2==0)
{
a = sendcount+area;
if (zpos!=(psizez-1))
{
MPI_Send(&Eta1[a],tempI,MPI_DOUBLE,(prank+psizey),204+(prank+psizey),MPI_COMM_WORLD);
}
}
else
{
MPI_Recv(Etabrcvdown,tempI,MPI_DOUBLE,(prank-psizey),204+prank,MPI_COMM_WORLD,&Status);
}

if(zpos%2!=0)
{
a = sendcount+area;
{
MPI_Send(&Eta1[a],tempI,MPI_DOUBLE,(prank+psizey),204+(prank+psizey),MPI_COMM_WORLD);
}
}

else
{
if(zpos!=0)
{
MPI_Recv(Etabrcvdown,tempI,MPI_DOUBLE,(prank-psizey),204+prank,MPI_COMM_WORLD,&Status);
}
}

//sending corners
//receiving top right corner
if (zpos%2!=0)
{if (ypos!=0){
MPI_Send(&Etatr[0],dim.xlength,MPI_DOUBLE,(prank-1-psizey), 2004 + (prank-1-psizey),MPI_COMM_WORLD);
}}

else
{
if ((zpos!=(psizez-1))&&(ypos!=(psizey-1)))
{
MPI_Recv(rhorup,dim.xlength,MPI_DOUBLE,(prank+psizey+1),2004+prank,MPI_COMM_WORLD,&Status);
}
}

if(zpos%2==0)
{
if ((zpos!=0)&&(ypos!=0))
{
MPI_Send(&Etatr[0],dim.xlength,MPI_DOUBLE,(prank-psizey-1),2004 + (prank-psizey-1),MPI_COMM_WORLD);
}
}

else
{
if(ypos!=(psizey-1)){
MPI_Recv(rhorup,dim.xlength,MPI_DOUBLE,(prank+psizey+1), 2004+prank,MPI_COMM_WORLD,&Status);
}}

//receiving top left corner
if (zpos%2!=0)
{
if (ypos!=(psizey-1)){
a = dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank-psizey+1), 2204 + (prank-psizey+1),MPI_COMM_WORLD);
}}

else
{
if ((zpos!=(psizez-1))&&(ypos!=0))
{
MPI_Recv(rholup,dim.xlength,MPI_DOUBLE,(prank+psizey-1),2204+prank,MPI_COMM_WORLD,&Status);
}
}

if(zpos%2==0)
{
if ((zpos!=0)&&(ypos!=(psizey-1)))
{
a = dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank-psizey+1),2204 + (prank-psizey+1),MPI_COMM_WORLD);
}
}

else
{if (ypos!=0){
MPI_Recv(rholup,dim.xlength,MPI_DOUBLE,(prank+psizey-1), 2204+prank,MPI_COMM_WORLD,&Status);
}}

//receiving bottom right corner
if(zpos%2==0)
{
a = 2*dim.xlength;
if((zpos!=(psizez-1))&&(ypos!=0))
{
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey-1), 2404+(prank+psizey-1),MPI_COMM_WORLD);
}
}

else
{if (ypos!=(psizey-1)){
MPI_Recv(rhordn,dim.xlength,MPI_DOUBLE,(prank-psizey+1),2404+prank,MPI_COMM_WORLD,&Status);
}}

if(zpos%2!=0)
{ if (ypos!=0){
a = 2*dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey-1),2404+(prank+psizey-1),MPI_COMM_WORLD);
}}

else
{
if((zpos!=0)&&(ypos!=(psizey-1)))
{
MPI_Recv(rhordn,dim.xlength,MPI_DOUBLE,(prank-psizey+1),2404+prank,MPI_COMM_WORLD,&Status);
}
}

//receiving bottom left corner
if(zpos%2==0)
{
if ((zpos!=(psizez-1))&&(ypos!=(psizey-1)))
{
a = 3*dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey+1), 2604+(prank+psizey+1),MPI_COMM_WORLD);
}
}

else
{ if (ypos!=0){
MPI_Recv(rholdn,dim.xlength,MPI_DOUBLE,(prank-psizey-1),2604+prank,MPI_COMM_WORLD,&Status);
}}

if(zpos%2!=0)
{if(ypos!=(psizey-1)){
a = 3*dim.xlength;
MPI_Send(&Etatr[a],dim.xlength,MPI_DOUBLE,(prank+psizey+1),2604+(prank+psizey+1),MPI_COMM_WORLD);
}}

else
{
if((zpos!=0)&&(ypos!=0))
{
MPI_Recv(rholdn,dim.xlength,MPI_DOUBLE,(prank-psizey-1),2604+prank,MPI_COMM_WORLD,&Status);
}
}
}

clear1D(Eta1);
clear1D(Etatr);

Eta1y = store1D(2*(sendcounty+areay));
a = 0;

for (int z = 0; z<zheight; z++)
{
    for (int x = 0; x<dim.xlength; x++)
    {
       int  k = x + (z*area);
        for (int g = 0; g<=dim.noOfparticles; g++)
        {
        Eta1y[a] = eta[k][g];
        a = a+1;
        }

    }
}

for (int z = 0; z<zheight; z++)
{
    for (int x = 0; x<dim.xlength; x++)
    {

    int k = x +dim.xlength+ (z*area);
    Eta1y[a] = eta[k][0];
    a = a+1;

    }
}

for (int z = 0; z<zheight; z++)
{
    for (int x = 0; x<dim.xlength; x++)
    {
       int  k = x + (z*area) + (ysize-1)*dim.xlength;
        for (int g = 0; g<=dim.noOfparticles; g++)
        {
        Eta1y[a] = eta[k][g];
        a = a+1;
        }

    }
}

for (int z = 0; z<zheight; z++)
{
    for (int x = 0; x<dim.xlength; x++)
    {

     int k = x + (z*area) + (ysize-2)*dim.xlength;
        Eta1y[a] = eta[k][0];
        a = a+1;

    }
}

//even number of y processors
if(psizey%2==0)
{

tempI = sendcounty+areay;
if (ypos%2!=0)
{
MPI_Send(&Eta1y[0],tempI,MPI_DOUBLE,(prank-1), 1004 + (prank-1),MPI_COMM_WORLD);
}

else
{
MPI_Recv(Etabrcvright,tempI,MPI_DOUBLE,(prank+1),1004+prank,MPI_COMM_WORLD,&Status);
}

if(ypos%2==0)
{

if(ypos!=0)
{
MPI_Send(&Eta1y[0],tempI,MPI_DOUBLE,(prank-1),1004 + (prank-1),MPI_COMM_WORLD);
}

}

else
{
if(ypos!=(psizey-1))
{
MPI_Recv(Etabrcvright,tempI,MPI_DOUBLE,(prank+1), 1004+prank,MPI_COMM_WORLD,&Status);
}
}

if(ypos%2==0)
{
a = tempI;
MPI_Send(&Eta1y[a],tempI,MPI_DOUBLE,(prank+1),1204+(prank+1),MPI_COMM_WORLD);
}
else
{
MPI_Recv(Etabrcvleft,tempI,MPI_DOUBLE,(prank-1),1204+prank,MPI_COMM_WORLD,&Status);
}

if(ypos%2!=0)
{
a = tempI;
if(ypos!=(psizey-1))
{
MPI_Send(&Eta1y[a],tempI,MPI_DOUBLE,(prank+1),1204+(prank+1),MPI_COMM_WORLD);
}
}

else
{
if(ypos!=0)
{
MPI_Recv(Etabrcvleft,tempI,MPI_DOUBLE,(prank-1),1204+prank,MPI_COMM_WORLD,&Status);
}
}
}


//odd number of y processors
else
{

tempI = sendcounty+areay;
//sending etabounds
if (ypos%2!=0)
{
MPI_Send(&Eta1y[0],tempI,MPI_DOUBLE,(prank-1), 1004 + (prank-1),MPI_COMM_WORLD);
}

else
{
if(ypos!=(psizey-1))
{
MPI_Recv(Etabrcvright,tempI,MPI_DOUBLE,(prank+1),1004+prank,MPI_COMM_WORLD,&Status);
}
}

if(ypos%2==0)
{
if(ypos!=0)
{
MPI_Send(&Eta1y[0],tempI,MPI_DOUBLE,(prank-1),1004 + (prank-1),MPI_COMM_WORLD);
}
}

else
{
MPI_Recv(Etabrcvright,tempI,MPI_DOUBLE,(prank+1), 1004+prank,MPI_COMM_WORLD,&Status);
}

if(ypos%2==0)
{
a = tempI;
if (ypos!=(psizey-1))
{
MPI_Send(&Eta1y[a],tempI,MPI_DOUBLE,(prank+1),1204+(prank+1),MPI_COMM_WORLD);
}
}
else
{
MPI_Recv(Etabrcvleft,tempI,MPI_DOUBLE,(prank-1),1204+prank,MPI_COMM_WORLD,&Status);
}

if(ypos%2!=0)
{
a = tempI;
{
MPI_Send(&Eta1y[a],tempI,MPI_DOUBLE,(prank+1),1204+(prank+1),MPI_COMM_WORLD);
}
}

else
{
if(ypos!=0)
{
MPI_Recv(Etabrcvleft,tempI,MPI_DOUBLE,(prank-1),1204+prank,MPI_COMM_WORLD,&Status);
}
}

}

clear1D(Eta1y);

//if (prank==chosenrank){cout<<"2. Done with Sending and Receiving from processors\n";}
//timeME(prank, chosenrank, start);

//-------------------------------------------------------------//
//>>end of sending and receiving data from other processors<<//
//>>Transferring data from buffers --> preparing for analysis<<//
//-------------------------------------------------------------//

Etaboundsup = store2D(area, dim.noOfparticles+1);
Etaboundsdown = store2D(area, dim.noOfparticles+1);
Etaboundsright = store2D(areay, dim.noOfparticles+1);
Etaboundsleft = store2D(areay, dim.noOfparticles+1);

int countb = 0;
if (zpos!=0)
{
for (int i = 0; i< area; i++)
{
    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    tempEN = Etabrcvdown[countb];
    Etaboundsdown[i][g] = tempEN;
    countb  = countb +1;
    }

    tempE = Etabrcvdown[sendcount+i];
    rhobounds2down[i] = tempE;
}
}

else
{
for (int i = 0; i< area; i++)
{
    Etaboundsdown[i][0] = dim.rhoVap;
    rhobounds2down[i] = dim.rhoVap;

    for (int g = 1; g<=dim.noOfparticles; g++)
    {
    Etaboundsdown[i][g] = dim.etaVap;
    }
}
}

countb = 0;
if(ypos!=0)
{
for (int i = 0; i<areay; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    tempEN = Etabrcvleft[countb];
    Etaboundsleft[i][g] = tempEN;
    countb = countb+1;
    }

    tempE = Etabrcvleft[sendcounty+i];
    rhobounds2left[i] = tempE;
}
}

else
{
for (int i = 0; i< areay; i++)
{
    Etaboundsleft[i][0] = dim.rhoVap;
    rhobounds2left[i] = dim.rhoVap;

    for (int g = 1; g<=dim.noOfparticles; g++)
    {
    Etaboundsleft[i][g] = dim.etaVap;
    }
}
}

countb = 0;
if (zpos!=(psizez-1))
{
for (int i = 0; i< area; i++)
{
    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    temp = Etabrcvup[countb];
    Etaboundsup[i][g] = temp;
    countb  = countb +1;
    }

    tempE = Etabrcvup[sendcount+i];
    rhobounds2up[i] = tempE;
}
}

else
{
for (int i = 0; i< area; i++)
{
    Etaboundsup[i][0] = dim.rhoVap;
    rhobounds2up[i] = dim.rhoVap;

    for (int g = 1; g<=dim.noOfparticles; g++)
    {
    Etaboundsup[i][g] = dim.etaVap;
    }
}
}

countb = 0;
if (ypos != (psizey-1))
{
for (int i = 0; i< areay; i++)
{
    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    temp = Etabrcvright[countb];
    Etaboundsright[i][g] = temp;
    countb  = countb +1;
    }

    tempE = Etabrcvright[sendcounty+i];
    rhobounds2right[i] = tempE;
}
}

else
{
for (int i = 0; i< areay; i++)
{
    Etaboundsright[i][0] = dim.rhoVap;
    rhobounds2right[i] = dim.rhoVap;

    for (int g = 1; g<=dim.noOfparticles; g++)
    {
    Etaboundsright[i][g] = dim.etaVap;
    }
}
}

if (zpos==0)
{
for (int x = 0; x< dim.xlength; x++)
    {
        rhordn[x] = dim.rhoVap;
        rholdn[x] = dim.rhoVap;
    }
}

if (zpos==(psizez-1))
{
for (int x = 0; x< dim.xlength; x++)
    {
        rhorup[x] = dim.rhoVap;
        rholup[x] = dim.rhoVap;
    }
}

if (ypos==0)
{
for (int x = 0; x< dim.xlength; x++)
    {
        rholup[x] = dim.rhoVap;
        rholdn[x] = dim.rhoVap;
    }
}

if (ypos==(psizey-1))
{
for (int x = 0; x< dim.xlength; x++)
    {
        rhorup[x] = dim.rhoVap;
        rhordn[x] = dim.rhoVap;
    }
}

//freeing memory in uneeded buffer

if (zpos!=0)
{
clear1D(Etabrcvdown);
}
if (zpos!=(psizez-1))
{
clear1D(Etabrcvup);
}
if (ypos!=0)
{
clear1D(Etabrcvleft);
}
if (ypos!=(psizey-1))
{
clear1D(Etabrcvright);
}

//----------------------------------------------------//
//>>Preparing y-z bounds on simulation for analysis<<//
//----------------------------------------------------//

//transferring rho and eta data into NEW array for faster analysis
Eta = store2D(dim.noOfparticles+1,procSize);

for (int i = 0; i<procSize; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    Eta[g][i] = eta[i][g];
    }

}

//allocating memory needed for filling bounds
rhot = store1D(area*3);
rhob = store1D(area*3);

// filling the top and bottom of the sample space (i.e side three z stacks)
for(int z = 0; z < 3; z++)
{ for(int y = 0; y < ysize; y++)
{ for(int x = 0; x < dim.xlength; x++)
{
int l = (area*z)+(y*dim.xlength) + x;
int i = (y*dim.xlength) + x;
int j = ((zheight-1)*area)+i;

if (z ==0)
{
tempn = Eta[0][i];
rhob[l] = tempn;

tempE = Eta[0][j];
rhot[l] = tempE;
}

else if (z ==1)
{
tempn = Etaboundsup[i][0];
rhot[l] = tempn;

tempE = Etaboundsdown[i][0];
rhob[l] = tempE;
}

else
{
tempn = rhobounds2up[i];
rhot[l] = tempn;

tempE = rhobounds2down[i];
rhob[l] = tempE;
}
}
}
}

//freeing uneeded memory buffers
clear1D(rhobounds2down);
clear1D(rhobounds2up);

//allocating memory needed for filling bounds
rhor = store1D(areay*3);
rhol = store1D(areay*3);


// filling the left and right of the sample space (i.e side three y stacks)
for(int z = 0; z < zheight; z++)
{ for(int y = 0; y < 3; y++)
{ for(int x = 0; x < dim.xlength; x++)
{
int l = (dim.xlength*3*z)+(y*dim.xlength) + x;
int i = (z*dim.xlength) + x;
int bl = x + (z*area);
int br = x + (z*area) + (ysize-1)*dim.xlength;

if (y ==0)
{
tempn = Eta[0][br];
rhor[l] = tempn;

tempE = Eta[0][bl];
rhol[l] = tempE;
}

else if (y ==1)
{
tempn = Etaboundsright[i][0];
rhor[l] = tempn;

tempE = Etaboundsleft[i][0];
rhol[l] = tempE;
}

else
{
tempn = rhobounds2right[i];
rhor[l] = tempn;

tempE = rhobounds2left[i];
rhol[l] = tempE;
}
}
}
}

//freeing uneeded memory buffers
clear1D(rhobounds2left);
clear1D(rhobounds2right);

//allocating memory needed for filling bounds
ht = store1D(area);
hb = store1D(area);
hr = store1D(areay);
hl = store1D(areay);
Dift = store1D(area);
Difb = store1D(area);
Difr = store1D(areay);
Difl = store1D(areay);

//ht -> first derivative term related to rho at the topmost edge of simulation box
for(int y = 0; y < ysize; y++)
{ for(int x = 0; x < dim.xlength; x++)
{
int a = (y*dim.xlength) + x;
int j = area + a;

etaSqSum = 0.0;
etaCubeSum = 0.0;
etaij = 0.0;

rhoSq = rhot[j]*rhot[j];
rhoQuad = rhoSq*rhoSq;

for (int g = 1; g<=dim.noOfparticles; g++)
{
//for (int gh = g+1; gh<=dim.noOfparticles; gh++) {etaij = etaij+Etaboundsup[a][g]*Etaboundsup[a][gh];}
etaSqSum = etaSqSum + (Etaboundsup[a][g]*Etaboundsup[a][g]);
etaCubeSum = etaCubeSum + (Etaboundsup[a][g]*Etaboundsup[a][g]*Etaboundsup[a][g]);
}

BFE_partialRhot = A2*rhot[j]*(1-rhot[j])*(1-2*rhot[j]) + con.B*((2*rhot[j])-(6*etaSqSum)+(4*etaCubeSum));

phi = rhoQuad*((7*rhoSq)-(18*rhot[j])+12);
Dift[a] = (con.voldif*phi)+(con.surfdif*rhoSq*(1-rhot[j])*(1-rhot[j]))+(con.gbdif*rhot[j]*(1-etaSqSum));

rho4 = rhot[N_proc[a][4]]; rho5 = rhot[N_proc[a][5]];
rho2 = rhot[N_proc[a][2]]; rho3 = rhot[N_proc[a][3]];
rho0 = rhot[N_proc[a][0]]; rho1 = rhot[N_proc[a][1]];

if (y == 0){rho3 = rholup[x];}
if (y == ysize-1){rho2 = rhorup[x];}
if (x == 0){rho1 = dim.rhoVap;}
if (x == dim.xlength-1){rho0 = dim.rhoVap;}


ht[a] = BFE_partialRhot - con.surfenergy*(rho0+rho1+rho2+rho3+rho4+rho5-(6*rhot[j]));

}
}

//hb -> first derivative term related to rho at the bottommost edge of simulation box

for(int y = 0; y < ysize; y++)
{ for(int x = 0; x < dim.xlength; x++)
{
int a = (y*dim.xlength) + x;
int j = area + a;

etaSqSum = 0.0;
etaCubeSum = 0.0;
etaij = 0.0;

rhoSq = rhob[j]*rhob[j];
rhoQuad = rhoSq*rhoSq;

for (int g = 1; g<=dim.noOfparticles; g++)
{
//for (int gh = g+1; gh<=dim.noOfparticles; gh++) {etaij = etaij+Etaboundsdown[a][g]*Etaboundsdown[a][gh];}
etaSqSum = etaSqSum + (Etaboundsdown[a][g]*Etaboundsdown[a][g]);
etaCubeSum = etaCubeSum + (Etaboundsdown[a][g]*Etaboundsdown[a][g]*Etaboundsdown[a][g]);
}

BFE_partialRhot = A2*rhob[j]*(1-rhob[j])*(1-2*rhob[j]) + con.B*((2*rhob[j])-(6*etaSqSum)+(4*etaCubeSum));

phi = rhoQuad*((7*rhoSq)-(18*rhob[j])+12);
Difb[a] = (con.voldif*phi)+(con.surfdif*rhoSq*(1-rhob[j])*(1-rhob[j]))+(con.gbdif*rhob[j]*(1-etaSqSum));

rho4 = rhob[N_proc[a][4]]; rho5 = rhob[N_proc[a][5]];
rho2 = rhob[N_proc[a][2]]; rho3 = rhob[N_proc[a][3]];
rho0 = rhob[N_proc[a][0]]; rho1 = rhob[N_proc[a][1]];

if (y == 0){rho3 = rholdn[x];}
if (y == ysize-1){rho2 = rhordn[x];}
if (x == 0){rho1 = dim.rhoVap;}
if (x == dim.xlength-1){rho0 = dim.rhoVap;}


hb[a] = BFE_partialRhot - con.surfenergy*(rho0+rho1+rho2+rho3+rho4+rho5-(6*rhob[j]));

}
}

//hr -> first derivative term related to rho at the rightmost edge of simulation box
for(int z = 0; z < zheight; z++)
{ for(int x = 0; x < dim.xlength; x++)
{
int a = (z*dim.xlength) + x;
int j = z*3*dim.xlength + dim.xlength + x;

etaSqSum = 0.0;
etaCubeSum = 0.0;
etaij = 0.0;
rhoSq = rhor[j]*rhor[j];
rhoQuad = rhoSq*rhoSq;

for (int g = 1; g<=dim.noOfparticles; g++)
{
//for (int gh = g+1; gh<=dim.noOfparticles; gh++) {etaij = etaij+Etaboundsright[a][g]*Etaboundsright[a][gh];}
etaSqSum = etaSqSum + (Etaboundsright[a][g]*Etaboundsright[a][g]);
etaCubeSum = etaCubeSum + (Etaboundsright[a][g]*Etaboundsright[a][g]*Etaboundsright[a][g]);
}

BFE_partialRhot = A2*rhor[j]*(1-rhor[j])*(1-2*rhor[j]) + con.B*((2*rhor[j])-(6*etaSqSum)+(4*etaCubeSum));

phi = rhoQuad*((7*rhoSq)-(18*rhor[j])+12);
Difr[a] = (con.voldif*phi)+(con.surfdif*rhoSq*(1-rhor[j])*(1-rhor[j]))+(con.gbdif*rhor[j]*(1-etaSqSum));

rho4 = rhor[Ny_proc[a][4]]; rho5 = rhor[Ny_proc[a][5]];
rho2 = rhor[Ny_proc[a][2]]; rho3 = rhor[Ny_proc[a][3]];
rho0 = rhor[Ny_proc[a][0]]; rho1 = rhor[Ny_proc[a][1]];

if (z == 0){rho5 = rhordn[x];}
if (z == zheight-1){rho4 = rhorup[x];}
if (x == 0){rho1 = dim.rhoVap;}
if (x == dim.xlength-1){rho0 = dim.rhoVap;}

hr[a] = BFE_partialRhot - con.surfenergy*(rho0+rho1+rho2+rho3+rho4+rho5-(6*rhor[j]));

}
}

//hl -> first derivative term related to rho at the leftmost edge of simulation box
for(int z = 0; z < zheight; z++)
{ for(int x = 0; x < dim.xlength; x++)
{
int a = (z*dim.xlength) + x;
int j = z*3*dim.xlength + dim.xlength + x;

etaSqSum = 0.0;
etaCubeSum = 0.0;
etaij = 0.0;

rhoSq = rhol[j]*rhol[j];
rhoQuad = rhoSq*rhoSq;

for (int g = 1; g<=dim.noOfparticles; g++)
{
//for (int gh = g+1; gh<=dim.noOfparticles; gh++) {etaij = etaij+Etaboundsleft[a][g]*Etaboundsleft[a][gh];}
etaSqSum = etaSqSum + (Etaboundsleft[a][g]*Etaboundsleft[a][g]);
etaCubeSum = etaCubeSum + (Etaboundsleft[a][g]*Etaboundsleft[a][g]*Etaboundsleft[a][g]);
}

BFE_partialRhot = A2*rhol[j]*(1-rhol[j])*(1-2*rhol[j]) + con.B*((2*rhol[j])-(6*etaSqSum)+(4*etaCubeSum));

phi = rhoQuad*((7*rhoSq)-(18*rhol[j])+12);
Difl[a] = (con.voldif*phi)+(con.surfdif*rhoSq*(1-rhol[j])*(1-rhol[j]))+(con.gbdif*rhol[j]*(1-etaSqSum));

rho4 = rhol[Ny_proc[a][4]]; rho5 = rhol[Ny_proc[a][5]];
rho2 = rhol[Ny_proc[a][2]]; rho3 = rhol[Ny_proc[a][3]];
rho0 = rhol[Ny_proc[a][0]]; rho1 = rhol[Ny_proc[a][1]];

if (z == 0){rho5 = rholdn[x];}
if (z == zheight-1){rho4 = rholup[x];}
if (x == 0){rho1 = dim.rhoVap;}
if (x == dim.xlength-1){rho0 = dim.rhoVap;}

hl[a] = BFE_partialRhot - con.surfenergy*(rho0+rho1+rho2+rho3+rho4+rho5-(6*rhol[j]));

}
}

//freeing uneeded memory buffers
clear1D(rhorup);
clear1D(rhordn);
clear1D(rholup);
clear1D(rholdn);
clear1D(rhot);
clear1D(rhob);
clear1D(rhor);
clear1D(rhol);


//if (prank==chosenrank){cout<<"3. Done with preparing the parts for the bounds\n";}
//timeME(prank, chosenrank, start);

//----------------------------------------------------//
//----------------------------------------------------//
//>>!!SINTERING OF SIMULATION BOX IS IN THIS LOOP!!<<//
//----------------------------------------------------//
//----------------------------------------------------//

//alocating memory needed for analysis
NewEta = store2D(dim.noOfparticles, procSize);
NewRho = store1D(procSize);
Dif = store1D(procSize);
h = store1D(procSize);
BFE_partialEta = store1D(dim.noOfparticles);

//calculating the h and Diffusion variables
for(int z = 0; z < zheight; z++)
{ for(int y = 0; y < ysize; y++)
{ for(int x = 0; x < dim.xlength; x++)
{
int a = (area*z)+(y*dim.xlength) + x;
int xy = (y*dim.xlength) + x;
int xz = (z*dim.xlength) + x;

etaSqSum = 0.0;
etaCubeSum = 0.0;
etaij = 0.0;

rhoSq = Eta[0][a]*Eta[0][a];
rhoQuad = rhoSq*rhoSq;

for (int g = 1; g<=dim.noOfparticles; g++)
{
//for (int gh =g+1; gh<=dim.noOfparticles; gh++) {etaij = etaij+Eta[g][a]*Eta[gh][a];}
etaSqSum = etaSqSum + (Eta[g][a]*Eta[g][a]);
etaCubeSum = etaCubeSum + (Eta[g][a]*Eta[g][a]*Eta[g][a]);
}

BFE_partialRho = A2*Eta[0][a]*(1-Eta[0][a])*(1-2*Eta[0][a]) + con.B*((2*Eta[0][a])-(6*etaSqSum)+(4*etaCubeSum));

phi = rhoQuad*((7*rhoSq)-(18*Eta[0][a])+12);
Dif[a] = (con.voldif*phi)+(con.surfdif*rhoSq*(1-Eta[0][a])*(1-Eta[0][a]))+(con.gbdif*Eta[0][a]*(1-etaSqSum));

rho4 = Eta[0][N[a][4]]; rho5 = Eta[0][N[a][5]]; //z
rho2 = Eta[0][N[a][2]]; rho3 = Eta[0][N[a][3]]; //y
rho0 = Eta[0][N[a][0]]; rho1 = Eta[0][N[a][1]]; //x

if (z == 0){rho5 = Etaboundsdown[xy][0];}
if (z == zheight-1){rho4 = Etaboundsup[xy][0];}
if (y == 0){rho3 = Etaboundsleft[xz][0];}
if (y == ysize-1){rho2 = Etaboundsright[xz][0];}
if (x == 0){rho1 = dim.rhoVap;}
if (x == dim.xlength-1){rho0 = dim.rhoVap;}

h[a] = BFE_partialRho - con.surfenergy*(rho0+rho1+rho2+rho3+rho4+rho5-(6*Eta[0][a]));

}
}
}

//analysis loop broken into: interior boxes and the edges of the simulation to find the newrho and eta values
for(int z = 0; z < zheight; z++)
{ for(int y = 0; y < ysize; y++)
{ for(int x = 0; x < dim.xlength; x++)
{
int a = (area*z)+(y*dim.xlength) + x;

int xy = (y*dim.xlength) + x;
int xz = (z*dim.xlength) + x;

etaSqSum = 0.0;

for (int g = 1; g<=dim.noOfparticles; g++)
{
etaSqSum = etaSqSum + (Eta[g][a]*Eta[g][a]);
}

for (int g = 1; g<=dim.noOfparticles; g++)
{
BFE_partialEta[g-1] = B12*Eta[g][a]*(1-Eta[0][a] - (2-Eta[0][a])*Eta[g][a] + etaSqSum);
}

Dif4 = Dif[N[a][4]]; Dif5 = Dif[N[a][5]];
Dif2 = Dif[N[a][2]]; Dif3 = Dif[N[a][3]];
Dif0 = Dif[N[a][0]]; Dif1 = Dif[N[a][1]]; 

h4 = h[N[a][4]]; h5 = h[N[a][5]];
h2 = h[N[a][2]]; h3 = h[N[a][3]];
h0 = h[N[a][0]]; h1 = h[N[a][1]];

if (z == 0)
{
    Dif5 = Difb[xy];
    h5 = hb[xy];
}
if (z == zheight-1)
{
    Dif4 = Dift[xy];
    h4 = ht[xy];
}
if (y == 0)
{
    Dif3 = Difl[xz];
    h3 = hl[xz];
}
if (y == ysize-1)
{
    Dif2 = Difr[xz];
    h2 = hr[xz];
}
if (x == 0)
{
    rho1 = dim.rhoVap; eta1 = dim.etaVap;

    rhoSq = rho1*rho1;
    rhoQuad = rhoSq*rhoSq;

    etaSqSum = dim.noOfparticles*eta1*eta1;
    etaCubeSum = etaSqSum*eta1;

    BFE_partialRho = A2*rho1*(1-rho1)*(1-2*rho1) + con.B*((2*rho1)-(6*etaSqSum)+(4*etaCubeSum));

    phi = rhoQuad*((7*rhoSq)-(18*rho1)+12);
    Dif1 = (con.voldif*phi)+(con.surfdif*rhoSq*(1-rho1)*(1-rho1))+(con.gbdif*rho1*(1-etaSqSum));

    h1 = BFE_partialRho;
}
if (x == dim.xlength-1)
{
    rho0 = dim.rhoVap; eta0 = dim.etaVap;

    rhoSq = rho0*rho0;
    rhoQuad = rhoSq*rhoSq;

    etaSqSum = dim.noOfparticles*eta0*eta0;
    etaCubeSum = etaSqSum*eta0;

    BFE_partialRho = A2*rho0*(1-rho0)*(1-2*rho0) + con.B*((2*rho0)-(6*etaSqSum)+(4*etaCubeSum));

    phi = rhoQuad*((7*rhoSq)-(18*rho0)+12);
    Dif0 = (con.voldif*phi)+(con.surfdif*rhoSq*(1-rho0)*(1-rho0))+(con.gbdif*rho0*(1-etaSqSum));

    h0 = BFE_partialRho;
}

//rho(t+1) after the time step the new variable is NewRho. Same goes for eta- eta(t+1)-NewEta
NewRho[a] = Eta[0][a] + con.deltaT*0.5*((Dif0*(h0-h[a]))-(Dif1*(h[a]-h1))+(Dif2*(h2-h[a]))-(Dif3*(h[a]-h3))+(Dif4*(h4-h[a]))-(Dif5*(h[a]-h5))+(Dif[a]*(h0+h1+h2+h3+h4+h5-(6*h[a]))));

for (int g = 1; g<=dim.noOfparticles; g++)
{
    rho4 = Eta[g][N[a][4]]; rho5 = Eta[g][N[a][5]];
    rho2 = Eta[g][N[a][2]]; rho3 = Eta[g][N[a][3]];
    rho0 = Eta[g][N[a][0]]; rho1 = Eta[g][N[a][1]];

    if (z == 0){rho5 = Etaboundsdown[xy][g];}
    if (z == zheight-1){rho4 = Etaboundsup[xy][g];}
    if (y == 0){rho3 = Etaboundsleft[xz][g];}
    if (y == ysize-1){rho2 = Etaboundsright[xz][g];}
    if (x == 0){rho1 = dim.etaVap;}
    if (x == dim.xlength-1){rho0 = dim.etaVap;}

    NewEta[g-1][a] = Eta[g][a] - con.gbmobility*con.deltaT*(BFE_partialEta[g-1]-(con.gbenergy*(rho0+rho1+rho2+rho3+rho4+rho5-(6*Eta[g][a]))));
}

}
}
}

//if (prank==chosenrank){cout<<"4. Done with sintering in main for loop\n";}
//timeME(prank, chosenrank, start);

//freeing uneeded memory
clear2D(Etaboundsup, area);
clear2D(Etaboundsdown, area);
clear2D(Etaboundsright, areay);
clear2D(Etaboundsleft, areay);
clear1D(ht);
clear1D(hb);
clear1D(hr);
clear1D(hl);
clear1D(Dift);
clear1D(Difb);
clear1D(Difr);
clear1D(Difl);

clear1D(Dif);
clear1D(h);
clear1D(BFE_partialEta);

//----------------------------------------------------//
////>>Two limit scheme ensuring density conservation<<//
////----------------------------------------------------//
//transfering data from old to new time step and restricting field variable values to a minimum of the cut off

for (int a = 0; a<procSize; a++)
{

    tempn = NewRho[a];

    if (tempn < cutoffR)
        {
        	Eta[0][a]=dim.rhoVap;
        }
    else
        {
        	Eta[0][a] = tempn;
        }

     for (int g = 1; g<=dim.noOfparticles; g++)
    {
	tempEN = NewEta[g-1][a];
        if (tempEN < cutoffR)
            {
            	Eta[g][a]=dim.etaVap;
            }
        else
            {
    	    	Eta[g][a] = tempEN;
            }
    }

}

//freeing uneeded memory
clear2D(NewEta,dim.noOfparticles);
clear1D(NewRho);

temprho = 0.0;

for (int i = 0; i<procSize; i++)
{
    temprho = temprho + Eta[0][i];
}

tempbv = 0;

for (int i = 0; i<procSize; i++)
{
    if (Eta[0][i] > limR)
    {
	tempbv = tempbv + 1;
    }
}

MPI_Reduce(&temprho,&sumrhotot,1,MPI_DOUBLE,MPI_SUM,pid,MPI_COMM_WORLD);
MPI_Reduce(&tempbv,&sizebvtot,1,MPI_INT,MPI_SUM,pid,MPI_COMM_WORLD);
MPI_Bcast(&sumrhotot,1,MPI_DOUBLE,pid,MPI_COMM_WORLD);
MPI_Bcast(&sizebvtot,1,MPI_INT,pid,MPI_COMM_WORLD);

intvrho = (sumrhototInit-sumrhotot)/sizebvtot;

if ((sumrhotot/sumrhototInit)>0 && (sumrhotot/sumrhototInit)<100)
{
simcheckV= 0;
}
else
{
simcheckV= 1;
}

/*
if (prank==pid)
{
cout<<sumrhototInit<<"\t"<<sumrhotot<<"\t"<<sizebvtot<<"\n";
}
*/

//this loop SHOULD NOT be commented out, this is the enforcement of conservaton laws
for (int i = 0; i<procSize; i++)
{
    if (Eta[0][i] > limR)
    {
    tempE = Eta[0][i];
    Eta[0][i] = tempE + intvrho;
    }
}

/*
sumrho = 0.0;

for (int i = 0; i<procSize; i++)
{
    sumrho = sumrho +Eta[0][i];
}
*/

//....................//
//>>Printing to file<<//
//....................//

evat = intvs + countt*intvs;
evatf = intvf + counttf*intvf;

if (noisy == 1 && prank == 0)
{
cout<<sumrhototInit<<"\t"<<sumrhotot<<"\t"<<sizebvtot<<"\n";
cout<<"time: "<<t<<"\n";
}

//Printing only rho or sum of rho and all eta variables (or both if desired)
if ((evat>t-con.deltaT)&&(evat<t+con.deltaT))
{

if (prank==0)
{
cout<<sumrhototInit<<"\t"<<sumrhotot<<"\t"<<sizebvtot<<"\n";
cout<<"time: "<<t<<"\n";
}

//MPIfileprintsum(t,evat,con.deltaT,Eta,"3p","0.0001.dat",100*evat,procSize,prank,dim,zpos,ypos,zheight,ysize);
MPIfileprintrho(t,evat,con.deltaT,Eta,"rho",".dat",100*evat,SNo,procSize,prank,dim,zpos,ypos,zheight,ysize);
countt = countt+1;

}

//Printing full data i.e rho and eta
if ((evatf>t-con.deltaT)&&(evatf<t+con.deltaT))
{
//MPIfileprintsum(t,evat,con.deltaT,Eta,"3p","0.0001.dat",100*evat,procSize,prank,dim,zpos,ypos,zheight,ysize);
MPIprintfullE(t,evatf,con.deltaT,Eta,"fullT",".dat",10*evatf,SNo,procSize,prank,dim,zpos,ypos,zheight,ysize);
counttf = counttf+1;
}
if (printEND==1)
{
if ((finalt-con.deltaT>t-con.deltaT/2)&&(finalt-con.deltaT<t+con.deltaT/2))
{
MPIprintfullE(t,finalt-con.deltaT,con.deltaT,Eta,"fullT",".dat",10*finalt-con.deltaT,SNo,procSize,prank,dim,zpos,ypos,zheight,ysize);
}
}

//checking simulation stability
//if any processor is unstable 
//all processors are sent to the final time
//and simulation closes

if (simcheckV > 0)
{
//cout<<prank<<"\t"<<sumrhototInit<<"\t"<<sumrhotot<<"\t"<<sizebvtot<<"\n";
if (prank==1)
{
cout<<sumrhototInit<<"\t"<<sumrhotot<<"\t"<<sizebvtot<<"\n";
cout<<"time: "<<t<<"\n-->Simulation unstable: Prematurely ended!<--\n";
}
t = finalt;
}


for (int i = 0; i< procSize; i++)
{
    for (int g = 0; g<=dim.noOfparticles; g++)
    {
    eta[i][g] = dim.etaVap;
    }

}

for (int i = 0; i< procSize; i++)
{

    for (int g = 0; g<=dim.noOfparticles; g++)
    {

    eta[i][g] = Eta[g][i];
    }

}

/*
if (prank==1)
{
cout<<"time: "<<t<<"\n";
}
*/

//freeing rho and eta memory
clear2D(Eta, dim.noOfparticles+1);

t = t + con.deltaT;

//if (prank==chosenrank){cout<<"5. Done with post processing. End of while loop\n";}
//timeME(prank, chosenrank, start);

}

clear2Dint(N, procSize);
clear2Dint(N_proc, area);
clear2Dint(Ny_proc, areay);

if (prank==chosenrank){cout<<"after\nwhile\n";}
MPI_Finalize ();

}

