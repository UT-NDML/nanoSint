#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <new>
#include <iomanip>
#include <ctime>

using namespace std;

const double err = 0.0000000001;
const double pi = 2*asin(1);

/*ARRAY MEMORY ALLOCATION*/
//memory allocation for 1D double array
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

//memory allocation for 1D integer array
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

//memory allocation for 2D double array
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

//memory allocation for 2D integer array
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

/*ARRAY MEMORY DEALLOCATION [CLEARING MEMORY]*/
//clear 1D double array
void clear1D(double* j1)
{
delete [] j1;
return;
}

//clear 1D integer array
void clear1Dint(int* j1)
{
delete [] j1;
return;
}

//clear 2D double array
void clear2D(double **j2, int a)
{
for (int i=0;i<a;i++)
{
delete [] j2 [i];
}
delete [] j2;
return;
}

//clear 2D integer array
void clear2Dint(int **j2, int a)
{
for (int i=0;i<a;i++)
{
delete [] j2 [i];
}
delete [] j2;
return;
}

//print 2d variable to a text file
void print2Dtotxtfile(string file1,string file2, double **E, int rowSize, int colSize, double timet)
{
    ostringstream oss;
    string var;
    char *filename;

    oss<<file1<<timet<<file2;

    var = oss.str();

    filename = new char [var.length()];

    strcpy(filename, var.c_str());

    ofstream myfile;

    myfile.open(filename);

    for (int ir = 0; ir<rowSize; ir++)
    {
        for (int ic = 0; ic<colSize; ic++)
        {
            myfile<<E[ir][ic]<<"\t";
        }
        myfile<<"\n";
    }
    myfile.close();
    return;
}

//print 1d integer variable to a text file
void print1Dinttotxtfile(char *filename, int *E, int rowSize, int colSize)
{
    ofstream myfile;

    myfile.open(filename);

    for (int ir = 0; ir<rowSize; ir++)
    {
        myfile<<E[ir]<<"\n";
    }
    myfile.close();
    return;
}

/*DEFINING AUXILARY FUNCTIONS*/
void crossprod(double *v1, double *v2, double *result)
{
    //Crossproduct assuming that inputs are two 3d vectors
    //result is an empty 3D array

    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v2[0]*v1[1];

    return;
}

double dotprod(double *v1, double *v2)
{
    //Dotproduct assuming that inputs are two 3d vectors
    //result is a single double value

    double result;

    result = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

    return result;
}

double magnit(double *v1)
{
    //Magnitude assuming that input is a 3d vector
    //result is a single double value

    double result;

    result = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

    return result;
}

void vectorsub(double *v1, double *v2, double *vsub)
{
    for (int i = 0; i<3; i++)
    {
        vsub[i] = v1[i] - v2[i];
    }
}

void vectoradd(double *v1, double *v2, double *vadd)
{
    for (int i = 0; i<3; i++)
    {
        vadd[i] = v1[i] + v2[i];
    }
}

void vectordiv(double *v1, double divs, double *vres)
{
    for (int i = 0; i<3; i++)
    {
        vres[i] = v1[i]/divs;
    }
}

void vectormult(double *v1, double mult, double *vres)
{
    for (int i = 0; i<3; i++)
    {
        vres[i] = v1[i]*mult;
    }
}

void vectorprint(double *v1, string vecName)
{
    cout<<vecName<<"\n";
    for (int b = 0; b<3; b++)
    {
        cout<<v1[b]<<"\t";
    }
    cout<<"\n";
}

void Initdob1D(double *vec1, int len)
{
	for (int a = 0; a<len; a++)
	{
		vec1[a] = 0.0;
	}
}

void Initint1D(int *intt1, int len)
{
	for (int a = 0; a<len; a++)
	{
		intt1[a] = 0.0;
	}
}

void Initdob2D(double **vec2, int len1, int len2)
{
	for (int a = 0; a<len1; a++)
	{
		for (int b = 0; b<len2; b++)
		{
			vec2[a][b] = 0.0;
		}
	}
}

//doesn't seem to work :(
void transferNewtoVec(vector <vector<double> > E, double **En, int rowSize, int colSize)
{
    for (int a = 0; a < rowSize; a++)
    {
        cout<<a<<"\t";
        for (int b = 0; b < colSize; b++)
        {
            E[a][b] = En[a][b];
            cout<<b<<"\t"<<En[a][b]<<"\t"<<E[a][b]<<"\n";
        }
    }

    return;
}

/*Find all Particles in Neighboring cells to compare for neighboring search algorithm method 2
ParticleList has to be an empty integer array to put the satisfactory Particles in*/
void partsinNeigh(double *boxSize, double *gridSize, double **Parts, int noParts, int **ParticleList)
{
    double xgrid, ygrid;

    xgrid = boxSize[0]/gridSize[0];
    ygrid = boxSize[1]/gridSize[1];

    for (int aa = 0; aa<noParts; aa++)
    {
	    //find grid location of the particle searched for, particleno aa
	    for (int ix = 0; ix<gridSize[0];ix++)
	    {
	        if (abs(Parts[aa][0]-ix*xgrid)<err)
	        {
	            ParticleList[aa][0] = ix;
	            break;
	        }
	        else if ((Parts[aa][0]>ix*xgrid) && (Parts[aa][0]<(ix+1)*xgrid))
	        {
	            ParticleList[aa][0] = ix;
	            break;
	        }

	    }

	        for (int iy = 0; iy<gridSize[1];iy++)
	    {
	        if (abs(Parts[aa][1]-iy*ygrid)<err)
	        {
	            ParticleList[aa][1] = iy;
	            break;
	        }
	        else if ((Parts[aa][1]>iy*ygrid) && (Parts[aa][1]<(iy+1)*ygrid))
	        {
	            ParticleList[aa][1] = iy;
	            break;
	        }

	    }


	    //cout<<aa<<"\t"<<ParticleList[aa][0]<<"\t"<<Parts[aa][0]<<"\t"<<ParticleList[aa][1]<<"\t"<<Parts[aa][1]<<"\n";
	}

	return;
}

/*NEIGHBOURING SEARCH ALGORITHM, finds particle neighbors, using one of two methods ::
particleno -> Number for particle to find neighbors of should be between 0 and noparts-1
Parts -> 2d array containing all particle information where the row corresponds to particle number
            and the columns correspond to the x,y,z location as well as the particle radius
noparts -> Number of particles in the system. Size of the 2d array
K -> Proportionality constant used for comparing radius to distance to determine if particles are neighbors
gridSize -> An array containing the number of grids to be created in each axis, array should have a size of 3
boxSize -> An array containing the (x,y,z) size of the simulation box
lev: method 1-> checks all the particles in the space, so does the test for all the particles in the space
lev: method 2-> breaks search up into grids, simulation space is broken into grids and only particles in the near
            grids are considered for the particle search
neighs -> an empty array where the list of neighboring particles will be stored
output of the function is an integer showing the number of neighboring particles*/

int NeighboringSearch(double **Parts, int noparts, double K, double *gridSize, double *boxSize, int lev, int **neighs)
{
    int countneighs = 0;
    double tmpmag, tmppart[3];
    int **ParticleList, ctir;

    if (lev == 1)
    {
	    for (int particleno = 0; particleno<noparts; particleno++)
	    {
	        for (int cti = particleno+1; cti<noparts; cti++)
	        {
	            if (cti!=particleno)
	            {
	                tmppart[0] = Parts[particleno][0] - Parts[cti][0];
	                tmppart[1] = Parts[particleno][1] - Parts[cti][1];
	                tmppart[2] = Parts[particleno][2] - Parts[cti][2];

	                tmpmag = magnit(tmppart);

	                if (tmpmag<K*(Parts[particleno][3]+Parts[cti][3]))
	                {
	                	neighs[countneighs][0] = particleno;
	                    neighs[countneighs][1] = cti;
	                    countneighs++;
	                }
	            }
	        }

	    }
    }

    else
    {
	    ParticleList = store2Dint(noparts,2);

		partsinNeigh(boxSize,gridSize,Parts,noparts,ParticleList);

	    for (int particleno = 0; particleno<noparts; particleno++)
		{
	        for (int ctir = particleno+1; ctir<noparts; ctir++)
	            {
	                if (ParticleList[ctir][0]>=ParticleList[particleno][0]-1 && ParticleList[ctir][0]<=ParticleList[particleno][0]+2 && ParticleList[ctir][1]>=ParticleList[particleno][1]-1 && ParticleList[ctir][1]<=ParticleList[particleno][1]+2)
	                {
        				//cout<<particleno<<"\t"<<ctir<<" =>\t"<<ParticleList[particleno][0]<<"\t"<<ParticleList[ctir][0]<<"\t"<<ParticleList[particleno][1]<<"\t"<<ParticleList[ctir][1]<<"\n";
		                tmppart[0] = Parts[particleno][0] - Parts[ctir][0];
		                tmppart[1] = Parts[particleno][1] - Parts[ctir][1];
		                tmppart[2] = Parts[particleno][2] - Parts[ctir][2];

		                tmpmag = magnit(tmppart);

		                if (tmpmag<K*(Parts[particleno][3]+Parts[ctir][3]))
		                {
		                	neighs[countneighs][0] = particleno;
		                    neighs[countneighs][1] = ctir;
		                    countneighs++;
		                }	                	
	                }

	            }
	    }

        clear2Dint(ParticleList,noparts);
	}
	
	//cout<<"Level: "<<lev<<"\tNumber of neighbours:"<<countneighs<<"\n";
    return countneighs;
}

/*WALLCOLLISIONCHECK Algorithm =>
Function determines wall collision entering integer values in an array corresponding to:
0 - No wall collision
1 - -x (x = 0) wall collision
2 - +x (x = xEnd) wall collision
3 - -y (y = 0) wall collision
4 - +y (y = yEnd) wall collision
5 - -z (z = 0) wall collision
6 - +z (z = zEnd) wall collision
Returns the number of wall collisions*/

int WallCollisionCheck(int particleno, double **Parts, int noparts, double Kw, double *boxSize, int *Walls)
{
    int countwallc = 0;
    double xl,xh,yl,yh,zl,zh;

    xl = abs(Parts[particleno][0]-0);
    xh = abs(Parts[particleno][0]-boxSize[0]);
    yl = abs(Parts[particleno][1]-0);
    yh = abs(Parts[particleno][1]-boxSize[1]);
    zl = abs(Parts[particleno][2]-0);
    zh = abs(Parts[particleno][2]-boxSize[2]);

    if (xl<Kw*Parts[particleno][3])
    {
        Walls[countwallc] = 1;
        countwallc++;
    }
    if (xh<Kw*Parts[particleno][3])
    {
        Walls[countwallc] = 2;
        countwallc++;
    }
    if (yl<Kw*Parts[particleno][3])
    {
        Walls[countwallc] = 3;
        countwallc++;
    }
    if (yh<Kw*Parts[particleno][3])
    {
        Walls[countwallc] = 4;
        countwallc++;
    }
    if (zl<Kw*Parts[particleno][3])
    {
        Walls[countwallc] = 5;
        countwallc++;
    }
    if (zh<Kw*Parts[particleno][3])
    {
        Walls[countwallc] = 6;
        countwallc++;
    }

    return countwallc;
}

void nearestPtWall(double *oldmet, double *Xi, int nofac, int *Walls, double *boxS, double *nearpt, double Ri, double *overlapPt)
{
    double *P0, *P1, *P2, *Base, box0, box1, box2, tmp;
    double *v1, *v2, *QP, *v1v2, bot, top, *QPv2, u, *Norm;
    int fil;

    P0 = store1D(3);
    Initdob1D(P0,3);
    P1 = store1D(3);
    Initdob1D(P1,3);
    P2 = store1D(3);
    Initdob1D(P2,3);
    v1 = store1D(3);
    Initdob1D(v1,3);
    v2 = store1D(3);
    Initdob1D(v2,3);
    QP = store1D(3);
    Initdob1D(QP,3);
    v1v2 = store1D(3);
    Initdob1D(v1v2,3);
    QPv2 = store1D(3);
    Initdob1D(QPv2,3);
    Norm = store1D(3);
    Initdob1D(Norm,3);

    Base = store1D(7);

    for (int bs = 0; bs<7; bs++)
    {
        Base[bs] = 0.0;
    }

    box0 = boxS[0];
    box1 = boxS[1];
    box2 = boxS[2];

    Base[2] = box0;
    Base[4] = box1;
    Base[6] = box2;

    for (int pt = 0; pt<3; pt++)
    {
        tmp = Xi[pt];
        P0[pt] = tmp;
        P1[pt] = tmp;
        P2[pt] = tmp;
    }

    for (int b = 0; b<nofac; b++)
    {
        if (b==0)
        {
            for (int cord = 1; cord<7; cord++)
            {
                if (Walls[b]==cord)
                {
                    fil = floor((cord-1)/2);
                    tmp = Base[cord];
                    P1[fil] = tmp;
                }                
            }

        }

        if (b==1)
        {
            for (int cord = 1; cord<7; cord++)
            {
                if (Walls[b]==cord)
                {
                    fil = floor((cord-1)/2);
                    tmp = Base[cord];
                    P2[fil] = tmp;
                }                
            }

        }

        if (b==2)
        {
            for (int cord = 1; cord<7; cord++)
            {
                if (Walls[b]==cord)
                {
                    fil = floor((cord-1)/2);
                    tmp = Base[cord];
                    P0[fil] = tmp;
                }                
            }

        }
       
    }

    if (nofac==2)
    {
        for (int k = 0; k<3; k++)
        {
            v1[k] = P2[k] - P1[k];
            v2[k] = Xi[k] - oldmet[k];
            QP[k] = oldmet[k] - P1[k];
        }

        crossprod(v1,v2,v1v2);
        crossprod(QP,v2,QPv2);
        top = dotprod(QPv2,v1v2);
        bot = dotprod(v1v2,v1v2);

        u = top/bot;

        for (int k = 0; k<3; k++)
        {
            nearpt[k] = P1[k] + (P2[k] - P1[k])*u;
        }
    }

    if (nofac==3)
    {
        for (int k = 0; k<3; k++)
        {
            v1[k] = P0[k] - P1[k];
            v2[k] = P0[k] - P2[k];
        }

        crossprod(v1,v2,Norm);

         for (int k = 0; k<3; k++)
        {
            v1[k] = Xi[k] - oldmet[k];
            v2[k] = oldmet[k] - P0[k];
        }

        top = dotprod(Norm,v2);
        bot = dotprod(Norm,v1);

        u = -top/bot;     

        for (int k = 0; k<3; k++)
        {
            nearpt[k] = oldmet[k] + (Xi[k] - oldmet[k])*u;
        }  
    }

    //cout<<nofac<<"\t"<<u<<"\n";
    //vectorprint(Xi,"Xi");
    //vectorprint(nearpt,"nearestpt");

    clear1D(P0);
    clear1D(P1);
    clear1D(P2);
    clear1D(v1);
    clear1D(v2);
    clear1D(QP);
    clear1D(v1v2);
    clear1D(QPv2);
    clear1D(Norm);
    clear1D(Base);

    for (int i = 0; i<3; i++)
    {

        if (abs(Xi[i]-oldmet[i])>err)
        {
            for (int wi = 0; wi<nofac; wi++)
            {
                if (floor((Walls[wi]-1)/2)==i)
                {
                    if ((Walls[wi]-1)%2==0)
                    {
                        overlapPt[i] = (Ri - Xi[i]);
                    }
                    else
                    {
                        overlapPt[i] = -(Ri + Xi[i] - oldmet[i]);
                    }
                }
            }
        }

        else
        {
            overlapPt[i] = 0.0;
        }

    }

    return;
}

void openPVD()
{
    ofstream myfilePVD;

    myfilePVD.open("PARA.pvd");

    myfilePVD<<"<?xml version=\"1.0\"?>\n";
    myfilePVD<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"BigEndian\">\n";
    myfilePVD<<"<Collection>\n";

    myfilePVD.close();
    return;
}

void printtoPVD(double timet, char *filename)
{
    ofstream myfilePVD;

    myfilePVD.open("PARA.pvd",ofstream::app);

    myfilePVD<<"<DataSet timestep=\""<<timet<<"\" group=\"\" part=\"0\" file=\""<<filename<<"\"/>\n";

    return;
}

void closePVD()
{
    ofstream myfilePVD;
    myfilePVD.open("PARA.pvd",ofstream::app);
    myfilePVD<<"</Collection>\n</VTKFile>\n";

    myfilePVD.close();

    return;
}

void paraviewWrite(string file1,string file2, double **E, int noParticles, double timet, double fileval)
{

    ostringstream oss;
    string var;
    char *filename;

    oss<<file1<<fileval<<file2;

    var = oss.str();

    filename = new char [var.length()];

    strcpy(filename, var.c_str());

    ofstream myfile;

    myfile.open(filename);

    myfile<<"<?xml version=\"1.0\"?>\n";
    myfile<<"<!--Time = "<<timet<<"s -->\n";
    myfile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    myfile<<"<PolyData>\n";
    myfile<<"<Piece NumberOfPoints=\""<<noParticles<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    myfile<<"<Points>\n";
    myfile<<"<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";


    for (int ir = 0; ir<noParticles; ir++)
    {
        myfile<<E[ir][0]<<"\t"<<E[ir][1]<<"\t"<<E[ir][2]<<"\t";
    }

    myfile<<"\n";

    myfile<<"</DataArray>\n</Points>\n";
    myfile<<"<PointData Scalars=\"Diameter\" Vectors=\"Velocity\">\n";
    myfile<<"<DataArray type=\"Float32\" Name=\"Diameter\" format=\"ascii\">\n";

    for (int ir = 0; ir<noParticles; ir++)
    {
        myfile<<2*E[ir][3]<<"\t";
    }

    myfile<<"\n";

    myfile<<"</DataArray>\n";
    myfile<<"<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (int ir = 0; ir<noParticles; ir++)
    {
        myfile<<E[ir][5]<<"\t"<<E[ir][6]<<"\t"<<E[ir][7]<<"\t";
    }

    myfile<<"\n";
    myfile<<"</DataArray>\n";
    myfile<<"<DataArray type=\"Float32\" Name=\"ID\" format=\"ascii\">\n";

    for (int ir = 0; ir<noParticles; ir++)
    {
        myfile<<ir+1<<"\t";
    }

    myfile<<"\n</DataArray>\n";
    myfile<<"</PointData>\n<CellData></CellData>\n<Verts></Verts>\n<Lines></Lines>\n";
    myfile<<"<Strips></Strips>\n<Polys></Polys>\n</Piece>\n</PolyData>\n</VTKFile>\n";

    myfile.close();

    printtoPVD(timet,filename);
    return;
}

int main (int argc, char *argv[])
{

    clock_t start;
    start = clock();

//cout<<setprecision(15);
//cout<<fixed;
    double **ParticleInfo,Xi0,Xi1,Xi2,box0,box1,box2,ForceTmag;
    double Xminmag, Di, Dj, ri, rj, Li, Lj;
    double kn, kt, bn, bt, m, I, g, rhoj, rhoi, Vndot, Vavg, Xavg, fcx, fcy, fcz;
    double knw, ktw, bnw, btw, intvpr, countpr, evatpr, transf, mew, meww;
    double dn, overlapn, mindnt, **ForceC, **TorQ, Vtijmag;
    double ForceNmag, overlaptmag, res, resw, meff, mj,countr, transfr,ftmd;
    int noParticles,  **Neighbs, nones, tmpvar, sizep, nb, wallc, Walls2Old;
    int inw, countneb, transfrint, countrwc, nonesOld, togodd, togeven, timctr, hertz, gridSearch;
    double startTime, endTime, deltaT,t, fx, fy, fz, Tx, Ty, Tz;

    double dtSdot, dttolddot,cross1mag;


    double *BoxSize,  *GridSize, *ForceN, *ForceT, **Tor, *Ttemp, **overlaptOld;
    double *nij, *Xi, *Xj, *Xmin, *Vij, *Vtij, *Vnij, *Grav, *tij, *Fcohp;
    double *Vi, *Vj, *Vmin, *wi, *wj, *Lwsum;
    double *WallXj,*Linij, *nearPT;
    double *Lwcn, *overlapt, **Fc, *overlap2NWall;
    double *pft_nb, **overlaptOldPtmp, **nijOldPtmp;

    int *Walls, **WallXjOld, **OldPairs;

    double *told, *tnew, *Sig, **nijOld, *cross1, **nijOldPairs, **overlaptOldPairs;
    double Fcoh, VDWin, VDWout, VDWinw, VDWoutw, asperity, conHamaker, conHamakerw, eqRad, phi, rsep; 
    double knl,ktl,knwl,ktwl,Ep,Ew,vp,vw,Gp,Gw;

    BoxSize = store1D(3);
    GridSize = store1D(2);
    Grav = store1D(3);

    stringstream strg(argv[1]);
    strg >> noParticles;

    //noParticles = 31;
    sizep = 11; //shouldn't change unless input file is changing and reading in differently

    vector <vector<double> > partcon;

    partcon.resize(noParticles);

    for (int i = 0; i < noParticles; i++)
    {
        partcon[i].resize(sizep);
    }

    double K = 1;
    double Kw = 1;

    hertz = 0; //0 for linear damping, 1 for hertzian
    gridSearch = 1; //1 for no gridding, 2 for gridding

    BoxSize[0] = 2;
    BoxSize[1] = 2;
    BoxSize[2] = 0.8;

    GridSize[0] = 15; //should be integers please
    GridSize[1] = 15; //same as above please

    deltaT = 0.00001;
    startTime = 0.0;
    endTime = 100;

    intvpr = 0.1;

    g = 9.8;//m/s^2 also um/ms^2...I think

    knl = 10e2;
    ktl = 2.0*knl/5.0;//7000;

    knwl = 10e2;
    ktwl = 2.0*knwl/5.0;//70000;

    Ep = 1.3e11;
    vp = 0.34;
    Ew = 1.6e11;
    vw = 0.27;
    Gp = Ep/(2*(1+vp));
    Gw = Ew/(2*(1+vw));

    res = 0.9; //restitution coefficients have to be between 0 and 1, not including 0
    resw = 0.9;

    mew = 0;
    meww = 0;

    VDWin = 0.04;
    VDWout = 0.5;
    VDWinw = 0.005;
    VDWoutw = 5;
    asperity = 0.5;
    conHamaker = 0.0;//28.4e-3; 
    conHamakerw = 0.0;//14e-5;

    t = startTime;

    countpr = t/intvpr;
    countr = 0;

    Grav[0] = 0;
    Grav[1] = 0;
    Grav[2] = -g;

    overlaptOld = store2D(noParticles,3);
    nijOld = store2D(noParticles,3);
    WallXjOld = store2Dint(noParticles,4);
    nijOldPairs = store2D(noParticles*noParticles,3);
    overlaptOldPairs = store2D(noParticles*noParticles,3);
	OldPairs = store2Dint(noParticles*noParticles,2);

    ParticleInfo = store2D(noParticles,sizep); //order of columns:::x,y,z,Radius,rho,vx,vy,vz

    ifstream partFile ("particle_input.txt");

    for (int i = 1; i <= noParticles; i++)
    {
        int j = i-1;
        for (int k = 1; k<=sizep; k++)
        {
            int l = k-1;
            partFile>>ParticleInfo[j][l];
        }
        //partFile>>ParticleInfo[j][0]>>ParticleInfo[j][1]>>ParticleInfo[j][2]>>ParticleInfo[j][3]>>ParticleInfo[j][4]>>ParticleInfo[j][5]>>ParticleInfo[j][6]>>ParticleInfo[j][7]>>ParticleInfo[j][8]>>ParticleInfo[j][9]>>ParticleInfo[j][10];
    }

    print2Dtotxtfile("particleOut",".txt",ParticleInfo,noParticles,sizep,0);

    openPVD();
    paraviewWrite("PARA",".vtp",ParticleInfo,noParticles,t,0);
    countr++;

    for (int z = 0; z < noParticles; z++)
    {
        for (int y = 0; y < sizep; y++)
        {
            partcon[z][y] = 0.0;
            partcon[z][y] = ParticleInfo[z][y];
        }


       for (int k = 0; k<3; k++)
		{
           overlaptOld[z][k] = 0.0;
           nijOld[z][k] = 0.0;
           WallXjOld[z][k] = 0;
		}

		WallXjOld[z][3] = 0;
    }

    for (int z = 0; z < noParticles*noParticles; z++)
    {

       for (int k = 0; k<3; k++)
		{
           overlaptOldPairs[z][k] = 0.0;
           nijOldPairs[z][k] = 0.0;
		}

		OldPairs[z][0] = 0;
		OldPairs[z][1] = 0;	
    }

    clear2D(ParticleInfo,noParticles);

    nonesOld = 0;
    timctr = 0;

    while (t<endTime)
    {
		ParticleInfo = store2D(noParticles,sizep);
		ForceC = store2D(noParticles,3);
		TorQ = store2D(noParticles,3);
		Tor = store2D(noParticles,3);
		Fc = store2D(noParticles,3);

		ForceN = store1D(3);
		Initdob1D(ForceN,3);
		ForceT = store1D(3);
		Initdob1D(ForceT,3);
		Fcohp = store1D(3);
		Initdob1D(Fcohp,3);
		Ttemp = store1D(3);
		Initdob1D(Ttemp,3);
		nij = store1D(3);
		Initdob1D(nij,3);
		Xi = store1D(3);
		Initdob1D(Xi,3);
		Xj = store1D(3);
		Initdob1D(Xj,3);
		Xmin = store1D(3);
		Initdob1D(Xmin,3);
		Vij = store1D(3);
		Initdob1D(Vij,3);
		Vtij = store1D(3);
		Initdob1D(Vtij,3);
		Vnij = store1D(3);
		Initdob1D(Vnij,3);
		tij = store1D(3);
		Initdob1D(tij,3);
		Vi = store1D(3);
		Initdob1D(Vi,3);
		Vj = store1D(3);
		Initdob1D(Vj,3);
		Vmin = store1D(3);
		Initdob1D(Vmin,3);
		wi = store1D(3);
		Initdob1D(wi,3);
		wj = store1D(3);
		Initdob1D(wj,3);
		Lwsum = store1D(3);
		Initdob1D(Lwsum,3);
		WallXj = store1D(3);
		Initdob1D(WallXj,3);
		Lwcn = store1D(3);
		Initdob1D(Lwcn,3);
		overlapt = store1D(3);
		Initdob1D(overlapt,3);

		pft_nb = store1D(3);
		Initdob1D(pft_nb,3);
		told = store1D(3);
		Initdob1D(told,3);
		tnew = store1D(3);
		Initdob1D(tnew,3);
		Sig = store1D(3);
		Initdob1D(Sig,3);
		cross1 = store1D(3);
		Initdob1D(cross1,3);
        nearPT = store1D(3);
        Initdob1D(nearPT,3);
		Walls = store1Dint(3);
		Initint1D(Walls,3);

		Linij = store1D(3);
		Initdob1D(Linij,3);

        overlap2NWall = store1D(3);
        Initdob1D(overlap2NWall,3);

		Vavg = 0.0;
		Xavg = 0.0;

		for (int z = 0; z < noParticles; z++)
		{

			for (int y = 0; y < sizep; y++)
			{
				ParticleInfo[z][y] = 0.0;
				ParticleInfo[z][y] = partcon[z][y];
			}

			for (int k = 0; k<3; k++)
			{
				ForceC[z][k] = 0.0;
				TorQ[z][k] = 0.0;
				Fc[z][k] = 0.0;
				Tor[z][k] = 0.0;
			}

		}

		for (int z = 0; z < noParticles; z++)
		{

		//cout<<"before: time "<<t<<" x "<<ParticleInfo[z][0]<<" y "<<ParticleInfo[z][1]<<"\n";

			if (ParticleInfo[z][0] <= 0)
			{
				ParticleInfo[z][0] = ParticleInfo[z][0] + BoxSize[0] - ParticleInfo[z][3];
			}
			else if (ParticleInfo[z][0] >= BoxSize[0])
			{
				ParticleInfo[z][0] = ParticleInfo[z][0] - BoxSize[0] + ParticleInfo[z][3];
			}


			if (ParticleInfo[z][1] <= 0)
			{
				ParticleInfo[z][1] = ParticleInfo[z][1] + BoxSize[1] - ParticleInfo[z][3];
			}
			else if (ParticleInfo[z][1] >= BoxSize[1])
			{
				ParticleInfo[z][1] = ParticleInfo[z][1] - BoxSize[1] + ParticleInfo[z][3];
			}
		//cout<<"after: time "<<t<<" x "<<ParticleInfo[z][0]<<" y "<<ParticleInfo[z][1]<<"\n";
		}

		//CHANGE TO PAIRS AND STORE THE OLD TANGENTIAL DISPLACEMENT HISTORY FOR THE PAIRS SPECIFICALLY
		//Particle - Particle Collision:

		Neighbs = store2Dint(noParticles*noParticles,2);

		for (int z = 0; z < noParticles*noParticles; z++)
		{
			Neighbs[z][0] = 0;
			Neighbs[z][1] = 0;			
		}

		nones = NeighboringSearch(ParticleInfo,noParticles,K,GridSize,BoxSize,gridSearch,Neighbs);

		if (nones>0)
		{
			nijOldPtmp = store2D(nones,3);
			overlaptOldPtmp = store2D(nones,3);

			for (int dd = 0; dd<nones; dd++)
			{
				for (int kk = 0; kk<3; kk++)
				{
					nijOldPtmp[dd][kk] = 0.0;
					overlaptOldPtmp[dd][kk] = 0.0;
				}

			}

			if (nonesOld>0)
			{
				for (int no = 0; no<nonesOld; no++)
				{
					for (int nbc = 0; nbc<nones; nbc++)
					{
						if (OldPairs[no][0]-Neighbs[nbc][0]==0 && OldPairs[no][1]-Neighbs[nbc][1]==0)
						{
							//cout<<OldPairs[no][0]<<" + "<<OldPairs[no][1]<<"...\t..."<<Neighbs[nbc][0]<<" + "<<Neighbs[nbc][1]<<"\n";
							for (int kk = 0; kk<3; kk++)
							{
								transf = nijOldPairs[no][kk];
								nijOldPtmp[nbc][kk] = transf;

								transf = overlaptOldPairs[no][kk];
								overlaptOldPtmp[nbc][kk] = transf;
							}
						}
					}
				}
			}

			countneb = 0;
			//cout<<a<<"\t#col"<<nones<<"\tcollides:\t";
			for (int nbc = 0; nbc < nones; nbc++)
			{
				int a = Neighbs[nbc][0];
				nb = Neighbs[nbc][1];

				for (int inm = 0; inm<3; inm++)
				{
					ForceT[inm] = 0.0;
					ForceN[inm] = 0.0;
					Fcohp[inm] = 0.0;
				}

				transfr = ParticleInfo[a][0];
				Xi[0] = transfr;
				transfr = ParticleInfo[a][1];
				Xi[1] = transfr;
				transfr = ParticleInfo[a][2];
				Xi[2] = transfr;
				transfr = ParticleInfo[a][3];
				Di = 2.0*transfr;
				transfr = ParticleInfo[a][4];
				rhoi = transfr;
				transfr = ParticleInfo[a][5];
				Vi[0] = transfr;
				transfr = ParticleInfo[a][6];
				Vi[1] = transfr;
				transfr = ParticleInfo[a][7];
				Vi[2] = transfr;
				transfr = ParticleInfo[a][8];
				wi[0] = transfr;
				transfr = ParticleInfo[a][9];
				wi[1] = transfr;
				transfr = ParticleInfo[a][10];
				wi[2] = transfr;
				ri = Di/2.0;

				m = rhoi*pi*Di*Di*Di/6.0;
				I = m*Di*Di/10.0;

				//cout<<nb<<"\t";
				transfr = ParticleInfo[nb][0];
				Xj[0] = transfr;
				transfr = ParticleInfo[nb][1];
				Xj[1] = transfr; 
				transfr = ParticleInfo[nb][2];
				Xj[2] = transfr;
				transfr = ParticleInfo[nb][3];
				Dj = 2.0*transfr;
				transfr = ParticleInfo[nb][4];
				rhoj = transfr;
				transfr = ParticleInfo[nb][5];
				Vj[0] = transfr;
				transfr = ParticleInfo[nb][6];
				Vj[1] = transfr;
				transfr = ParticleInfo[nb][7];
				Vj[2] = transfr;
				transfr = ParticleInfo[nb][8];
				wj[0] = transfr;
				transfr = ParticleInfo[nb][9];
				wj[1] = transfr;
				transfr = ParticleInfo[nb][10];
				wj[2] = transfr;
				rj = Dj/2.0;
				mj = rhoj*pi*Dj*Dj*Dj/6.0;

				meff = m*mj/(m+mj);

				//cout<<"Part-part b=> "<<bn<<"\t"<<bt<<"\n";
				//cout<<"Collision => "<<a<<"   +   "<<nb<<"\n";
				//vectorprint(Xi,"Xi");
				//vectorprint(Xj,"Xj");
				//vectorprint(Vi,"Vi");
				//vectorprint(Vj,"Vj");
				//vectorprint(wi,"wi");
				//vectorprint(wj,"wj");


				vectorsub(Xj,Xi,Xmin);
				Xminmag = magnit(Xmin);

				overlapn = 0.5*(Di+Dj)-Xminmag;
				vectordiv(Xmin,Xminmag,nij);

				vectorsub(Vi,Vj,Vmin);

				Li = ((Xminmag*Xminmag)+(ri*ri)-(rj*rj))/(2*Xminmag);
				Lj = Xminmag-Li;

				if (hertz==0)
				{
					kn = knl;
					kt = ktl;
				}
				
				else
				{
                    eqRad = ri*rj/(ri+rj);
					kn = 4*Ep*pow(abs(eqRad*overlapn),0.5)/(6*(1-vp*vp));
					kt = 16*Gp*pow(abs(eqRad*overlapn),0.5)/(6*(2-vp));
				}

				bn = 2.0*pow(meff*kn,0.5)*abs(log(res))/(pow((pi*pi+log(res)*log(res)),0.5));
				bt = bn/2.0;

				for (int cd = 0; cd<3; cd++)
				{
					Lwsum[cd] = (Li*wi[cd])+(Lj*wj[cd]);
				}

				crossprod(Lwsum,nij,Lwcn);
				vectoradd(Vmin,Lwcn,Vij);

				Vndot = dotprod(Vij,nij);
				vectormult(nij,Vndot,Vnij);

				for (int inm = 0; inm<3; inm++)
				{
					Vtij[inm] = Vij[inm] - Vnij[inm];
				}

				if (Vndot>err)
				{
					dn = abs(overlapn/Vndot);
					mindnt = min(dn,deltaT);
				}
				else
				{
					mindnt = deltaT;
				}

				vectormult(Vtij,mindnt,overlapt);

                if (Xminmag<ri+rj+VDWout)
                {
                    eqRad = 2*ri*rj/(ri+rj);

                    if (Xminmag>ri+rj+VDWin)
                    {
                        rsep = Xminmag-(ri+rj);
                        phi = conHamaker/(12.0*rsep*rsep);
                        Fcoh = phi*eqRad*((asperity/(asperity+eqRad))+(1.0/((1.0+asperity/rsep)*(1.0+asperity/rsep))));
                    }

                    else
                    {
                        phi = conHamaker/(12.0*VDWin*VDWin);
                        Fcoh = phi*eqRad*((asperity/(asperity+eqRad))+(1.0/((1.0+asperity/VDWin)*(1.0+asperity/VDWin))));                        
                    }

                    for (int inm = 0; inm<3; inm++)
                    {
                        Fcohp[inm] = (Xj[inm]-Xi[inm])*Fcoh/Xminmag;
                    }
                }                  

				if (magnit(nijOldPtmp[nbc])!=0 && magnit(overlaptOldPtmp[nbc])!=0)
				{
					crossprod(nijOldPtmp[nbc],nij,cross1);
					cross1mag = magnit(cross1);

					if (cross1mag>err)
					{
						vectordiv(cross1,cross1mag,Sig);
						crossprod(Sig,nijOldPtmp[nbc],cross1);
						cross1mag = magnit(cross1);
						vectordiv(cross1,cross1mag,told);
						crossprod(Sig,nij,cross1);
						cross1mag = magnit(cross1);
						vectordiv(cross1,cross1mag,tnew);

						dtSdot = dotprod(overlaptOldPtmp[nbc],Sig);
						dttolddot = dotprod(overlaptOldPtmp[nbc],told);

						for (int inm = 0; inm<3; inm++)
						{
							pft_nb[inm] = dtSdot*Sig[inm] + dttolddot*tnew[inm];
							overlapt[inm] = overlapt[inm] + pft_nb[inm];
						}
					}

					else
					{
						for (int inm = 0; inm<3; inm++)
						{
							transfr = overlaptOldPtmp[nbc][inm];
							pft_nb[inm] = transfr;
							overlapt[inm] = overlapt[inm] + pft_nb[inm];
						}
					}
				}

				for (int inm = 0; inm<3; inm++)
				{
					//transfr = nij[inm];
					ForceN[inm] = -(kn*overlapn*nij[inm] + bn*Vnij[inm]);
					//cout<<-kn*overlapn*nij[inm]<<"\t";
				}
				//cout<<"\n";
				for (int inm = 0; inm<3; inm++)
				{
					ForceT[inm] = -(kt*overlapt[inm] + bt*Vtij[inm]);
				}
				//Calculating tangential force when sliding occurs between particles

				ForceTmag = magnit(ForceT);
				ForceNmag = mew*magnit(ForceN);
				overlaptmag = magnit(overlapt);

				if (ForceTmag>ForceNmag)
				{
					Vtijmag = magnit(Vtij);
					vectordiv(Vtij,Vtijmag,tij);

					if (magnit(tij) > err)
					{
						for (int inm = 0; inm<3; inm++)
						{
							ForceT[inm] = -ForceNmag*tij[inm];
							overlapt[inm] = ForceNmag*tij[inm]/kt;
						}
					}

					else if ((magnit(tij) < err) && (overlaptmag>err))
					{
						for (int inm = 0; inm<3; inm++)
						{
							ForceT[inm] = -ForceNmag*overlapt[inm]/overlaptmag;
							overlapt[inm] = ForceNmag*overlapt[inm]/(kt*overlaptmag);
						}
					}

					else
					{
						for (int inm = 0; inm<3; inm++)
						{
						ForceT[inm] = 0.0;
						overlapt[inm] = 0.0;
						}
					}


					for (int inm = 0; inm<3; inm++)
					{
						ForceT[inm] = ForceT[inm]-bt*Vtij[inm];
					}


				}
				crossprod(nij,ForceT,Ttemp);

                //cout<<kn<<"\n"<<overlapn<<"\n";
                //cout<<Vndot<<"\n";
                //vectorprint(overlapt,"overlapt");
                //vectorprint(nijOldPtmp[nbc],"nijOld");    
                //vectorprint(overlaptOldPtmp[nbc],"overlaptOld");
                //vectorprint(Vij,"Vij");
                //vectorprint(Vnij,"Vnij");    
                //vectorprint(Vtij,"Vtij");
                //vectorprint(nij,"nij");

				/*if (a==18||nb==18){
                if (t>=17.46 && t<=19)
                {        
				cout<<"Collision => "<<a<<"   +   "<<nb<<"\n";      
		        vectorprint(Xi,"Xi");
		        vectorprint(Vi,"Vi");            
		        vectorprint(wi,"wi");
		        vectorprint(ForceT,"ForceT");
		        vectorprint(ForceN,"ForceN");
                }}*/

                //cout<<"SumOfForce==>\t";
				for (int inm = 0; inm<3; inm++)
				{
                    //cout<<ForceN[inm]+ForceT[inm]<<"\t";
					Fc[a][inm] = Fc[a][inm]+(ForceN[inm]+ForceT[inm]+Fcohp[inm]);
					Tor[a][inm] = Tor[a][inm]+Li*Ttemp[inm];

					ForceC[nb][inm] = ForceC[nb][inm]-(ForceN[inm]+ForceT[inm]+Fcohp[inm]);
					TorQ[nb][inm] = TorQ[nb][inm]+Lj*Ttemp[inm];

					transfr = nij[inm];
					nijOldPairs[nbc][inm] = transfr;

					transfr = overlapt[inm];
					overlaptOldPairs[nbc][inm] = transfr;
				}
	            //cout<<"\n";
				transfrint = Neighbs[nbc][0];
				OldPairs[nbc][0] = transfrint;
				transfrint = Neighbs[nbc][1];
				OldPairs[nbc][1] = transfrint;

				countneb++;
			}

			nonesOld = nones;

			clear2D(nijOldPtmp,nones);
			clear2D(overlaptOldPtmp,nones);
		}

		//vectorprint(ForceT,"ForceT");
		//vectorprint(overlapt,"overlapt");
		//vectorprint(Vtij,"Vtij");
		//vectorprint(nij,"nij");
		//vectorprint(Vij,"Vij");

		//cout<<"\n";
		//cout<<"p: "<<a<<" colliding p => "<<nones<<"\n";
		//vectorprint(ForceC[a],"ForceC");
		//vectorprint(Fc,"Fc");
		//vectorprint(TorQ[a],"TOrQ");
		//vectorprint(Tor,"Tor");
		//cout<<"\n\n";

		clear2Dint(Neighbs,noParticles*noParticles);

		//vectorprint(ForceN,"ForceN");
		//vectorprint(Vij,"Vij");
		//vectorprint(Xmin,"Xmin");
		//vectorprint(Lwcn,"Lwcn");

        for (int a = 0; a < noParticles; a++)
        {
            transfr = ParticleInfo[a][0];
            Xi[0] = transfr;
            transfr = ParticleInfo[a][1];
            Xi[1] = transfr;
            transfr = ParticleInfo[a][2];
            Xi[2] = transfr;
            transfr = ParticleInfo[a][3];
            Di = 2.0*transfr;
            transfr = ParticleInfo[a][4];
            rhoi = transfr;
            transfr = ParticleInfo[a][5];
            Vi[0] = transfr;
            transfr = ParticleInfo[a][6];
            Vi[1] = transfr;
            transfr = ParticleInfo[a][7];
            Vi[2] = transfr;
            transfr = ParticleInfo[a][8];
            wi[0] = transfr;
            transfr = ParticleInfo[a][9];
            wi[1] = transfr;
            transfr = ParticleInfo[a][10];
            wi[2] = transfr;
            ri = Di/2.0;

            //cout<<"particle=>"<<a<<"\n";
            //vectorprint(Vi,"Vi");
            //vectorprint(Xi,"Xi");
            //vectorprint(wi,"wi");

            m = rhoi*pi*Di*Di*Di/6.0;
            I = m*Di*Di/10.0;

            fx = 0.0;
            fy = 0.0;
            fz = 0.0;

            Tx = 0.0;
            Ty = 0.0;
            Tz = 0.0;

            for (int inm = 0; inm<3; inm++)
            {
                pft_nb[inm] = 0.0;
                ForceT[inm] = 0.0;
                ForceN[inm] = 0.0;
            }

            /*for (int inm = 0; inm<3; inm++)
            {
                cout<<"ParticleStart "<<a<<"\t"<<Fc[a][inm]+ForceC[a][inm]<<"\t"<<Tor[a][inm]+TorQ[a][inm]<<"\n";
            }*/

            /*Particle - Wall Collision:*/
            wallc = WallCollisionCheck(a,ParticleInfo,noParticles,Kw,BoxSize,Walls);
            //cout<<"\n#Wallcollisions = "<<wallc<<"\n";

            if (wallc==0)
            {

                for (int z = 0; z < noParticles; z++)
                {
                   for (int k = 0; k<3; k++)
                    {
                       overlaptOld[z][k] = 0.0;
                       nijOld[z][k] = 0.0;
                       WallXjOld[z][k] = 0;
                    }

                    WallXjOld[z][3] = 0;
                }                
            }

            if (wallc!=0)
            {

                for (int inm = 0; inm<3; inm++)
                {
                    nijOld[a][inm] = 0.0;
                    overlaptOld[a][inm] = 0.0;
                }                       

                if (WallXjOld[a][0]==wallc)
                {
                    countrwc = 0;

                    for (int wc = 0; wc<wallc; wc++)
                    {
                        countrwc = countrwc + (WallXjOld[a][wc+1]-Walls[wc]);
                    }

                    if (countrwc==0)
                    {
                        for (int inm = 0; inm<3; inm++)
                        {
                            transfr = nijOld[a][inm];
                            nijOld[a][inm] = transfr;
                            transfr = overlaptOld[a][inm];
                            overlaptOld[a][inm] = transfr;
                        }
                    }
                }
                
                //cout<<"Wall Collision"<<a<<"\n";
                box0 = BoxSize[0];
                box1 = BoxSize[1];
                box2 = BoxSize[2];
                Xi0 = Xi[0];
                Xi1 = Xi[1];
                Xi2 = Xi[2];

                if (wallc==1)
                {
                    inw = Walls[wallc-1];

                    switch(inw)
                    {
                        case 1:
                            WallXj[0] = 0;
                            WallXj[1] = Xi1;
                            WallXj[2] = Xi2;
                        break;

                        case 2:
                            WallXj[0] = box0;
                            WallXj[1] = Xi1;
                            WallXj[2] = Xi2;
                        break;

                        case 3:
                            WallXj[0] = Xi0;
                            WallXj[1] = 0;
                            WallXj[2] = Xi2;
                        break;

                        case 4:
                            WallXj[0] = Xi0;
                            WallXj[1] = box1;
                            WallXj[2] = Xi2;
                        break;

                        case 5:
                            WallXj[0] = Xi0;
                            WallXj[1] = Xi1;
                            WallXj[2] = 0;
                        break;

                        case 6:
                            WallXj[0] = Xi0;
                            WallXj[1] = Xi1;
                            WallXj[2] = box2;
                        break;
                    }

                    if(timctr%2==1)
                    {
                        togodd = Walls[0];
                    }
                    else
                    {
                        togeven = Walls[0];
                    }
                }

                else if (wallc==2)
                {
                    inw = Walls[0] + Walls[wallc-1];

                    switch(inw)
                    {
                        case 4: //xlow & ylow
                            WallXj[0] = 0;
                            WallXj[1] = 0;
                            WallXj[2] = Xi2;
                        break;

                        case 5:
                            if (Walls[0]==1) //xlow and yhigh
                            {
                                WallXj[0] = 0;
                                WallXj[1] = box1;
                                WallXj[2] = Xi2;
                            }
                            if (Walls[0]==2) //xhigh and ylow
                            {
                                WallXj[0] = box0;
                                WallXj[1] = 0;
                                WallXj[2] = Xi2;
                            }
                        break;

                        case 6:
                            if (Walls[0]==2) //xhigh and yhigh
                            {
                                WallXj[0] = box0;
                                WallXj[1] = box1;
                                WallXj[2] = Xi2;
                            }
                            if (Walls[0]==1) //xlow and zlow
                            {
                                WallXj[0] = 0;
                                WallXj[1] = Xi1;
                                WallXj[2] = 0;
                            }
                        break;

                        case 7:
                            if (Walls[0]==1) //xlow and zhigh
                            {
                                WallXj[0] = 0;
                                WallXj[1] = Xi1;
                                WallXj[2] = box2;
                            }
                            if (Walls[0]==2) //xhigh and zlow
                            {
                                WallXj[0] = box0;
                                WallXj[1] = Xi1;
                                WallXj[2] = 0;
                            }
                        break;

                        case 8:
                            if (Walls[0]==2) //xhigh and zhigh
                            {
                                WallXj[0] = box0;
                                WallXj[1] = Xi1;
                                WallXj[2] = box2;
                            }
                            if (Walls[0]==3) //ylow and zlow
                            {
                                WallXj[0] = Xi0;
                                WallXj[1] = 0;
                                WallXj[2] = 0;
                            }
                        break;

                        case 9:
                            if (Walls[0]==3) //ylow and zhigh
                            {
                                WallXj[0] = Xi0;
                                WallXj[1] = 0;
                                WallXj[2] = box2;
                            }
                            if (Walls[0]==4) //yhigh and zlow
                            {
                                WallXj[0] = Xi0;
                                WallXj[1] = box1;
                                WallXj[2] = 0;
                            }
                        break;

                        case 10:
                            WallXj[0] = Xi0;
                            WallXj[1] = box1;
                            WallXj[2] = box2;
                        break;

                    }

                }

                else if (wallc==3)
                {
                    inw = Walls[0] + Walls[1] + Walls[wallc-1];

                    switch(inw)
                    {
                        case 9: //xlow,ylow,zlow
                            WallXj[0] = 0;
                            WallXj[1] = 0;
                            WallXj[2] = 0;
                        break;

                        case 10:
                            if ((Walls[0]==1)&&(Walls[1]==4)&&(Walls[2]==5)) //xlow,yhigh, zlow
                            {
                                WallXj[0] = 0;
                                WallXj[1] = box1;
                                WallXj[2] = 0;
                            }
                            else if ((Walls[0]==1)&&(Walls[1]==3)&&(Walls[2]==6)) //xlow,ylow,zhigh
                            {
                                WallXj[0] = 0;
                                WallXj[1] = 0;
                                WallXj[2] = box2;
                            }
                            else if ((Walls[0]==2)&&(Walls[1]==3)&&(Walls[2]==5)) //xhigh,ylow,zlow
                            {
                                WallXj[0] = box0;
                                WallXj[1] = 0;
                                WallXj[2] = 0;
                            }
                        break;

                        case 11:
                            if ((Walls[0]==1)&&(Walls[1]==4)&&(Walls[2]==6)) //xlow,yigh,zhigh
                            {
                                WallXj[0] = 0;
                                WallXj[1] = box1;
                                WallXj[2] = box2;
                            }
                            else if ((Walls[0]==2)&&(Walls[1]==4)&&(Walls[2]==5)) //xhigh,yhigh,zlow
                            {
                                WallXj[0] = box0;
                                WallXj[1] = box1;
                                WallXj[2] = 0;
                            }
                            else if ((Walls[0]==2)&&(Walls[1]==3)&&(Walls[2]==6)) //xhigh,ylow,zhigh
                            {
                                WallXj[0] = box0;
                                WallXj[1] = 0;
                                WallXj[2] = box2;
                            }
                        break;

                        case 12: //xhigh, yhigh, zhigh
                            WallXj[0] = box0;
                            WallXj[1] = box1;
                            WallXj[2] = box2;
                        break;
                    }
                }

                if (wallc>1)
                {
                    nearestPtWall(WallXj,Xi,wallc,Walls,BoxSize,nearPT,ri,overlap2NWall);

                    for (int ik = 0; ik<3; ik++)
                    {
                        transfr = nearPT[ik];
                        WallXj[ik] = transfr;
                    }

                }

                /*
                cout<<inw<<"[";

                for (int ch = 0; ch <3; ch++)
                {
                cout<<WallXj[ch]<<" ";
                }
                cout<<"]\n";

                cout<<"X-[";

                for (int ch = 0; ch <3; ch++)
                {
                cout<<Xi[ch]<<" ";
                }
                cout<<"]\n";
                */
                //vectorprint(WallXj,"WallXj");
                meff = m;

                //cout<<"Part-wall b=> "<<bnw<<"\t"<<btw<<"\n";

                vectorsub(WallXj,Xi,Xmin);
                Xminmag = magnit(Xmin);
                overlapn = ri-Xminmag;

                vectordiv(Xmin,Xminmag,nij);

				if (hertz==0)
				{
					knw = knwl;
					ktw = ktwl;
				}
						
				else
				{
                	eqRad = ri;
					knw = 4*Ep*Ew*pow(abs(eqRad*overlapn),0.5)/(3*(Ep*(1-vp*vp)+Ew*(1-vw*vw)));
					ktw = 16*Gp*Gw*pow(abs(eqRad*overlapn),0.5)/(3*(Gp*(2-vp)+Gw*(2-vw)));
				}

                bnw = 2.0*pow(meff*knw,0.5)*abs(log(resw))/(pow((pi*pi+log(resw)*log(resw)),0.5));
                btw = bnw/2.0;

                Li = Xminmag;//((Xminmag*Xminmag)+(ri*ri))/(2*Xminmag);

                for (int cd = 0; cd<3; cd++)
                {
                    Lwsum[cd] = Li*wi[cd];
                }

                crossprod(Lwsum,nij,Lwcn);
                vectoradd(Vi,Lwcn,Vij);

                Vndot = dotprod(Vij,nij);

                vectormult(nij,Vndot,Vnij);

                //vectorprint(Vi,"Vi");
                //cout<<"Li "<<Li<<"\n";
                //vectorprint(wi,"wi");
                //vectorprint(Vij,"Vij");
                //vectorprint(nij,"nij");

                for (int inm = 0; inm<3; inm++)
                {
                    Vtij[inm] = Vij[inm] - Vnij[inm];
                }


                if (Vndot>err)
                {
                    dn = overlapn/Vndot;
                    mindnt = min(dn,deltaT);
                }
                else
                {
                    mindnt = deltaT;
                }

                vectormult(Vtij,mindnt,overlapt);

		
				if(Walls[wallc-1]==5)
				{
	                if (Xminmag<ri+VDWoutw)
	                {
	                    eqRad = ri;

	                    if (Xminmag>ri+VDWinw)
	                    {
	                        rsep = Xminmag-ri;
	                        phi = conHamakerw/(12.0*rsep*rsep);
	                        Fcoh = phi*eqRad*((asperity/(asperity+eqRad))+(1.0/((1.0+asperity/rsep)*(1.0+asperity/rsep))));
	                    }

	                    else
	                    {
	                        phi = conHamakerw/(12.0*VDWinw*VDWinw);
	                        Fcoh = phi*eqRad*((asperity/(asperity+eqRad))+(1.0/((1.0+asperity/VDWinw)*(1.0+asperity/VDWinw))));                        
	                    }

		            	Fc[a][2] = Fc[a][2] + (WallXj[2]-Xi[2])*Fcoh/abs(WallXj[2]-Xi[2]);
	                }    
				}

                if (magnit(nijOld[a])!=0 && magnit(overlaptOld[a])!=0)
                {
                    crossprod(nijOld[a],nij,cross1);
                    cross1mag = magnit(cross1);

                    if (cross1mag>err)
                    {
                        vectordiv(cross1,cross1mag,Sig);
                        crossprod(Sig,nijOld[a],cross1);
                        cross1mag = magnit(cross1);
                        vectordiv(cross1,cross1mag,told);
                        crossprod(Sig,nij,cross1);
                        cross1mag = magnit(cross1);
                        vectordiv(cross1,cross1mag,tnew);

                        dtSdot = dotprod(overlaptOld[a],Sig);
                        dttolddot = dotprod(overlaptOld[a],told);

                        for (int inm = 0; inm<3; inm++)
                        {
                            pft_nb[inm] = dtSdot*Sig[inm] + dttolddot*tnew[inm];
                            overlapt[inm] = overlapt[inm] + pft_nb[inm];
                        }
                    }

                    else
                    {
                        for (int inm = 0; inm<3; inm++)
                        {
                            transfr = overlaptOld[a][inm];
                            pft_nb[inm] = transfr;
                            overlapt[inm] = overlapt[inm] + pft_nb[inm];
                        }
                    }

                }

                if (wallc>1)
                {
                    for (int inm = 0; inm<3; inm++)
                    {
                        ForceN[inm] = (knw*overlap2NWall[inm]);
                    }                   
                }
                else
                {
                    for (int inm = 0; inm<3; inm++)
                    {
                        ForceN[inm] = -(knw*overlapn*nij[inm]);
                    }                    
                }

				/*if (a==18){
                if (t>=17.46 && t<=19)
                {

		    	cout<<"particle=>"<<a<<" just spring\n";
		        vectorprint(ForceN,"ForceN");
                }}*/


                if (wallc>1)
                {
                    for (int inm = 0; inm<3; inm++)
                    {
                        ForceN[inm] = (knw*overlap2NWall[inm] - bnw*Vnij[inm]);
                    }                   
                }
                else
                {
                    for (int inm = 0; inm<3; inm++)
                    {
                        ForceN[inm] = -(knw*overlapn*nij[inm] + bnw*Vnij[inm]);
                    }                    
                }

                for (int inm = 0; inm<3; inm++)
                {
                    ForceT[inm] = -(ktw*overlapt[inm] + btw*Vtij[inm]);
                }

                //Calculating tangential force when sliding occurs between particles

                ForceTmag = magnit(ForceT);
                ForceNmag = meww*magnit(ForceN);
                overlaptmag = magnit(overlapt);

                if (ForceTmag>ForceNmag)
                {
                    Vtijmag = magnit(Vtij);
                    vectordiv(Vtij,Vtijmag,tij);

                    if (magnit(tij) > err)
                    {
                        for (int inm = 0; inm<3; inm++)
                        {
                            ForceT[inm] = -ForceNmag*tij[inm];
                            overlapt[inm] = ForceNmag*tij[inm]/ktw;
                        }
                    }

                    else if ((magnit(tij) < err) && (overlaptmag>err))
                    {
                        for (int inm = 0; inm<3; inm++)
                        {
                            ForceT[inm] = -ForceNmag*overlapt[inm]/overlaptmag;
                            overlapt[inm] = ForceNmag*overlapt[inm]/(ktw*overlaptmag);
                        }
                    }

                    else
                    {
                        for (int inm = 0; inm<3; inm++)
                        {
                            ForceT[inm] = 0.0;
                            overlapt[inm] = 0.0;
                        }
                    }

                    for (int inm = 0; inm<3; inm++)
                    {
                        ForceT[inm] = ForceT[inm]-btw*Vtij[inm];
                    }

                }

                for (int inm = 0; inm<3; inm++)
                {
                    Linij[inm] = Li*nij[inm];
                }

                crossprod(Linij,ForceT,Ttemp);

				/*if (a==18){
                if (t>=17.46 && t<=19)
                {
		    	cout<<"particle=>"<<a<<"\n";
		        //cout<<"overlapn "<<overlapn<<"\n";
		        //cout<<"Xminmag "<<Xminmag<<"\n";                
		        vectorprint(Xi,"Xi");
		        vectorprint(Vi,"Vi");
		        //vectorprint(nij,"nij");               
		        vectorprint(wi,"wi");
		        //vectorprint(Xmin,"Xmin");
		        //vectorprint(overlap2NWall,"overlap2NWall");
		        vectorprint(Vnij,"Vnij");
		        //vectorprint(overlapt,"overlapt"); 
		        //vectorprint(WallXj,"WallXj");
		        vectorprint(ForceT,"ForceT");
		        vectorprint(ForceN,"ForceN");
                }}*/

                for (int inm = 0; inm<3; inm++)
                {
                    Fc[a][inm] = Fc[a][inm]+ForceN[inm]+ForceT[inm];
                    Tor[a][inm] = Tor[a][inm]+Li*Ttemp[inm];

                    transfr = nij[inm];
                    nijOld[a][inm] = transfr;
                    transfr = overlapt[inm];
                    overlaptOld[a][inm] = transfr;
                }

                WallXjOld[a][0] = wallc;

                for (int wc = 0; wc<wallc; wc++)
                {
                    transfrint = Walls[wc];
                    WallXjOld[a][wc+1] = transfrint;
                }

            }

            /*for (int inm = 0; inm<3; inm++)
            {
                cout<<"ParticleAfterWall "<<a<<"\t"<<Fc[a][inm]+ForceC[a][inm]<<"\t"<<Tor[a][inm]+TorQ[a][inm]<<"\n";
            }*/

        }


		for (int a = 0; a<noParticles; a++)
		{

			transfr = ParticleInfo[a][0];
			Xi[0] = transfr;
			transfr = ParticleInfo[a][1];
			Xi[1] = transfr;
			transfr = ParticleInfo[a][2];
			Xi[2] = transfr;
			transfr = ParticleInfo[a][3];
			Di = 2.0*transfr;
			transfr = ParticleInfo[a][4];
			rhoi = transfr;
			transfr = ParticleInfo[a][5];
			Vi[0] = transfr;
			transfr = ParticleInfo[a][6];
			Vi[1] = transfr;
			transfr = ParticleInfo[a][7];
			Vi[2] = transfr;
			transfr = ParticleInfo[a][8];
			wi[0] = transfr;
			transfr = ParticleInfo[a][9];
			wi[1] = transfr;
			transfr = ParticleInfo[a][10];
			wi[2] = transfr;

			m = rhoi*pi*Di*Di*Di/6.0;
			I = m*Di*Di/10.0;

            //vectorprint(Xi,"Xi");
            //vectorprint(Vi,"Vi");

			for (int inm = 0; inm<3; inm++)
			{
				//cout<<"ParticleEnd "<<a<<"\t"<<Fc[a][inm]<<"\t"<<Fc[a][inm]+ForceC[a][inm]<<"\t"<<Tor[a][inm]<<"\t"<<Tor[a][inm]+TorQ[a][inm]<<"\n";
                //cout<<m<<"\n"<<Grav[inm]<<"\n";
                //cout<<a<<" "<<Fc[a][inm] + ForceC[a][inm]<<"\n";
				Vi[inm] = Vi[inm] + (Grav[inm] + (Fc[a][inm] + ForceC[a][inm])/m)*deltaT;
				Xi[inm] = Xi[inm] + Vi[inm]*deltaT;
				wi[inm] = wi[inm] + (Tor[a][inm]+TorQ[a][inm])*deltaT/I;
			}

            //vectorprint(Xi,"Xi");
            //vectorprint(Vi,"Vi");

			transfr = Xi[0];
			ParticleInfo[a][0] = transfr;
			transfr = Xi[1];
			ParticleInfo[a][1] = transfr;
			transfr = Xi[2];
			ParticleInfo[a][2] = transfr;
			transfr = Vi[0];
			ParticleInfo[a][5] = transfr;
			transfr = Vi[1];
			ParticleInfo[a][6] = transfr;
			transfr = Vi[2];
			ParticleInfo[a][7] = transfr;
			transfr = wi[0];
			ParticleInfo[a][8] = transfr;
			transfr = wi[1];
			ParticleInfo[a][9] = transfr;
			transfr = wi[2];
			ParticleInfo[a][10] = transfr;

		}


		for (int h = 0; h<noParticles; h++)
		{
			Xavg = Xavg+(ParticleInfo[h][0]+ParticleInfo[h][1]+ParticleInfo[h][2])/3;
			Vavg = Vavg+(ParticleInfo[h][5]+ParticleInfo[h][6]+ParticleInfo[h][7])/3;
		}

		evatpr = intvpr + countpr*intvpr;

		if ((evatpr>t-deltaT)&&(evatpr<t+deltaT))
		{
			cout<<"Time: "<<t<<"\tAverage Position:"<<Xavg<<"\tAverage Velocity: "<<Vavg<<"\n";

			print2Dtotxtfile("particleOut",".txt",ParticleInfo,noParticles,sizep,countr);
			countpr++;

			paraviewWrite("PARA",".vtp",ParticleInfo,noParticles,t,countr);
			countr++;

			//cout<<"Time so far: "<<(clock()-start)/((double)CLOCKS_PER_SEC*60)<<"min\n";
		}

		/*if ((t+deltaT>17.4)&&(t-deltaT<17.5))
		{
			print2Dtotxtfile("particleOut",".txt",ParticleInfo,noParticles,sizep,countr);
			paraviewWrite("PARA",".vtp",ParticleInfo,noParticles,t,countr);
			countr++;
		}*/

	

		for (int z = 0; z < noParticles; z++)
		{/*
			if (ParticleInfo[z][2]<0)
			{
				transf = ParticleInfo[z][2];
				ParticleInfo[z][2] = -transf;
			}

			else if (ParticleInfo[z][2]>BoxSize[2])
			{
				transf = ParticleInfo[z][2] - BoxSize[2];
				ParticleInfo[z][2] = BoxSize[2]-transf;
			}*/

			for (int y = 0; y < sizep; y++)
			{
				partcon[z][y] = 0;
				partcon[z][y] = ParticleInfo[z][y];
			}
		}

		//cout<<"Time: "<<t<<"\tAverage Position:"<<Xavg<<"\tAverage Velocity: "<<Vavg<<"\n";

		clear2D(ParticleInfo,noParticles);
		clear2D(ForceC,noParticles);
		clear2D(TorQ,noParticles);
		clear1D(ForceN);
		clear1D(ForceT);
		clear2D(Tor,noParticles);
		clear1D(Ttemp);
		clear1D(nij);
		clear1D(Xi);
		clear1D(Xj);
		clear1D(Xmin);
		clear1D(Vij);
		clear1D(Vtij);
		clear1D(Vnij);
		clear1D(tij);
		clear1D(Vi);
		clear1D(Vj);
		clear1D(Vmin);
		clear1D(wi);
		clear1D(wj);
		clear1D(Lwsum);
		clear1D(WallXj);
		clear1D(Lwcn);
		clear1D(overlapt);
		clear2D(Fc,noParticles);
		clear1D(pft_nb);
		clear1D(told);
		clear1D(tnew);
		clear1D(Sig);
		clear1D(cross1);
		clear1D(Linij);
		clear1D(nearPT);
		clear1D(Fcohp);
        clear1D(overlap2NWall);
        clear1Dint(Walls);

        timctr++;
		t = t + deltaT;
    }

    clear2D(overlaptOld,noParticles);
    clear2D(nijOld,noParticles);
    clear2D(nijOldPairs,noParticles*noParticles);
	clear2D(overlaptOldPairs,noParticles*noParticles);

    clear2Dint(WallXjOld,noParticles);
    clear2Dint(OldPairs,noParticles*noParticles);

    clear1D(BoxSize);
    clear1D(GridSize);
    clear1D(Grav);

    closePVD();

    cout<<"Total time: "<<(clock()-start)/((double)CLOCKS_PER_SEC*60)<<"min\n";
    return 0;
}
