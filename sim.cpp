//Copyright 2019 Thomas M. Poulsen

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string.h> 
using namespace std;

struct randomwalkstruct
{
	int nTimeStepsToReach;
	int nTimeStepsToInitTrans;
	double dPercentOfmRNAInitTrans;
	double dEntraped;	
	double dPercentReach;
	double dMLR;
	double dAlpha;
	double pMitoImport;
};


double GetRandomDouble(double minimum, double maximum)
{
     double r = ( (double)rand() / (double)(RAND_MAX) ); 
	 return r*(maximum-minimum) + minimum;
}

void PrintHelp()
{
	cout << "\n"; 
	cout << "How to run simulation:\n";
	cout << "---------------------- \n";
	cout << "sim [output file] (optional arguments)\n";
	cout << "\n";	

	cout << "Optional arguments:\n";
	cout << "-------------------\n";
	cout << "-a" << "   " << "Anchoring probability                                        " << "   " << "[calculated]" << "\n";
	cout << "-e" << "   " << "Overall probability of entrapment                            " << "   " << "[0.4]" << "\n";
	cout << "-m" << "   " << "Number of mitochondria                                       " << "   " << "[5]" << "\n";	
	cout << "-x" << "   " << "Fold increase of entrapment prob. upon translation initiation" << "   " << "[5]" << "\n";	
	cout << "\n";
	cout << "-r" << "   " << "Number of simulations to run                                 " << "   " << "[1]" << "\n";
	cout << "-w" << "   " << "Number of random walks for each mRNA molecule                " << "   " << "[1000]" << "\n";
	cout << "-s" << "   " << "Simulation time (seconds)                                    " << "   " << "[660]" << "\n";
	cout << "\n";
	cout << "-n" << "   " << "Diffusion constant for non-translating mRNA                  " << "   " << "[3.7*10^-9]" << "\n";
	cout << "-t" << "   " << "Diffusion constant for translating mRNA                      " << "   " << "[0.9*10^-9]" << "\n";
	cout << "-d" << "   " << "One-axis or three-axis simulation (= 1 or 3)                 " << "   " << "[1]" << "\n";
	cout << "\n";
	cout << "-g" << "   " << "Gene information file                                        " << "   " << "[geneinfo.csv]" << "\n";	
	cout << "-i" << "   " << "Translation initiation times input file name                 " << "   " << "[inittimes.csv]" << "\n";	
		
	cout << "\n";
	cout << "One/three-axis and diffusion constants:\n";
	cout << "---------------------------------------\n";
	cout << "The diffusion constants reflect the values for a three-axis simulation, while\nrunning with a one-axis simulation effectively results in constants that are\n1/3 of the values given. Use -b 1 for one-axis and -b 3 for three-axis\nsimulations (other input values not possible). Difussion constants need to be\nwritten in their full format, for example: 0.0000000037.\n\n";
	cout << "Input files:\n";
	cout << "------------\n";
	cout << "The number of entries in the gene information file should be the same as the\ntranslation initiation times file. Note: the gene information file has a header\nwhile the initiation times file doesn't (see geneinfo.csv and inittimes.csv).\n";
	
	cout << "\n";
}


randomwalkstruct RandomWalk3d(double dTimeToInitTranslation, int nmRNA, int nTotalMito, double dtotaltime, double dt, 
					int nsteps, double Dnotrans, double Dtrans, int nSimulationAxis,
					double pEntrapmentNoInitTrans, double pMitoImport, double pMitoReach, 
					double dInitTransEntrapIncrease, double cellradius, double dmitominradius, double dmitomaxradius, int nStepsToMito, 
					double dPercentOfmRNAInitTrans, double dAlpha)
{
	int i;
	int j;
	int k;
	int ndirections;	
	int bMove;
	int bMitoContact;
	int nbuf;
	int nmitobuf;
	int bRanOk;
	int nCyto;
	int nMito;	
	int n;
	int n1;
	int n2;
	int nExpectedStepsForInitTrans;
	int nCheckTransInitEveryNSteps;
	int nMitoLocalize;
	int nReach;
	int nEntraped;
	int nNumberOfmRNAInitTrans;
	int nTimeStepsToReach;
	int nTimeStepsToInitTrans;
	
	double r;
	double p;
	double dx;
	double D;
	double x;
	double dxtrans;
	double dxnotrans;
	double dtemp;
	double xnew;
	double ynew;
	double znew;
	double xmito;
	double ymito;
	double zmito;
	double rmito;
	double xmitoprev;
	double ymitoprev;
	double zmitoprev;
	double rmitoprev;

	double x0;
	double y0;
	double z0;

	double ddist;
	double dMLR;
	double pEntrapment;
	double pMitoLocalize;

	double pStepsToInitTranslation;
	double pInitTranslation;
	double dAvgInitTransSteps;
	double dCurrentTime;
	double pInitiatedTrans;
	double dp;
		
	randomwalkstruct rs;			
	
	const int nMaxmRNA = 100000;
	
	double xmitos[nTotalMito];
	double ymitos[nTotalMito];
	double zmitos[nTotalMito];
	double rmitos[nTotalMito];		
	
	int bInit[nMaxmRNA];
	double xs[nMaxmRNA];
	double ys[nMaxmRNA];
	double zs[nMaxmRNA];
	int bInitTranslations[nMaxmRNA];	
	int nTranslationSteps[nMaxmRNA];	
	
	ndirections= 6;
	double pmove = 1.0/(double)ndirections;
	
	nReach = 0;
	
	x0 = 0;
	y0 = 0;
	z0 = 0;

	pStepsToInitTranslation = (int)(dTimeToInitTranslation/dt);
	pInitTranslation = 1.0/pStepsToInitTranslation;
	
	nEntraped = 0;
	
	nbuf = 10;         //Minimum dx steps that a mitochondria must be away from mRNA starting point x0,y0,z0 and from other mitochondria
	nmitobuf = 0;      //If a mRNA is nmitobuf*dx steps near a mitochondria it can be imported for localization
		
			
	dxnotrans = sqrt(2*Dnotrans*dt);		
	dxtrans = sqrt(2*Dtrans*dt); 
	
	nTimeStepsToReach = 0;
	nTimeStepsToInitTrans = 0;
	nNumberOfmRNAInitTrans = 0;
		
	for(i = 0; i < nTotalMito; i++) 
	{
		bRanOk = 0;
		
		while(bRanOk == 0)
		{
			rmito = GetRandomDouble(dmitominradius,dmitomaxradius);
			xmito = GetRandomDouble(-cellradius,cellradius);
			ymito = GetRandomDouble(-cellradius,cellradius);
			zmito = GetRandomDouble(-cellradius,cellradius);
			
			bRanOk = 1;
		
			ddist = sqrt((xmito*xmito) + (ymito*ymito));
			
			if(ddist+rmito >= cellradius) bRanOk = 0;
			if(ddist < (rmito+(nbuf*dx))) bRanOk = 0;
			
			
			for(j = 0; j < i; j++)
			{
				xmitoprev = xmitos[j];
				ymitoprev = ymitos[j];
				zmitoprev = zmitos[j];
				rmitoprev = rmitos[j];								
				
				ddist = sqrt( ((xmito-xmitoprev)*(xmito-xmitoprev)) + ((ymito-ymitoprev)*(ymito-ymitoprev)) + ((zmito-zmitoprev)*(zmito-zmitoprev)) );
				
				if(ddist <= (rmito+rmitoprev)) bRanOk = 0;
			}
			
			if(bRanOk == 0)
			{
				rmito = GetRandomDouble(dmitominradius,dmitomaxradius);
				xmito = GetRandomDouble(-cellradius,cellradius);
				ymito = GetRandomDouble(-cellradius,cellradius);
				zmito = GetRandomDouble(-cellradius,cellradius);
			}
		}
		
		rmitos[i] = rmito;
		xmitos[i] = xmito;
		ymitos[i] = ymito;
		zmitos[i] = zmito;				
	}		
	
	
	for(i = 0; i < nMaxmRNA; i++) 
	{
		bInit[i] = 0;
		bInitTranslations[i] = 0;
	}
	

	for(i = 0; i < nsteps; i++)
	{
		for(j = 0; j < nmRNA; j++)
		{
			if(bInit[j] == 0)
			{
				xs[j] = x0;
				ys[j] = y0;
				zs[j] = z0;
				bInit[j] = 1;													
				
			}
			else if(bInit[j] == 1)
			{			
				if(bInitTranslations[j] == 0)
				{
					r = GetRandomDouble(0,1.0);
									
					if(r < pInitTranslation) 
					{										
						bInitTranslations[j] = 1;
						
						nTimeStepsToInitTrans += i;
											
						nNumberOfmRNAInitTrans++;
						
						nTranslationSteps[j] = 0;
					}
				}				
								
				xnew = xs[j];
				ynew = ys[j];
				znew = zs[j];
							
				if(bInitTranslations[j] == 0) dx = dxnotrans; 
				else dx = dxtrans;
					
				if(nSimulationAxis == 1)
				{
					r = GetRandomDouble(0,1.0);
					
					if(r < pmove)        xnew += dx; 
					else if(r < 2*pmove) xnew -= dx;
					else if(r < 3*pmove) ynew += dx;
					else if(r < 4*pmove) ynew -= dx;
					else if(r < 5*pmove) znew += dx;
					else if(r < 6*pmove) znew -= dx;					
				}
				else
				{			
					r = GetRandomDouble(0,1.0);										
					if(r < 0.5) xnew += dx; 
					else 		xnew -= dx;
				
					r = GetRandomDouble(0,1.0);
					if(r < 0.5) ynew += dx; 
					else 		ynew -= dx;

					r = GetRandomDouble(0,1.0);
					if(r < 0.5) znew += dx; 
					else 		znew -= dx;
				}
				
				bMitoContact = 0;
				bMove = 1;
				
				for(k = 0; k < nTotalMito; k++)
				{
					xmito = xmitos[k];
					ymito = ymitos[k];
					zmito = zmitos[k];
					rmito = rmitos[k];
					
					ddist = sqrt( (xnew-xmito)*(xnew-xmito)+(ynew-ymito)*(ynew-ymito)+(znew-zmito)*(znew-zmito));
					
					if(ddist <= (rmito+(nmitobuf*dx))) bMitoContact = 1;	
										
					if(ddist <= rmito) 
					{
						bMove = 0;
						bMitoContact = 1;
					}				
				}
				
				ddist = sqrt( (xnew*xnew) + (ynew*ynew) + (znew*znew) );
				
				if(ddist >= cellradius) bMove = 0;
				
				if(bMove == 1)
				{
					xs[j] = xnew;
					ys[j] = ynew;
					zs[j] = znew;					
				}
				
				if(bMitoContact == 1) //mRNA near mitochondria, see if it is imported this step
				{									
					r = GetRandomDouble(0,1.0);
						
					if(r < pMitoImport)
					{
						bInit[j] = -1;
					}		
					
									
				}
				else //mRNA not close to mitocondria, see if it becomes entraped this step
				{
					r = GetRandomDouble(0,1.0);
										
					if(bInitTranslations[j] == 1) pEntrapment = pEntrapmentNoInitTrans*dInitTransEntrapIncrease;
					else pEntrapment = pEntrapmentNoInitTrans;
															
					if(r < pEntrapment) 
					{					
						bInit[j] = -2;		
						nEntraped++;		
					}
				}												
			}			
		}

	}	
	

	for(i = 0; i < nmRNA; i++) 
	{
		if(bInit[i] != -1) 
		{
			nTimeStepsToReach += nsteps;
		}
	}
	
	nTimeStepsToReach = (int)((double)nTimeStepsToReach/(double)nmRNA);
	
	if(nNumberOfmRNAInitTrans > 0)
	{
		nTimeStepsToInitTrans = (int)((double)nTimeStepsToInitTrans/(double)nNumberOfmRNAInitTrans);
	}
		
	nCyto = 0;
	nMito = 0;
	
	for(i = 0; i < nmRNA; i++)
	{
		if(bInit[i] == -1) nMito++;
		else nCyto++;	
	}

	dMLR = (log((double)(nMito+1))-log((double)(nCyto+1)))/log((double)(nmRNA+2));
	
	
	rs.dPercentOfmRNAInitTrans = (double)nNumberOfmRNAInitTrans/(double)nmRNA;		
	rs.nTimeStepsToReach = nTimeStepsToReach;	
	rs.nTimeStepsToInitTrans = nTimeStepsToInitTrans;	
	rs.dEntraped = (double)nEntraped/(double)nmRNA;	
	rs.pMitoImport = pMitoImport;
	
	if(pMitoReach == -1) 
	{
		rs.dPercentReach = ((double)nMito)/((double)(nMito+nCyto));
	}
	else rs.dMLR = dMLR;
	
	return rs;
	
}


			   
void RunMitoSim(string sOutFile, string sInitTimesFile, string sInfoFile, double dProbOfMitoImport, double dProbOfOverallEntrapment, 
int nNumberOfRandomWalksPermRNA, double dFoldIncreaseOfEntrapProbOnTrans, int nSimulations,
double Dnotrans, double Dtrans, int nSimulationAxis, int nTotalMito, double dTotaltime)
{
	int i;	
	int j;
	int ncount;	
	int nStepsToMito;
	int nStepsToInitTrans;
	int nStepsToInitTransAvg;
	int nStepsToMitoAvg;
	int nMitoLocalize;	
	int nSimulationSteps;		
	int nmRNAForEstimate;
	int nFactorStep;
	int nFactor;
	int nFactorOld;	
	
	double d;
	double dstep;
	double dp;
	double dtemp;
	double dtotal;
	double dmlr;
	double davgmlr;	
	double dminmlr;
	double dmaxmlr;
	double stdev;
	double dMin;
	double dMax;	
	double pOverallEntrapment;
	double dt;
	double pInitiatedTrans;
	double pMitoImport;
	double pMitoReach;
	double pMitoLocalize;
	double pInitTranslation;
	double dInitTransEntrapIncrease;
	double cellradius;
	double dmitominradius;
	double dmitomaxradius;
	double dPercentOfmRNAInitTrans;
	double dAlpha;
	double dEntraped;
	double dInitTranslationTime;
	double pAvgInitTranslation;	
	double pEntrapmentNoInitTrans;
	double dPercentEntraped;	
	double dStep;
	double minEstimatePercentEntrapted;
	double dStart;
	
	string line;
	string lineinfo;
	string soutmlrsfile;
	string soutinfofile;
	
	randomwalkstruct rs;	
	
	
	const int nmaxcount = 10000;	
	double dmlrs[nmaxcount];
	
	ifstream fin;
	ifstream finfo;
	
	
	finfo.open(sInfoFile.c_str());
	
	if(!fin)
	{
		cout << "Failed to open " << sInitTimesFile << "\n";
		return;
	}
	
	if(!finfo)
	{
		cout << "Failed to open " << sInfoFile << "\n";
		return;
	}
	
	finfo.close();
	
	
	soutmlrsfile = sOutFile;
	soutmlrsfile += "_mlrs";
	
	soutinfofile = sOutFile;
	soutinfofile += "_info";
	
	ncount = 0;
	davgmlr = 0;
	dMin = 1000;
	dMax = -1000;
	
	srand( (unsigned)time( NULL ) );
	
	cellradius     = 0.00025;	//2.5um radius
	dmitominradius = 0.000025;  //0.25um radius
	dmitomaxradius = 0.000050;  //0.5um  radius
	

	//For parameter grid search
	dPercentEntraped = 1;
	nmRNAForEstimate = 100;
	pEntrapmentNoInitTrans = 0.1;
	minEstimatePercentEntrapted = 0.05;
	nFactorStep = 10;
		
	pOverallEntrapment = dProbOfOverallEntrapment;
	dInitTransEntrapIncrease = dFoldIncreaseOfEntrapProbOnTrans;
		
	dt = 0.01;	   
	pMitoImport = dProbOfMitoImport;
	pMitoReach = -1;
	
	dtotal = 0;		
	nStepsToMito = 0;
	nStepsToMitoAvg = 0;
	nStepsToInitTrans = 0;
	 
	
	dPercentOfmRNAInitTrans = 0;
	
	pMitoReach = -1;
	
	nSimulationSteps = (int)(dTotaltime/dt);		
	pMitoReach = 1;
	
	cout << "Deriving parameter search space (this may take some time)...\n";
		
	while(dPercentEntraped > minEstimatePercentEntrapted)
	{
		dEntraped = 0;				
		pEntrapmentNoInitTrans = pEntrapmentNoInitTrans*0.1;
		
		fin.open(sInitTimesFile.c_str());	
		ncount = 0;
		while(fin.eof() == 0)
		{
			getline(fin,line);
			if(line.length() > 0)
			{
				dInitTranslationTime = atof(line.c_str());
			
				rs = RandomWalk3d(dInitTranslationTime, nmRNAForEstimate, nTotalMito, dTotaltime, dt, nSimulationSteps, Dnotrans, Dtrans, nSimulationAxis, pEntrapmentNoInitTrans, pMitoImport, pMitoReach, dInitTransEntrapIncrease, cellradius, dmitominradius, dmitomaxradius, nStepsToMito, dPercentOfmRNAInitTrans, dAlpha);
			
				dEntraped += rs.dEntraped;		
				ncount++;				
			}
		}
		fin.close();
		dPercentEntraped = dEntraped/(double)ncount;
	}
	
	dStep = pEntrapmentNoInitTrans;
	
	
	
	nFactor = nFactorStep;
	nFactorOld = 1;
	dStart = dStep*nFactor;
	while(dPercentEntraped < (dProbOfOverallEntrapment-0.01))
	{		
		pEntrapmentNoInitTrans = dStart;
		dEntraped = 0;
				
		fin.open(sInitTimesFile.c_str());	
		ncount = 0;
		while(fin.eof() == 0)
		{
			getline(fin,line);
			if(line.length() > 0)
			{
				dInitTranslationTime = atof(line.c_str());
			
				rs = RandomWalk3d(dInitTranslationTime, nmRNAForEstimate, nTotalMito, dTotaltime, dt, nSimulationSteps, Dnotrans, Dtrans, nSimulationAxis, pEntrapmentNoInitTrans, pMitoImport, pMitoReach, dInitTransEntrapIncrease, cellradius,dmitominradius,dmitomaxradius,nStepsToMito, dPercentOfmRNAInitTrans, dAlpha);
			
				dEntraped += rs.dEntraped;		
				ncount++;				
			}
		}
		fin.close();
		dPercentEntraped = dEntraped/(double)ncount;
		
		
		if(dPercentEntraped < (dProbOfOverallEntrapment-0.01))
		{
			dStart = dStep*nFactor;
			nFactorOld = nFactor;
			nFactor += nFactorStep;
		}
	}
	
	dStart = dStep*nFactorOld;
	
	
	cout << "Performing parameter grid search...\n";	
	
	dPercentEntraped = 0;
	pEntrapmentNoInitTrans = dStart;
	
	while(dPercentEntraped < dProbOfOverallEntrapment)
	{
		dEntraped = 0;	
		pEntrapmentNoInitTrans += dStep;	
				
		fin.open(sInitTimesFile.c_str());	
		ncount = 0;
		while(fin.eof() == 0)
		{
			getline(fin,line);
			if(line.length() > 0)
			{
				dInitTranslationTime = atof(line.c_str());
			
				rs = RandomWalk3d(dInitTranslationTime, nmRNAForEstimate, nTotalMito, dTotaltime, dt, nSimulationSteps, Dnotrans, Dtrans, nSimulationAxis, pEntrapmentNoInitTrans, pMitoImport, pMitoReach, dInitTransEntrapIncrease, cellradius, dmitominradius, dmitomaxradius, nStepsToMito, dPercentOfmRNAInitTrans, dAlpha);
			
				dEntraped += rs.dEntraped;		
				ncount++;				
			}
		}
		fin.close();
		dPercentEntraped = dEntraped/(double)ncount;
	}
	
	if(dPercentEntraped > dProbOfOverallEntrapment) pEntrapmentNoInitTrans -= dStep;
	

	cout << "Running simulation...\n\n";

	///////// Perform simulation

	string sOutFileRun;
	string sOutMlrsFileRun;
	string sOutInfoFileRun;
	
	ofstream fout;
	ofstream foutmlrs;
	ofstream foutinfo;
		
	for(i = 0; i < nSimulations; i++)
	{			
		ncount = 0;
		dEntraped = 0;				
			
		sOutFileRun = sOutFile;
		sOutMlrsFileRun = soutmlrsfile;
		sOutInfoFileRun = soutinfofile;
		
		if(nSimulations > 1)
		{
			std::ostringstream os;
			os << i+1;
		
			sOutFileRun     += "_";
			sOutMlrsFileRun += "_";
			sOutInfoFileRun += "_";
				
			sOutFileRun     += os.str();
			sOutMlrsFileRun += os.str();
			sOutInfoFileRun += os.str();
		}
		
		
		fin.open(sInitTimesFile.c_str());	
		finfo.open(sInfoFile.c_str());
		
		fout.open(sOutFileRun.c_str());
		foutmlrs.open(sOutMlrsFileRun.c_str());
		foutinfo.open(sOutInfoFileRun.c_str());
	
		foutmlrs << "simmlrs\n";
	
		//Remove end of line character if present	
		getline(finfo,lineinfo);			
		if(lineinfo.length() > 0 && (int)lineinfo[lineinfo.length()-1] == 13) lineinfo.erase(lineinfo.length()-1,1); 
		fout << lineinfo << ",simmlrs\n";
				
		//Run simulation for each gene and its translation initiation 
		//time (for each gene we run nNumberOfRandomWalksPermRNA times) 
		while(fin.eof() == 0)
		{
			getline(fin,line);
			getline(finfo,lineinfo);
			if(lineinfo.length() > 0 && (int)lineinfo[lineinfo.length()-1] == 13) lineinfo.erase(lineinfo.length()-1,1); 
		
			if(line.length() > 0)
			{
				dInitTranslationTime = atof(line.c_str());
			
				rs = RandomWalk3d(dInitTranslationTime, nNumberOfRandomWalksPermRNA, nTotalMito, dTotaltime, dt, nSimulationSteps,
				Dnotrans, Dtrans, nSimulationAxis, pEntrapmentNoInitTrans, pMitoImport, pMitoReach, dInitTransEntrapIncrease, cellradius,dmitominradius, dmitomaxradius, nStepsToMito, dPercentOfmRNAInitTrans, dAlpha);
			
				dmlr = rs.dMLR;			
				dEntraped += rs.dEntraped;						
				dmlrs[ncount] = dmlr;			
				davgmlr += dmlr;
				ncount++;
			
				if(dmlr > dMax) dMax = dmlr;
				if(dmlr < dMin) dMin = dmlr;
			
				fout << lineinfo << "," << dmlr << "\n";
				foutmlrs << dmlr << "\n";			
			}
		}
		
		fin.close();
		finfo.close();
		
				
		dEntraped = dEntraped/ncount;
		davgmlr = davgmlr/(double)ncount;
			
		stdev = 0;
	
		for(j = 0; j < ncount; j++)
		{	
			d = dmlrs[i]-davgmlr;
			stdev += d*d;		
		}	
		stdev = sqrt(stdev/(double)(ncount-1));

		foutinfo << "\n";
		foutinfo << "Ration entraped = " << dEntraped << "\n\n";
	
		foutinfo << "Avgmlr = " << davgmlr << "\n";
		foutinfo << "Stdmlr = " << stdev << "\n";
		foutinfo << "Minmlr = " << dMin << "\n";
		foutinfo << "Maxmlr = " << dMax << "\n";	
						
		fout.close();
		foutmlrs.close();
		foutinfo.close();					
		
	}

	
	
}

int main(int argc,char *argv[])
{
	int i;
	int n;
	int nNumberOfRandomWalksPermRNA;
	int nSimulations; 
	int nTotalMito;
	int nSimulationAxis;
	
	double d;
	double dProbOfOverallEntrapment;
	double dProbOfMitoImport;  //Anchoring probability
	double dFoldIncreaseOfEntrapProbOnTrans;
	double dTotaltime; 
	double Dnotrans;
	double Dtrans;
	
	string s1;
	string s2;
	string sOutFile;
	string sInfoFile;
	string sInitTimesFile;
	
	ifstream fin;
	
	
	dProbOfMitoImport = 1;
	dFoldIncreaseOfEntrapProbOnTrans = 5;
	dProbOfOverallEntrapment = 0.40;
	nNumberOfRandomWalksPermRNA = 1000;
	
	nSimulations = 1;
	nTotalMito = 5;
	dTotaltime = 660; //sec

	//Diffusion constants - running with the 3-axis option results in the below values while the 
	//1-axis simulation effectivly results in difussion constants that are 1/3 of the two 
	//values listed below:	
	Dnotrans = 3.7*pow(10,-9);  
	Dtrans = 0.9*pow(10,-9); 
	nSimulationAxis = 1;
		
	sInfoFile = "geneinfo.csv";
	sInitTimesFile = "inittimes.csv";
		
	if(argc < 2)
	{
		PrintHelp();
		return 0;
	}	
	
	sOutFile = argv[1];	
	
	i = 2;	
	while(i < argc-1)
	{
		s1 = argv[i];
		s2 = argv[i+1];
		
		if(s1.length() > 1 && s1[0] == '-')
		{
			switch(s1[1])
			{
				case 'a':
					d = atof(s2.c_str());
					if(d <= 0 || d > 1) 
					{
						cout << "Please input anchor probability that is > 0 and <= 1\n";
						return 0;
					}					
					dProbOfMitoImport = d;
				break;
				
				case 'e':
				
					d = atof(s2.c_str());
					if(d < 0 || d > 1) 
					{
						cout << "Please input entrapment probability that is >= 0 and <= 1\n";
						return 0;
					}		
					dProbOfOverallEntrapment = d;
				break;
				
				case 'g':
					fin.open(s2.c_str());
					if(!fin)
					{
						cout << "Gene information file not found\n";
						return 0;
					}
					fin.close();
					sInfoFile = s2;
				break;
				
				case 'i':
					fin.open(s2.c_str());
					if(!fin)
					{
						cout << "Translation times input file not found\n";
						return 0;
					}
					fin.close();
					sInitTimesFile = s2;
				break;
				
				case 'n':
					Dnotrans = atof(s2.c_str());
				break;
					
				case 't':
					Dtrans = atof(s2.c_str());
				break;
					
				case 'd':
					nSimulationAxis = atoi(s2.c_str());
					if(nSimulationAxis != 1 && nSimulationAxis != 3)
					{
						cout << "One/three axis value should be 1 or 3\n";
						return 0;
					}
				break;

	
				case 'w':
					nNumberOfRandomWalksPermRNA = atoi(s2.c_str());
				break;
				
				case 'r':
					nSimulations = atoi(s2.c_str());
				break;
								
				case 'm':
					nTotalMito = atoi(s2.c_str());
				break;
				
				case 's':
					dTotaltime = atof(s2.c_str());
				break;
					
				case 'x':
					dFoldIncreaseOfEntrapProbOnTrans = atoi(s2.c_str());
				break;
				
				case 'h':
					PrintHelp();
				break;
				
			}
		}
		
		i += 2;			
	}
	

	RunMitoSim(sOutFile, sInitTimesFile, sInfoFile, dProbOfMitoImport, dProbOfOverallEntrapment, 
			   nNumberOfRandomWalksPermRNA, dFoldIncreaseOfEntrapProbOnTrans, nSimulations, 
			   Dnotrans, Dtrans, nSimulationAxis, nTotalMito, dTotaltime);	
	
		
	return 0;
}




	
