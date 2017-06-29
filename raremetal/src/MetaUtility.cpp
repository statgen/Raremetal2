#include "MetaUtility.h"
#include <math.h>
#include <complex>
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include "StringBasics.h"
#include "qfc.h"
#include <algorithm>
#include "math.h"
#define MATHLIB_STANDALONE 
#include <Rmath.h>
//#include <gsl/gsl_integration.h> // integration for Dajiang's method

bool SetIfilePosition( IFILE & sfile, Tabix & myTabix, String Chr, int pos )
{
	const char* chr = Chr.c_str();
	uint64_t fstart = 0;
    bool st = myTabix.getStartPos( chr, pos, fstart );
    if ( !st ) {
		printf("Warning:[SetIfilePosition] Unable to locate %s:%d\n", chr, pos);
        return false;
    }
// seek position now with fstart
//      if ( fstart != (uint64_t)iftell(sfile) ) {
    if ( ifseek(sfile,fstart,SEEK_SET) != true ) {
        printf("Unable to seek position %s:%d\n", chr, pos);
        exit(1);
    }
//      }


    if ( ifeof(sfile) ) {
        printf("Warning:[SetIfilePosition] Reached end!\n");
        return false;
    }
// then move on to that position
    String buffer;
    StringArray tmp;
    bool found = false;
    while( !ifeof(sfile) ) {
        buffer.ReadLine(sfile);
        tmp.Clear();
        tmp.AddTokens( buffer, "\t" );
        if ( tmp[1].AsInteger() >= pos ) {
            found = true;
            break;
        }
    }
    return found;
}

// RMW: adjust = 0;
// RVtest: adjust = 1;
void tellRvOrRmw( String & buffer, bool & adjust, int marker_col, int cov_col )
{
    if(buffer.Find("RareMetalWorker")==-1) // rvt
        adjust =1;
    else  //rmw
        adjust = 0;
    setFromRvOrRmwAdjust( adjust, marker_col, cov_col );
}       

bool setFromRvOrRmwAdjust( bool adjust, int marker_col, int cov_col )
{
    if (adjust) { // rv test
        marker_col = 4;
        cov_col = 5;
    }
    else { // rmw
        marker_col = 2;
        cov_col = 3;        
    }
}


// open result file as prefix.meta.method.results
void openMetaResultFile( String & prefix, String & filename, IFILE & output, String & method )
{
    if(prefix =="")
        filename = "meta."+method + ".results";
    else if(prefix.Last()=='.' || prefix.Last()=='/')
        filename = prefix +  "meta."+method +".results";
    else
        filename = prefix + ".meta."+method +".results";
    output=ifopen(filename,"w",InputFile::UNCOMPRESSED);
}

/*
double GetGenomicControlFromChisq(Vector& study_chisqs)
{
    if (study_chisqs.Length()<1) {
        printf("Warning: No p value available for calculating Genomic Control in this study. Check if markers were all filtered out!\n\n");
        return 0;
    }
    study_chisqs.Sort();
    double chisq_median = study_chisqs[0.5];
    double GC = chisq_median / 0.456; // qchisq(0.5,1)
    return GC;
}
*/

double GetGenomicControlFromPvalue(Vector& pvalues)
{
    if (pvalues.Length()<1) {
        printf("Warning: No p value available for calculating Genomic Control in this study. Check if markers were all filtered out!\n\n");
        return 0;
    }
    pvalues.Sort();
    double p_median = pvalues[0.5];
    double GC = qchisq(p_median,1,0,0);
    GC /= 0.456;
    return GC;
}


double GetBetaDensity(double a, double b, double x)
{
    double density;
    //density = exp(gammln(a+b)-gammln(a)-gammln(b)+(a-1.0)*log(x)+(b-1.0)*log(1.0-x));
    density = dbeta(x,a,b,0);
    return density;
}


double CalculateCorrCoef(Vector & a,Vector & b)
{
    double coef;
    int n = a.Length();
    double sum_a = a.Sum();
    double sum_b = b.Sum();
    coef =  n*a.InnerProduct(b) - sum_a*sum_b;
    coef /= sqrt((n*a.SumSquares()-sum_a*sum_a)*(n*b.SumSquares()-sum_b*sum_b));
    return coef;
}

void RevertAllele(String SNP, String & newSNP)
{
    StringArray tmp;
    tmp.AddTokens(SNP,":");
    newSNP = tmp[0]+":"+tmp[1]+":"+tmp[3]+":"+tmp[2];
}

void removeChrFromString( String & chr_token )
{
    if (chr_token.Find("chr")!=-1)
        chr_token = chr_token.SubStr(3);
}

char getDirection(String & chr_pos,double effsize,bool flip)
{
    char direction = '+';
    if(flip) {
        if(effsize>0)
            direction = '-';
    }
    else {
        if(effsize<0)
            direction = '-';
    }
    return direction;
}

double MixChidist(double lambda [], int n, double q_obs,String method)
{
    double * nc1 = new double [n];
    int * n1 = new int [n];
    double * trace = new double [7];
    double sigma=0.0;
    int lim1 = 10000;
    double  acc = 0.0001;
    int  ifault = 0;
    double  res = 0.0;

    for(int i=0;i<n;i++) {
        nc1[i] = 0.0; //noncentral parameters
        n1[i] = 1; //h in R function
    }
    for(int i=0;i<7;i++)
        trace[i] = 0.0;

    if(method=="Davies") {
        //probability will be saved in res.
        qfc(lambda,nc1,n1,& n,& sigma, & q_obs,& lim1,& acc,trace,& ifault,& res);
        res = 1.0-res;
    }

    //From Lee's SKAT package
    if(method=="Liu") {
        double  c1=0.0,c2=0.0,c3=0.0,c4=0.0;
        for(int i=0;i<n;i++) {
            c1 += lambda[i];
            c2 += lambda[i]*lambda[i];
            c3 += lambda[i]*lambda[i]*lambda[i];
            c4 += lambda[i]*lambda[i]*lambda[i]*lambda[i];
        }

        //printf("%g,%g,%g,%g\n",c1,c2,c3,c4);
        double s1 = c3/sqrt(c2*c2*c2);
        double s2 = c4/(c2*c2);
        double muQ = c1;
        double sigmaQ = sqrt(2.0*c2);

        double tstar = (q_obs-muQ)/sigmaQ;
        double delta,l,a;

        if (s1*s1 > s2)  {
            a = 1.0/(s1-sqrt(s1*s1-s2));
            delta = s1*a*a*a-a*a;
            l = a*a-2.0*delta;
        }
        else {
            a = 1.0/s1;
            delta = 0.0;
            l = c2*c2*c2/(c3*c3);
        }

        double muX = l+delta;
        double sigmaX = sqrt(2.0)*a;
        double q_new = tstar*sigmaX+muX;

        //printf("q_new=%g, df=%g, ncp=%g\n",q_new,l,delta);
        double Qq;
        if(delta==0)
        Qq = pchisq(q_new,l,0,0);
        else 
        Qq = pnchisq(q_new,l,delta,0,0);

        res = Qq;
    }
    return res;
}

/*
a = log(y.binary.mean / (1-y.binary.mean))
f = exp(a+e) / (1+exp(a+e)^2 * dnorm(e))
m = integrate: f(-500,500)
*/
/*
// density function for normal distribution N(0,sigma)
double dnorm_mean0( double x )
{
    double p;
    double cst = 0.3989423;
    p = cst * exp(-x*x/2);
    return p;
}

// integration of c
double c_inte( double x, void * params )
{
    double a = *(double*) params;
    double f = exp(a+x) / ( 1+exp(a+x)*exp(a+x) * dnorm_mean0(x) );
    return f;
}

double getCorrectC( double ymean,double sigma)
{
    double alpha = std::log( ymean / (1-ymean));

    gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
//    F.function = &f;
    F.params = &alpha;
    double result, error;
    gsl_integration_qags( &F,0,1,0,1e-7,1000,w,&result,&error );
    double c = result;
    gsl_integration_workspace_free(w);
    
    return c;
}
*/

