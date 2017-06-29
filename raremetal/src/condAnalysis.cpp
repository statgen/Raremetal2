#include "condAnalysis.h"
#include "MetaUtility.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "MathSVD.h"

bool condAnalysis::exclude_missing = false;
bool condAnalysis::missing_as_zero = false;

void condAnalysis::Initialize(int slen)
{
	study_number = slen;
	Us = new Vector [slen];
	Vs = new Vector [slen];
	Covs = new Matrix [slen];
	marker_number = 0;
	Ns.resize(slen,0);
}

void condAnalysis::SetCondMarkers(String& cond_file_name)
{
	IFILE condFile = ifopen(cond_file_name,"r");
	if(condFile==NULL)
		error("Can not open file %s.\n",cond_file_name.c_str());
	int idx = 0;
	while (!ifeof(condFile)) {
		String buffer;
		buffer.ReadLine(condFile);
		if(buffer.FindChar('#')!=-1 || buffer.FindChar(':')==-1)
			continue;
		condMarkers[buffer] = idx;
		idx++;
	}
	ifclose(condFile);
	cond_marker_number = idx;
	marker_number = cond_marker_number;
	// array initialization
	for(int s=0;s<study_number;s++) {
		Us[s].Dimension(cond_marker_number,0);
		Vs[s].Dimension(cond_marker_number,0);
	}
	// result storage initialization
	condUs.resize(marker_number,0);
	Ns.resize(marker_number,0);
	condVs.resize(marker_number,0);
}


void condAnalysis::SetSampleSize(std::vector<int>& ss)
{
	for(int s=0;s<study_number;s++)
		Ns[s] = ss[s];
}

// set U and N at the same time
void condAnalysis::SetCondStats(int s, StringArray& tokens, double u, double v, int N, bool flip)
{
	String marker = tokens[0] + ":" + tokens[1] + ":" + tokens[2] + ":" + tokens[3];
	if (condMarkers.find(marker)==condMarkers.end()) {
		condMarkers[marker] = marker_number;
		Us[s].Push(u);
		Vs[s].Push(v);
		marker_number++;
	}
	else {
		int idx = condMarkers[marker];
		Us[s][idx] = u;
		Vs[s][idx] = v;
	}
}

// derive GX and XX_inv
// calculate condU and condV
void condAnalysis::Run()
{
	// sanity check
	runSanityCheck();

	// derive XX_inv and betaHat
	for(int s=0;s<study_number;s++) {
		SVD svd;
		CovInvs[s] = Covs[s];
		svd.InvertInPlace(CovInvs[s]);
		for(int i=0;i<study_number;i++)
			betaHats[s][i] = CovInvs[s][i].InnerProduct(Us[s]);
	}

	// calculate u and v
	for(int m=0;m<marker_number;m++) {
		double cond_u = 0;
		double cond_v = 0;
		for(int s=0;s<study_number;s++) {
			// derive GX
			Vector GX = Covs[s][m];
			double u = Us[s][m] - GX.InnerProduct(betaHats[s]);
			cond_u += u;
			// calculate condU and V for each marker
			Vector tmp;
			for(int i=0;i<GX.dim;i++)
				tmp.Push(GX.InnerProduct(CovInvs[s][i]));
			double v = Vs[s][m]*Vs[s][m] - tmp.InnerProduct(GX);
			cond_v += v;
		}
		condUs[m] = cond_u;
		condVs[m] = cond_v;
	}
}

void condAnalysis::runSanityCheck()
{
	// constants
	if (study_number<=0)
		error("Wrong study number of %d\n\n",study_number);
	if (study_number==1)
		printf("Warning: only 1 study for conditional analysis...is meta-analysis necessary?\n\n");
	if (cond_marker_number<=0)
		error("Wrong number of markers to be conditioned on: %d\n\n",cond_marker_number);
	if (marker_number<=cond_marker_number)
		error("Wrong total marker number of %d\n\n",marker_number);

	for(int i=0;i<Ns.size();i++) {
		if (Ns[i]<0)
			error("Study #%d has abnormal sample size!\n\n",Ns[i]+1);
		if (Ns[i]==0)
			printf("Warning: study #%d is empty. Skipped for conditional analysis!\n\n",Ns[i]+1);
	}

	// fill missing based on parameters
	if (missing_as_zero) {
		for(int s=0;s<study_number;s++) {
			for(int m=0;m<marker_number;m++) {
				if (Us[s][m]==_NAN_) {
					Us[s][m] = 0;
					Vs[s][m] = 0;
					for(int c=0;c<cond_marker_number;c++)
						CovInvs[s][m][c] = 0;
				}
			}
		}
	}
}

void condAnalysis::PrintCondHeader(IFILE& f)
{
	ifprintf(f,"\tCOND_EFFSIZE\tCOND_EFFSIZE_SD\tCOND_H2\tCOND_PVALUE");
}

void condAnalysis::PrintSingleLine(int i, IFILE& f)
{
	double cond_eff = condUs[i] / condVs[i];
	double cond_sd = 1/sqrt(condVs[i]);
	double cond_h2 = condVs[i]*cond_eff*cond_eff/Ns[i];
	double chisq = condUs[i]*condUs[i]/condVs[i];
	double cond_p = pchisq(chisq,1,0,0);
	bool disect = false;
	while(cond_p==0) {
		disect = true;
		chisq *= 0.999;
		cond_p = pchisq(chisq,1,0,0);
	}
	ifprintf(f,"\t%g\t%g\t%g\t%g",cond_eff,cond_sd,cond_h2,cond_p);
}

