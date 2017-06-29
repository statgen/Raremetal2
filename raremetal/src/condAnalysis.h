#ifndef CONDANALYSIS_H
#define CONDANALYSIS_H

#include "MathMatrix.h"
#include "MathVector.h"
#include "StringArray.h"
#include "StringHash.h"
#include "InputFile.h"
#include <vector>
#include <map>

// for doing conditional condAnalysis

class condAnalysis
{
 public:
 	static bool exclude_missing; // exclude studies with missing data AT the conditioned variant
 	static bool missing_as_zero; // mark studies with missing data AT the conditioned variant as zero
 	
 	void Initialize(int slen);
 	void SetCondMarkers( String& cond_file_name );
 	void SetCondStats(int s, StringArray& tokens, double u, double v, int N, bool flip);
 	void SetSampleSize(std::vector<int>& ss);
	Matrix* Covs; // loading for orginal cov values
	void Run();
	void PrintCondHeader(IFILE& f);
	void PrintSingleLine( int idx,IFILE& f);

 private:
 	void runSanityCheck();
 	int study_number;
 	int marker_number; // total marker number (including #markers to be conditioned)
 	int cond_marker_number;
 	std::map<String,int> condMarkers;
 	// cond markers are at the beginning of each container
 	Vector* Us; // study -> u
 	Vector* Vs;
 	Matrix* CovInvs; // study -> (-1)covariance of (row) all markers (include the conditioned markers) ~ (col) conditioned markers
 	Vector* betaHats;
	// results
	std::vector<double> condUs;
	std::vector<int> Ns;
	std::vector<double> condVs;
};

#endif

