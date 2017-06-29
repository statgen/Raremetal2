#ifndef METAELEMENT_H
#define METAELEMENT_H

#include <vector>
#include "StringArray.h"

typedef struct {
	String ref;
	String alt;
	double U;
	double V2;
	int N;
	int AC;
	double disect; // indicate if pvalue < something
	double pvalue;
// 		double pool_maf;
	std::vector<char> directions;
	std::vector<int> study_N;
	std::vector<double> study_mafs;
	std::vector<double> ref_mafs; // maf in 1000G
	bool skip_pop_correct;
} metaElement;

void initializeMetaElement( metaElement & me, int n);

double getMinorMAF(metaElement& me);

#endif
