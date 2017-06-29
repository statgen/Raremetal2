#include "metaElement.h"

void initializeMetaElement( metaElement & me, int n)
{
	me.U = 0;
	me.V2 = 0;
	me.N = 0;
	me.AC = 0;
	me.disect = false;
	me.pvalue = -1;
	me.directions.resize(n,'?');
	me.study_N.resize(n,0);
	me.study_mafs.resize(n,0);
	me.skip_pop_correct = false;
}

// minor allele frequency
double getMinorMAF(metaElement& me)
{
	double maf = (double)me.AC/me.N/2;
	if (maf>0.5)
		maf = 1-maf;
	return maf;
}

