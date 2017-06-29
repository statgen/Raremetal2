#ifndef AF1KG_H
#define AF1KG_H

#include <vector>
#include <map>
#include "metaElement.h"
#include "StringArray.h"
#include "MathVector.h"

typedef struct {
	String ref;
	String alt;
	std::vector<double> mafs;	
} OneKGelement;

class af1KG
{
  public:
  	int GetPopCount();
	void loadPopList(String& pop_list_name);
	void loadPosIndex(String& filename,String& Chr, int Start, int End);
	void loadPopVcf( String& pop_vcf_name, bool NoPreloading);
	void setLinearRegressionCoefficients( Vector& gamma,Matrix& X, Vector& Y);
	bool setVariantPopMAF( metaElement& me, String& chr, int position );

  private:
 	int nPop;
 	std::map<int, bool> posIndex; // record the positions that show in summary files
	StringArray popList;
	std::map<String, int> popMap;
	std::map<String, std::map<int, std::vector<OneKGelement> > > refMAFs; // chr->pos->1KG element
};

#endif
