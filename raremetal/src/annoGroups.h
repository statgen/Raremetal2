#ifndef ANNOGROUPS_H
#define ANNOGROUPS_H

#include <map>
#include <vector>
#include "StringArray.h"
#include "MathMatrix.h"

class annoGroups
{
  public:
 	void LoadGroupFile( String& group_file_name );
	void LoadCovFile(String& cov_file_name,int _adjust, int _marker_col, int _cov_col,int _sampleSize);
	void LoadCondFile(String& cond_file_name);
	void FlipAllele(int g,int i);
	void SetU(int g, int i, double u);
	void SetV(int g, int i, double v2);
	void SetMaf(int g, int i, double maf);
	String GetMarkerName(int g, int i);
	int VariantCountInGenes(int g);
	int GetGeneNumber();
	String GetChr(int g,int i);
	int GetPosition(int g,int i);
	String GetRef(int g,int i);
	String GetAlt(int g,int i);
	double GetOneCov(int g,int m1,int m2);
	void UpdateCovValue(double new_val,int g,int m1,int m2);
	void LoadSingleMarker(int g,String& chr,int position,String& ref,String& alt);

	void MakeFullMatrix();
	void RemoveMissingData();

	void ExportCov( Matrix& cov, int g);
	void GroupTest(String& method, String& prefix);

// 	void GetFourParts( int index, String& chr, int& position, String& ref, String& alt );

  private:
  	int cond_number; // #markers to be conditioned. Covariance will only be store between markers and these cond_markers. Always the first #cond_number of anno*
  	int marker_number; // total marker number to be analyzed. include cond_number
 	std::vector<String> annoGenes;
	std::vector<std::vector<String> > annoChrs; // group -> element of groups
	std::vector<std::vector<int> > annoPositions;
	std::vector<std::vector<String> > annoRefs;
	std::vector<std::vector<String> > annoAlts;
	void initializeMarkerIndex();
	bool isInMarkerIndex( StringArray& tokens, bool& flip_status );
	void loadCovStrings(String& cov_file_name);
	void setMarkersInCovs();
	void clearLoadData();
	void initializeDataStorage();
	bool addMarkerGroupToAnno( int offset, StringArray& tmp );
	void addNewFormatCov( int mexp,String& cov_str,StringArray& covs);

	void printGroupResult(int g, IFILE& f);
	void runVt();
	void runSKAT();
	void setWeight(String& method, int g);
	void runBurdenTest(String& method);
		
//	std::map<String, std::map<int,std::vector<String> > > annoVariantList; // to check if cov this variant need to be stored: marker->exp,
//	std::vector<bool> annoCondStatus; // if true, it's conditional analysis instead of group test

	// need to clear before load each cov file
	std::map< String, std::map<String, std::pair<int,int> > > markerIndex; // chr -> pos:ref:alt -> in-genome index, vector index in RM
	std::map<String, bool> flip_allele_map; // for old format,record the flipped markers due to MAF in annoGroups::FlipAllele. key: chr:pos
	std::vector<bool> markersFlip; // indicate if the marker is flipped or not
	std::vector<int> markersExp; // for new format, exp of markers
	StringArray markersInWindow; // for old format
	std::vector<String> markersCov; // string of covs in the window

	// data storage (group->variants)
	Vector* groupUs; // here u is meta-analyzed single variant U
	Vector* groupVs;
	Matrix* Covs;

	// for performing group test
	Vector* Mafs;
	Vector Weight;
	Vector metaUs; // each group has one U
	Vector metaVs;

	// constants&stats
	bool is_cov_half_matrix = true;
	bool newFormat;
	int sampleSize;
	int adjust;
	int marker_col;
	int cov_col;
};

#endif