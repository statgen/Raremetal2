#ifndef __INITIAL_H__
#define __INITIAL_H__

#include "VcfRecord.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "MathMatrix.h"
#include "StringHash.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "SummaryFileReader.h"
#include "WritePDF.h"
#include "metaElement.h"
#include "af1KG.h"
#include "annoGroups.h"
#include "condAnalysis.h"

#include <map>
#include <vector>

class Meta
{
  public:
	Meta( FILE * plog );
	void Prepare();
	void ConditionalAnalysis();
	void SingleVariantMetaAnalysis();
	void CalculateMetaPvalues();
	void ExactNormPop();

	void GroupMetaAnalysis();

	void ErrorToLog(const char* msg);
        
    //Input/Output options  
	static String summaryFiles; // name of summary file list
	static String covFiles; // name of cov file list
	static String groupFile;
	static String vcfInput;
	static String prefix; // output prefix
	static String cond; // conditional analysis or not
	static bool correctGC;
	static bool outvcf;
	static bool Burden;
	static bool MB; // MB test
	static bool MAB; // MAF as weight burden test
	static bool BBeta; // truncated beta as weight burden test
	static bool SKAT; // SKAT
//	static bool STO;
	static bool VT; // Variant Threshold test
	static bool report;
	static bool fullResult;
//	static bool founderAF;
	static bool dosage;
	static double CALLRATE; // call rate cut off
	static double report_pvalue_cutoff;
	static double MAF_cutoff;
	static double HWE;
	static int marker_col; // index of column of u
	static int cov_col; // index of column of cov
	static double maf_threshold;
	static bool altMAF; // exclude size of studies that do not contain that variant
	static bool RegionStatus; // restrict gene-based test to the specified region
	static String Region; // raw region option
	static String Chr;
	static int Start;
	static int End; // 3 variables to define region
	FILE * log; // log file

	static bool useExactMetaMethod; // Jingjing's exact method for unbalanced studies
	static bool normPop; // normalize population for Jingjing's method
	static String pop_vcf_name; // population vcf for the purpose of pop correction in exact method
	static bool matchOnly; // only pop-correct matched variants
	static double matchDist; // only pop-correct variant with matching distance smaller than this
	static double minMatchMAF;
	static double maxMatchMAF;
	static bool relateBinary; // Dajing's method for related binary samples
	static bool debug; // print out debug info
	static bool NoPlot;
//	static bool popVar;
	static bool NoPreloading;

	static String dosageOptionFile;
	static String pop_list_name;

	//saved single variant information from pooled studies
	StringArray scorefile;
	StringArray covfile;
 	String pdf_filename;
 	
 	std::map<String, std::map<int, std::vector<metaElement> > > variantMap; // chr-> pos->info
 	af1KG AF1KG; // for exact test

  private:
	void openMetaFiles();
	bool poolSingleRecord( int study, double& current_chisq, int& duplicateSNP, bool adjust, String & buffer);
 	bool calculateSinglePvalue( int position, metaElement & me);
 	bool adjustStatsForExact( metaElement& me );
 	void printSingleMetaHeader( String & filename, IFILE & output );
 	void printOutVcfHeader( String & vcf_filename, IFILE & vcfout );
 	void printSingleMetaVariant(String& chr, int position, metaElement& me, IFILE& output, IFILE& vcfout );
 	bool isGroupTestRequired();
 	void plotSingleMetaGC( IFILE & output, bool calc_gc );

 	bool updateYstat( int study );
 	void updateExactCov( int s, int g, int m1, int m2);
 	
 	std::vector< std::vector<metaElement*> > groupInfos; // pointer to metaElement for single var in each group
 	void LoadGroupUandInfos();
	metaElement* getPointerMetaElement( String chr, int position, String ref, String alt);

 	void runCondAnalysis();

 	// general options
 	std::vector<bool> dosageOptions;
 	IntArray FormatAdjust;

 	// some house keeping numerics
 //	int Nvariants;
	int Nsamples;
	int flip_count;
	int skip_count;
	Vector pvalueAll,pvalue1,pvalue5,pvalueAll_cond,pvalue1_cond,pvalue5_cond; // store p-value for plot
 	// /net/fantasia/home/yjingj/METAL/1KG/MAF_1KG.txt
	// for exact method
	Vector Ydelta; // yk - ymean
	Vector Ysigma2; // sigmak(variance), divided by nk-1
	Matrix pGamma; // coefficients for population in each study
	int nPop; // #populations

	// for related binary
	Vector Const_binary; // Dajiang's C for binary
	Vector Ysigma2g;
	Vector Ysigma2e;
	Vector Ymean;

	// record variants in group
	std::map<String, bool> markersExclude; // study:chr:pos:ref:alt, record SNPs excluded in single variant meta-analysis

	std::vector<int> SampleSize;
	Vector GCbyStudy;
	PDF pdf;
	WritePDF writepdf;

	// group test
	annoGroups AnnoGroups;

	// conditional analysis
	condAnalysis CondAna;
	Matrix UsByStudy; // for the transfer purpose to CondAna
	Matrix VsByStudy;
};

/*
	// for annotation
	String target_chr;
	int target_pos;
	int target;
	double target_pvalue;


	for group test, store variant (to test) in each group


	// for cov load
	StringIntHash markerPosHash; // old format

	// temporary structure for load
	//std::map< String, std::map<int,std::bool> > markerPosInGroup; // chr -> position -> bool

	// store group test stats
	std::vector<metaElement*> groupInfos; // pointer to metaElement for single var in each group
	std::vector< Vector > groupUs; // U of each group
	std::vector< Matrix > groupVs; // V of each groupquora 

 old structures
	std::map< String, std::vector<double> > variant_fk; //fk of each variant for each study. will turn into fk-tilda if --normPop. last one is total maf
	std::map< String, std::vector<double> > variant_nk; // nk of each variant for each study. last one is new
//	StringDoubleHash V2; // sum(4*nk*fk^2)
//	StringDoubleHash residual_adj; // sigma tilda outside	
//	StringDoubleHash NkDeltak; // nk * deltak
//	StringDoubleHash regressedTotalAF; // for pop stratification, update f-tilda
	std::map< String, std::vector<double> > af1KG; // 1000g pop af for all markers...


// let's update cond analysis later

	//these are for conditional analysis
	StringArray commonVar,common_chr,common_ref,common_alt;
	StringIntHash conditionVar;
	IntArray common_pos,FormatAdjust;
	Matrix  * XX_inv;
	IntArray * commonVar_study;
	Vector * commonVar_effSize;
	Vector * commonVar_V;
	Vector * commonVar_U;
	Vector * commonVar_betaHat;
	IntArray ** commonVar_markers_in_window;
	Vector ** commonVar_marker_cov;
	bool * cond_status;

	Vector * maf;
	Vector * stats;
	Vector * cond_stats;
	Vector * singlePvalue;
	Vector * singleEffSize;
	Matrix * cov;
	Matrix * cond_cov;

	StringArray chr_plot,geneLabel;
	Vector pos_plot,chisq_before_GC;

	void Prepare();
	void PoolSummaryStat(GroupFromAnnotation & group);
	void Run(GroupFromAnnotation & group);
	void CalculateGenotypeCov(SummaryFileReader & covReader,String chr, int pos,int study,Vector & result);
	double GrabGenotypeCov(SummaryFileReader & covReader,int study,String chr1,String pos1,String chr2,String pos2,String & SNP1, String & SNP2);
	void CalculateXXCov(int study,Matrix & result);

	void setMetaStatics();
	void openMetaFiles();
	void prepareConditionalAnalysis();
	bool setCondMarkers();

	bool updateYstat( int study );
	void UpdateDirection(int & direction_idx,int study,char marker_direction,String & chr_pos,bool exclude);
	void UpdateExcludedMarker(int & study, String & chr_pos, int filter,String markername);
	void UpdateStrIntHash(String & chr_pos, int val, StringIntHash & sihash);
	void UpdateStrDoubleHash(String & chr_pos, double val, StringDoubleHash & sdhash);
	void UpdateACInfo(String & chr_pos,double AC);
	void UpdateStats(int study, String & markerName,double stat,double vstat,bool flip);
	char GetDirection(String & chr_pos,double effsize,bool flip);
	int MatchTwoAlleles(String refalt_current,int & marker_idx,String & chr_pos);
	int MatchOneAllele(String refalt_current, int & marker_idx);
	
	bool poolSingleRecord( int study, double & current_chisq, int & duplicateSNP, bool adjust, String & buffer, SummaryFileReader & covReader );

	bool isDupMarker( String & chr_str, String & chr_pos );
	void setRefAltHashKey( String & refalt_current, StringArray & tokens, int c1, int c2 );
	void setPolyMatch( int & match, String & chr_pos, String & refalt_current, StringArray & tokens, int marker_idx );
	void updateSNPcond( int study, bool flip, int adjust, String & chr_pos, StringArray & tokens, SummaryFileReader & covReader );
	void setPooledAF();
	void printSingleMetaHeader( String & filename, IFILE & output );
	void printOutVcfHeader( String & vcf_filename, IFILE & vcfout );
	void printSingleMetaVariant( GroupFromAnnotation & group, int i, IFILE & output, IFILE & vcfout );
	void annotateSingleVariantToGene( GroupFromAnnotation & group, double pvalue, double cond_pvalue, StringArray & tmp );
	void plotSingleMetaGC( IFILE & output, bool calc_gc );

	// group test functions
	void loadSingleStatsInGroup( GroupFromAnnotation & group );
	void loadSingleCovInGroup( GroupFromAnnotation & group );
	void BurdenAssoc(String method,GroupFromAnnotation & group,Vector *& maf,Vector *& stats,Vector *& con_stats,Matrix *& cov,Matrix *& cond_cov,Vector *& singleEffSize,Vector *& singlePvalue);
	void VTassoc( GroupFromAnnotation & group );
	double VTassocSingle(GroupFromAnnotation & group, Vector & maf_cutoff, IFILE reportOutput, IFILE output, int & g,bool condition,String & method);
	void SKATassoc( GroupFromAnnotation & group );
	void CalculateLambda(Matrix & cov, Vector& weight, double * lambda);

	// pop correction
	void addToMapStrVec( std::map<String, std::vector<double> >& variant, int study, String & markername, double fk, int vsize);
	void SetWeight( String & method, Vector & weight, Vector& maf);
	void FitPgamma();
	void fitpGammaForSingleStudy( int study);
	bool add1kgMafToX( Matrix & X, String & markername, int current_index );
	void setGammaFromRegression( int study, Matrix & X, Vector & Y );
	void load1kgPopMAF();
	void updateRegressedTotalAF( String & markername, double total );
	double getAFtilda( String & markername, double raw_af, int study );
	void addToMapStrVec( std::map<String, Vector>& variant, int study, String & markername, double fk);
	void ErrorToLog(const char* msg);

	// cov matrix
	void updateGroupStats(GroupFromAnnotation& group,int study,int g,bool newFormat);
	void updateSingleVariantGroupStats(GroupFromAnnotation& group,int study,int g,Matrix& cov_i,Matrix& GX,StringArray& chr,StringArray&pos,int m,int gvar_count,bool newFormat);
	void updateSingleVariantGroupStatsOldFormat(GroupFromAnnotation& group,int study,int g,Matrix& cov_i,Matrix& GX,StringArray& chr,StringArray&pos,int loc,int m,int gvar_count,double multiplyFactor);
	void updateSingleVariantGroupStatsNewFormat(GroupFromAnnotation& group,int study,int g,Matrix& cov_i,Matrix& GX,StringArray& chr,StringArray&pos,int loc,int m,int gvar_count,double multiplyFactor);
};
*/

#endif
