#include "Error.h"
#include "MixChidist.h"
#include "MathSVD.h"
#include "Meta.h"
#include "MetaUtility.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "Calculate_mvt_pvalue.h"
#include "My_mvt.h"
#include "InputFile.h"
#include "SummaryFileReader.h"
#include "QuickIndex.h"
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "QuadProg.h"
#include <iterator>
#include <math.h> // pow
#include <stdio.h>
#include <iostream>

String Meta::summaryFiles="";
String Meta::covFiles = "";
String Meta::pop_vcf_name = "";
String Meta::prefix = "";
bool Meta::correctGC = false;
bool Meta::Burden = false;
bool Meta::MB = false;
bool Meta::MAB = false;
bool Meta::BBeta = false;
bool Meta::SKAT = false;
bool Meta::VT = false;
bool Meta::outvcf =false;
bool Meta::fullResult =false;
//bool Meta::report = false;
double  Meta::report_pvalue_cutoff = 1e-06;
//bool Meta::founderAF = false;
double Meta::HWE = 1e-05;
double Meta::CALLRATE = 0.95;
double Meta::MAF_cutoff = 0.05;
int Meta::marker_col = 2;
int Meta::cov_col = 3;
bool Meta::altMAF = false;
//bool Meta::tabix = false;
bool Meta::dosage = false;
bool Meta::RegionStatus = false;
String Meta::Region = "";
String Meta::cond = "";
String Meta::Chr = "";
int Meta::Start = -1;
int Meta::End = -1;
bool Meta::useExactMetaMethod = false; // use Jingjing's exact method
bool Meta::normPop = false; // correct for population stratification
bool Meta::relateBinary = false;
bool Meta::debug = false; // show debug info & debug output
bool Meta::matchOnly = true; // only use matched SNP in 1000G
double Meta::matchDist = 0; // in --normPop, exclude variants with dist > matchDist
//bool Meta::noAdjustUnmatch = false; // do not adjust f~tilda for variants that do not match 1000G

Meta::Meta( FILE * plog )
{
	log = plog;
	Nsamples = -1;
}

//This function read all study names from a file
//and save the names in the StringArray files
void Meta::Prepare()
{
	setMetaStatics();
	openMetaFiles();

	if (dosageOptionFile != "") { // read dosage
		IFILE file = ifopen(dosageOptionFile,"r");
		while(!ifeof(file)) {
			String buffer;
			buffer.ReadLine(file);
			bool status;
			if (buffer=="dosage")
				status = true;
			else if (buffer=="genotype")
				status = false;
			else
				error("Invalid option in dosageOptionFile! Allowed option is dosage or genotype!\n");
			dosageOptions.push_back(status);
		}
		ifclose(file);
	}
	if (dosageOptionFile != ""&& dosageOptions.size() != scorefile.Length())
		error("DosageOptionFile should have the same line number as score files. Please check your input files!\n");

	if (normPop) {
		AF1KG.loadPopList(pop_list_name); // names of populations
		AF1KG.loadPopVcf( pop_vcf_name ); // read 1000g vcf
	}
}

void Meta::ConditionalAnalysis()
{
	setMetaStatics();
	openMetaFiles();

	CondAna.Initialize( scorefile.Length() );
	CondAna.SetCondMarkers( cond );

	// loadUs
	loadCondUs();

	// load Vs
	loadCondVs();

	// analysis and print
	runCondAnalysis();
}

void Meta::loadCondUs()
{
	// load Us
	for(int s=0;s<scorefile.Length();s++) {
		printf("Pooling summary statistics from study %d ...\n",s+1);
		IFILE file;
		file = ifopen(scorefile[s],"r");
		if(file == NULL)
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",scorefile[s].c_str());
		
		String buffer;
		buffer.ReadLine(file);
		bool adjust;
		tellRvOrRmw( buffer, adjust, marker_col, cov_col );

		ifclose(file);
		FormatAdjust.Push(adjust);

		// now read through tabix
		int duplicateSNP = 0;
		SummaryFileReader tbx;
		String tabix_name = scorefile[s] + ".tbi";
		tbx.ReadTabix(tabix_name);
		for(int i=0;i<AnnoGroups.GetGeneNumber();i++) {
			String chr = AnnoGroups.GetChr(0,i);
			int position = AnnoGroups.GetPosition(0,i);
			bool status = tbx.ReadRecord( chr,position);
			if (!status)
				error("cannot get marker at study!\n");
			double current_pvalue;
			status = poolSingleRecord(s, current_pvalue, duplicateSNP, adjust, tbx.buffer);
			if (!status)
				error("pool status");
		}
	}
	CalculateMetaPvalues();
}

// run and print
void Meta::runCondAnalysis()
{
	CondAna.Run();
	// print them out
	IFILE output;
	String filename;
	printSingleMetaHeader( filename, output );
	CondAna.PrintCondHeader(output);
	for(std::map<String, std::map<int, std::vector<metaElement> > >::iterator p1=variantMap.begin(); p1!=variantMap.end(); p1++) {
		String chr = p1->first;
		for(std::map<int, std::vector<metaElement> >::iterator p2=p1->second.begin(); p2!=p1->second.end();p2++) {
			for(int i=0;i<p2->second.size();i++) {
				printSingleMetaVariant(chr,p2->first,p2->second[i],output,output);
				CondAna.PrintSingleLine(i,output);
			}
		}
	}
	ifclose(output);
}

//this function will read through summary statistics of each study and pool the information.
//At the end, single variant meta-analysis will be completed.
// information will be stored for gene-based analysis
void Meta::SingleVariantMetaAnalysis()
{	
	// pool summary stats by reading
	for(int study=0;study<scorefile.Length();study++) {
		printf("Pooling summary statistics from study %d ...\n",study+1);
		String filename = scorefile[study];
		IFILE file;
		file = ifopen(filename,"r");
		if(file == NULL)
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",scorefile[study].c_str());
		
		String buffer;
		buffer.ReadLine(file);
		bool adjust;
		tellRvOrRmw( buffer, adjust, marker_col, cov_col );

		ifclose(file);
		FormatAdjust.Push(adjust);
		
		int duplicateSNP = 0;
		Vector study_pvalues;
		file = ifopen(filename,"r");
		bool pass_header = 0;
		while (!ifeof(file)) {
			buffer.ReadLine(file);
			if (!pass_header) {
				if (SampleSize[study]==-1) {
					if(buffer.Find("##AnalyzedSamples")!=-1) {
						StringArray tokens;
						tokens.AddTokens(buffer,"=");
						SampleSize[study] = tokens[1].AsInteger();
						continue;
					}
				}
				if(buffer.FindChar('#') ==-1 && buffer.Find("CHROM")==-1)
					pass_header = 1;				
			}
			if (!pass_header)
				continue;
			double current_pvalue;
			bool status = poolSingleRecord( study, current_pvalue, duplicateSNP, adjust, buffer);
			if (status)
				study_pvalues.Push(current_pvalue);
		}
		ifclose(file);
		// now calculate gc by study in standard method
		if (study_pvalues.Length()==0) {
			printf("\nWarning: no GC calculated in study #%d!\n",study+1);
			GCbyStudy[study] = 0;
		}
		else
			GCbyStudy[study] = GetGenomicControlFromPvalue(study_pvalues);
	}

	// now calculate meta p value
	CalculateMetaPvalues();

	// print them out
	IFILE output;
	String filename;	
	IFILE vcfout;
	String vcf_filename;
	printSingleMetaHeader( filename, output );
	if (outvcf)
		printOutVcfHeader( vcf_filename, vcfout );		
	for(std::map<String, std::map<int, std::vector<metaElement> > >::iterator p1=variantMap.begin(); p1!=variantMap.end(); p1++) {
		String chr = p1->first;
		for(std::map<int, std::vector<metaElement> >::iterator p2=p1->second.begin(); p2!=p1->second.end();p2++) {
			for(int i=0;i<p2->second.size();i++)
				printSingleMetaVariant(chr,p2->first,p2->second[i],output,vcfout);
		}
	}
	ifclose(vcfout);
	ifclose(output);
}

// calculate meta p value using stored summary stats
void Meta::CalculateMetaPvalues()
{
	if (useExactMetaMethod && normPop)
		ExactNormPop();

	StringArray plot_chrs;
	Vector plot_positions;
	Vector pvalue_all;
	Vector pvalue_001;
	Vector pvalue_005;

	for(std::map<String, std::map<int, std::vector<metaElement> > >::iterator p1=variantMap.begin();p1!=variantMap.end();p1++) {
		for(std::map<int, std::vector<metaElement> >::iterator p2=p1->second.begin();p2!=p1->second.end();p2++) {
			for(int i=0;i<p2->second.size();i++) {
				double maf = (double)p2->second[i].AC/2/p2->second[i].N;
				bool status = calculateSinglePvalue( p2->second[i]);
				if (status) {
					plot_chrs.Push(p1->first);
					plot_positions.Push(p2->first);
					pvalue_all.Push(p2->second[i].pvalue);
					if (maf<0.05)
						pvalue_005.Push(p2->second[i].pvalue);
					if (maf<0.01)
						pvalue_001.Push(p2->second[i].pvalue);
				}
			}
		}
	}

	// now do qq plot
	String title = "Single variant analysis";
	String demo_all = "GC=";
	String demo_005 = demo_all;
	String demo_001 = demo_all;
	if (pvalue_all.Length()==0)
		error("No available p values in single variant meta-analysis. Something goes wrong?\n");
	double gc_all = GetGenomicControlFromPvalue(pvalue_all);
	demo_all += gc_all;
	double gc_005,gc_001;
	if (pvalue_001.Length()==0) {
		gc_001 = 0;
		printf("\nWarning: no MAF<0.01 variants. This group GC = 0!\n");
	}
	else
		gc_001 = GetGenomicControlFromPvalue(pvalue_001);
	if (pvalue_005.Length()==0) {
		gc_005 = 0;
		printf("\nWarning: no MAF<0.05 variants. This group GC = 0!\n");
	}
	else
		gc_005 = GetGenomicControlFromPvalue(pvalue_005);
	demo_005 += gc_005;
	demo_001 += gc_001;
	StringArray geneLabel; // need implementation
	writepdf.Draw(pdf,geneLabel,pvalue_all,pvalue_001,pvalue_005,plot_chrs,plot_positions,title,demo_all,demo_005,demo_001,true);
}

// the normalization step in exact method
// match variants & calculate pop coefficients based on af1KG and metaElement map
void Meta::ExactNormPop()
{
	printf("Match and calculating population correction coefficients...\n\n");
	// compute
	for(int s=0;s<scorefile.Length();s++) {
		String match_id_file_name = prefix + ".debug.match_id_study";
		match_id_file_name += s;
		match_id_file_name += ".txt";
		IFILE match_id_file;
		if (debug) {
			match_id_file = ifopen(match_id_file_name,"w",InputFile::UNCOMPRESSED);
			ifprintf(match_id_file,"#ID MAF Match1KG_MAF\n");
		}
		Matrix X;
		Vector Y;
		int matchonly_skipped = 0;
		int dist_skipped = 0;
		// match variants
		for(std::map<String, std::map<int, std::vector<metaElement> > >::iterator p1=variantMap.begin();p1!=variantMap.end();p1++) {
			for(std::map<int, std::vector<metaElement> >::iterator p2=p1->second.begin();p2!=p1->second.end();p2++) {
				String chr = p1->first;
				for(int i=0;i<p2->second.size();i++) {
					bool status = AF1KG.setVariantPopMAF(p2->second[i],chr,p2->first);
					if (status) { // use this record for regression
						double dist = 0;
						double maf = (double)p2->second[i].AC/p2->second[i].N/2;
						if (debug) {
							String markername = p1->first + ":" + String(p2->first) + ":" + p2->second[i].ref + ":" + p2->second[i].alt;
							ifprintf(match_id_file,"%s %g",markername.c_str(),maf);
							for(int k=0; k<nPop; k++)
								ifprintf(match_id_file," %g",X[i][k]);
							ifprintf(match_id_file,"\n");
						}
						for(int k=0;k<nPop;k++)
							dist += (maf-p2->second[i].ref_mafs[k])*(maf-p2->second[i].ref_mafs[k]);
						dist = sqrt(dist/nPop);
						if (matchDist>0 && dist>matchDist) {
							dist_skipped++;
							continue;
						}
						int n = X.rows;
						X.Dimension(n,nPop);
						for(int j=0;j<nPop;j++) {
							X[n][j] = p2->second[i].ref_mafs[j];
						}
						Y.Push(p2->second[i].study_mafs[s]);
					}
					else
						matchonly_skipped++;
				}
			}
		}
		// compute gamma
		Vector gamma;
		AF1KG.setLinearRegressionCoefficients( gamma,X,Y );
		for(int i=0;i<nPop;i++)
			pGamma[s][i] = gamma[i];
		if (debug)
			ifclose(match_id_file);
		if (matchOnly)
			printf("In pop correction, total %d variants that do not have a match in 1000G were skipped.\n",matchonly_skipped);
		if (matchDist > 0)
			printf("In pop correction, total %d variants with dist>%g were skipped.\n",dist_skipped,matchDist);

	}
}

bool Meta::calculateSinglePvalue( metaElement & me)
{
	// monomorphic
	if (me.AC==0 || me.AC==me.N)
		return false;

	// exact method
	if (useExactMetaMethod)
		bool status = adjustStatsForExact(me);

	// now compute
	double chisq = me.U*me.U/me.V2;
	me.pvalue = pchisq(chisq,1,0,0);
	bool disect;
	while(me.pvalue==0.0) {
		disect=true;
		chisq *= 0.999;
		me.pvalue = pchisq(chisq,1,0,0);
	}
}

bool Meta::adjustStatsForExact( metaElement& me )
{
	bool status = true;

	double exact_maf;
	vector<int>::iterator pnk = me.study_N.begin();
	vector<double>::iterator pfk;

	if (normPop) {
		std::vector<double> fk_tilda;
		if (me.ref_mafs.empty()) { // do not adjust
			fk_tilda.resize(scorefile.Length(),0);
			exact_maf = 0;
			status = false;
		}
		else {
			double f = 0;
			for(int k=0;k<scorefile.Length();k++) {
				double fk = me.study_mafs[k];
				for(int i=0;i<nPop;i++)
					fk -= pGamma[k][i] * me.ref_mafs[i];
				fk_tilda[k] = fk;
				f += fk/me.N*me.study_N[k];
			}
			exact_maf = f;
		}
		pfk = fk_tilda.begin();
	}
	else { // use original maf
		exact_maf = me.AC / me.N;
		pfk = me.study_mafs.begin();
	}

	double nkdeltak = 0;
	double nkdeltakfk = 0;
	double v2 = 0;
	double new_r = 0;
	for(int i=0;i<scorefile.Length();i++) {
		nkdeltak += (*pnk) * Ydelta[i];
		double ac = (*pfk) * me.study_N[i];
		nkdeltakfk += ac * Ydelta[i];
		v2 += ac * (*pfk);
		new_r += (*pnk - 1) * Ysigma2[i] + (*pnk) * Ydelta[i]*Ydelta[i];
		pfk++;
		pnk++;
	}
	new_r /= (me.N-1);
	double maf = (double)me.AC/me.N/2;
	me.U -= 2*maf*nkdeltak + 2*nkdeltakfk;
	me.V2 = new_r * (me.V2 + 4*v2 - 4 *me.N*exact_maf*exact_maf);

	return status;
}


bool Meta::isGroupTestRequired()
{
	bool status = false;
	if (MB||MAB||BBeta||VT||SKAT)
		status = true;
	return status;
}

// perform group test
void Meta::GroupMetaAnalysis()
{
	if(outvcf)
		return;
	
	if(!isGroupTestRequired()) {
		printf("\nWarning: none of the gene-level tests was requested; only single variant meta-analysis was done.\n");
		return;
	}

	AnnoGroups.LoadGroupFile( groupFile );
	
	// load Vs
	for(int s=0;s<scorefile.Length();s++)
		AnnoGroups.LoadCovFile(scorefile[s]);

	// load Us
	LoadGroupUandInfos();
/*
	// perform group test
	String method = "";
	if(Burden) {
		method = "burden";
		BurdenAssoc(method);
	}
	if(MB) {
		method = "MB";
		BurdenAssoc(method);
	}
	if(MAB) {
		method = "MAB";
		BurdenAssoc(method);
	}
	if(BBeta) {
		method = "BBeta";
		BurdenAssoc(method);
	}

	if(VT)
		VTassoc();
	
	if(SKAT)
		SKATassoc();

*/		
}


// retrieve info from variantMap
// U stored in Vector* groupUs
// pointer to metaElement stored in vector<metaElement*> groupInfos
void Meta::LoadGroupUandInfos()
{
	// initialize
	int n_genes = AnnoGroups.GetGeneNumber();
	groupInfos.resize(n_genes);
	for(int g=0;g<n_genes;g++)
		groupInfos[g].resize( AnnoGroups.VariantCountInGenes(g), NULL );

	// set pointers & groupUs
	for(int g=0;g<n_genes;g++) {
		for(int i=0;i<AnnoGroups.VariantCountInGenes(g);i++) {
			String marker = AnnoGroups.GetMarkerName(g,i);
			metaElement* p = getPointerMetaElement( marker );
			std::map<String, std::map<int, std::vector<metaElement> > > variantMap;
//			if (p==NULL)
//				error("%s:%d:%s:%s cannot get metaElement pointer. Something is wrong!\n",annoChrs[g][i].c_str(),annoPositions[g][i],annoRefs[g][i].c_str(),annoAlts[g][i].c_str());
//			groupInfos[g][i] = p;
			AnnoGroups.SetU(g,i,p->U);
		}
	}
}

metaElement* Meta::getPointerMetaElement( String& marker )
{
	StringArray tokens;
	tokens.AddTokens(marker,":");
	metaElement* p = &(variantMap[tokens[0]][tokens[1]][0]);
	// implement multi alleles here
}

/*
void Meta::BurdenAssoc( String& method )
{
	printf("Performing %s tests ...\n",method.c_str());
	
	Vector pvalue_burden,pvalue_burden_cond;
	
	IFILE output;
	String filename;
	openMetaResultFile( prefix, filename, output, method );
	
	IFILE reportOutput;
	String reportFile;
	if(report)
	{
		if(prefix =="")
			reportFile = "meta.tophits."+method+".tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
			reportFile = prefix +  "meta.tophits."+method+".tbl";
		else
			reportFile = prefix + ".meta.tophits."+method+".tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		if(cond!="")
			ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\tCOND_EFFSIZES\tCOND_PVALUES\n");
		else
			ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
	}
	ifprintf(output,"##Method=Burden\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(cond!="")
		if(fullResult)
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
		else
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
	else
		if(fullResult)
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
		else
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
	StringArray chr_plot,geneLabels;
	Vector pos_plot;
	for(int g=0;g<annoGroups.size();g++)
	{
		int n = annoPositions[g].size();
		Vector weight;
		weight.Dimension(n);
		Vector mafs;
		fillGroupMAF(g,mafs);
		SetWeight(method,weight,mafs);
		// for burden test, need to 1/w
		for(int w=0;w<weight.Length();w++)
			weight[w] = 1/weight[w];
		
		// caluclate group p value
		wU = weight.InnerProduct( groupUs[g] );
		Vector wV;
		wV.Dimension(n);	
		for(int i=0;i<n;i++)
			wV[i] = weight.InnerProduct(groupVs[g][i]);
		double wVw = wV.InnerProduct(weight);
		
		if(denominator==0) // unable to perform group test for this gene
			continue;
		double chisq = wU * wU / wVw;
		double pvalue = pchisq(chisq,1,0,0);
		double effSize = numerator/denominator;
		double cond_chisq=_NAN_,cond_effSize=_NAN_,cond_pvalue=_NAN_;

		bool disect=false;
		while(pvalue==0.0)
		{
			disect=true;
			chisq *= 0.999;
			pvalue = pchisq(chisq,1,0,0);
		}
		bool cond_disect =false;
		if(cond!="")
		{
			if(cond_denom==0)
			{
				cond_effSize =0.0;
				cond_pvalue =1.0;
			}
			else
			{
				cond_chisq = cond_num*cond_num/cond_denom;
				cond_effSize = cond_num/cond_denom;
				cond_pvalue = pchisq(cond_chisq,1,0,0);
				while(cond_pvalue==0.0)
				{
					cond_disect=true;
					cond_chisq *= 0.999;
					cond_pvalue = pchisq(cond_chisq,1,0,0);
				}
			}
			pvalue_burden_cond.Push(cond_pvalue);
		}
		
		if(fullResult)
		{
			ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
			
			for(int i=0;i<maf[g].Length()-1;i++)
				ifprintf(output,"%g,",maf[g][i]);
			ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
			
			for(int i=0;i<singleEffSize[g].Length()-1;i++)
				ifprintf(output,"%g,",singleEffSize[g][i]);
			ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
			
			for(int i=0;i<singlePvalue[g].Length()-1;i++)
				ifprintf(output,"%g,",singlePvalue[g][i]);
			ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
			
			if(cond!="")
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
			else
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
		}
		else
		{
			if(cond!="")
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
		}
		if(pvalue <report_pvalue_cutoff && report)
		{
			StringArray variants;
			variants.AddTokens(var,";");
			if(cond!="")
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,cond_disect?"<":"",cond_pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			else
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
		}
		pvalue_burden.Push(pvalue);
		geneLabels.Push(group.annoGroups[g]);
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
	}
	String name = method;
	name += " (maf<";
	name += MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo;
	if(pvalue_burden.Length()>0)
	{
		//Calculate genomic control
		double GC = GetGenomicControlFromPvalue(pvalue_burden);
		demo="GC = ";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,pvalue_burden,chr_plot,pos_plot,name,extraname,demo,true);
		if(cond!="")
		{
			name += "conditional analysis";
			double GC = GetGenomicControlFromPvalue(pvalue_burden_cond);
			demo="GC = ";
			demo += GC;
			writepdf.Draw(pdf,geneLabels,pvalue_burden_cond,chr_plot,pos_plot,name,extraname,demo,true);
		}
	}
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("  done.\n\n");
}
*/

/**** for Meta::Prepare() ******/

// set static stat in meta class
void Meta::setMetaStatics()
{
	// region
	if ( Region != "" ) {
		RegionStatus = true;
		printf("Restrict analysis to region %s!\n", Region.c_str());
		StringArray tf;
		tf.AddTokens(Region, ":-");
		Chr = tf[0];
		Start = tf[1].AsInteger();
		End = tf[2].AsInteger();
	}

	// exact method
	flip_count =0;
	skip_count=0;
	
	SampleSize.Dimension( scorefile.Length() );
	for(int s=0; s<scorefile.Length(); s++)
		SampleSize[s] = -1;

	if (relateBinary) {
		Const_binary.Dimension(scorefile.Length());
		for(int s=0;s<scorefile.Length();s++) {
			double sigma = Ysigma2g[s] + Ysigma2e[s];
		}
	}

	if (useExactMetaMethod) {
		Ydelta.Dimension( scorefile.Length() );
		Ysigma2.Dimension( scorefile.Length() );
		for(int s=0; s<scorefile.Length(); s++) // no need to initialize ydelta
			Ysigma2[s] = -1;
		for(int s=0;s<scorefile.Length();s++) {
			bool y_status = updateYstat(s);
			if (!y_status)
				error("cannot update y stat at study #%d\n", s+1);
		}
		// caculate overall mean & adjust
		int n = 0;
		for(int s=0; s<scorefile.Length(); s++)
			n += SampleSize[s];
		double ymean = 0;
		for(int s=0; s<scorefile.Length(); s++)
			ymean += (double)SampleSize[s] / (double)n * Ydelta[s];
		for(int s=0; s<scorefile.Length(); s++)
			Ydelta[s] -= ymean;
	}
}


// open pdf
// load list of summary files
// load list of cov files if needed
void Meta::openMetaFiles()
{
	// pdf file name
	if(prefix.Last()=='.' || prefix.Last()=='/')
		pdf_filename = prefix + "meta.plots.pdf";
	else
		pdf_filename = prefix + ".meta.plots.pdf";
	pdf.OpenFile(pdf_filename);

	// summary file
	if(summaryFiles!="") {
		IFILE inFile = ifopen(summaryFiles,"r");
		if(inFile==NULL)
			error("FATAL ERROR! Please check file name for --summaryFiles  %s\n",summaryFiles.c_str());
		String buffer;
		while (!ifeof(inFile)) {
			buffer.ReadLine(inFile);
			if(buffer.FindChar('#')!=-1)
				continue;
			scorefile.AddTokens(buffer, "\n");
		}
		ifclose(inFile);
	}
	else
		error("FATAL ERROR! --studyName can not be empty! Usage: --summaryFiles your.list.of.summary.files\n");
	
	// cov file
	if(covFiles!="") {
		IFILE inFile = ifopen(covFiles,"r");
		if(inFile==NULL)
			error("Cannot open file %s\n",covFiles.c_str());
		String buffer;
		while (!ifeof(inFile)) {
			buffer.ReadLine(inFile);
			if(buffer.FindChar('#')!=-1)
				continue;
			covfile.AddTokens(buffer, "\n");
		}
		ifclose(inFile);
		//check if summary files and cov files have the same length
		if(scorefile.Length()!=covfile.Length())
			error("There are %d summary files and %d covariance files. Please check to make sure the same number of files are included in the list.\n");
	}
	else if(!Burden || !MB || !VT || !SKAT)
		error("Covariance files are essential to do gene-level tests. Pleaes ues --covFiles your.list.of.cov.files option.\n");

}


// load cond file
void Meta::loadCondVs()
{
	printf("Loading covariance for conditional analysis...\n");
	AnnoGroups.LoadCondFile(cond);
	for(int s=0;s<scorefile.Length();s++) {
		AnnoGroups.LoadCovFile(scorefile[s]);
		AnnoGroups.ExportCov( CondAna.Covs[s] );
	}
	printf("done.\n\n");
}

/* load into covLociMap
void setCovLociMapFromCondFile()
{
	IFILE condFile = ifopen(cond,"r");
	if(condFile==NULL)
		error("Can not open file %s.\n",cond.c_str());

	String buffer;
	while (!ifeof(condFile)) {
		buffer.ReadLine(condFile);
		if(buffer.FindChar('#')!=-1 || buffer.FindChar(':')==-1)
			continue;
		StringArray tmpMarker;
		tmpMarker.AddTokens(buffer, " \t");
		for(int i=0;i<tmpMarker.Length();i++) {
			StringArray tokens;
			tokens.AddTokens(tmpMarker[i],":");
			int position = tokens[1].AsInteger();
			std::map<String, std::map<int,vector<covInfo> > >::iterator pc = covLociMap.find(tokens[0]);
			if (pc==covLociMap.end()) {
				covLociMap[tokens[0]];
				pc = covLociMap.find(tokens[0]);
			}
			std::map<int,covInfo> >::iterator pi = pc->second.find(position);
			if (pi==pc->second.end()) { // new entry
				pc->second[position];
				pi = pc->second.find(position);
				pi->second.resize(1);
				pi->second[0].ref = tokens[2];
				pi->second[0].alt = tokens[3];
			}
		}
	}
}

// add info to covLociMap from reading group file
void setCovLociMapFromGroupFile()
{
	FILE * file = fopen(groupFile,"r");
	if(file==NULL) {
	  printf("ERROR! Cannot open group file %s.\n",groupFile.c_str());
	  error("ERROR! Cannot open group file %s.\n",groupFile.c_str());
	}
	String buffer;
	while (!feof(file)) {
		buffer.ReadLine(file);
		tokens.AddTokens(buffer,SEPARATORS);
		if (tokens.Length()<=1)
			error("ERROR! Group file has no variant at line: %s\n",buffer.c_str());
		for(int i=1;i<tokens.Length();i++) {
			int position = tokens[1].AsInteger();
			std::map<String, std::map<int,std::vector<covInfo> > >::iterator pc = covLociMap.find(tokens[0]);
			if (pc==covLociMap.end()) {
				covLociMap[tokens[0]];
				pc = covLociMap.find(tokens[0]);
			}
			std::map<int,covInfo> >::iterator pi = pc->second.find(position);
			if (pi==pc->second.end()) { // new entry
				pc->second[position];
				pi = pc->second.find(position);
				pi->second.resize(1);
				pi->second[0].ref = tokens[2];
				pi->second[0].alt = tokens[3];
			}			
		}
	}
}


// read cond marker info from score & cov file
// if no marker info was load, return false
bool Meta::setCondMarkers()
{
	for(int s=0;s<scorefile.Length();s++) {
		SummaryFileReader statReader,covReader;
		String filename = scorefile[s];
		String covFilename = covfile[s];
			
		if(!statReader.ReadTabix(filename))
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",scorefile[s].c_str());
		
		bool adjust;
		String buffer;
		IFILE dup = ifopen(filename,"r");
		buffer.ReadLine(dup);
		tellRvOrRmw( buffer, adjust, marker_col, cov_col );
		ifclose(dup);
		
		if (SampleSize[s]==-1) {
			dup=ifopen(filename,"r");
			while (!ifeof(dup)) {
				buffer.ReadLine(dup);
				//get sample size
				if(buffer.Find("##AnalyzedSamples")!=-1) {
					StringArray tokens;
					tokens.AddTokens(buffer,"=");
					SampleSize[s] = tokens[1].AsInteger();
					break;
				}
			}
			ifclose(dup);
		}
			
		if(!covReader.ReadTabix(covFilename))
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.covariance.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.covariance.file.gz.\n",covfile[s].c_str());
			
		for(int i=0;i<commonVar.Length();i++) {
			//if this variant is genotyped in this study
			bool cov_status = covReader.ReadRecord(common_chr[i],common_pos[i]);
			if(!cov_status)
				continue;
			if(!statReader.ReadRecord(common_chr[i],common_pos[i]))
				continue;
			StringArray record;
			record.AddTokens(statReader.buffer,"\t");
				
			if((record[2]==common_ref[i] && record[3]==common_alt[i]) || (record[3]==common_ref[i] && record[2]==common_alt[i])) {
				double v = record[14-adjust].AsDouble();
				if(v>0) {
					commonVar_study[s].Push(i);
					commonVar_V[s].Push(v*v);
					commonVar_effSize[s].Push(record[15-adjust].AsDouble());
					commonVar_U[s].Push(record[13-adjust].AsDouble());
				}
			}
		}
		
		int dim = commonVar_study[s].Length();
		if(dim!=0)
			cond_status[s] =true;
		else {
			cond_status[s] = false;
			continue;
		}
		commonVar_markers_in_window[s] = new IntArray [dim];
		commonVar_marker_cov[s] = new Vector [dim];
			
		bool cov_status = covReader.ReadTabix(covFilename);
		if ( !cov_status )
			error("Cannot open cov file: %s\n", covFilename.c_str());

		for(int i=0;i<dim;i++) {
			int idx = commonVar_study[s][i];
			if(!covReader.ReadRecord(common_chr[idx],common_pos[idx]))
				continue;
			StringArray tmp_markers;
			tmp_markers.AddTokens(covReader.marker_nearby,",");
			for(int j=0;j<tmp_markers.Length();j++)
				commonVar_markers_in_window[s][i].Push(tmp_markers[j].AsInteger());
			tmp_markers.Clear();
			tmp_markers.AddTokens(covReader.marker_cov,",");
			for(int j=0;j<tmp_markers.Length();j++)
				commonVar_marker_cov[s][i].Push(tmp_markers[j].AsDouble());
		}
		XX_inv[s].Dimension(commonVar_study[s].Length(),commonVar_study[s].Length(),0.0);
		CalculateXXCov(s,XX_inv[s]);
		//calculate beta hat for common variants
		commonVar_betaHat[s].Dimension(dim);
		for(int i=0;i<commonVar_betaHat[s].dim;i++)
			commonVar_betaHat[s][i] = XX_inv[s][i].InnerProduct(commonVar_U[s]);
	}

	// check status of cond markers
	bool status = false;
	for(int i=0;i<scorefile.Length();i++) {
		if(cond_status[i]) {
			status=true;
			break;
		}
	}
	return status;	
}
*/

// read RMW output header that records y
// update y and y-mean
bool Meta::updateYstat( int study )
{
	bool found_y = 0;
	int status = 0;
	String fname = scorefile[study];
	IFILE file = ifopen(fname, "r");
	if (file==NULL)
		error("Cannot open file: %s!\n",fname.c_str());
	while( !ifeof(file) ) {
		String buffer;
		buffer.ReadLine(file);
		if (SampleSize[study]==-1) {
			if (buffer.Find("##AnalyzedSamples")!=-1) {
				StringArray tokens;
				tokens.AddTokens(buffer, "=");
				SampleSize[study] = tokens[1].AsInteger();
			}
		}
/*		if (buffer.Find("#AnalyzedTrait") != -1) {
			StringArray tokens;
			tokens.AddTokens(buffer, "\t ");
			Ysigma2[study] = tokens[7].AsDouble();
			status++;
			if (status==2)
				break;
		}*/
		if (buffer.Find("#")==-1)
			break;
		if (found_y) { // name min 25% median 75th max mean(7th) variance(8th)
			StringArray tokens;
			tokens.AddTokens(buffer, "\t ");
			Ydelta[study] = tokens[6].AsDouble();
			Ysigma2[study] = tokens[7].AsDouble();
			status++;
			found_y = 0;
			if (status==1)
				break;
		}
		if (buffer.Find("##TraitSummaries") != -1) {
			found_y = 1;
			continue;
		}
	}
	if (status ==1)
		return 1;
	else
		return 0;
}

bool matchAlleles(StringArray & tokens, metaElement & me)
{
	if (me.ref==tokens[2] && me.alt==tokens[3])
		return true;
	else
		return false;
}

bool isFlipRecord(StringArray & tokens, metaElement & me)
{
	if (me.ref==tokens[3] && me.alt==tokens[2])
		return true;
	else
		return false;
}


bool matchRefAllele(StringArray& tokens, metaElement& me, bool& flip)
{
	if (me.ref == tokens[2]) {
		flip = false;
		return true;
	}
	if (me.alt == tokens[3]) {
		flip = true;
		return true;
	}
	return false;
}


/**************************************************/


/*** pool summary stat related ****/

// pool single record from score file
// do not calculate p at this time
bool Meta::poolSingleRecord( int study, double & current_pvalue, int & duplicateSNP, bool adjust, String & buffer)
{
	StringArray tokens;
	tokens.AddTokens(buffer, SEPARATORS);
	if (tokens[0].Find("#")!=-1)
		return 0;

	if(tokens[0].Find("chr")!=-1)
		tokens[0] = tokens[0].SubStr(3);
		
	String chr_pos = tokens[0] + ":" + tokens[1];	
	//POOLING STEP1: if fail HWE or CALLRATE then skip this record
	//if a variant has a missing allele but not mornomorphic then exclude this variant without updating the total sample size
	//if(((tokens[2]=="." || tokens[3]==".") && (c1+c2!=0 && c2+c3!=0)) || (tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))

	if (tokens[2]=="0")
		tokens[2] = ".";
	if (tokens[3]=="0")
		tokens[3] = ".";

	//check allele counts to see if the site is monomorphic
	int c1,c2,c3; // genotype counts
	if(dosage)
	{
		c3 = tokens[4].AsDouble()*tokens[6].AsDouble()*tokens[6].AsDouble();
		c2 = tokens[4].AsDouble()*2.0*tokens[6].AsDouble()*(1.0-tokens[6].AsDouble());
		c1 = tokens[4].AsDouble()*(1.0-tokens[6].AsDouble())*(1.0-tokens[6].AsDouble());
	}
	else
	{
		c1 = tokens[10-adjust].AsDouble();
		c2 = tokens[11-adjust].AsDouble();
		c3 = tokens[12-adjust].AsDouble();
	}

	double current_AC;
	int  current_N;
	current_N = c1+c2+c3;
	current_AC = 2*c3+c2;
	double raw_af = (double)current_AC/current_N/2;

	bool is_fail = 0;
	int filter_type = 1;
	if ((tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))
		is_fail = 1;
	if(tokens[8-adjust].AsDouble()<CALLRATE || tokens[9-adjust].AsDouble()<HWE) {
		is_fail = 1;
		filter_type = 0;
	}

	if (is_fail) {
		if(filter_type==0)
			fprintf(log,"Warning: variant %s:%s:%s has at least one allele missing but is polymorphic; and is excluded from meta-analysis.\n",chr_pos.c_str(),tokens[2].c_str(),tokens[3].c_str(),study);
		else if(filter_type==1)
			fprintf(log,"Warning: variant %s:%s:%s from study %d failed to pass HWE or CALLRATE filter and it is excluded from meta analysis.\n",chr_pos.c_str(),tokens[2].c_str(),tokens[3].c_str(),study);
		return false;
	}

	//STEP2: add info to variantMap if not skipped
	bool flip=false;
	int mindex = 0; // multi-allelic index
	if (variantMap.find(tokens[0])==variantMap.end())
		variantMap[tokens[0]];
	int position = tokens[1].AsInteger();
	std::map<int, std::vector<metaElement> >::iterator vp = variantMap[tokens[0]].find(position);
	if (vp==variantMap[tokens[0]].end()) { // this position is not hashed before
		variantMap[tokens[0]][position].resize(1);
		vp = variantMap[tokens[0]].find(position);
	}
	else { // this position is hashed before
		bool is_match = false;
		for(int i=0;i<=vp->second.size();i++) {
			is_match = matchAlleles(tokens,vp->second[i]);
			if (!is_match) {
				flip = isFlipRecord(tokens,vp->second[i]);
				if (flip)
					is_match = true;
			}
			if (is_match)
				break;
		}
		if (is_match) { // check duplicates
			if (vp->second[mindex].study_N[study]>0) {
				duplicateSNP++;
				printf("Duplicate SNP found at %s:%s:%s from study %d. Skipped!\n", chr_pos.c_str(),tokens[2].c_str(),tokens[3].c_str(),study);
				return false;
			}
		}
		else {
			mindex = vp->second.size();
			vp->second.resize(mindex+1);
		}
	}
	// initialize
	initializeMetaElement( vp->second[mindex],scorefile.Length());
	vp->second[mindex].ref = tokens[2];
	vp->second[mindex].alt = tokens[3];
	vp->second[mindex].study_N[study] = current_N;
	if (flip)
		vp->second[mindex].study_mafs[study] = 1-raw_af;
	else
		vp->second[mindex].study_mafs[study] = raw_af;
	// stats
	double u = tokens[13-adjust].AsDouble();
	double v = tokens[14-adjust].AsDouble();
	if(c2+c3==0 || c1+c2==0) { // if monomorphic
		vp->second[mindex].N += current_N; // ac+=0
		vp->second[mindex].directions[study] = '?';
	}
	else { // update normally
		char direct = getDirection(chr_pos,tokens[15-adjust].AsDouble(),flip);
		vp->second[mindex].directions[study] = direct;
		vp->second[mindex].N += current_N;
		if (flip)
			vp->second[mindex].AC += current_N*2-current_AC;
		else
			vp->second[mindex].AC += current_AC;
		double flip_factor = 1;
		if (flip)
			flip_factor = -1;
		double new_u = u * flip_factor;
		double new_v2 = v*v;
		if (useExactMetaMethod) {
			new_u *= Ysigma2[study];
			new_v2 *= Ysigma2[study];
		}
		vp->second[mindex].U += new_u;
		vp->second[mindex].V2 += new_v2;
	}
	
	if(flip)
		flip_count++;
	//if a variant is monomorphic then update the count of sample size and generate warning
	if(c1+c2==0 || c2+c3==0)
		return false;

	if(tokens[14-adjust].AsDouble()>0) {
		double chisq = u*u/v*v;
		current_pvalue = pchisq(chisq,1,0,0);
	}

	// conditional analysis
	if (cond!="" && tokens[14-adjust].AsDouble()>0.0) {
		CondAna.SetCondStats(study,tokens,u,v,current_N,flip);
	}
	return true;
}


// after reading & pooling stat from score file, set overall af of each marker
// results stored in SNPmaf
/*
void Meta::setPooledAF()
{
	StringArray chr_AC,unique_chr,SNPname_AC;
	IntArray pos_AC;		
	//get the unique chromosomes
	for(int i=0;i<directionByChrPos.Capacity();i++) {
		if(!directionByChrPos.SlotInUse(i))
			continue;
		String SNPname = directionByChrPos[i];
		int idx = directionByChrPos.Integer(SNPname);
		StringArray tmp;
		tmp.AddTokens(SNPname,":");
		chr_AC.Push(tmp[0]);
		pos_AC.Push(tmp[1].AsInteger());
		if(directions[idx].Length()<scorefile.Length()) {
			for(int l=directions[idx].Length();l<scorefile.Length();l++)
				directions[idx] += '?';
		}
		SNPname += ':';
		SNPname += refalt[idx];
		SNPname_AC.Push(SNPname);
		if(unique_chr.Find(tmp[0])==-1)
			unique_chr.Push(tmp[0]);
	}
	QuickIndex chr_AC_idx(chr_AC);
	unique_chr.Sort();
	StringArray chr_cp,character_chr;

	// pool
	for(int i=0;i<unique_chr.Length();i++) {
		if(unique_chr[i].AsInteger()<=22 && unique_chr[i].AsInteger()>=1)
			chr_cp.Push(unique_chr[i]);
		else
			character_chr.Push(unique_chr[i]);
	}
	for(int i=0;i<character_chr.Length();i++)
		chr_cp.Push(character_chr[i]);
	unique_chr = chr_cp; //now unique_chr are sorted as 1,2,...,22,X,Y,M...
	chr_cp.Clear();
	character_chr.Clear();
	for(int i=0;i<unique_chr.Length();i++)
	{
		IntArray pos_i;
		StringArray SNPname_i;
		for(int j=0;j<chr_AC.Length();j++) {
			if(chr_AC[chr_AC_idx[j]] == unique_chr[i]) {
				pos_i.Push(pos_AC[chr_AC_idx[j]]);
				SNPname_i.Push(SNPname_AC[chr_AC_idx[j]]);
			}
		}
		QuickIndex pos_i_idx(pos_i);
		for(int j=0;j<pos_i.Length();j++) {
			StringArray tmp;
			tmp.AddTokens(SNPname_i[pos_i_idx[j]],":");
			double AC = usefulAC.Double(tmp[0]+":"+tmp[1]);
				
			// usefulSize has the # of samples to be excluded
			// recSize has # of samples truly added from vcf/ped. Use this when altMAF is toggled
			int N;
			if ( this->altMAF ) {
				N = recSize.Integer( tmp[0]+":"+tmp[1] );
				if ( N==-1 )
					N = 0;
			}
			else { // use default
				N = usefulSize.Integer(tmp[0]+":"+tmp[1]);
				if(N!=-1)
					N = total_N-N;
				else
					N = total_N;
			}
				
			double maf;
			if(founderAF)
				maf = AC/(2.0*N);
			else
				maf = AC/(2.0*N);
				
			int idx = directionByChrPos.Integer(tmp[0]+":"+tmp[1]);
			if(directions[idx].FindChar('+')==-1 && directions[idx].FindChar('-')==-1)
				maf = 0.0;
				
			SNPmaf_maf.Push(maf);
			SNPmaf_name.Push(SNPname_i[pos_i_idx[j]]);
			SNP_effect_N.Push(N);
			SNPmaf.SetInteger(SNPname_i[pos_i_idx[j]],SNPmaf_maf.Length()-1);
		}
	}
}
*/

void Meta::printSingleMetaHeader( String & filename, IFILE & output )
{
	if(prefix =="")
		filename = "meta.singlevar.results";
	else if(prefix.Last()=='.' || prefix.Last()=='/')
		filename = prefix +  "meta.singlevar.results";
	else
		filename = prefix + ".meta.singlevar.results";

	output=ifopen(filename,"w",InputFile::UNCOMPRESSED);	
	ifprintf(output,"##Method=SinglevarScore\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",Nsamples);
	ifprintf(output,"##Genomic_Control=");
	for(int i=0;i<GCbyStudy.Length();i++)
		ifprintf(output," %g",GCbyStudy[i]);
	ifprintf(output,"\n");
	if(cond=="")
		ifprintf(output,"#CHROM\tPOS\tREF\tALT\tN\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tEFFECT_SIZE_SD\tH2\tPVALUE\n");
	else
		ifprintf(output,"#CHROM\tPOS\tREF\tALT\tN\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tEFFECT_SIZE_SD\tH2\tPVALUE\tCOND_EFFSIZE\tCOND_EFFSIZE_SD\tCOND_H2\tCOND_PVALUE\n");
}

void Meta::printOutVcfHeader( String & vcf_filename, IFILE & vcfout )
{
	if(prefix=="")
		vcf_filename = "pooled.variants.vcf";
	else if(prefix.Last()=='.' || prefix.Last()=='/')
		vcf_filename = prefix + "pooled.variants.vcf";
	else
		vcf_filename = prefix + ".pooled.variants.vcf";
		
	vcfout = ifopen(vcf_filename,"w",InputFile::UNCOMPRESSED);
	ifprintf(vcfout,"#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}


// print ith record  of single var result
void Meta::printSingleMetaVariant(String& chr, int position, metaElement& me, IFILE& output, IFILE& vcfout )
{	
	// monomorphic
	double maf = (double)me.AC/me.N/2;
	if (maf==1||maf==0) {
		ifprintf(output,"%s\t%d\t%s\t%s\t%d\t%g\t",chr.c_str(),position,me.ref.c_str(),me.alt.c_str(),me.N,maf);
		for(int i=0;i<me.directions.size();i++)
			ifprintf(output,"?");
		ifprintf(output,"\tNA\tNA\tNA\tNA");
		if(cond!="")
			ifprintf(output,"\tNA\tNA\tNA\tNA\n");
		else
			ifprintf(output,"\n");
		return;
	}

	double effSize = me.U/me.V2;
	double h2 = me.V2*effSize*effSize/me.N;
	double effSize_se = 1.0/sqrt(me.V2);
	ifprintf(output,"%s\t%d\t%s\t%s\t%d\t%g\t",chr.c_str(),position,me.ref.c_str(),me.alt.c_str(),me.N,maf);
	for(int i=0;i<me.directions.size();i++)
		ifprintf(output,"%c",me.directions[i]);
	ifprintf(output,"\t%g\t%g\t%g\t",effSize,effSize_se,h2);
	if (me.disect)
		ifprintf(output,"<");
	ifprintf(output,"%g\n",me.pvalue);
	if(outvcf)
		ifprintf(vcfout,"%s\t%d\t%s\t%s\t%s\t.\t.\tALT_AF=;\n",chr.c_str(),position,me.ref.c_str(),me.alt.c_str(),maf);

	// annotation if needed
//	annotateSingleVariantToGene( group, pvalue, cond_pvalue, tmp );
}

/*Annotate single variants to gene
void Meta::annotateSingleVariantToGene( GroupFromAnnotation & group, double pvalue, double cond_pvalue, StringArray & tmp )
{
	// initialize
	if(geneLabel.Length()==0) {
		target_chr="";
		target_pos=0;
		target=0;
		target_pvalue = _NAN_;
	}	

	bool skip = 0;
	if(pvalue > 0.05/SNPmaf_maf.Length()) {
		if (cond != "") {
			if ( cond_pvalue > 0.05/SNPmaf_maf.Length() )
				skip = 1;
		}
		else
			skip = 1;
	}
	if (skip) {
		geneLabel.Push("");
		return;
	}

	// only annotate those significant signals
	if(target_chr==tmp[0] && target_pos>tmp[1].AsInteger()-1000000) { //annotated already
		if(pvalue<target_pvalue)
		{ //if this is the higher peak then update the annotation of this locus
			String current_anno = group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger());
			int distance = 1000;
			while(current_anno=="" && distance<=1000000)
			{
				current_anno = group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger()+distance);
				if(current_anno=="")
					current_anno = group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger()-distance);
				distance += 1000;
			}
			if(geneLabel[target]!="" && current_anno!="" && geneLabel[target].Find(current_anno)==-1)
				current_anno = geneLabel[target] + "/" + current_anno;
			if(geneLabel[target]!="" && current_anno=="")
				current_anno = geneLabel[target];				
			geneLabel.Push(current_anno);
			geneLabel[target]="";
			target = geneLabel.Length()-1;
			target_pvalue = pvalue;
			target_pos = tmp[1].AsInteger();
			target_chr = tmp[0];
		}
		else //otherwise, leave the original annotation of this locus
			geneLabel.Push("");
	}
	else
	{
		geneLabel.Push(group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger()));
		target_chr = tmp[0];
		target_pos = tmp[1].AsInteger();
		target_pvalue = pvalue;
		target = geneLabel.Length()-1;
	}
}
*/

// single variant GC plot by MAF
void Meta::plotSingleMetaGC( IFILE & output, bool calc_gc )
{
	String title,demo1,demo2,demo3;
	title = "single variant analysis";
	if (pvalueAll.size <= 0)
		error("No available p values. Something goes wrong?\n");
	double GC1 = GetGenomicControlFromPvalue(pvalueAll);
	demo1="GC=";
	demo1 += GC1;
	double GC2 = 0;
	double GC3 = 0;
	if (pvalue1.size <= 0)
		printf("\nWarning: no MAF<0.01 variants. This group GC = 0!\n");
	else
		GC2 = GetGenomicControlFromPvalue(pvalue1);
	demo2="GC=";
	demo2 += GC2;
	if (pvalue5.size <= 0)
		printf("\nWarning: no MAF<0.05 variants. This group GC = 0!\n");
	else
		GC3 = GetGenomicControlFromPvalue(pvalue5);
	demo3="GC=";
	demo3 += GC3;
	StringArray geneLabel, chr_plot, pos_plot;
	writepdf.Draw(pdf,geneLabel,pvalueAll,pvalue1,pvalue5,chr_plot,pos_plot,title,demo1,demo2,demo3,true);
		
	//Calculate genomic control
	if (calc_gc) {
		ifprintf(output,"#Genomic Control for pooled sample is: %g\n",GC1);
			printf("  Genomic Control for all studies are listed below:\n");
		for(int s=0;s<scorefile.Length();s++) {
			printf("  %g\t",GCbyStudy[s]);
			ifprintf(output,"#Genomic Control for study %d is: %g\n",s,GCbyStudy[s]);
		}
	}
}

/****************************************************************************************/
/*

//for Meta::Run
void Meta::loadSingleStatsInGroup( GroupFromAnnotation & group )
{
	for(int g=0;g<group.annoGroups.Length();g++) {
		IntArray del; //del has the SNPs to be deleted from SNPlist
		maf[g].Dimension(group.SNPlist[g].Length());
		stats[g].Dimension(group.SNPlist[g].Length());
		singlePvalue[g].Dimension(group.SNPlist[g].Length());
		singleEffSize[g].Dimension(group.SNPlist[g].Length());
		
		for(int m=0;m<group.SNPlist[g].Length();m++) {
			String newSNP;
			bool flipStatus = false;
			double af=0.0;
			double singleP = getSinglePvalue( g,m );
			
			if(singleP==_NAN_) {
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reasons:monomorphic or failed QC.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}

			af = getSingleMaf(g,m);
			double singleEff = getSingleEff(g,m);
			
			if(af==0.0 || af==1.0) {
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reasons: monomorphic.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}
			
			if(min(af,1-maf) > MAF_cutoff) {
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reason: maf>cutoff.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}
			
			double u = getSingleU(g,m); // be aware of flip

			stats[g][m] = (u);
			maf[g][m] = (af);
			singlePvalue[g][m] = singleP;
			singleEffSize[g][m] = singleEff;
		}
		//delete the SNPs that are not genotyped in any of the studies,
		//or not polymorphic, or SNPs with maf>cutoff.
		if(cond!="") {
			cond_stats[g].Dimension(group.SNPlist[g].Length(),0.0);
			cond_cov[g].Dimension(group.SNPlist[g].Length(),group.SNPlist[g].Length(),0.0);
		}
	}
}

//loop through cov matrices of all studies and update cov
void Meta::loadSingleCovInGroup( GroupFromAnnotation & group )
{
	for(int study=0;study<covfile.Length();study++)
	{
		SummaryFileReader covReader;
		String covFilename = covfile[study];
		setFromRvOrRmwAdjust( FormatAdjust[study], marker_col, cov_col );
		printf("Reading cov matrix from study %d ...\n",study+1);
		String filename = covfile[study];
		IFILE covfile_;
		covfile_  = ifopen(filename,"r");
		if(covfile_ == NULL)
			error("ERROR! Cannot open file: %s! Input cov file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.cov.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.cov.txt.gz\n",filename.c_str());
		Tabix covtabix;
		String tabix_name = filename + ".tbi";
		StatGenStatus::Status libstatus = covtabix.readIndex( tabix_name.c_str() );
		if ( RegionStatus ) {
			if ( libstatus != StatGenStatus::SUCCESS )
				error("Cannot open tabix file %s!\n", tabix_name.c_str());
			bool status = SetIfilePosition( covfile_, covtabix, Chr, Start );
			if ( !status )
				error( "Cannot find position %s:%d-%d in cov file %s!\n", Chr.c_str(), Start, End, filename.c_str() );
		}

		int m=0;
		bool pass_header = 0;
		int newFormat = 0;
		while (!ifeof(covfile_)) {
			String buffer;
			buffer.ReadLine(covfile_);
			if (!pass_header) {
				if (buffer.Find("CHROM")==-1)
					continue;
				// now check new or old format
				StringArray tokens;
				tokens.AddTokens(buffer,"\t ");
				if (tokens[2]=="MARKERS_IN_WINDOW" && tokens[3]=="COV_MATRICES")
					newFormat = 0;
				else if (tokens[2]=="EXP" && tokens[3]=="COV_MATRICES")
					newFormat = 1;
				else if (tokens[2]=="REF" && tokens[3]=="ALT" && tokens[4]=="EXP" && tokens[5]=="COV_MATRICES")
					newFormat = 2;
				else
					error("Covariance matrix is neither new or old format...are you using the right file?\n");
				pass_header = 1;
				continue;
			}
			StringArray tokens;
			tokens.AddTokens(buffer,"\t ");
			if ( RegionStatus ) {
				if ( tokens[1].AsInteger() > End || tokens[0] != Chr ) // out of this region or into another chromosome
				break;
			}
			if (FormatAdjust[study]) // rvtest
				readCovOldFormatLine(study,tokens,m);
			else {
				if (newFormat)
					readCovNewFormatLine(study,tokens,m,newFormat);
				else
					readCovOldFormatLine(study,tokens,m,newFormat);
			}
		}
		ifclose(covfile_);
		printf("done\n");

		//   printf("Updating group stats ...\n");
		//update group statistics
		for(int g=0;g<group.annoGroups.Length();g++)
			updateGroupStats(study,g,newFormat);

		// clear after filling up each study
		markerPosHash.Clear();
		markersExp.Clear(); // for new format
		markersInWindow.Clear(); // for old format
		markersCov.Clear();			
	}	
}
*/
/*
void Meta::BurdenAssoc(String method,GroupFromAnnotation & group,Vector *& maf,Vector *& stats,Vector *& cond_stats,Matrix *& cov,Matrix *& cond_cov,Vector *& singleEffSize,Vector *& singlePvalue)
{
	printf("Performing %s tests ...\n",method.c_str());
	//calculate final results here
	
	Vector pvalue_burden,pvalue_burden_cond;
	
	IFILE output;
	String filename;
	openMetaResultFile( prefix, filename, output, method );
	
	String method_out = method;
	method_out +="_";
	method_out += MAF_cutoff;
	
	IFILE reportOutput;
	String reportFile;
	if(report)
	{
		if(prefix =="")
			reportFile = "meta.tophits."+method+".tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
			reportFile = prefix +  "meta.tophits."+method+".tbl";
		else
			reportFile = prefix + ".meta.tophits."+method+".tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		if(cond!="")
			ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\tCOND_EFFSIZES\tCOND_PVALUES\n");
		else
			ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
	}
	ifprintf(output,"##Method=Burden\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(cond!="")
		if(fullResult)
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
		else
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
	else
		if(fullResult)
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
		else
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
	double numerator=_NAN_,denominator=_NAN_,chisq=0.0,pvalue=_NAN_,cond_num=_NAN_,cond_denom=_NAN_;
	StringArray chr_plot,geneLabels;
	Vector pos_plot;
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		if(maf[g].Length()==0)
			continue;
		
		double average_af = maf[g].Average();
		double min_af = maf[g].Min();
		double max_af = maf[g].Max();
		
		String var;
		for(int i=0;i<maf[g].Length()-1;i++)
			var += group.SNPlist[g][i] + ";";
		var += group.SNPlist[g][maf[g].Length()-1];
		
		Vector weight;
		weight.Dimension(maf[g].Length());
		SetWeight(method,weight,maf[g]);
		// for burden test, need to 1/w
		for(int w=0;w<weight.Length();w++)
			weight[w] = 1/weight[w];
		
		numerator  = weight.InnerProduct(stats[g]);
		Vector tmp;
		tmp.Dimension(group.SNPlist[g].Length());
		
		for(int i=0;i<tmp.Length();i++)
			tmp[i] = weight.InnerProduct(cov[g][i]);
		denominator = tmp.InnerProduct(weight);
			
		if(cond!="") {
			cond_num = weight.InnerProduct(cond_stats[g]);
			for(int i=0;i<tmp.Length();i++)
				tmp[i] = weight.InnerProduct(cond_cov[g][i]);
			cond_denom = tmp.InnerProduct(weight);
		}
		
		if(denominator==0.0)
			continue;
		
		chisq = numerator*numerator/denominator;
		pvalue = pchisq(chisq,1,0,0);
		double effSize = numerator/denominator;
		double cond_chisq=_NAN_,cond_effSize=_NAN_,cond_pvalue=_NAN_;

		bool disect=false;
		while(pvalue==0.0)
		{
			disect=true;
			chisq *= 0.999;
			pvalue = pchisq(chisq,1,0,0);
		}
		bool cond_disect =false;
		if(cond!="")
		{
			if(cond_denom==0)
			{
				cond_effSize =0.0;
				cond_pvalue =1.0;
			}
			else
			{
				cond_chisq = cond_num*cond_num/cond_denom;
				cond_effSize = cond_num/cond_denom;
				cond_pvalue = pchisq(cond_chisq,1,0,0);
				while(cond_pvalue==0.0)
				{
					cond_disect=true;
					cond_chisq *= 0.999;
					cond_pvalue = pchisq(cond_chisq,1,0,0);
				}
			}
			pvalue_burden_cond.Push(cond_pvalue);
		}
		
		if(fullResult)
		{
			ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
			
			for(int i=0;i<maf[g].Length()-1;i++)
				ifprintf(output,"%g,",maf[g][i]);
			ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
			
			for(int i=0;i<singleEffSize[g].Length()-1;i++)
				ifprintf(output,"%g,",singleEffSize[g][i]);
			ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
			
			for(int i=0;i<singlePvalue[g].Length()-1;i++)
				ifprintf(output,"%g,",singlePvalue[g][i]);
			ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
			
			if(cond!="")
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
			else
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
		}
		else
		{
			if(cond!="")
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
		}
		if(pvalue <report_pvalue_cutoff && report)
		{
			StringArray variants;
			variants.AddTokens(var,";");
			if(cond!="")
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,cond_disect?"<":"",cond_pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			else
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
		}
		pvalue_burden.Push(pvalue);
		geneLabels.Push(group.annoGroups[g]);
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
	}
	String name = method;
	name += " (maf<";
	name += MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo;
	if(pvalue_burden.Length()>0)
	{
		//Calculate genomic control
		double GC = GetGenomicControlFromPvalue(pvalue_burden);
		demo="GC = ";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,pvalue_burden,chr_plot,pos_plot,name,extraname,demo,true);
		if(cond!="")
		{
			name += "conditional analysis";
			double GC = GetGenomicControlFromPvalue(pvalue_burden_cond);
			demo="GC = ";
			demo += GC;
			writepdf.Draw(pdf,geneLabels,pvalue_burden_cond,chr_plot,pos_plot,name,extraname,demo,true);
		}
	}
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("  done.\n\n");
}

void Meta::VTassoc( GroupFromAnnotation & group )
{
	printf("Performing Variable Threshold tests ...\n");
	//calculate final results here
	Vector pvalue_VT,pos_plot,cond_pvalue_VT;
	StringArray chr_plot,geneLabels;
	
	IFILE output;
	String filename;
	String method = "VT_";
	openMetaResultFile( prefix, filename, output, method );
	
	method += MAF_cutoff;
	IFILE reportOutput;
	if(report)
	{
		String reportFile;
		if(prefix =="")
		reportFile = "meta.tophits.VT.tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
		reportFile = prefix +  "meta.tophits.VT.tbl";
		else
			reportFile = prefix + ".meta.tophits.VT.tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
	}
	
	ifprintf(output,"##Method=VT\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(fullResult)
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
	else
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
	
	if(cond!="")
		ifprintf(output,"EFFECT_SIZE\tMAF_CUTOFF\tCOND_PVALUE\n");
	else
		ifprintf(output,"\n");
	
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		if(g>1 && g%1000==1)
		printf("Finished analyzing %d genes.\n",g-1);
		
		if(maf[g].Length()==0)
		continue;
		
		if(maf[g].Length()==1)
		{
			if(fullResult)
				ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],singleEffSize[g][0],singlePvalue[g][0],maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
			else
				ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
			pvalue_VT.Push(singlePvalue[g][0]);
			
			if(cond!="")
			{
				String SNPname_noallele;
				StringArray tmp;
				tmp.AddTokens(group.SNPlist[g][0],":");
				SNPname_noallele=tmp[0]+":"+tmp[1];
				double cond_pvalue_=_NAN_,cond_U=_NAN_,cond_V=_NAN_,chisq=_NAN_,cond_effSize_=_NAN_;
				bool disect=false;
				if(conditionVar.Integer(SNPname_noallele)==-1)
				{
					cond_U = SNPstat_cond.Double(SNPname_noallele);
					cond_V = SNP_Vstat_cond.Double(SNPname_noallele);
					chisq = cond_U*cond_U/cond_V;
					cond_pvalue_ = pchisq(chisq,1,0,0);
					cond_effSize_ = cond_U/cond_V;
					while(cond_pvalue_==0.0)
					{
						disect=true;
						chisq *= 0.999;
						cond_pvalue_ = pchisq(chisq,1,0,0);
					}
				}
				else {
					cond_effSize_=0.0;
					cond_pvalue_ = 1.0;
				}
				cond_pvalue_VT.Push(cond_pvalue_);
				ifprintf(output,"%g\t%g\t%s%g",cond_effSize_,maf[g][0],disect?"<":"",cond_pvalue_);
			}
			ifprintf(output,"\n");
			StringArray tmp_SNPname;
			tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
			chr_plot.Push(tmp_SNPname[0]);
			pos_plot.Push(tmp_SNPname[1].AsDouble());
			geneLabels.Push(group.annoGroups[g]);
			continue;
		}
		
		//STEP1: sort maf[g] to find # of cutoffs
		Vector cp_maf;
		cp_maf.Copy(maf[g]);
		cp_maf.Sort();
		Vector maf_cutoff;
		maf_cutoff.Push(cp_maf[0]);
		
		for(int i=1;i<cp_maf.Length();i++)
		{
			if(cp_maf[i]>maf_cutoff[maf_cutoff.Length()-1])
				maf_cutoff.Push(cp_maf[i]);
		} //now unique maf cutoffs are saved in maf_cutoff.
		double pvalue = _NAN_,cond_pvalue = _NAN_;
		bool condition_status = false;
		pvalue = VTassocSingle(group,maf_cutoff,reportOutput,output,g,condition_status,method);
		if(cond!="") {
			condition_status=true;
			cond_pvalue = VTassocSingle(group,maf_cutoff,reportOutput,output,g,condition_status,method);
		}
		pvalue_VT.Push(pvalue);
		if(cond!="")
			cond_pvalue_VT.Push(cond_pvalue);
		
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
		geneLabels.Push(group.annoGroups[g]);
	}
	
	String name = "VT (maf<";
	name +=  MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo="";
	double GC = GetGenomicControlFromPvalue(pvalue_VT);
	demo="GC=";
	demo += GC;
	writepdf.Draw(pdf,geneLabels,pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);
	if(cond!="")
	{
		name += " Conditional Analysis";
		double GC = GetGenomicControlFromPvalue(cond_pvalue_VT);
		demo="GC=";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,cond_pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);
	}
	
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("Done.\n\n");
}

double Meta::VTassocSingle(GroupFromAnnotation & group, Vector & maf_cutoff, IFILE reportOutput, IFILE output, int & g,bool condition,String & method)
{
	double pvalue=_NAN_,chosen_cutoff=_NAN_,chosen_effSize=_NAN_;
	double numerator=0.0,denominator=0.0,t_max=_NAN_;
	Vector weight,tmp,chosen_weight,score;
	Matrix cov_weight;
	weight.Dimension(maf[g].Length());
	tmp.Dimension(group.SNPlist[g].Length());
	cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
	score.Dimension(maf_cutoff.Length());
	for(int i=0;i<maf_cutoff.Length();i++)
	{
		for(int w=0;w<weight.Length();w++)
		{
			if(maf[g][w]<=maf_cutoff[i])
				weight[w]=1.0;
			else
				weight[w]=0.0;
			cov_weight[i][w]=weight[w];
		}
		if(condition)
			numerator = weight.InnerProduct(cond_stats[g]);
		else
			numerator = weight.InnerProduct(stats[g]);
		for(int d=0;d<tmp.Length();d++)
		{
			if(condition)
				tmp[d] = weight.InnerProduct(cond_cov[g][d]);
			else
				tmp[d] = weight.InnerProduct(cov[g][d]);
		}
		denominator = tmp.InnerProduct(weight);
		
		if(denominator != 0.0)
		{
			double t_stat = fabs(numerator/sqrt(denominator));
			score[i]=t_stat;
			if(t_max==_NAN_)
			{
				t_max = t_stat;
				chosen_cutoff = maf_cutoff[i];
				chosen_weight.Copy(weight);
				chosen_effSize = numerator/denominator;
			}
			else
			{
				if(t_stat>t_max)
				{
					t_max = t_stat;
					chosen_cutoff = maf_cutoff[i];
					chosen_weight.Copy(weight);
					chosen_effSize = numerator/denominator;
				}
			}
		}
		else
			score[i]=0.0;
	}
	if(score.Max()==0.0)
	{
		printf("Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
		fprintf(log,"Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
		return pvalue;
	}
	Vector tmp_maf,tmp_eff,tmp_pvalue;
	for(int i=0;i<maf[g].Length();i++)
	{
		if(chosen_weight[i]==1.0)
			tmp_maf.Push(maf[g][i]);
	}
	for(int i=0;i<maf[g].Length();i++)
	{
		if(chosen_weight[i]==1.0)
		{
			tmp_eff.Push(singleEffSize[g][i]);
			tmp_pvalue.Push(singlePvalue[g][i]);
		}
	}
	
	double average_af = tmp_maf.Average();
	double min_af = tmp_maf.Min();
	double max_af = tmp_maf.Max();
	String var;
	for(int i=0;i<maf[g].Length()-1;i++)
	{
		if(chosen_weight[i]==1.0)
			var += group.SNPlist[g][i] + ";";
	}
	if(chosen_weight[maf[g].Length()-1]==1.0)
		var += group.SNPlist[g][maf[g].Length()-1];
	
	//STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
	Matrix cov_U,cov_U_tmp;
	if(condition)
		cov_U_tmp.Product(cov_weight,cond_cov[g]);
	else
		cov_U_tmp.Product(cov_weight,cov[g]);
	Matrix cov_weight_trans;
	cov_weight_trans.Transpose(cov_weight);
	cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
	//Calculate covariance matrix for (T_1 ... T_#cutoff)
	Matrix cov_T;
	cov_T.Dimension(cov_U.rows,cov_U.cols);
	cov2cor(cov_U,cov_T);
	//STEP4: calculate VT pvalue and report.
	int cutoff = maf_cutoff.Length();
	double * lower = new double [cutoff];
	double * upper = new double [cutoff];
	double * mean = new double [cutoff];
	
	for(int i=0;i<cutoff;i++)
	{
		mean[i] = 0.0;
		lower[i] = -t_max;
		upper[i] = t_max;
	}
	
	//Use pmvnorm to calculate the asymptotic p-value
	Vector result;
	pmvnorm(lower,upper,mean,cov_T,false,result);
	
	if(result[0]==-1.0)
	{
		if(!condition)
		{
			if(cond!="")
			{
				if(fullResult)
				{
					ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());
					
					for(int i=0;i<tmp_maf.Length()-1;i++)
						ifprintf(output,"%g,",tmp_maf[i]);
					ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);
					
					for(int i=0;i<tmp_eff.Length()-1;i++)
						ifprintf(output,"%g,",tmp_eff[i]);
					ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);
					
					for(int i=0;i<tmp_pvalue.Length()-1;i++)
						ifprintf(output,"%g,",tmp_pvalue[i]);
					ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);
					
					ifprintf(output,"%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\t",average_af,min_af,max_af,chosen_effSize,chosen_cutoff);
				}
				else
					ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff);
			}
			else
			{
				ifprintf(output,"\n");
			}
		}
		else
		{
			ifprintf(output,"%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\n",chosen_effSize,chosen_cutoff);
		}
	}
	else
	{
		if(1.0-result[0]==0.0)
		{
			//           printf("gene %s has result %g\n",group.annoGroups[g].c_str(),1.0-result[0]);
			printf("Using Shuang's algorithm to calculate MVN pvalue for gene %s ... ",group.annoGroups[g].c_str());
				if(maf_cutoff.Length()>20)
			{
				while(maf_cutoff.Length()>20)
					maf_cutoff.Delete(0);
				
				double numerator,denominator,t_max=_NAN_;
				Vector weight,tmp,chosen_weight,score;
				Matrix cov_weight;
				weight.Dimension(maf[g].Length());
				tmp.Dimension(group.SNPlist[g].Length());
				cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
				for(int i=0;i<maf_cutoff.Length();i++)
				{
					for(int w=0;w<weight.Length();w++)
					{
						if(maf[g][w]<=maf_cutoff[i])
							weight[w]=1.0;
						else
							weight[w]=0.0;
						cov_weight[i][w]=weight[w];
					}
					if(condition)
						numerator = weight.InnerProduct(cond_stats[g]);
					else
						numerator = weight.InnerProduct(stats[g]);
					
					for(int d=0;d<tmp.Length();d++)
					{
						if(condition)
							tmp[d] = weight.InnerProduct(cond_cov[g][d]);
						else
							tmp[d] = weight.InnerProduct(cov[g][d]);
					}
					denominator = tmp.InnerProduct(weight);
					if(denominator != 0)
					{
						double t_stat = fabs(numerator/sqrt(denominator));
						score.Push(t_stat);
						if(t_max==_NAN_)
						{
							t_max = t_stat;
							chosen_cutoff = maf_cutoff[i];
							chosen_weight.Copy(weight);
							chosen_effSize = numerator/denominator;
						}
						else
						{
							if(t_stat>t_max)
							{
								t_max = t_stat;
								chosen_cutoff = maf_cutoff[i];
								chosen_weight.Copy(weight);
								chosen_effSize = numerator/denominator;
							}
						}
					}
					else
						score.Push(0.0);
				}
				if(score.Max()==0.0)
				{
					printf("Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
					fprintf(log,"Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
					return pvalue;
					printf("completed!\n");
				}
				Vector tmp_maf,tmp_eff,tmp_pvalue;
				for(int i=0;i<maf[g].Length();i++)
				{
					if(chosen_weight[i]==1.0)
						tmp_maf.Push(maf[g][i]);
				}
				
				for(int i=0;i<maf[g].Length();i++)
				{
					if(chosen_weight[i]==1.0)
					{
						tmp_eff.Push(singleEffSize[g][i]);
						tmp_pvalue.Push(singlePvalue[g][i]);
					}
				}
				average_af = tmp_maf.Average();
				min_af = tmp_maf.Min();
				max_af = tmp_maf.Max();
				
				String var;
				for(int i=0;i<maf[g].Length()-1;i++)
				{
					if(chosen_weight[i]==1.0)
						var += group.SNPlist[g][i] + ";";
				}
				if(chosen_weight[maf[g].Length()-1]==1.0)
					var += group.SNPlist[g][maf[g].Length()-1];
				//STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
				Matrix cov_U,cov_U_tmp;
				if(condition)
					cov_U_tmp.Product(cov_weight,cond_cov[g]);
				else
					cov_U_tmp.Product(cov_weight,cov[g]);
				Matrix cov_weight_trans;
				cov_weight_trans.Transpose(cov_weight);
				cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
				//Calculate covariance matrix for (T_1 ... T_#cutoff)
				Matrix cov_T;
				cov_T.Dimension(cov_U.rows,cov_U.cols);
				cov2cor(cov_U,cov_T);
				
				pvalue = CalculateMVTPvalue(score,cov_T,t_max);
				printf("completed!\n");
			}
			else
			{
				pvalue = CalculateMVTPvalue(score,cov_T,t_max);
				printf("completed!\n");
			}
		}
		else
		{
			pvalue = 1.0 - result[0];
		}
		
		if((condition && cond!="") || cond=="")
		{
			if(pvalue <report_pvalue_cutoff && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<tmp_maf.Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),pvalue,MAF_cutoff,chosen_cutoff,variants[v].c_str(),tmp_maf[v],tmp_eff[v],tmp_pvalue[v]);
			}
		}
		
		if(cond=="" || (!condition && cond!=""))
		{
			if(fullResult)
			{
				ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());
				
				for(int i=0;i<tmp_maf.Length()-1;i++)
					ifprintf(output,"%g,",tmp_maf[i]);
				ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);
				
				for(int i=0;i<tmp_eff.Length()-1;i++)
					ifprintf(output,"%g,",tmp_eff[i]);
				ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);
				
				for(int i=0;i<tmp_pvalue.Length()-1;i++)
					ifprintf(output,"%g,",tmp_pvalue[i]);
				ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);
				
				ifprintf(output,"%g\t%g\t%g\t%g\t%g\t%g\t",average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
			}
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
			if(cond=="")
				ifprintf(output,"\n");
		}
		
		if(cond!="" && condition)
			ifprintf(output,"%g\t%g\t%g\n",chosen_effSize,chosen_cutoff,pvalue);
		
		if(pvalue>1.0)
		pvalue = 1.0;
	}
	if(lower) delete [] lower;
	if(upper) delete [] upper;
	if(mean) delete [] mean;
	return pvalue;
}

void Meta::SKATassoc( GroupFromAnnotation & group )
{
	printf("Performing SKAT ...\n");
	//calculate Q statistics here
	Vector pvalue_SKAT,pos_plot,cond_pvalue_SKAT;
	StringArray chr_plot,geneLabels;
	IFILE output;
	String filename;
	String method = "SKAT_";
	openMetaResultFile( prefix, filename, output, method );
	
	method += MAF_cutoff;
	IFILE reportOutput;
	if(report)
	{
		String reportFile;
		if(prefix =="")
			reportFile = "meta.tophits.SKAT.tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
			reportFile = prefix +  "meta.tophits.SKAT.tbl";
		else
			reportFile = prefix + ".meta.tophits.SKAT.tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE_DAVIES\tGENE_PVALUE_LIU\tMAF_CUTOFF\tACTUAL_CUTOFF\tVAR\tMAF\tEFFSIZE\tPVALUE\n");
	}
	
	ifprintf(output,"##Method=SKAT\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(fullResult)
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
	else
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
	
	double Qstat,pvalue,pvalue_liu;
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		//      printf("Now working on group %s\n",group.annoGroups[g].c_str());
		if(g%1000==1 && g>1000)
			printf("Finished analyzing %d genes.\n",((int)g/1000)*1000);
		int n = maf[g].Length();
		Vector weight;
		weight.Dimension(n,0);
		if(maf[g].Length()==0)
			continue;
		
		double average_af = maf[g].Average();
		double min_af = maf[g].Min();
		double max_af = maf[g].Max();
		
		String var;
		for(int i=0;i<maf[g].Length()-1;i++)
			var += group.SNPlist[g][i] + ";";
		var += group.SNPlist[g][maf[g].Length()-1];
		
		//get weight based on maf
//		double alpha = 1.0;
//		double beta=25.0;
		String skat_method = "BBeta";
		SetWeight( skat_method, weight,maf[g] );
		Vector tmp,cond_tmp;
		tmp.Dimension(n);
		if(cond!="")
			cond_tmp.Dimension(n);
		for(int i=0;i<n;i++)
		{
			tmp[i] = weight[i]*stats[g][i];
			if(cond!="")
				cond_tmp[i] = weight[i]*cond_stats[g][i];
		}
		Qstat = tmp.InnerProduct(stats[g]);

		double cond_Qstat = _NAN_;
		if(cond!="")
			cond_Qstat = cond_tmp.InnerProduct(cond_stats[g]);
		double * lambda = new double [n];
		CalculateLambda(cov[g],weight,lambda);
		// check lambda for dead loop
		double lambda_sum=0;
		for( int i=0; i<n; i++ )
			lambda_sum += lambda[i];
		if ( lambda_sum < 0.0000000001 ) {
			fprintf(log,"Gene %s lambda sum is zero. Skipped!\n",group.annoGroups[g].c_str());
			continue;
		}

		double Qstat_dav = Qstat;
		double Qstat_liu = Qstat;
		double cond_pvalue=_NAN_,cond_pvalue_liu=_NAN_;
		
		pvalue = MixChidist(lambda, n, Qstat,"Davies");
		
		bool disect_davies=false,disect_liu = false,cond_disect_davies=false,cond_disect_liu=false;
		int disect_itr=0;
		while( (pvalue<=0.0 ||pvalue==2.0 || std::isnan(pvalue)) && disect_itr<10000)
		{
			disect_davies=true;
			Qstat_dav*=0.9999;
			pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
			disect_itr++;
		}
		while((pvalue<=0.0 ||pvalue==2.0 || std::isnan(pvalue)))
		{
			Qstat_dav*=0.99;
			pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
		}
		pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
		disect_itr=0;
		while( (pvalue_liu<=0.0 ||pvalue_liu==2.0 || std::isnan(pvalue_liu)) && disect_itr<10000)
		{
			disect_liu=true;
			Qstat_liu*=0.9999;
			pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
			disect_itr++;
		}
		while((pvalue_liu<=0.0 ||pvalue_liu==2.0 || std::isnan(pvalue_liu)))
		{
			Qstat_liu*=0.99;
			pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
		}
		
		if(cond!="")
		{
			double * lambda = new double [n];
			CalculateLambda(cond_cov[g],weight,lambda);
			Qstat_dav = cond_Qstat;
			Qstat_liu = cond_Qstat;
			
			cond_pvalue = MixChidist(lambda, n, cond_Qstat,"Davies");
			
			int disect_itr=0;
			while( (cond_pvalue<=0.0 ||cond_pvalue==2.0 || std::isnan(cond_pvalue)) && disect_itr<10000)
			{
				cond_disect_davies=true;
				Qstat_dav*=0.9999;
				pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
				disect_itr++;
			}
			while((cond_pvalue<=0.0 ||cond_pvalue==2.0 || std::isnan(cond_pvalue)))
			{
				Qstat_dav*=0.99;
				cond_pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
			}
			cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
			disect_itr=0;
			while( (cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || std::isnan(cond_pvalue_liu)) && disect_itr<10000)
			{
				cond_disect_liu=true;
				Qstat_liu*=0.9999;
				cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
				disect_itr++;
			}
			while((cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || std::isnan(cond_pvalue_liu)))
			{
				Qstat_liu*=0.99;
				cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
			}
			if(lambda) delete [] lambda;
		}
		if(std::isnan(pvalue_liu) || std::isnan(pvalue))
		{
			if(fullResult)
			{
				ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
				
				for(int i=0;i<maf[g].Length()-1;i++)
					ifprintf(output,"%g,",maf[g][i]);
				ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
				
				for(int i=0;i<singleEffSize[g].Length()-1;i++)
					ifprintf(output,"%g,",singleEffSize[g][i]);
				ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
				
				for(int i=0;i<singlePvalue[g].Length()-1;i++)
					ifprintf(output,"%g,",singlePvalue[g][i]);
				ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
				if(cond=="")
					ifprintf(output,"%g\t%g\t%g\t%g\t-\t-\n",average_af,min_af,max_af,Qstat);
				else
					ifprintf(output,"%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n",average_af,min_af,max_af,Qstat,cond_Qstat);
			}
			else {
				if(cond=="")
					ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat);
				else
					ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,cond_Qstat);
			}
			continue;
		}
		//tabulate top results
		if(cond!="")
		{
			if((cond_pvalue <report_pvalue_cutoff && cond_pvalue_liu<report_pvalue_cutoff)  && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			}
			else
				if((pvalue <report_pvalue_cutoff || pvalue_liu<report_pvalue_cutoff)  && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			}
		}
		if(fullResult)
		{
			ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
			
			for(int i=0;i<maf[g].Length()-1;i++)
				ifprintf(output,"%g,",maf[g][i]);
			ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
			
			for(int i=0;i<singleEffSize[g].Length()-1;i++)
				ifprintf(output,"%g,",singleEffSize[g][i]);
			ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
			
			for(int i=0;i<singlePvalue[g].Length()-1;i++)
				ifprintf(output,"%g,",singlePvalue[g][i]);
			ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
			if(cond=="")
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%s%g\n",average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu);
			else
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,cond_Qstat,cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu);
		}
		else
		{
			if(cond=="")
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu);
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,cond_Qstat,cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu);
		}
		
		if(lambda) delete [] lambda;
		
		pvalue_SKAT.Push(pvalue_liu);
		if(cond!="")
			cond_pvalue_SKAT.Push(cond_pvalue_liu);
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
		geneLabels.Push(group.annoGroups[g]);
	}
	
	String name = "SKAT (maf<";
	name +=  MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo = "";
	double GC = GetGenomicControlFromPvalue(pvalue_SKAT);
	demo="GC=";
	demo += GC;
	writepdf.Draw(pdf,geneLabels,pvalue_SKAT,chr_plot,pos_plot,name,extraname,demo,true);
	if(cond!="")
	{
		name += " conditional analysis";
		double GC = GetGenomicControlFromPvalue(cond_pvalue_SKAT);
		demo="GC=";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,cond_pvalue_SKAT,chr_plot,pos_plot,name,extraname,demo,true);
	}
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("Done.\n\n");
}

/*
void Meta::SKAToptimized( GroupFromAnnotation & group)
{
	double rho;
}
*/
/*
void Meta::CalculateLambda(Matrix & cov,Vector& weight, double * lambda)
{
	int n = cov.rows;
	//calculat sqrt(V)
	Eigen::MatrixXd cov_eigen(n,n);
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
		cov_eigen(i,j) = cov[i][j];
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_cov(cov_eigen, Eigen::ComputeThinU);
	Eigen::MatrixXd final_eigen(n,n);
	Eigen::MatrixXd final_eigen_rhs(n,n);
	Eigen::MatrixXd final_eigen_lhs(n,n);
	Eigen::VectorXd tmp(n);
	for(int i=0;i<n;i++)
		tmp[i]=sqrt(svd_cov.singularValues()[i]);
	Eigen::VectorXd weight_eigen(n);
	for(int i=0;i<n;i++)
		weight_eigen[i] = weight[i];
	final_eigen_rhs = svd_cov.matrixU()*tmp.asDiagonal()*svd_cov.matrixU().transpose();
	final_eigen_lhs = final_eigen_rhs*weight_eigen.asDiagonal();
	final_eigen = final_eigen_lhs*final_eigen_rhs;
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(final_eigen, Eigen::ComputeThinU);
	const Eigen::VectorXd S = svd.singularValues();
	for(int i=0;i<n;i++)
		lambda[i] = fabs(S[i]);
}
*/
// print error to both log and std err
void Meta::ErrorToLog( const char* msg )
{
	fprintf(log, "Error [Meta.cpp]: ");
	fprintf(log, msg);
	fprintf(log, "\n");
	error(msg);
}


void Meta::SetWeight( String & method, Vector & weight, Vector& maf )
{
	if(method=="burden")// equal weight
		for(int w=0;w<weight.Length();w++)
				weight[w] = 1.0;
	else if(method=="MB") // weight by 1/sqrt( maf*(1-maf) )
		for(int w=0;w<weight.Length();w++)
			weight[w] = sqrt(maf[w]*(1.0-maf[w]));
	else if (method="MAB") {
		for(int w=0;w<weight.Length();w++)
			weight[w] = maf[w];
	}
	else if (method="BBeta") { // truncated beta
		double alpha = 0.5;
		double beta = 0.5;
		double f0 = 2 / Nsamples; // truncate at 4 alleles
		for(int w=0;w<weight.Length();w++) {
			double xmaf = maf[w];
			if (xmaf>0 && xmaf<f0)
				xmaf = f0;
			if (xmaf>0 && xmaf>(1-f0))
				xmaf = 1-f0;
			double beta_density = GetBetaDensity(alpha,beta,xmaf);
			weight[w] = (beta_density*beta_density);
		}
	}
	else
		error("Invalid weight %s!\n",method.c_str());
}

// update v in cov matrix in exact
void Meta::updateExactCov( int study, int m, int s, StringArray& chr,StringArray& pos,Matrix& cov_i)
{
	String markername = chr[m] + ":" + pos[m];
	String markername2 = chr[s] + ":" + pos[s];
	std::map< String, std::vector<double> >::iterator pnk = variant_nk.find(markername);
	/*if (p == variant_nk.end()) {
		markername = markerName + ":" + alt_allele[s] + ":"+ref_allele[s];
		p = variant_nk.find(markername);
	}*/
	if( pnk==variant_nk.end())
		error("error2137\n");
	std::map< String, std::vector<double> >::iterator pfk = variant_fk.find(markername);
	std::map< String, std::vector<double> >::iterator pnk2;
	std::map< String, std::vector<double> >::iterator pfk2;
	if (m==s) {
		pnk2 = pnk;
		pfk2 = pfk;
	}
	else {
		pnk2 = variant_nk.find(markername2);
		pfk2 = variant_fk.find(markername2);
	}
	int nk = pnk->second[study];
	double nkfk2 = nk * pfk->second[study] * pfk2->second[study];
	double new_r = pnk->second[scorefile.Length()];
	double ff = pfk->second[scorefile.Length()]*pfk2->second[scorefile.Length()];
	// Vexact=Vrmw*sigma4
	cov_i[m][s] = new_r*(cov_i[m][s]*Ysigma2[study] - 4*nk*ff + 4*nkfk2);
}


/******* update covariance matrix for group test ****/
void Meta::updateGroupStats( GroupFromAnnotation& group,int study,int g,bool newFormat)
{
	//printf("doing group %d\n",g);
	int gvar_count = group.SNPlist[g].Length();
	StringArray chr,pos;
	for(int i=0;i<gvar_count;i++) {
		StringArray tmp;
		tmp.AddTokens(group.SNPNoAllele[g][i],":");
		chr.Push(tmp[0]);
		pos.Push(tmp[1]);
	}
	//now pos has all the positions of markers in group g.
	Matrix cov_i,GX;
	cov_i.Dimension(gvar_count,gvar_count,0.0);
	if(cond!="")
		GX.Dimension(gvar_count,XX_inv[study].cols,0.0);

	for(int m=0;m<gvar_count;m++)
		updateSingleVariantGroupStats( group,study,g,cov_i,GX,chr,pos,m,gvar_count,newFormat);

	cov[g].Add(cov_i);
	if(cond!="") {
		cond_cov[g].Add(cov_i);
		Matrix GX_trans,tmp,extra_cov_i;
		GX_trans.Transpose(GX);
		tmp.Product(GX,XX_inv[study]);
		extra_cov_i.Product(tmp,GX_trans);
		extra_cov_i.Multiply(-1.0);
		cond_cov[g].Add(extra_cov_i);
	}
}

void Meta::updateSingleVariantGroupStats( GroupFromAnnotation& group,int study,int g,Matrix& cov_i,Matrix& GX,StringArray& chr,StringArray&pos,int m,int gvar_count,bool newFormat)
{
	int loc = markerPosHash.Integer(group.SNPNoAllele[g][m]);
	if(loc==-1)
		return;
	int skip = SNPexclude.Integer(study+":"+group.SNPNoAllele[g][m]);
	if(skip!=-1)
		return;			
	int flip = flipSNP.Integer(study+":"+group.SNPNoAllele[g][m]);
	double multiplyFactor=1.0;
	if(flip!=-1)
		multiplyFactor=-1.0;
	//read through markersInWindow and find the selected markers
	if (!newFormat)
		updateSingleVariantGroupStatsOldFormat(group,study,g,cov_i,GX,chr,pos,loc,m,gvar_count,multiplyFactor);
	else
		updateSingleVariantGroupStatsNewFormat(group,study,g,cov_i,GX,chr,pos,loc,m,gvar_count,multiplyFactor);
}

void Meta::updateSingleVariantGroupStatsOldFormat(GroupFromAnnotation& group,int study,int g,Matrix& cov_i,Matrix& GX,StringArray& chr,StringArray&pos,int loc,int m,int gvar_count,double multiplyFactor)
{
	StringArray markers,markerscov;
	markers.AddTokens(markersInWindow[loc-1],",");
	markerscov.AddTokens(markersCov[loc-1],",");
	for(int s=m;s<gvar_count;s++) {
		int p = markers.Find(pos[s]);
		if(p==-1)
			return;
		String markerName = study+":"+chr[s]+":"+pos[s];
		//String markername = markerName + ":"+ref_allele[s] + ":"+alt_allele[s];
		//if the marker in window is supposed to be excluded due to non-consistent ref/alt allele
		//then skip
		int skip = SNPexclude.Integer(markerName);
		if(skip!=-1)
			return;
		int flip = flipSNP.Integer(markerName);
		double factor=1.0;
		if(flip!=-1)
			factor=-1.0;
		cov_i[m][s]= multiplyFactor*factor*markerscov[p].AsDouble()*SampleSize[study];
		if (useExactMetaMethod)
			updateExactCov( study,m,s,chr,pos,cov_i );
	}
	// fill in GX
	if (cond!="") {
		StringArray name;
		name.AddTokens(group.SNPlist[g][m],":");
		cond_stats[g][m] = SNPstat_cond.Double(name[0]+":"+name[1]);
		String pos_str;
		for(int s=0;s<commonVar_study[study].Length();s++) {
			if(pos[m].AsInteger()<common_pos[commonVar_study[study][s]]) { // direct fill in
				pos_str = common_pos[commonVar_study[study][s]];
				int p = markers.Find(pos_str);
				if(p==-1)
					continue;
				GX[m][s]= markerscov[p].AsDouble()*SampleSize[study];
			}
			else { // need to put the smaller one first
				int loc = markerPosHash.Integer(common_chr[commonVar_study[study][s]]+":"+common_pos[commonVar_study[study][s]]);
				//If this SNP is not genotpyed in this study, then skip
				if(loc==-1)
					continue;
				//read through markersInWindow and find the selected markers
				pos_str = pos[m];
				StringArray new_markers, new_markerscov;
				new_markers.AddTokens(markersInWindow[loc-1],",");
				new_markerscov.AddTokens(markersCov[loc-1],",");
				int p = new_markers.Find(pos_str);
				if (p==-1)
					continue;
				GX[m][s]= new_markerscov[p].AsDouble()*SampleSize[study];
			}
		}
	}
}

void Meta::updateSingleVariantGroupStatsNewFormat(GroupFromAnnotation& group,int study,int g,Matrix& cov_i,Matrix& GX,StringArray& chr,StringArray&pos,int loc,int m,int gvar_count,double multiplyFactor)
{
	Vector markerscov;
	addNewFormatCov( markersExp[loc-1], markersCov[loc-1],markerscov);
	for(int s=m;s<gvar_count;s++) {
		String mkname = chr[s] + ":" + pos[s];
		int p = markerPosHash.Integer( mkname );
		if(p==-1)
			return;
		p--;
		p -= m;
		String markerName = study+":"+chr[s]+":"+pos[s];
		int skip = SNPexclude.Integer(markerName);
		if(skip!=-1)
			return;
		int flip = flipSNP.Integer(markerName);
		double factor=1.0;
		if(flip!=-1)
			factor=-1.0;
		if (p >= markerscov.Length()) // zeros in the tail
			cov_i[m][s]= 0;
		else
			cov_i[m][s]= multiplyFactor*factor*markerscov[p]*SampleSize[study];
		if (useExactMetaMethod)
			updateExactCov(study,m,s,chr,pos,cov_i);
	}
	// fill in GX
	if (cond!="") {
		String pos_str;
		StringArray name;
		name.AddTokens(group.SNPlist[g][m],":");
		cond_stats[g][m] = SNPstat_cond.Double(name[0]+":"+name[1]);
		for(int s=0;s<commonVar_study[study].Length();s++) {
			if(pos[m].AsInteger()<common_pos[commonVar_study[study][s]]) { // direct fill in
				pos_str = common_pos[commonVar_study[study][s]];
				String mkname = chr[s] + ":" + pos[s];
				int p = markerPosHash.Integer(mkname);
				if(p==-1)
					continue;
				p-=m;
				GX[m][s]= markerscov[p]*SampleSize[study];
			}
			else { // need to put the smaller one first
				int loc = markerPosHash.Integer(common_chr[commonVar_study[study][s]]+":"+common_pos[commonVar_study[study][s]]);
				//If this SNP is not genotpyed in this study, then skip
				if(loc==-1)
					continue;
				//read through markersInWindow and find the selected markers
				pos_str = pos[m];
				StringArray new_markers;
				Vector new_markerscov;
				String mkname = chr[s] + ":" + pos[s];
				int p = markerPosHash.Integer( mkname );
				p -= m;
				addNewFormatCov( markersExp[loc-1],markersCov[loc-1],new_markerscov);
				if (p==-1)
					continue;
				GX[m][s]= new_markerscov[p]*SampleSize[study];
			}
		}
	}	
}


