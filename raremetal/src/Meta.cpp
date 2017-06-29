#include "Error.h"
#include "MathSVD.h"
#include "Meta.h"
#include "MetaUtility.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "Calculate_mvt_pvalue.h"
#include "My_mvt.h"
#include "WritePDF.h"
#include "InputFile.h"
#include "SummaryFileReader.h"
#include "QuickIndex.h"
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <iterator>
#include <math.h> // pow
#include <stdio.h>
#include <iostream>

String Meta::summaryFiles="";
String Meta::covFiles = "";
String Meta::groupFile = "";
String Meta::pop_vcf_name = "";
String Meta::pop_list_name = "";
String Meta::dosageOptionFile = "";
bool Meta::report = false;
String Meta::vcfInput = "";
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
int Meta::marker_col = 2;
int Meta::cov_col = 3;
double Meta::maf_threshold = 0.5;
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
bool Meta::NoPlot = false; // do not generate Manhattan plot
//bool Meta::matchOnly = true; // only use matched SNP in 1000G
double Meta::matchDist = 0.5; // in --normPop, exclude variants with dist > matchDist
double Meta::minMatchMAF = 0;
double Meta::maxMatchMAF = 1;
//bool Meta::popVar = false; // load summary file position first then load only these positions in 1000G
bool Meta::NoPreloading = false;

Meta::Meta( FILE * plog )
{
	log = plog;
	Nsamples = 0;
}

//This function read all study names from a file
//and save the names in the StringArray files
void Meta::Prepare()
{
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
		if (!NoPreloading) {
			printf("Preloading variant positions to speed up loading population vcfs...\n");
			for(int s=0;s<scorefile.Length();s++)
				AF1KG.loadPosIndex(scorefile[s],Chr,Start,End);
			printf("  done.\n\n");
		}
		AF1KG.loadPopVcf( pop_vcf_name,NoPreloading); // read 1000g vcf
		nPop = AF1KG.GetPopCount();
	}
}

void Meta::ConditionalAnalysis()
{

	// load covariance
	AnnoGroups.LoadCondFile(cond);
	for(int s=0;s<scorefile.Length();s++) {
		AnnoGroups.LoadCovFile(scorefile[s],FormatAdjust[s],marker_col,cov_col,SampleSize[s]);
		AnnoGroups.ExportCov( CondAna.Covs[s],0);
	}
	printf("done.\n\n");

	// analysis and print
	runCondAnalysis();
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
	if (cond!="") {
		CondAna.Initialize( scorefile.Length() );
		CondAna.SetCondMarkers( cond );
	}
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
		Vector study_chisqs;
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
			double current_chisq;
			bool status = poolSingleRecord( study, current_chisq, duplicateSNP, adjust, buffer);
			if (status)
				study_chisqs.Push(current_chisq);
		}
		ifclose(file);
		if (duplicateSNP>0)
			printf("In study #%d, %d SNPs are duplicated and the duplicates were skipped!\n\n",study+1,duplicateSNP);
		// now calculate gc by study in standard method
		if (study_chisqs.Length()==0) {
			printf("\nWarning: no GC calculated in study #%d!\n",study+1);
			GCbyStudy[study] = 0;
		}
		else {
			study_chisqs.Sort();
			GCbyStudy[study] = study_chisqs[0.5] / 0.456;
		}
	}

	for(int s=0;s<SampleSize.size();s++)
		Nsamples += SampleSize[s];

	if (cond!="")
		CondAna.SetSampleSize(SampleSize);

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
	if (outvcf)
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
				double minor_maf = getMinorMAF(p2->second[i]);
				bool status = calculateSinglePvalue( p2->first, p2->second[i]);
				if (status) {
					plot_chrs.Push(p1->first);
					plot_positions.Push(p2->first);
					pvalue_all.Push(p2->second[i].pvalue);
					if (minor_maf<0.05)
						pvalue_005.Push(p2->second[i].pvalue);
					if (minor_maf<0.01)
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
	StringArray geneLabel;
	if (!NoPlot)
		writepdf.Draw(pdf,geneLabel,pvalue_all,pvalue_001,pvalue_005,plot_chrs,plot_positions,title,demo_all,demo_005,demo_001,false);
}

// the normalization step in exact method
// match variants & calculate pop coefficients based on af1KG and metaElement map
void Meta::ExactNormPop()
{
	printf("Match and calculating population correction coefficients...\n\n");
	pGamma.Dimension(scorefile.Length(),nPop,0);
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
		int maf_skipped = 0;;
		int matchonly_skipped = 0;
		int dist_skipped = 0;
		// match variants
		for(std::map<String, std::map<int, std::vector<metaElement> > >::iterator p1=variantMap.begin();p1!=variantMap.end();p1++) {
			for(std::map<int, std::vector<metaElement> >::iterator p2=p1->second.begin();p2!=p1->second.end();p2++) {
				String chr = p1->first;
				for(int i=0;i<p2->second.size();i++) {
					if (p2->second[i].N==0) {
						p2->second[i].skip_pop_correct = true;
						continue;
					}
					double maf = (double)p2->second[i].AC/p2->second[i].N/2;
					double minor_maf = maf;
					if (maf>0.5)
						minor_maf = 1-maf;
					if (minor_maf>maxMatchMAF||minor_maf<minMatchMAF) {
						p2->second[i].skip_pop_correct = true;
						maf_skipped++;
						continue;
					}
					bool status = AF1KG.setVariantPopMAF(p2->second[i],chr,p2->first);
					if (status) { // use this record for regression
						double dist = 0;
						for(int k=0;k<nPop;k++)
							dist += (maf-p2->second[i].ref_mafs[k])*(maf-p2->second[i].ref_mafs[k]);
						dist = sqrt(dist/nPop);
						if (matchDist>0 && dist>matchDist) {
							p2->second[i].skip_pop_correct = true;
							dist_skipped++;
							continue;
						}
						int n = X.rows;
						X.Dimension(n+1,nPop);
						for(int j=0;j<nPop;j++) {
							X[n][j] = p2->second[i].ref_mafs[j];
						}
						Y.Push(p2->second[i].study_mafs[s]);
						if (debug) {
							String markername = p1->first + ":" + p2->first + ":" + p2->second[i].ref + ":" + p2->second[i].alt;
							ifprintf(match_id_file,"%s %g",markername.c_str(),maf);
							for(int k=0; k<nPop; k++)
								ifprintf(match_id_file," %g",X[n][k]);
							ifprintf(match_id_file,"\n");
						}
					}
					else
						matchonly_skipped++;
				}
			}
		}
		if (minMatchMAF!=0||maxMatchMAF!=1)
			printf("In pop correction of study #%d, total %d variants with maf not between %g-%g were not corrected.\n",s+1,maf_skipped,minMatchMAF,maxMatchMAF);
		printf("In pop correction of study #%d, total %d variants that do not have a match in 1000G were not corrected.\n",s+1,matchonly_skipped);
		if (matchDist > 0)
			printf("In pop correction of study #%d, total %d variants with dist>%g were not corrected.\n",s+1,dist_skipped,matchDist);
		// compute gamma
		Vector gamma;
		AF1KG.setLinearRegressionCoefficients( gamma,X,Y );
		for(int i=0;i<nPop;i++)
			pGamma[s][i] = gamma[i];
		if (debug)
			ifclose(match_id_file);
	}
}

bool Meta::calculateSinglePvalue( int position, metaElement & me)
{
	// monomorphic
	if (me.AC==0 || me.AC==me.N*2)
		return false;

	// exact method
	if (useExactMetaMethod)
		bool status = adjustStatsForExact(me);

	// now compute
	if (me.V2==0)
		return false;
/*	if (me.V2==0) {
		printf("Warning: maf!=0 but metaV==0. Something is wrong!\n\n");
		me.pvalue = 0;
		return false;
	}*/
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
	// skip those unnecessary ones
	if (normPop&&!me.skip_pop_correct)
		return false;

	bool status = true;

	double exact_maf;
//	vector<int>::iterator pnk = me.study_N.begin();
	vector<double>::iterator pfk;

	std::vector<double> fk_tilda;
	if (normPop) {
		fk_tilda.resize(scorefile.Length(),0);
		if (me.ref_mafs.empty()) { // do not adjust
			exact_maf = 0;
			status = false;
		}
		else {
			double f = 0;
			for(int k=0;k<scorefile.Length();k++) {
				double fk = me.study_mafs[k];
				for(int i=0;i<nPop;i++)
					fk -= (double)pGamma[k][i] * me.ref_mafs[i];
				fk_tilda[k] = fk;
				f += fk/me.N*me.study_N[k];
			}
			exact_maf = f;
		}
		pfk = fk_tilda.begin();
	}
	else { // use original maf
		exact_maf = (double)me.AC / me.N/2;
		pfk = me.study_mafs.begin();
	}

	double nkdeltak = 0;
	double nkdeltakfk = 0;
	double v2 = 0;
	double new_r = 0;
	for(int i=0;i<scorefile.Length();i++) {
		nkdeltak += (double)me.study_N[i] * Ydelta[i];
		double ac = (*pfk) * (double)me.study_N[i];
		nkdeltakfk += ac * Ydelta[i];
		v2 += ac * (*pfk);
		new_r += (me.study_N[i] - 1) * Ysigma2[i] + me.study_N[i] * Ydelta[i]*Ydelta[i];
		pfk++;
	}
	new_r /= (me.N-1);
	me.U = me.U - 2*exact_maf*nkdeltak + 2*nkdeltakfk;
	me.V2 = new_r * (me.V2 + 4*v2 - 4 *me.N*exact_maf*exact_maf);

	return status;
}


bool Meta::isGroupTestRequired()
{
	bool status = false;
	if (Burden||MB||MAB||BBeta||VT||SKAT)
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
	// load Us & maf
	LoadGroupUandInfos();
	
	// load Vs
	for(int s=0;s<scorefile.Length();s++)
		AnnoGroups.LoadCovFile(covfile[s],FormatAdjust[s], marker_col, cov_col,SampleSize[s]);

	AnnoGroups.MakeFullMatrix();
	AnnoGroups.RemoveMissingData();

	if (useExactMetaMethod) {
		for(int s=0;s<scorefile.Length();s++)
			for(int g=0;g<AnnoGroups.GetGeneNumber();g++)
				for(int m1=0;m1<AnnoGroups.VariantCountInGenes(g);m1++)
					for(int m2=0;m2<AnnoGroups.VariantCountInGenes(g);m2++)
						updateExactCov(s,g,m1,m2);
	}

	// perform group test
	String method = "";
	if(Burden) {
		method = "burden";
		AnnoGroups.GroupTest(method, prefix);
	}
	if(MB) {
		method = "MB";
		AnnoGroups.GroupTest(method, prefix);
	}
	if(MAB) {
		method = "MAB";
		AnnoGroups.GroupTest(method, prefix);
	}
	if(BBeta) {
		method = "BBeta";
		AnnoGroups.GroupTest(method, prefix);
	}

	if(VT) {
		method = "VT";
		AnnoGroups.GroupTest(method, prefix);
	}
	
	if(SKAT) {
		method = "SKAT";
		AnnoGroups.GroupTest(method, prefix);
	}
}


// retrieve info from variantMap
// U stored in Vector* groupUs
// pointer to metaElement stored in vector<metaElement*> groupInfos
void Meta::LoadGroupUandInfos()
{
	// initialize
	int n_genes = AnnoGroups.GetGeneNumber();
//	groupInfos.resize(n_genes);
//	for(int g=0;g<n_genes;g++)
//		groupInfos[g].resize( AnnoGroups.VariantCountInGenes(g), NULL );

	// set pointers & groupUs
	for(int g=0;g<n_genes;g++) {
		for(int i=0;i<AnnoGroups.VariantCountInGenes(g);i++) {
			metaElement* p = getPointerMetaElement( AnnoGroups.GetChr(g,i),AnnoGroups.GetPosition(g,i),AnnoGroups.GetRef(g,i),AnnoGroups.GetAlt(g,i));
			if (p==NULL)
				continue;
//			if (p==NULL)
//				error("%s:%d:%s:%s cannot get metaElement pointer. Something is wrong!\n",annoChrs[g][i].c_str(),annoPositions[g][i],annoRefs[g][i].c_str(),annoAlts[g][i].c_str());
//			groupInfos[g][i] = p;
			double maf = p->N==0 ? 0 : (double)p->AC/p->N/2;
			if (maf==0||maf==1)
				continue;
			bool flip = false;
			double minor_maf = maf;
			if (maf>0.5) {
				minor_maf = 1 - maf;
				flip = true;
			}
			if (maf_threshold != 0.5) {
				if (minor_maf>maf_threshold)
					continue;
			}
			if (flip) {
				AnnoGroups.FlipAllele(g,i);
				AnnoGroups.SetU(g,i,-(p->U));
				AnnoGroups.SetV(g,i,-(p->V2));
			}
			else {
				AnnoGroups.SetU(g,i,p->U);
				AnnoGroups.SetV(g,i,p->V2);
			}
			AnnoGroups.SetMaf(g,i,minor_maf);
		}
	}
}

metaElement* Meta::getPointerMetaElement(String chr, int position, String ref, String alt)
{
	metaElement* p = NULL;
	if (variantMap.find(chr)==variantMap.end())
		return p;
	if (variantMap[chr].find(position)==variantMap[chr].end())
		return p;
	for(int i=0;i<variantMap[chr][position].size();i++) {
		if (variantMap[chr][position][i].ref==ref&&variantMap[chr][position][i].alt==alt) {
			p = &(variantMap[chr][position][i]);
			break;
		}
	}
	return p;
}


/**** for Meta::Prepare() ******/

// open pdf
// load list of summary files
// load list of cov files if needed
void Meta::openMetaFiles()
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
	

	SampleSize.resize( scorefile.Length(),-1 );
	GCbyStudy.Dimension(scorefile.Length(),0);

	// exact method
	flip_count =0;
	skip_count=0;
	
/*	if (relateBinary) {
		Const_binary.Dimension(scorefile.Length());
		for(int s=0;s<scorefile.Length();s++) {
			double sigma = Ysigma2g[s] + Ysigma2e[s];
		}
	}
*/
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
bool Meta::poolSingleRecord( int study, double& current_chisq, int& duplicateSNP, bool adjust, String & buffer)
{
	StringArray tokens;
	tokens.AddTokens(buffer, SEPARATORS);
	if (tokens[0].Find("#")!=-1)
		return 0;

	if (RegionStatus&&tokens[0]!=Chr)
		return false;
	int position = tokens[1].AsInteger();
	if (RegionStatus&&(position<Start||position>End))
		return false;

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

	bool is_fail = false;
	int filter_type = 1;
	if ((tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))
		is_fail = 1;
	if(tokens[8-adjust].AsDouble()<CALLRATE || tokens[9-adjust].AsDouble()<HWE) {
		is_fail = 1;
		filter_type = 0;
	}

	//STEP2: add info to variantMap if not skipped
	bool flip=false;
	int mindex = 0; // multi-allelic index
	if (variantMap.find(tokens[0])==variantMap.end())
		variantMap[tokens[0]];
	std::map<int, std::vector<metaElement> >::iterator vp = variantMap[tokens[0]].find(position);
	bool is_new_site = false;
	if (vp==variantMap[tokens[0]].end()) { // this position is not hashed before
		variantMap[tokens[0]][position].resize(1);
		vp = variantMap[tokens[0]].find(position);
		is_new_site = true;
	}
	else { // this position is hashed before
		bool is_match = false;
		for(int i=0;i<vp->second.size();i++) {
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
				if (debug)
					printf("Duplicate SNP found at %s:%s:%s from study #%d. Skipped!\n", chr_pos.c_str(),tokens[2].c_str(),tokens[3].c_str(),study+1);
				return false;
			}
		}
		else {
			mindex = vp->second.size();
			vp->second.resize(mindex+1);
			is_new_site = true;
		}
	}
	// initialize
	if (is_new_site) {
		initializeMetaElement( vp->second[mindex],scorefile.Length());
		vp->second[mindex].ref = tokens[2];
		vp->second[mindex].alt = tokens[3];
	}

	// now mark failed sites
	if (is_fail) {
		if(filter_type==0)
			fprintf(log,"Warning: variant %s:%s:%s from study %d has at least one allele missing but is polymorphic; and is excluded from meta-analysis.\n",chr_pos.c_str(),tokens[2].c_str(),tokens[3].c_str(),study);
		else if(filter_type==1)
			fprintf(log,"Warning: variant %s:%s:%s from study %d failed to pass HWE or CALLRATE filter and it is excluded from meta analysis.\n",chr_pos.c_str(),tokens[2].c_str(),tokens[3].c_str(),study);
		return false;
	}

	vp->second[mindex].study_N[study] = current_N;
	if (flip)
		vp->second[mindex].study_mafs[study] = 1-raw_af;
	else
		vp->second[mindex].study_mafs[study] = raw_af;
	// stats
	double u = tokens[13-adjust].AsDouble();
	double v = tokens[14-adjust].AsDouble();
	if (c2+c3==0) {
		vp->second[mindex].N += current_N; // ac+=0
		vp->second[mindex].directions[study] = '?';
	}
	else if (c1+c2==0) {
		vp->second[mindex].AC += current_N*2;
		vp->second[mindex].N += current_N;
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
		current_chisq = u*u/(v*v);
//		current_pvalue = pchisq(chisq,1,0,0);
	}

	// conditional analysis
	if (cond!="" && tokens[14-adjust].AsDouble()>0.0) {
		CondAna.SetCondStats(study,tokens,u,v,current_N,flip);
	}
	return true;
}

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
	double maf = me.AC==0 ? 0: (double)me.AC/me.N/2;
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
	StringArray geneLabel, chr_plot;
	Vector pos_plot;
	if (!NoPlot)
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
	fprintf(log, "%s\n",msg);
	fprintf(log, "\n");
	error(msg);
}

// update v in cov matrix in exact
void Meta::updateExactCov( int s, int g, int m1, int m2)
{
	metaElement* p1 = getPointerMetaElement( AnnoGroups.GetChr(g,m1),AnnoGroups.GetPosition(g,m1),AnnoGroups.GetRef(g,m1),AnnoGroups.GetAlt(g,m1));
	metaElement* p2 = getPointerMetaElement( AnnoGroups.GetChr(g,m2),AnnoGroups.GetPosition(g,m2),AnnoGroups.GetRef(g,m2),AnnoGroups.GetAlt(g,m2));

	if (p1==NULL) {
		String chr = AnnoGroups.GetChr(g,m1);
		String ref = AnnoGroups.GetRef(g,m1);
		String alt = AnnoGroups.GetAlt(g,m1);
		error("[Meta::updateExactCov] Cannot get information of %s:%d:%s:%s!\n",chr.c_str(),AnnoGroups.GetPosition(g,m1),ref.c_str(),alt.c_str());
	}
	if (p2==NULL) {
		String chr = AnnoGroups.GetChr(g,m2);
		String ref = AnnoGroups.GetRef(g,m2);
		String alt = AnnoGroups.GetAlt(g,m2);
		error("[Meta::updateExactCov] Cannot get information of %s:%d:%s:%s!\n",chr.c_str(),AnnoGroups.GetPosition(g,m2),ref.c_str(),alt.c_str());
	}

	int nk = p1->study_N[s];
	double nkfk2 = nk * p1->study_mafs[s] * p2->study_mafs[s];
	double new_r = p1->N;
	double maf1 = (double)p1->N / p1->AC;
	double maf2 = (double)p2->N / p2->AC;
	double ff = maf1*maf2;
	double old_cov = AnnoGroups.GetOneCov(g,m1,m2);
	double new_val = new_r*(old_cov*Ysigma2[s]-4*nk*ff + 4*nkfk2);
	AnnoGroups.UpdateCovValue(new_val,g,m1,m2);
}

