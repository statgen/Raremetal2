#include "af1KG.h"
#include "Error.h"
#include "MetaUtility.h"
#include "MathSVD.h"
#include "InputFile.h"

// return number of populations
int af1KG::GetPopCount()
{
	if (popList.Length()<0)
		error("Cannot get population count before loading popList!\n\n");
	return popList.Length();
}

// read pop info
void af1KG::loadPopList( String& pop_list_name )
{
	int index = 0;
	IFILE f_pop_list = ifopen(pop_list_name,"r");
	if (f_pop_list==NULL)
		error("Please provide a list of populations you want to include for correction in exact method using option --popList!\n");
	while(!(ifeof(f_pop_list))) {
		String buffer;
		buffer.ReadLine(f_pop_list);
		popList.Push(buffer);
		popMap[buffer] = index;
		index++;
	}
	ifclose(f_pop_list);
	if (popList.Length()==0)
		error("--popList file is empty!\n");
	nPop = popList.Length();
	printf("Loaded --popList. %d populations will be used for correction:\n",nPop);
	for(int i=0;i<popList.Length();i++)
		printf(" %s",popList[i].c_str());
	printf("\n\n");
}

void af1KG::loadPosIndex(String& filename, String& Chr, int Start, int End)
{
	bool region_status = false;
	if (Start!=-1 && End!=-1)
		region_status = true;
	IFILE file;
	file = ifopen(filename,"r");
	if (file==NULL)
		error("[af1KG::loadPosIndex] Cannot open file:%s\n",filename.c_str());
	while(!ifeof(file)) {
		String buffer;
		buffer.ReadLine(file);
		if (buffer.FindChar('#')!=-1)
			continue;
		StringArray tokens;
		tokens.AddTokens(buffer, '\t');
		int position = tokens[1].AsInteger();
		if (region_status) {
			if (tokens[0]!=Chr)
				continue;
			if (position<Start||position>End)
				continue;
		}
		posIndex[position] = true;
	}
	ifclose(file);
}


// load from pop vcf instead of plain text file
void af1KG::loadPopVcf( String& pop_vcf_name, bool NoPreloading)
{
	int n_outrage = 0;
	int max_outrage = 50;
	printf("Loading population vcf. This will take a while...\n");
	IFILE pop_vcf = ifopen(pop_vcf_name,"r");
	if (pop_vcf==NULL)
		error("Cannot open --popFile\n");
	bool header_line = true;
	int partial_added = 0;
	int lc = 0;
	int lc_interval = 1000000;
	while(!(ifeof(pop_vcf))) {
		String buffer;
		buffer.ReadLine(pop_vcf);
		if (header_line) {
			if (buffer.Find("#")!=-1)
				continue;
			else
				header_line = false;
		}
		lc++;
		int a = lc/lc_interval;
		if (a*lc_interval==lc)
			printf("Processed %d million lines...\n",a);
		StringArray tokens;
		tokens.AddTokens(buffer, '\t');
		int position = tokens[1].AsInteger();
		if (!NoPreloading) {
			if (posIndex.find(position)==posIndex.end())
				continue;
		}
		String chr = tokens[0];
		removeChrFromString(chr);
		String ref = tokens[3];
		StringArray alts;
		alts.AddTokens(tokens[4],',');
		// now check for alleles and allocate new element
		if(refMAFs.find(chr)==refMAFs.end())
			refMAFs[chr];
		if(refMAFs[chr].find(position)!=refMAFs[chr].end()) {
			for(int a=0;a<alts.Length();a++) {
				bool match = false;
				for(int a2=0;a2<refMAFs[chr][position].size();a2++) {
					if(refMAFs[chr][position][a2].ref==ref && refMAFs[chr][position][a2].alt==alts[a]) {
						warning("Skip duplicate record at vcf:\n %s\n",buffer.c_str());
						match = true;
						break;
					}
				}
				if (!match) {
					int old_n = refMAFs[chr][position].size();
					refMAFs[chr][position].resize(old_n+1);
					refMAFs[chr][position][old_n].ref = ref;
					refMAFs[chr][position][old_n].alt = alts[a];
					refMAFs[chr][position][old_n].mafs.resize(nPop,0);
				}
			}
		}
		else { // new
			refMAFs[chr][position].resize(alts.Length());
			for(int a=0;a<alts.Length();a++) {
				refMAFs[chr][position][a].ref = ref;
				refMAFs[chr][position][a].alt = alts[a];
				refMAFs[chr][position][a].mafs.resize(nPop,0);
			}
		}
		// now add mafs
		StringArray info_tokens;
		info_tokens.AddTokens(tokens[7],';');
		int added_pops = 0;
		for(int i=0;i<info_tokens.Length();i++) {
			StringArray etokens;
			etokens.AddTokens(info_tokens[i],'=');
			if (etokens.Length()!=2)
				continue;
			if (popMap.find(etokens[0])==popMap.end())
				continue;
			int idx = popMap[etokens[0]];
			StringArray maf_fields;
			maf_fields.AddTokens(etokens[1],',');
			if(maf_fields.Length() != alts.Length())
				error("#ALT and #MAFs do not match at line:\n %s\n",buffer.c_str());
			for(int a=0;a<alts.Length();a++) {
				double f = maf_fields[a].AsDouble();
				refMAFs[chr][position][a].mafs[idx] = f;
			}
			added_pops++;
		}
		if (added_pops<nPop) {
			partial_added++;
			printf("One or more pop maf is missing at line:\n%s\n",buffer.c_str());
		}
		if (partial_added>max_outrage)
			error(">%d lines has one or more missing pop MAF. Please check your if the vcf match the pop list file!\n",max_outrage);
	}
	ifclose(pop_vcf);
	printf("   done.\n\n");
}


bool af1KG::setVariantPopMAF( metaElement& me, String& chr, int position )
{
	if (refMAFs.find(chr)==refMAFs.end())
		return false;
	if (refMAFs[chr].find(position)==refMAFs[chr].end())
		return false;
	bool match = false;
	bool flip = false;
	int index;
	for(int i=0;i<refMAFs[chr][position].size();i++) {
		if (refMAFs[chr][position][i].ref==me.ref && refMAFs[chr][position][i].alt==me.alt) {
			match = true;
		}
		else if (refMAFs[chr][position][i].alt==me.ref && refMAFs[chr][position][i].ref==me.alt) {
			match = true;
			flip = true;
		}
		if (match) {
			index = i;
			break;
		}
	}
	if (!match)
		return false;

	me.ref_mafs.resize(nPop);
	for(int i=0;i<nPop;i++) {
		if (flip)
			me.ref_mafs[i] = 1 - refMAFs[chr][position][index].mafs[i];
		else
			me.ref_mafs[i] = refMAFs[chr][position][index].mafs[i];
	}
	return true;
}


// simple linear regression
void af1KG::setLinearRegressionCoefficients( Vector& gamma,Matrix& X, Vector& Y)
{
	if (X.rows != Y.Length())
		error("[getLinearRegressionCoefficients] X.row=%d, Y length=%d. Something is wrong...\n",X.rows,Y.Length());
	if (Y.Length()==0)
		error("No available non-zero maf in study. Please check your input file!\n");
	Matrix transX;
	transX.Transpose(X);
	Matrix inv;
	inv.Product(transX,X);
	SVD svd;
	svd.InvertInPlace(inv);
	// X'Y
	Vector XY;
	XY.Dimension(nPop);
	for(int i=0;i<nPop; i++)
		XY[i] = Y.InnerProduct(transX[i]);
//	printf("study=%d",study);
	gamma.Dimension(nPop);
	for(int i=0;i<nPop; i++)
		gamma[i] = inv[i].InnerProduct(XY);	

}

/* do regression, set coefficient
// fit linear model with beta >=0
void af1KG::setGammaFromRegression( Matrix& pgamma, int study, Matrix & X, Vector & Y )
{
	if (X.rows != Y.Length())
		error("[Meta::setGammaFromRegression] X.row=%d, Y length=%d. Something is wrong...\n",X.rows,Y.Length());
	if (Y.Length()==0)
		error("No available non-zero maf in study #%d. Please check your input file!\n",study+1);

	// (X'X)^-1
	Matrix transX;
	transX.Transpose(X);
	Matrix inv;
	inv.Product(transX,X);
	SVD svd;
	svd.InvertInPlace(inv);
	// X'Y
	Vector XY;
	XY.Dimension(nPop);
	for(int i=0;i<nPop; i++)
		XY[i] = Y.InnerProduct(transX[i]);
//	printf("study=%d",study);
	for(int i=0;i<nPop; i++) {
		pgamma[study][i] = inv[i].InnerProduct(XY);	
//printf(",%g",pgamma[study][i]);
	}	
//printf("\n\n");
}*/

