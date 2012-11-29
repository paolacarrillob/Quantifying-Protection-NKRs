/*
 * Genes.cpp
 *
 *  Created on: Apr 20, 2011
 *      Author: paola
 */

#include "Genes.h"


/* FUNCTIONS OF CLASS GenePool*/
/* constructs a pool with HLA C Alleles (14) */
void MHCGenePool::FillMHCGenePool(int size)
{
	poolSize = size;
	for(int i=0; i<poolSize; i++)
	{
		MHCGene dummy;
		//int geneID = dummy.Bit2Int();
		int geneID = dummy.GetGeneID();
		if(!GenePool::GeneAlreadyInPool(geneID))
			GenePool::genes.push_back(geneID);
	}
}

KIRGenePool :: KIRGenePool(GenePool& mhcPool, bool dist, int specificity)
{
	genePoolFile.open("GenePool.data", ios::out);
	genePoolFile << "# MHC \t KIR \n";
	for(int i=0; i<10; i++)
	{
		for(unsigned int j=0; j<mhcPool.GetPoolSize(); j++)
		{
			Bitstring mhc_molecule(mhcPool.GetGenes().at(j));
			//Bitstring mhc_molecule(mhcPool.RandomlyPickGene(dist));
			BitstringKIR kir_molecule(mhc_molecule,specificity);
			int kir_id = kir_molecule.Bit2Int();
			if(!GenePool::GeneAlreadyInPool(kir_id))
			{
				GenePool::genes.push_back(kir_id);
				genePoolFile << mhc_molecule.Bit2Int() << "\t" << kir_id <<endl;
				//cout << mhc_molecule.Bit2Int() << "\t" ;
				mhc_molecule.PrintBits();
				//cout << "\t"<< kir_id << "\t";
				kir_molecule.PrintBits();
				//cout<< "\t"<<endl;
			}
		}

	}
	genePoolFile.close();
}


void Map::FillMap(GenePool& mhcpool, KIRGene& kirGene)
{
	int M_id = 0;
	//calculate the value of M_id to determine whether the gene is pseudogene or not
	for(unsigned int i = 0; i < mhcpool.GetPoolSize(); i++)
	{
		Gene mhcGene;
		mhcGene.SetGeneID(mhcpool.GetGenes().at(i));
		int L = kirGene.BindMolecule(mhcGene);
		if(L >= kirGene.GetGeneSpecificity())
			M_id += (1<<i);
	}

	//set the Gene pseudo
	kirGene.SetPseudogene(M_id);

	//save the infos in pairs into a map
	if(!IsGeneInMap(kirGene))
	{
		pair <int, int> genePhenotype;
		genePhenotype = make_pair(kirGene.GetGeneSpecificity(), M_id);
		mapGenes.insert(make_pair(kirGene.GetGeneID(),genePhenotype)); //!!! i don't completely
																// understand why I need to make_pair of (int and pair)!
	}

	/*
	 * for accessing the entries of the map!
	map< int, pair <int, int> > ::iterator it;
	for(it = mapGenes.begin(); it != mapGenes.end(); it ++)
	 {
		cout << (*it).first << "|" << (*it).second.first <<"|"<< (*it).second.second <<endl;

	 }*/
}

/*This function determines if the map contains a particular gene.*/
bool Map ::IsGeneInMap(Gene& gene)
{
	int geneNumber = gene.GetGeneID();
	map<int, pair <int, int> > :: iterator it = mapGenes.find(geneNumber);
	if(it != mapGenes.end()) //if it finds the gene
		return true;
	else
		return false;

}

//restoring the map from the file...
string Map :: RestoreMap(stringstream& smline)
{
	Gene dummy;

	int gene_id;
	int m_id;
	int L;
	smline >> gene_id;
	smline >> L;
	smline >> m_id;
	pair <int, int> genePhenotype;
	dummy.SetGeneID(gene_id);
	genePhenotype = make_pair(L, m_id);
	if(!IsGeneInMap(dummy))
		mapGenes.insert(make_pair(gene_id,genePhenotype));
	//cout << gene_id << "|" << L << "|" << m_id << "\n";
	string mstring = smline.str();
	return mstring;
}

//saving the Map into the backup file
void Map :: SaveBackupMap(fstream&  backupFile)
{
	//cout << "saving Map..."<<endl;
	map< int, pair <int, int> > ::iterator it;
	backupFile << mapGenes.size() <<"\t";
	for(it = mapGenes.begin(); it != mapGenes.end(); it ++)
	 {
		backupFile <<(*it).first << "\t" << (*it).second.first <<"\t"<< (*it).second.second <<"\t";

	 }
}
/* This function assures that every single gene in the pool is unique*/
bool GenePool:: GeneAlreadyInPool(int geneID)
{
	for(unsigned int i=0; i<genes.size(); i++)
	{
		int number = genes.at(i);
		if(number == geneID)
			return true;
	}
	return false;
}


/* This function sets each MHC allele with a predefined frequency. (adopted from dbMHC Project)*/
float GenePool :: GetAlleleFreq(int alleleIndex)
{
	static float allele_freq=0;
	allele_freq = 0.2*exp(-0.22*alleleIndex); // probability distribution for hla-c alleles in european population
	return allele_freq;
}


/*This function returns an allele "randomly". i.e. either according to the HLA -C distribution or with equal probability*/
int GenePool :: RandomlyPickGene(bool distribution)
{
	int j;
	if(distribution==true)
	{
		//double size = genes.size();
		static double probabilities[50];
		for (unsigned int i=0; i<genes.size(); i++)
		{
			probabilities[i] = GetAlleleFreq(i);
			//cout << probabilities[i] << " ";

		}
		j = RandomNumber(probabilities, genes.size()-1);
		//cout << endl;
	}
	else
	{
		j = RandomNumber(0,genes.size()-1);
	}
	return genes.at(j);
}

/* FUNCTIONS OF CLASS Gene*/
/*constructs default genes*/
Gene::Gene ()
{
	functional = true;
	isExpressed = true;
	geneID = 0;
	specificity = 0;
	pseudoGene = false;
	mID = 0;

}

/*This function checks whether two genes are equal*/
bool Gene :: operator == (Gene& rhs)
{
	if((geneID == rhs.geneID)&& (specificity == rhs.specificity))
		return true;
	else
		return false;
}

void Gene :: SetGeneFunctionality(bool functionality)
{
	functional = functionality;
}

void Gene:: SetGeneExpression(bool expression)
{
	isExpressed = expression;
}

void Gene:: SetGeneSpecificity()
{
	specificity = RandomNumber(1,16);

}

void Gene:: SetPseudogene(bool expression)
{
	pseudoGene = expression;

}

void Gene :: SetPseudogene(int M_id)
{

	if(M_id > 0)
	{
		pseudoGene = false;
		mID = M_id;
	}
	else
		pseudoGene = true;
}


/*This function performs the different mutational operators for a Gene*/
void Gene :: Mutate(double mutationRate)
{
	Bitstring molecule(geneID);
	SetGeneSpecificity();
	//molecule.PointMutation(mutationRate);
}

/*This function flips one random bit of the bit string*/
void Gene :: PointMutation()
{
	int bit = RandomNumber(1,MOLECULE_LENGTH-1);
	geneID ^= (1 << bit); //XOR function will flip the bit
}

void Gene::PrintBits()
{
	for(int l = 1; l < pow(2.0,16); l = l*2)
	{
		if(geneID & l)
			cout << "1";
		else
			cout << "0";
	}
	cout <<endl;
}

/*This function mutates only the specificity... it's a more gradual change*/
void Gene ::MutateSpecificity()
{
	if(specificity > 1)
	{
		if(RandomNumberDouble() < 0.5)
			specificity += 1;
		else
			specificity -=1;
	}
	else
		specificity +=1;
}

/*This function counts unique genes within a given pool*/
bool Gene :: IsGeneUnique(vector<Gene>& genePool, int counter)
{
	//cout <<"this is the counter : "<<counter << "\n";
	for (int i = 0; i<counter; i++)
	{
		Gene currentGene = genePool.at(i);
		if(geneID == currentGene.geneID)
			return false;
	}
	return true;
}

int Gene::BindMolecule(Gene& anotherMolecule)
{

	/* old function, when I was using the Bitstring class... now it's not necessary anymore!
	 * Bitstring anotherBitstring(anotherMolecule.geneID);
	Bitstring currentBitstring(geneID);

	int bindingStrength = currentBitstring.AdjacentMatch(anotherBitstring);
	return bindingStrength;*/

	int counter = 0;
	int max = 0;
	for( int l = 1; l < pow(2.0,16); l = l*2)
	{
		if((geneID & l) ^ (anotherMolecule.geneID & l)) //check with XOR operator whether the bits are  complementary
		{
			counter ++;
			if(counter > max)
				max = counter;
		}
		else
			counter = 0;
	}
	if(counter > max)
		max = counter;
	return max;
}

void Gene::SaveBackupGenes(fstream& backupFile)
{
	backupFile << geneID << "\t" << functional << "\t"<<isExpressed <<"\t"<< specificity << "\t" << pseudoGene  << "\t"<< mID <<"\t";
	//cout << "saving genes..."<<endl;
}

void Gene::SaveGenes(fstream& outfile)
{
	outfile << geneID << "_" << specificity <<"_";
	if(pseudoGene)
		outfile << "x\t";
	else
		outfile <<mID<<"\t";

}

Gene& Gene::Copy(Gene& rhsGene)
{
	this->functional = rhsGene.functional;
	this->geneID = rhsGene.geneID;
	this->isExpressed = rhsGene.isExpressed;
	this->specificity = rhsGene.specificity;
	//this->genePhenotype = rhsGene.genePhenotype;
	this->pseudoGene = rhsGene.pseudoGene;
	this->mID = rhsGene.mID;
	return *this;
}

string Gene::RestoreGenes(stringstream& sgline)
{
	sgline >> geneID;
	sgline >> functional;
	sgline >> isExpressed;
	sgline >> specificity;
	sgline >> pseudoGene;
	sgline >> mID;
	//cout <<geneID <<"\t" <<functional << "\t";
	string gstring = sgline.str();
	return gstring;
}

MHCGene::MHCGene() //constructs an MHC gene with L = 0
{
	functional = true;
	isExpressed = true;
	geneID = 0;
	specificity = 0;
	pseudoGene = false;
	mID = 0;
	int random_bit = 0;
	for(int i = 0; i< MOLECULE_LENGTH; i++)
	{
		random_bit = RandomNumberDouble()<0.5;
		geneID += random_bit * (1<<i);
		//cout << random_bit;
	}
	//genePhenotype.first = geneID;
	//genePhenotype.second = specificity;
	//cout << genePhenotype.first << " " <<genePhenotype.second <<endl;
}


KIRGene::KIRGene(int L) //constructs a KIR gene with a random L from 1-16
{
	functional = true;
	isExpressed = true;
	geneID = 0;
	specificity = RandomNumber(L-1,L+1);
	pseudoGene = false;
	mID = 0;
	if(RandomNumberDouble()<0.5)
		geneType = inhibitory;
	else
		geneType = activating;

	int random_bit = 0;
	for(int i = 0; i< MOLECULE_LENGTH; i++)
	{
		random_bit = RandomNumberDouble()<0.5;
		geneID += random_bit * (1<<i);
		//cout << random_bit;
	}
	//genePhenotype.first = geneID;
	//genePhenotype.second = specificity;
	//cout << genePhenotype.first << " " <<genePhenotype.second <<endl;
}

void KIRGene ::SetGeneType(int type)
{
	switch(type)
	{
	case 0:{geneType = inhibitory;}break;
	case 1:{geneType = activating;}break;
	}
}


void KIRGene::MutateReceptorType()
{
	if(geneType == inhibitory)
	{
		geneType = activating;
		return;
	}
	if(geneType == activating)
		geneType = inhibitory;
}

bool KIRGene :: operator == (KIRGene& rhs)
{
	if(Gene::operator ==(rhs) && geneType == rhs.geneType)
		return true;
	else
		return false;
}

bool KIRGene :: IsInhibitory()
{
	if(geneType == inhibitory)
		return true;
	else
		return false;
}

bool KIRGene :: IsActivating()
{
	if(geneType == activating)
		return true;
	else
		return false;
}

KIRGene& KIRGene::Copy(KIRGene &rhsGene)
{
	Gene::Copy(rhsGene);
	this->geneType = rhsGene.geneType;
	//cout << this->geneType<< endl;
	return *this;
}

string KIRGene::RestoreGenes(stringstream& sgline)
{
	Gene::RestoreGenes(sgline);
	int gene_type;
	sgline >> gene_type;
	switch(gene_type)
	{
	case 0:{geneType = inhibitory;}break;
	case 1:{geneType = activating;}break;
	}
	//cout <<geneID <<"\t" <<functional << "\t";
	string gstring = sgline.str();
	return gstring;
}

void KIRGene::SaveGenes(fstream &outfile)
{
	if(geneType == activating)
	{
		outfile << "a_";
	}
	if(geneType == inhibitory)
	{
		outfile << "i_";
	}
	Gene::SaveGenes(outfile);

}

void KIRGene::SaveBackupGenes(fstream &backupFile)
{
	Gene::SaveBackupGenes(backupFile);
	backupFile << geneType << "\t";
	//cout <<"saving kir genes......" <<endl;
}


