/*
 * Host.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: paola
 */

#include "Host.h"
/* FUNTIONS OF CLASS VIRUS
 * 	Constructs the virus*/
Virus::Virus()
{
	SetViralParameters(0.0, 0.0, 0.0, 0.0, 0, false);
}

void Virus::SetViralParameters(double _downregulation, double _decoy, double _viralLoad, double _lifeTimeVirus, int virus_type, bool only_acute)
{
	mutationRateDownregulation = _downregulation;
	mutationRateDecoy = _decoy;
	viralLoad = _viralLoad;
	originalViralLoad = _viralLoad;
	lifeTimeVirus = _lifeTimeVirus;
	switch(virus_type)
	{
	case 0:{virusType = wildType;};break;
	case 1:{virusType = downregulating;}; break;
	case 2:{virusType = decoy;};break;
	}
	mhcDecoy.SetGeneID(0);
	onlyAcute = only_acute;
}

/*This function sets the virus to downregulate the MHC expression*/
void Virus ::DownregulateMHC()
{
	if(virusType == wildType)
	{
		if(RandomNumberDouble()<mutationRateDownregulation)
		{
			virusType = downregulating;
			return;
		}
	}
}

/*This function sets the virus to steal a decoy from the host*/
void Virus::BuildDecoy(int mhcID)
{
	//if(virusType == downregulating || virusType == decoy)
	//{
		if(RandomNumberDouble() <mutationRateDecoy)
		{
			mhcDecoy.SetGeneID(mhcID);
			virusType = decoy;
		}
	//}
}

/*This function performs a deep copy*/
Virus & Virus::Copy(Virus& rhsVirus)
{
	// checking if it is a self-assignment
	if(this == &rhsVirus)
		return *this;
	//copying member variables
	this->mutationRateDownregulation = rhsVirus.mutationRateDownregulation;
	this->mutationRateDecoy = rhsVirus.mutationRateDecoy;
	this->lifeTimeVirus = rhsVirus.lifeTimeVirus;
	this->viralLoad = rhsVirus.viralLoad;
	this->originalViralLoad = rhsVirus.originalViralLoad;
	this->virusType = rhsVirus.virusType;
	//this->mhcDecoy.SetGeneID(rhsVirus.mhcDecoy.GetGeneID());
	this->mhcDecoy.Copy(rhsVirus.mhcDecoy);
	this->onlyAcute = rhsVirus.onlyAcute;
	return *this;
}

bool Virus :: operator == (Virus& rhs)
{
	if(this->virusType == rhs.virusType)
		return true;
	else
		return false;
}

bool Virus :: IsDownregulatingMHC()
{
	if(virusType == downregulating)
		return true;
	else
		return false;
}

bool Virus :: IsStealingDecoy()
{
	if(virusType == decoy)
		return true;
	else
		return false;
}

bool Virus :: IsWildType()
{
	if(virusType == wildType)
		return true;
	else
		return false;
}

void Virus::SaveBackupVirus(fstream& file)
{
	file << mutationRateDownregulation << "\t" << mutationRateDecoy << "\t" << lifeTimeVirus << "\t" << viralLoad <<"\t" <<originalViralLoad<<"\t"<< virusType <<"\t" << mhcDecoy.GetGeneID()<< "\t"<<onlyAcute <<"\t";
}

void Virus::PrintParametersVirus()
{
	cout << mutationRateDownregulation << "\t" << mutationRateDecoy << "\t" << lifeTimeVirus << "\t" << viralLoad <<"\t" <<originalViralLoad<<"\t"<< virusType <<"\t" << mhcDecoy.GetGeneID()<< "\t"<<onlyAcute <<endl;
}

string Virus::RestoreVirus(stringstream& svline)
{
	int type;
	int id;
	svline >> mutationRateDownregulation;
	svline >> mutationRateDecoy;
	svline >> lifeTimeVirus;
	svline >> viralLoad;
	svline >> originalViralLoad;
	svline >> type;
	switch(type)
	{
	case 0:{virusType = wildType;};break;
	case 1:{virusType = downregulating;}; break;
	case 2:{virusType = decoy;};break;
	}
	svline >> id;
	mhcDecoy.SetGeneID(id);
	svline >> onlyAcute;
	string vstring = svline.str();
	return vstring;
}

void Virus :: SaveParametersVirus(fstream& outfile)
{
	outfile << virusType << "\t" << viralLoad <<"\t" <<originalViralLoad << "\t"<< mhcDecoy.GetGeneID() <<"\t"<<onlyAcute<<"\t";
}


/*FUNCTIONS OF CLASS HOST
 *Constructs a host for the initialization of the population: it fills the MHC genes with a randomly picked allele from the population
 * and creates KIRs that match their own MHC according to the specificity*/
Host::Host(int loci_kir, int loci_mhc, double _mutationRate, bool _tuning, int numberOfExtraKirs,Map& kirMap, MHCGenePool& mhcPool, bool hla, int gene_type) {

	InitializeHostParameters(_mutationRate,_tuning, loci_kir, loci_mhc);
	//fill the mhc Genes
	for(int i = 0; i <LOCI_MHC*TWO; i ++)
	{
		Gene firstGene;
		int mhc1 = mhcPool.RandomlyPickGene(hla);
		firstGene.SetGeneID(mhc1);
		mhcGenes.push_back(firstGene);
	}

/*
	Gene secondGene;
	secondGene.SetGeneID(mhc2);
	mhcGenes.push_back(secondGene);*/

	//create random KIRs with random specificity for ONE haplotype
	// (at the beginning of the simulation, the size of the map should equal the LOCI_NUMBER)
	map< int, pair <int, int> > ::iterator it;
	for(it = kirMap.GetMap().begin(); it != kirMap.GetMap().end(); it ++)
	{
		int id = it->first;
		int L = it->second.first;
		int pseudo = it->second.second;
		KIRGene kir;
		kir.SetGeneSpecificity(L);
		kir.SetGeneID(id);
		kir.SetPseudogene(pseudo);
		kir.SetGeneFunctionality(false);
		kir.SetGeneExpression(false);
		kir.SetGeneType(gene_type); // 0 for inhibitory receptor; 1 for activating receptors!
		kirGenes.push_back(kir);

	}
	//copy the 5 first kirs into the other haplotype
	for(int i=0; i<LOCI_KIR; i++)
	{
		KIRGene kir = kirGenes.at(i);
		kirGenes.push_back(kir);

	}

	if(tuning == true)
		EducateKIRs();
	ExpressKIRs(numberOfExtraKirs);
	CountFunctionalKIRs();
	age = RandomNumber(1,70); //population initialized with a random age between 1 and 70
}

/*Constructs a baby host out of two parents*/
Host::Host(int loci_kir, int loci_mhc, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, vector<KIRGene>& kirGenesMother, vector<KIRGene>& kirGenesFather,double _mutationRate, bool _tuning, int numberOfExtraKirs, Map& kirMap, int mutationType, int gene_type)
//Host::Host(int loci_kir, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, GenePool& kirPool, vector<Gene>& kirGenesMother, vector<Gene>& kirGenesFather, int specificity, double _mutationRate, bool _tuning, int numberOfExtraKirs)
{	//to create a NEW host: the haplotypes of KIR of BOTH parents are needed. Besides one MHC haplotype of one parent plus one of the pool

	InitializeHostParameters(_mutationRate, _tuning, loci_kir, loci_mhc);
	 //cout<< "do I get stuck here?? | Host constructor!"<<endl;

	int KIR_init_mum = 0;
	int KIR_end_mum = 0;
	int MHC_init_mum = 0;
	int MHC_end_mum = 0;

	int KIR_init_dad = 0;
	int KIR_end_dad = 0;

	//pick haplotype 1/0 of each parent for the KIRs!
	int hap_mum=(RandomNumberDouble()<0.5);
	int hap_dad= (RandomNumberDouble()<0.5);
	//int k= 0;
	if(hap_mum)
	{
		KIR_init_mum = LOCI_KIR;
		KIR_end_mum = LOCI_KIR*TWO;
	}
	else
	{
		KIR_init_mum = 0;
 		KIR_end_mum = LOCI_KIR;
	}

	if(hap_dad)
	{
		KIR_init_dad = LOCI_KIR;
		KIR_end_dad = LOCI_KIR*TWO;
	}

	else
	{
		KIR_init_dad = 0;
		KIR_end_dad = LOCI_KIR;
	}

	//pick haplotype 1/0 of the parent for the MHCs!
	int hap_mhc = (RandomNumberDouble()<0.5);
	if(hap_mhc)
	{
		MHC_init_mum = LOCI_MHC;
		MHC_end_mum = LOCI_MHC*TWO;
	}

	else
	{
		MHC_init_mum = 0;
		MHC_end_mum = LOCI_MHC;
	}

	//copy the KIR haplotype into the new host (mutation occurs!)
	int hap_child = (RandomNumberDouble()<0.5);
	if(hap_child)
	{
		for(int i=KIR_init_dad; i<KIR_end_dad; i++)
		{
			KIRGene kir_hap2;
			kir_hap2.Copy(kirGenesFather.at(i));
			//MUTATION!
			if(RandomNumberDouble() < mutationRateHost)
			{
				MutateGenes(mutationType, kir_hap2, kirMap, mhcPool, gene_type);
			}
			kirGenes.push_back(kir_hap2);
		}
		for(int j=KIR_init_mum; j<KIR_end_mum; j++)
		{
			KIRGene kir_hap1;
			kir_hap1.Copy(kirGenesMother.at(j));
			if(RandomNumberDouble() < mutationRateHost)
			{
				MutateGenes(mutationType, kir_hap1, kirMap, mhcPool, gene_type);
			}
			kirGenes.push_back(kir_hap1);
		}

		//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
		{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			//cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
		}

		//generate new mhc genes (from the pool) for the new born
		for(int mm = 0; mm<LOCI_MHC; mm++)
		{
			Gene mhc2;
			int mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
			mhc2.SetGeneID(mhcFromTheGenePool);
			mhcGenes.push_back(mhc2);
		}
	}
	else
	{
		for(int i=KIR_init_mum; i<KIR_end_mum; i++)
		{
			KIRGene kir_hap1;
			kir_hap1.Copy(kirGenesMother.at(i));
			if(RandomNumberDouble() < mutationRateHost)
			{
				MutateGenes(mutationType, kir_hap1, kirMap, mhcPool, gene_type);
			}
			kirGenes.push_back(kir_hap1);
		}
		for(int j=KIR_init_dad; j<KIR_end_dad; j++)
		{
			KIRGene kir_hap2;
			kir_hap2.Copy(kirGenesFather.at(j));
			if(RandomNumberDouble() < mutationRateHost)
			{
				MutateGenes(mutationType, kir_hap2, kirMap, mhcPool, gene_type);
			}
			kirGenes.push_back(kir_hap2);
		}
		//copy the MHC haplotype into the new host
		//generate new mhc genes for the new born

		//generate new mhc genes (from the pool) for the new born
		for(int mm = 0; mm<LOCI_MHC; mm++)
		{
			Gene mhc2;
			int mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
			mhc2.SetGeneID(mhcFromTheGenePool);
			mhcGenes.push_back(mhc2);
		}
		//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
		{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			//cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
		}

	}

	if(tuning == true)
	{
		EducateKIRs();
	}
	ExpressKIRs(numberOfExtraKirs);
	CountFunctionalKIRs();
	age = 1.0; //newborns are given the age of 15
}

void Host::InitializeHostParameters(double mutationRate, bool _tuning, int loci_kir, int loci_mhc)
{
	dead = false;
	mutationRateHost = mutationRate;
	tuning = _tuning;
	mainInfectionType = susceptible;
	LOCI_KIR = loci_kir;
	LOCI_MHC = loci_mhc;
	viralDeathRate = 0.0;
	age = 0.0;
	inhibitoryKIRs = 0;
	activatingKIRs = 0;
	//MHCsearchinPool = 0;

}

void Host :: MutateGenes(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool, int gene_type)
{
	if(mutationType == 1)//pick another molecules as mutation
	{
		KIRGene newGene(RandomNumber(2,16));
		if(!kirMap.IsGeneInMap(newGene))
		{
			kirMap.FillMap(mhcPool, newGene);
			kir_hap2.SetGeneType(gene_type);////to test just what happens if i ONLY have activating receptors!
			kir_hap2.Copy(newGene);
		}

		/*if(RandomNumberDouble()<0.5) //change it to switch receptor types
		{
			kir_hap2.MutateReceptorType();
		}*/
	}

	if(mutationType == 2)//point mutation + L
	{
		//cout<< "original: ";
		//cout << kir_hap2.GetGeneSpecificity() << "|" << kir_hap2.GetGeneMid()<< "|";
		//kir_hap2.PrintBits();
		if(RandomNumberDouble()<0.8)
		{
			kir_hap2.PointMutation();
			kirMap.FillMap(mhcPool, kir_hap2);
			//cout<< "point:    ";
			//cout << kir_hap2.GetGeneSpecificity() << "|" << kir_hap2.GetGeneMid()<< "|";
			//kir_hap2.PrintBits();
		}

		if(RandomNumberDouble()<0.2)
		{
			kir_hap2.MutateSpecificity();
			kirMap.FillMap(mhcPool, kir_hap2);
			/*if(RandomNumberDouble()<0.5)
			{
				kir_hap2.MutateReceptorType();
			}*/
			//cout<< "spec:     ";
			//cout << kir_hap2.GetGeneSpecificity() << "|" << kir_hap2.GetGeneMid()<< "|";
			//kir_hap2.PrintBits();

		}
		/*if(RandomNumberDouble()<0.2)
		{
			kir_hap2.MutateReceptorType();
		}*/
	}
}


/*This functions tunes the KIR repertoire according to the self MHC repertoire and whether they are inhibiting or activating*/
void Host :: EducateKIRs()
{
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;

	for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
	{
		int education_signal_activating = 0; // to keep track of the activating signal per receptor
		for(mhcIt = mhcGenes.begin(); mhcIt !=mhcGenes.end(); mhcIt ++)
		{
			int bindingStrength = kirIt->BindMolecule(*mhcIt);
			//mhcIt->PrintBits();
			//kirIt->PrintBits();

			if(kirIt->IsInhibitory()) // if KIRs are inhibitory
			{
				if(bindingStrength>=kirIt->GetGeneSpecificity()) //AND it binds to the MHC
				{
					//cout << "it binds!"<<endl;
					kirIt->SetGeneFunctionality(true); //it is licensed!
					kirIt->SetGeneExpression(true);
					//cout<< "inhibitory |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
				}

				else // but if it doesn't bind
				{
					//cout << "it doesn't bind!"<<endl;
					kirIt->SetGeneFunctionality(false); // it should not be licensed
					kirIt->SetGeneExpression(false);
					//cout<< "inhibitory |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
				}
			}

			if(kirIt->IsActivating()) // if KIRs are activating
			{
				if(bindingStrength>=kirIt->GetGeneSpecificity()) // AND it binds to the MHC
				{
					//cout << "it binds!"<<endl;
					kirIt->SetGeneFunctionality(false); //this gene is unlicensed
					kirIt->SetGeneExpression(false);
					//cout<< "activating |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
				}
				else //BUT if it deosn't bind
				{
					//cout << "it doesn't bind!"<<endl;
					education_signal_activating ++; //keep track of how many MHC it doesn't bind
					//cout << education_signal_activating <<endl;
					//cout<< "activating |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
				}
			}
		}
		if(education_signal_activating == mhcGenes.size()) //if the activating KIR doesn't recognize ANY of the MHC
		{
			//cout << education_signal_activating << "|" << mhcGenes.size() << ": ";
			kirIt->SetGeneFunctionality(true); //the gene is licensed
			kirIt->SetGeneExpression(true);
			//cout<< "activating |"<<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
		}
	}
}

void Host :: ExpressKIRs(int numberOfExtraKirs)
{
	int counter = 0;
	vector<KIRGene>::iterator kirIt;
	if(tuning == false)
	{
		for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
			kirIt->SetGeneExpression(true);
		// if there is no education, all KIRs should be expressed!
	}
	else //otherwise express more KIRs besides those which are already functional
	{
		for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
		{
			if(kirIt->IsFunctional()) //KIRs that are functional are already expressed (see EducateKIRs!)
				continue;
			else
			{
				if(counter<numberOfExtraKirs)
				{
					kirIt->SetGeneExpression(true);
					counter++;
				}
				else
					return;
			}
		}
	}
}

/* This function counts the UNIQUE KIRs within one host that are functional, i.e. that are able to recognize their own MHC*/
void Host::CountFunctionalKIRs()
{
	for(unsigned int i= 0; i<kirGenes.size(); i++)
	{
		if(!kirGenes.at(i).IsFunctional()) //ignore non-functional KIRs
			continue;
		if(!kirGenes.at(i).IsActivating())
		{
			int flagNotUnique = 1;
			for(unsigned int j= 0; j<i; j++)
			{
				if(kirGenes.at(i).GetGeneID()==kirGenes.at(j).GetGeneID())
					flagNotUnique++;
			}
			if (flagNotUnique == 1)
				inhibitoryKIRs ++;
		}
		else //if(!kirGenes.at(i).IsInhibitory)
		{
			int flagNotUnique = 1;
			for(unsigned int j= 0; j<i; j++)
			{
				if(kirGenes.at(i).GetGeneID()==kirGenes.at(j).GetGeneID())
					flagNotUnique++;
			}
			if (flagNotUnique == 1)
				activatingKIRs ++;
		}
	}
	//cout << inhibitoryKIRs << "|" <<activatingKIRs <<endl;
}

/* This function counts ALL KIRs within one host that are expressed*/
int Host::CountExpressedKIRs()
{
	int kirsThatAreExpressed = 0;
	vector<KIRGene>::iterator kirIt;
	for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
	{
		if(kirIt->IsExpressed())
		{
				kirsThatAreExpressed ++;
		}
	}

	return kirsThatAreExpressed;
}
/*This function performs a deep copy*/
Host& Host:: Copy(Host& rightHandSideHost)
{
	// checking if it is a self-assignment
	if(this == &rightHandSideHost)
		return *this;

	//copying member variables
	this->dead = rightHandSideHost.dead;
	this->tuning = rightHandSideHost.tuning;
	this->age = rightHandSideHost.age;
	this->LOCI_KIR = rightHandSideHost.LOCI_KIR;
	this->LOCI_MHC = rightHandSideHost.LOCI_MHC;
	this->mainInfectionType = rightHandSideHost.mainInfectionType;
	this->mutationRateHost = rightHandSideHost.mutationRateHost;
	this->viralDeathRate = rightHandSideHost.viralDeathRate;
	this->inhibitoryKIRs = rightHandSideHost.inhibitoryKIRs;
	this->activatingKIRs = rightHandSideHost.activatingKIRs;

	//copy genes
	for(unsigned int i=0; i<rightHandSideHost.mhcGenes.size(); i++)
	{
		this->mhcGenes.push_back(rightHandSideHost.mhcGenes.at(i));
		//cout << "mhc genes: \t"<<mhcGenes.at(i).GetGeneID() << "\t" << rightHandSideHost.mhcGenes.at(i).GetGeneID() <<endl;
	}
	for(unsigned int i=0; i<rightHandSideHost.kirGenes.size(); i++)
	{
		this->kirGenes.push_back(rightHandSideHost.kirGenes.at(i));
	}
	
	list<Infection>::iterator it;
	for (it = rightHandSideHost.infections.begin() ; it !=rightHandSideHost.infections.end(); it++)
	{
		Infection dummy;
		dummy.Copy((*it));
		this->infections.push_back(dummy);
	}
	//check whether it works! do i need a copy function for the infection class???copy infections TO DO!!!!!!!!!!!!!!!!!!1



	return *this; // returns self-reference so cascaded assignment works
}

#include <boost/math/special_functions/round.hpp>
/*
double Host :: GetIntrinsicDeathRate(const vector<double>& rates)
{
	int a = round(age*52.0);
	if(a >= rates.size())
		cout << a << "|" << rates.size() << "age" << "ERROR in the deathrate!!!!" <<endl;
	intrinsicDeathRate = rates.at(a) + viralDeathRate;
// 	cout << "a, age, rate at, intrinsicdeathRates: "<< a << ", " << age << " , " << rates.at(a) << " , "<<intrinsicDeathRate << endl;
	return intrinsicDeathRate;
}

double Host :: GetAgeDependentBirthRate(vector<double>& rates)
{
	int a = round(age*52.0);
	if(a >= rates.size())
		cout << a << "|" << rates.size() << "age" << "ERROR in the birth rate!!!!" <<endl;
	ageDependentBirthrate = rates.at(a);
	return ageDependentBirthrate;
}
*/
/*this function sets the vector of infections*/
void Host::InfectWith(Virus& newVirus, double simulationTime, int maxNumberInfections) //function of the host receiving a new virus
{
	//cout << "do i get stuck here??? InfectWith()"<<endl;
	Infection newInfection;
	list<Infection>::iterator it;
	list<Infection>::iterator end = infections.end();

	if(!IsInfected()) // if the host is NOT infected yet
	{
		newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
		infections.push_back(newInfection);
	}
	else //but if the host IS infected with some viruses
	{
		if(infections.size()<maxNumberInfections) //check first whether it is already at its maximum
		{
			int howManyInfections = 0;
			for (it = infections.begin(); it != end; it++)
			{
				if(!it->IsPathogenNew(newVirus)) // if the virus is already there (i.e. as incubating, acute, chronic or immune) ignore it
					continue;
				else //but if it's not the same, keep track of how many infections are different from the new one
					howManyInfections++;
			}
			
			if(howManyInfections == infections.size()) // if the new virus is different from ALL infections present in that host
			{
				newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
				infections.push_back(newInfection);
			}
		}
	}
}

int Host :: IsInfected() //this functions basically counts whether there is more than zero infections!
{
	return infections.size();
}

/*This functions sets each host with an infection type according to the "main "state of all his infections*/
int Host::GetMainInfectionType()
{
	list<Infection>::iterator it;
	int chronicInfections = 0;
	int recoveredInfections = 0;

	for(it = infections.begin(); it!= infections.end(); it++)
	{
		int type = it->GetInfectionType();
		switch(type)
		{
		case 0:{mainInfectionType = susceptible;}break; //is incubating... so is not really infectious yet!
		case 1:
		{
			mainInfectionType = acute;
			return mainInfectionType;//if only one of them is acute return it as acute and leave the loop!
		}break;
		case 2:{chronicInfections++;} break; // is chronic
		case 3:{recoveredInfections++;} break; //is memory
		case 4:{mainInfectionType = susceptible;}break;//is cleared
		}
	}
	if(chronicInfections > 0)
		mainInfectionType = chronic;
	else
	{
		if(recoveredInfections > 0)
			mainInfectionType = immune;
	}
	return mainInfectionType;
}


void Host::ClearInfection(double simulationTime, Infection& _infection)
//void Host::ClearInfection(double simulationTime, int specificity)
{
	int virusType = _infection.pathogen.GetVirusType();
	switch (virusType)
	{
		case 0: //if it's a wild-type virus
		{
			if(RandomNumberDouble()<0.85)
			{
				_infection.ResetInfection(simulationTime);
				///cout << "clearing infection"<<endl;
			}
				
		}break;
		case 1://if it's an mHC downregulating virus
		{

			if(inhibitoryKIRs) //only if there is at least one functional INHIBITORY KIR, clear the infection
			{
				if(RandomNumberDouble()<0.5)
				{
					_infection.ResetInfection(simulationTime);
					//cout << "clearing infection"<<endl;
				}
					
			}
		}break;
		case 2: //if it's a decoy virus
		{
			int inhibiting_kirs_recognizing_decoy = 0;
			int activating_kirs_recognizing_decoy = 0;
			vector<KIRGene>::iterator kirIt;
			//cout <<"decoy | ";
			for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
			{
				if(!kirIt->IsExpressed()) //ignore KIRs that are not functional / expressed!
					continue;
				else
				{
					int score = kirIt->BindMolecule(_infection.pathogen.mhcDecoy); //check if they bind to the decoy!
					if(score>=kirIt->GetGeneSpecificity()) //receptors binds to decoy
					{
						//now check what kind of receptors we have
						if(kirIt->IsInhibitory()) //if it is inhibiting
						{
							inhibiting_kirs_recognizing_decoy++;//keep track of how many KIRs are recognizing decoys
							//cout <<"inhibiting| "<<inhibiting_kirs_recognizing_decoy;
						}
						
						if(kirIt->IsActivating()) //if it is activating
						{
							activating_kirs_recognizing_decoy ++;//keep track of how many KIRs are recognizing decoys
							//cout <<"activating| "<<activating_kirs_recognizing_decoy;
						}
					}
				}
			}
			
			ClearDecoyWithActivatingAndInhibitory(inhibiting_kirs_recognizing_decoy, activating_kirs_recognizing_decoy, simulationTime, _infection);
			
			/*This is with different of protection depending which receptors we have. Above is the
			 *same... I don't completely agree, but let's see what happens! -> try out the new one!
			 *
			 * if(inhibitoryKIRs && !activatingKIRs) //if that host has ONLY inhibitory receptors
			 {
			 ClearDecoyWithInhibitoryOnly(inhibiting_kirs_recognizing_decoy, simulationTime);
			 }
			 
			 if(!inhibitoryKIRs && activatingKIRs) //if that host has ONLY activating receptors
			 {
			 ClearDecoyWithActivatingOnly(activating_kirs_recognizing_decoy, simulationTime);
			 }
			 
			 if(inhibitoryKIRs && activatingKIRs) //if that host has both types of receptors
			 {
			 ClearDecoyWithActivatingAndInhibitory(inhibiting_kirs_recognizing_decoy, activating_kirs_recognizing_decoy, simulationTime);
			 }*/
		}break;
	}
}

void Host :: ClearDecoyWithInhibitoryOnly(int inhibiting_signal, double simulationTime, Infection& _infection)
{
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	if(!inhibiting_signal) //if there is no inhibiting signal, receptor didn't bind to the decoy: protection like MHC down
	{
		if(RandomNumberDouble()<0.5)
		{
			_infection.ResetInfection(simulationTime);
			cout <<"clearing infection..."<<endl;
		}
			
	}
}

void Host :: ClearDecoyWithActivatingOnly(int activating_signal, double simulationTime, Infection& _infection)
{
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	if(activating_signal) //if activating receptor recognizes decoy, but there are no inhibitory receptors, the virus still escapes response of the T-cells
	{                       // protection as with an MHC-downregulating one
		if(RandomNumberDouble()<0.5)
			_infection.ResetInfection(simulationTime);
	}
	
}

void Host ::ClearDecoyWithActivatingAndInhibitory(int inhibiting_signal, int activating_signal, double simulationTime, Infection& _infection)
{
	//if functional KIRs recognize MHC down
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	if(inhibitoryKIRs > 0 && activating_signal) //and there is enough activating signal
	{
		if(!inhibiting_signal) //best protection
		{
			if(RandomNumberDouble()<0.5)
			{
				//cout <<"best protection!" <<endl;
				_infection.ResetInfection(simulationTime);
				return;
			}
		}
		if(inhibiting_signal) //good protection
		{
			if(RandomNumberDouble()<0.5)
			{
				//cout <<"good protection!" <<endl;
				_infection.ResetInfection(simulationTime);
				return;
			}
		}
	}
	//if functional KIRs recognize MHC down
	if(inhibitoryKIRs > 0 && !activating_signal)	//but there is not enough activation OR only inhibitory receptors!
	{
		//check whether the host is fooled or not
		if(!inhibiting_signal) // if no inhibitory receptor recognizes the decoy:
		{                                             //same as MHC downregulation p= 0.5
			if(RandomNumberDouble()<0.5)
			{
				//cout <<"same as mhc!" <<endl;
				_infection.ResetInfection(simulationTime);
				return;
			}
		}
		//cout <<"screwed!!!!!!!!!" <<endl;
	}
}

/*This functions returns the Virus of the acute infection*/
Virus& Host :: GetAcuteInfection()
{
	list<Infection>:: iterator it;
	vector<Virus> randomAcuteInfections;;
	int inf_type = 0;
	for(it = infections.begin(); it!= infections.end(); it++)
	{
		inf_type = it->GetInfectionType();
		if(inf_type != 1) //if it's not an acute infection
			continue;
		else //if it is, return the acute virus
			randomAcuteInfections.push_back(it->pathogen);
	}
	
	int random = RandomNumber(0,randomAcuteInfections.size()-1);
	//cout <<"do i get stuck here? GetAcuteInfection!!! "<<endl;
	return randomAcuteInfections.at(random); // in case the list is empty, return an empty virus... shoudln't happen anyway! if the list is empty, this function should not be called!
}

/*This functions returns a ramdonly chose chronic virus infection  !!!!!! WHAT IS HAPPENING TO THIS VECTOR???? SHALL I DELETE IT SOMEHOW????*/
Virus& Host :: GetChronicInfection()
{
	vector<Virus> randomChronicInfections;
	list<Infection>:: iterator it;
	int inf_type = 0;
	for(it = infections.begin(); it!= infections.end(); it++)
	{
		inf_type = it->GetInfectionType();
		//cout << inf_type << "infection type! "<<endl;
		if(inf_type!=2) //if it's not chronic
			continue;
		else
			randomChronicInfections.push_back(it->pathogen);
	}
	int random = RandomNumber(0,randomChronicInfections.size()-1);
	//cout << randomChronicInfections.size() << "   size....." <<endl;
	return randomChronicInfections.at(random);
}

void Host :: UpdateParameters(double timeStep, double simulationTime)
{
	age += (timeStep/YEAR);
	
	//update (get) all the information for each infection
	double vl = 0.0;
	list<Infection>::iterator it = infections.begin();
	while(it != infections.end())
	{
		it->SetInfectionType(simulationTime);
		if(it->IsCured())
		{//if the infection has been cleared again,remove it from the infections list!
			it = infections.erase(it);
			continue;
		}
		if(it->GetViralLoad() > vl)
			vl = it->GetViralLoad();
		if(it->GetDeadFlag()==true)
			SetDead();
		
		it++;
	}
	viralDeathRate = vl;
}

//functions for SAVING PARAMETERS, BACKUP, etc
/*this function saves each locus -> to keep track of KIR and MHC diversity*/
void Host::SaveGenes(fstream& outfile)
{
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	for(mhcIt = mhcGenes.begin(); mhcIt!=mhcGenes.end(); mhcIt++)
	{
		mhcIt->SaveGenes(outfile);
		//outfile << mhcGenes.at(i).GetGeneID() << "\t";
	}
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveGenes(outfile);
		//outfile << kirGenes.at(i).GetGeneID() << "\t";
	}
	int functional_kirs = inhibitoryKIRs + activatingKIRs;
	outfile << functional_kirs << "\t"<< CountExpressedKIRs() << "\t" << Get_Host_ID()<<"\n";
}

void Host ::SaveParameters(fstream& outfile)
{
	outfile <<age << "\t"<<infections.size()<< "\t";
	list<Infection>::iterator infIt;
	for(infIt = infections.begin(); infIt!=infections.end(); infIt++)
	{
		 infIt->SaveParametersInfection(outfile);
	}
	outfile << Get_Host_ID() <<  "\n";
}

void Host :: PrintParametersHost()
{
	cout <<age << "\t"<<infections.size()<< "\t"<< acute << chronic << immune << "\t";
	list<Infection>::iterator infIt;
	for(infIt = infections.begin(); infIt!=infections.end(); infIt++)
	{
		infIt->PrintParameters();
	}
	cout <<  "\n";
}

void Host :: SaveAgeDyingHost(fstream& outfile)
{
	outfile << age << "\t";
	vector<KIRGene>::iterator kirIt;
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveGenes(outfile);
	}
	//pathogen.SaveParametersVirus(outfile);
	//outfile << infectionType << "\t"<< ageInfection << "\t"<< ageClearance<<"\t";
}

void Host::SaveBackupHost(fstream& backupFile)
{
	//cout << "saving hosts..."<<endl;
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	list<Infection>::iterator infIt;

	backupFile << LOCI_KIR << "\t"<<LOCI_MHC<<"\t"<< age << "\t"<< tuning << "\t"<< dead << "\t"<< mutationRateHost << "\t"<<inhibitoryKIRs<<"\t"<<activatingKIRs<<"\t"<< infections.size()<<"\t";

	for(infIt = infections.begin(); infIt != infections.end(); infIt++)
	{
		infIt->SaveBackupInfection(backupFile);
	}

	for(mhcIt = mhcGenes.begin(); mhcIt!=mhcGenes.end(); mhcIt++)
	{
		mhcIt->SaveBackupGenes(backupFile);
	}
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveBackupGenes(backupFile);
	}
	backupFile << "\n";
}

void Host::RestoreHost(const string& sline)
{
	stringstream ssline(sline);
	ssline >> LOCI_KIR;
	ssline >> LOCI_MHC;
	ssline >>age;
	ssline >> tuning;
	ssline >> dead;
	ssline >> mutationRateHost;
	ssline >> inhibitoryKIRs;
	ssline >> activatingKIRs;
	int totalInfections;
	ssline >> totalInfections;
//	cout << LOCI_KIR<<"\t"<<age << "\t"<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<< tuning << "\t"<< dead << "\t"<< mutationRateHost << "\t"<< viralDeathRate << "\t"<<endl;


	string infectionString;
	for(int i=0; i<totalInfections; i++)
	{
		Infection dummyInfection;
		infectionString = dummyInfection.RestoreInfection(ssline);
		ssline.str() = infectionString;
	}
	string geneString = infectionString;
	ssline.str() = geneString;

	for(int i=0; i<LOCI_MHC*TWO; i++)
	{
		Gene mhc;
		geneString = mhc.RestoreGenes(ssline);
		mhcGenes.push_back(mhc);
		ssline.str() = geneString;
	}

	for(int i=0; i<LOCI_KIR*TWO; i++)
	{
		KIRGene kir;
		geneString = kir.RestoreGenes(ssline);
		kirGenes.push_back(kir);
//		cout <<"id:\t"<<kir.GetGeneID()<<" func:\t"<<kir.IsFunctional()<<" kirSize: "<< kirGenes.size()<<endl;
//		cout <<ssline.str()<<endl;
		ssline.str() = geneString;
	}
}


//FUNCTIONS OF CLASS INFECTION
Infection::Infection()
{
	SetInfectionParameters(cleared);
}

Infection::Infection(Virus& _virus)
{
	pathogen.Copy(_virus);
	SetInfectionParameters(cleared);

}

void Infection :: SetInfectionParameters(state _type)
{
	infectionType = _type;
	infectionTime = 0.0;
	immunityTime = 0.0;
	clearanceTime = 0.0;
	deadFlagForHost = false;
	//cout << incubating  <<acute<< chronic<< memory<< cleared <<endl;
	//ageInfection = 0.0;
	//ageClearance = 0.0;
}

/*This function sets the pathogen parameters as infectious*/
void Infection::TransmitInfection(Virus& nastyInfection, double simulationTime)
{
	pathogen.Copy(nastyInfection); //infect it with the pathogen
	pathogen.SetViralLoad(pathogen.GetOriginalViralLoad());//get the original one for a new infection being transmission!
	infectionTime = simulationTime;
	if(pathogen.IsOnlyAcute())
		immunityTime = 80;
	else
		immunityTime = 10;
}

bool Infection :: IsPathogenNew(Virus& _newVirus)
{
	if(pathogen == _newVirus) //if the host is already infected with the same virus
		return false; //don't infect it
	else //but if it is a different one
		return true;
}

bool Infection::IsCured()
{
	if(infectionType == cleared)
		return true;
	else
		return false;
}

bool Infection :: IsAcute()
{
	if(infectionType == acute)
		return true;
	else
		return false;
}

void Infection::ResetInfection(double simulationTime)
{
	//set the viral load to zero, but keep the viral type to keep track of which infections is still present in the host!
	pathogen.SetViralLoad(0.0);
	infectionTime = 0.0;
	clearanceTime = simulationTime;
	deadFlagForHost = false;
}


/*This function changes the infection type of the host according to the time of infection*/
void Infection::SetInfectionType(double simulationTime)
{
	//infection has started but has not been cleared yet:
	if(infectionTime > 0.0 && clearanceTime == 0.0)
	{
		//check first if the host is in the infection period (i.e. between 1-4 weeks)
		if((simulationTime-infectionTime)<=1.0*WEEK)
			infectionType = incubating;

		if((simulationTime-infectionTime)>1.0*WEEK && (simulationTime-infectionTime)<1.0*WEEK+4.0*WEEK*pathogen.GetLifeTimeVirus())
		{
			//viralDeathRate = pathogen.GetViralLoad();STILL NEED TO IMPLEMENT THIS IN THE HOST!!!!!!!!
			infectionType = acute;
		}
		if((simulationTime -infectionTime) > 1.0*WEEK+ 4.0*WEEK*pathogen.GetLifeTimeVirus())
		{
			if(pathogen.IsOnlyAcute())
				deadFlagForHost = true;
			else
			{
				pathogen.SetViralLoad(0.6*pathogen.GetOriginalViralLoad());
				infectionType = chronic;
			}
		}
	}

	if(infectionTime == 0.0 && clearanceTime > 0.0)
	{

		if(simulationTime - clearanceTime < ImmunityTime(immunityTime,0.5)*YEAR)
		{
			infectionType = memory;
		}
		else
		{
			infectionType = cleared;
			clearanceTime = 0.0;
		}
	}

	if(infectionTime == 0.0 && clearanceTime == 0.0)
	{
		Virus deadVirus;
		pathogen.Copy(deadVirus);
		infectionType = cleared;
	}
}
void Infection::PrintParameters()
{
	cout <<"infectionType infectionTime clearanceTime immunityTime and other viral parameters"<<endl;
	cout << infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<<immunityTime << "\t";
	pathogen.PrintParametersVirus();
}

void Infection ::SaveParametersInfection(fstream& outfile)
{
	outfile << infectionTime/YEAR <<"\t" << infectionType << "\t";
	pathogen.SaveParametersVirus(outfile);

}

void Infection::SaveBackupInfection(fstream& backupFile)
{
	backupFile<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<<immunityTime << "\t";
	pathogen.SaveBackupVirus(backupFile);
}


string Infection::RestoreInfection(stringstream& siline)
{
	int inf;
	siline >> inf;
	//cout <<"inf:  "<<inf <<endl;
	switch(inf)
	{
	case 0:{infectionType = cleared;}break;
	case 1:(infectionType = incubating);break;
	case 2:{infectionType = acute;}break;
	case 3:{infectionType = chronic;}break;
	case 4:{infectionType = memory;}break;
	}
	//cout << "infection Type: "<<infectionType <<endl;
	siline >> infectionTime;
	siline >> clearanceTime;
	siline >> immunityTime;

	string istring = siline.str();
	return istring;
}


/*This function performs a deep copy*/
Infection & Infection::Copy(Infection& rhsInfection)
{
	// checking if it is a self-assignment
	if(this == &rhsInfection)
		return *this;
	//copying member variables
	this->infectionType = rhsInfection.infectionType;
	this->infectionTime = rhsInfection.infectionTime;
	this->clearanceTime = rhsInfection.clearanceTime;
	this->immunityTime = rhsInfection.immunityTime;
	this->pathogen.Copy(rhsInfection.pathogen);
	return *this;
}


/*
 * to do:
 * write down in the file how many infections the host  has
 */
