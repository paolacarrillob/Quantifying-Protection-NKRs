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
	this->virusType = rhsVirus.virusType;
	//this->mhcDecoy.SetGeneID(rhsVirus.mhcDecoy.GetGeneID());
	this->mhcDecoy.Copy(rhsVirus.mhcDecoy);
	this->onlyAcute = rhsVirus.onlyAcute;
	return *this;
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
	file << mutationRateDownregulation << "\t" << mutationRateDecoy << "\t" << lifeTimeVirus << "\t" << viralLoad <<"\t" << virusType <<"\t" << mhcDecoy.GetGeneID()<< "\t"<<onlyAcute <<"\t";
}

string Virus::RestoreVirus(stringstream& svline)
{
	int type;
	int id;
	svline >> mutationRateDownregulation;
	svline >> mutationRateDecoy;
	svline >> lifeTimeVirus;
	svline >> viralLoad;
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
	outfile << virusType << "\t" << viralLoad <<"\t" << mhcDecoy.GetGeneID() <<"\t"<<onlyAcute<<"\t";
}


/*FUNCTIONS OF CLASS HOST
 *Constructs a host for the initialization of the population: it fills the MHC genes with a randomly picked allele from the population
 * and creates KIRs that match their own MHC according to the specificity*/
Host::Host(int loci_kir, int loci_mhc, double _mutationRate, bool _tuning, int numberOfExtraKirs,Map& kirMap, MHCGenePool& mhcPool, bool hla) {

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
		kir.SetGeneType(0); // 0 for inhibitory receptor!
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
Host::Host(int loci_kir, int loci_mhc, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, vector<KIRGene>& kirGenesMother, vector<KIRGene>& kirGenesFather,double _mutationRate, bool _tuning, int numberOfExtraKirs, Map& kirMap, int mutationType)
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
				MutateGenes(mutationType, kir_hap2, kirMap, mhcPool);
			}
			kirGenes.push_back(kir_hap2);
		}
		for(int j=KIR_init_mum; j<KIR_end_mum; j++)
		{
			KIRGene kir_hap1;
			kir_hap1.Copy(kirGenesMother.at(j));
			if(RandomNumberDouble() < mutationRateHost)
			{
				MutateGenes(mutationType, kir_hap1, kirMap, mhcPool);
			}
			kirGenes.push_back(kir_hap1);
		}

		//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
		{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
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
				MutateGenes(mutationType, kir_hap1, kirMap, mhcPool);
			}
			kirGenes.push_back(kir_hap1);
		}
		for(int j=KIR_init_dad; j<KIR_end_dad; j++)
		{
			KIRGene kir_hap2;
			kir_hap2.Copy(kirGenesFather.at(j));
			if(RandomNumberDouble() < mutationRateHost)
			{
				MutateGenes(mutationType, kir_hap2, kirMap, mhcPool);
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
			cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
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
	LOCI_KIR = loci_kir;
	LOCI_MHC = loci_mhc;
	infectionType = susceptible;
	infectionTime = 0.0;
	immunityTime = 0.0;
	viralDeathRate = 0.0;
	clearanceTime = 0.0;
	ageInfection = 0.0;
	ageClearance = 0.0;
	age = 0.0;
	inhibitoryKIRs = 0;
	activatingKIRs = 0;
	//MHCsearchinPool = 0;

}

void Host :: MutateGenes(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool)
{
	if(mutationType == 1)//pick another molecules as mutation
	{
		KIRGene newGene(RandomNumber(2,16));
		if(!kirMap.IsGeneInMap(newGene))
		{
			kirMap.FillMap(mhcPool, newGene);
			//kir_hap2.SetGeneType(1);////to test just what happens if i ONLY have activating receptors!
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
		if(RandomNumberDouble()<0.2)
		{
			kir_hap2.MutateReceptorType();
		}
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
	this->infectionType = rightHandSideHost.infectionType;
	this->infectionTime = rightHandSideHost.infectionTime;
	this->immunityTime = rightHandSideHost.immunityTime;
	this->mutationRateHost = rightHandSideHost.mutationRateHost;
	this->viralDeathRate = rightHandSideHost.viralDeathRate;
	this->clearanceTime = rightHandSideHost.clearanceTime;
	this->ageInfection = rightHandSideHost.ageInfection;
	this->ageClearance = rightHandSideHost.ageClearance;
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

	return *this; // returns self-reference so cascaded assignment works
}

double Host :: GetIntrinsicDeathRate(vector<double>& rates)
{
	int a = int(ceil(age));
	//cout << a<<"|"<<age <<"\t";
	if(a > rates.size())
		cout << a << "|" << rates.size() << "age" << "ERROR in the deathrate!!!!" <<endl;

	intrinsicDeathRate = rates.at(a) + viralDeathRate;
	return intrinsicDeathRate;
}

double Host :: GetAgeDependentBirthRate(vector<double>& rates)
{
	int a = int(ceil(age));
	//cout << a<<"|"<<age <<"\t";
	if(a > rates.size())
		cout << a << "|" << rates.size() << "age" << "ERROR in the birth rate!!!!" <<endl;
	ageDependentBirthrate = rates.at(a);
	return ageDependentBirthrate;
}

/*This functions return the infection type of the host*/
bool Host :: IsSusceptible()
{
	if (infectionType == susceptible)
		return true;
	else
		return false;
	//return (infection == susceptible);
}
bool Host :: IsImmune()
{
	if (infectionType == immune)
		return true;
	else
		return false;
}
bool Host:: IsIncubating()
{
	if(infectionType == incubating)
		return true;
	else
		return false;
}
bool Host :: IsAcuteInfected()
{
	if (infectionType == acute)
		return true;
	else
		return false;
}
bool Host :: IsChronicInfected()
{
	if(infectionType == chronic)
		return true;
	else
		return false;
}

/*This function sets the pathogen parameters as infectious*/
void Host :: InfectWith(Virus& nastyVirus, double simulationTime)
{
	if(infectionType == susceptible||infectionType == incubating)
	{
		pathogen.Copy(nastyVirus);
		infectionTime = simulationTime;
		ageInfection = age;
		if(pathogen.IsOnlyAcute())
			immunityTime = 80;
		else
			immunityTime = 10;
	}
}

/*This function changes the infection type of the host according to the time of infection*/
void Host::SetInfectionType(double simulationTime, double lifeTimeVirus)
{
	//infection has started but has not been cleared yet:
	if(infectionTime > 0.0 && clearanceTime == 0.0)
	{
		//check first if the host is in the infection period (i.e. between 1-4 weeks)
		if((simulationTime-infectionTime)<=1.0*WEEK)
			infectionType = incubating;

		if((simulationTime-infectionTime)>1.0*WEEK && (simulationTime-infectionTime)<1.0*WEEK+4.0*lifeTimeVirus)
		{
			viralDeathRate = pathogen.GetViralLoad();
			infectionType = acute;
		}
		if((simulationTime -infectionTime) > 1.0*WEEK+ 4.0*lifeTimeVirus)
		{
			if(pathogen.IsOnlyAcute())
				dead = true;
			else
			{
				viralDeathRate = 0.6*pathogen.GetViralLoad();
				infectionType = chronic;
			}
		}
	}

	if(infectionTime == 0.0 && clearanceTime > 0.0)
	{

		if(simulationTime - clearanceTime<ImmunityTime(immunityTime,0.5)*YEAR)
		{
			infectionType = immune;
		}
		else
		{
			infectionType = susceptible;
			clearanceTime = 0.0;
		}
	}

	if(infectionTime == 0.0 && clearanceTime == 0.0)
		infectionType = susceptible;
}

void Host :: UpdateParameters(double timeStep, double simulationTime, double lifeTimeVirus)
{
	age += (timeStep/YEAR);
	SetInfectionType(simulationTime, lifeTimeVirus);

}

void Host::ClearInfection(double simulationTime)
//void Host::ClearInfection(double simulationTime, int specificity)
{
	if(pathogen.IsWildType()) //wildtype virus, regular infection
	{
		if(RandomNumberDouble()<0.85)
		{
			//cout <<"wild-type" <<endl;
			ResetInfection(simulationTime);
			return;
		}
	}
	if (pathogen.IsDownregulatingMHC()) //missing self
	{
		//change the model and make functional KIRs important also for the recognition of MHC downregulation

		if(inhibitoryKIRs) //only if there is at least one functional INHIBITORY KIR, clear the infection
		{
			if(RandomNumberDouble()<0.5)
			{
				//cout <<"mhc down" <<endl;
				ResetInfection(simulationTime);
				return;
			}
		}
	}
	if(pathogen.IsStealingDecoy()) //decoy+missing self
	{
		int inhibiting_kirs_recognizing_decoy = 0;
		int activating_kirs_recognizing_decoy = 0;
		vector<KIRGene>::iterator kirIt;
		//cout <<"decoy | ";
		for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
		{
			if(/*!kirGenes.at(i).IsFunctional() || */!kirIt->IsExpressed()) //ignore KIRs that are not functional / expressed!
				continue;
			else
			{
				int score = kirIt->BindMolecule(pathogen.mhcDecoy); //check if they bind to the decoy!
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

		/*if(inhibitoryKIRs && !activatingKIRs) //if that host has ONLY inhibitory receptors
		{
			ClearDecoyWithInhibitoryOnly(inhibiting_kirs_recognizing_decoy, simulationTime);
		}*/

		if(!inhibitoryKIRs && activatingKIRs) //if that host has ONLY activating receptors
		{
			ClearDecoyWithActivatingOnly(activating_kirs_recognizing_decoy, simulationTime);
		}

		if(inhibitoryKIRs && activatingKIRs) //if that host has both types of receptors
		{
			ClearDecoyWithActivatingAndInhibitory(inhibiting_kirs_recognizing_decoy, activating_kirs_recognizing_decoy, simulationTime);
		}
	}
}

void Host::ResetInfection(double simulationTime)
{
	Virus deadVirus;
	viralDeathRate = 0.0;
	pathogen.Copy(deadVirus);
	infectionTime = 0.0;
	clearanceTime = simulationTime;
	ageClearance = age;
}

void Host :: ClearDecoyWithInhibitoryOnly(int inhibiting_signal, double simulationTime)
{
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	if(!inhibiting_signal) //if there is no inhibiting signal, receptor didn't bind to the decoy: protection like MHC down
	{
		if(RandomNumberDouble()<0.5)
			ResetInfection(simulationTime);
	}
}

void Host :: ClearDecoyWithActivatingOnly(int activating_signal, double simulationTime)
{
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	if(activating_signal) //if activating receptor recognizes decoy, but there are no inhibitory receptors, the virus still escapes response of the T-cells
	{                       // protection as with an MHC-downregulating one
		if(RandomNumberDouble()<0.5)
			ResetInfection(simulationTime);
	}

}

void Host ::ClearDecoyWithActivatingAndInhibitory(int inhibiting_signal, int activating_signal, double simulationTime)
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
				ResetInfection(simulationTime);
				return;
			}
		}
		if(inhibiting_signal) //good protection
		{
			if(RandomNumberDouble()<0.5)
			{
				//cout <<"good protection!" <<endl;
				ResetInfection(simulationTime);
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
				ResetInfection(simulationTime);
				return;
			}
		}
		//cout <<"screwed!!!!!!!!!" <<endl;
	}
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
	outfile << functional_kirs << "\t"<< CountExpressedKIRs() <<"\n";
}

void Host ::SaveParameters(fstream& outfile)
{
	outfile <<age << "\t" << infectionTime/YEAR <<"\t" << infectionType << "\t"<< viralDeathRate<<"\t"<<ageInfection << "\t"<< ageClearance<<"\t";
	pathogen.SaveParametersVirus(outfile);
	outfile <<  "\n";
}


void Host :: SaveAgeDyingHost(fstream& outfile)
{
	outfile << age << "\t";
	vector<KIRGene>::iterator kirIt;
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveGenes(outfile);
	}
	pathogen.SaveParametersVirus(outfile);
	outfile << infectionType << "\t"<< ageInfection << "\t"<< ageClearance<<"\t";
}


void Host::SaveBackupHost(fstream& backupFile)
{
	//cout << "saving hosts..."<<endl;
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	backupFile << LOCI_KIR << "\t"<<LOCI_MHC<<"\t"<< age << "\t"<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<<immunityTime << "\t"<< tuning << "\t"
	<< dead << "\t"<< mutationRateHost << "\t"<< viralDeathRate << "\t"<<ageInfection <<"\t"<<ageClearance <<"\t"<<inhibitoryKIRs<<"\t"<<activatingKIRs<<"\t";

	pathogen.SaveBackupVirus(backupFile);

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
	//Ouss ssline >> infectionType;
	int inf;
	ssline >> inf;
	//cout <<"inf:  "<<inf <<endl;
	switch(inf)
	{
	case 0:{infectionType = susceptible;}break;
	case 1:(infectionType = incubating);break;
	case 2:{infectionType = acute;}break;
	case 3:{infectionType = chronic;}break;
	case 4:{infectionType = immune;}break;
	}
	//cout << "infection Type: "<<infectionType <<endl;
	ssline >> infectionTime;
	ssline >> clearanceTime;
	ssline >> immunityTime;
	ssline >> tuning;
	ssline >> dead;
	ssline >> mutationRateHost;
	ssline >> viralDeathRate;
	ssline >> ageInfection;
	ssline >> ageClearance;
	ssline >> inhibitoryKIRs;
	ssline >> activatingKIRs;
//	cout << LOCI_KIR<<"\t"<<age << "\t"<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<< tuning << "\t"<< dead << "\t"<< mutationRateHost << "\t"<< viralDeathRate << "\t"<<endl;
	string geneString = pathogen.RestoreVirus(ssline);
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
//	cout << "halloooooooooooooooooooooooooooo"<<endl;

}
