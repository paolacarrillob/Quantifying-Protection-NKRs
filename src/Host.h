/*
 * Host.h
 *
 *  Created on: Apr 21, 2011
 *      Author: paola
 */

#ifndef HOST_H_
#define HOST_H_
#include "Genes.h"
#define TWO 2
//#define LOCI_KIR 5
#define LOCI_MHC 1
const double YEAR = 31207680.0;
const double MONTH = 2600640.0;
const double WEEK = 604800.0;

class Virus{
public:
	Virus();
	virtual ~Virus(){};

	void SetViralParameters(double _downregulation, double _decoy, double _viralLoad, double _lifeTimeVirus, int virus_type, bool only_acute);
	enum type{wildType, downregulating, decoy};

	Virus& Copy(Virus& rhsVirus);//works
	double GetLifeTimeVirus()const{return lifeTimeVirus;}
	double GetViralLoad()const{return viralLoad;}
	double GetRateDownregulation()const{return mutationRateDownregulation;}
	double GetRateDecoy()const{return mutationRateDecoy;}
	type GetVirusType()const{return virusType;}
//	void BuildDecoy(Gene *mhcGene);
	void BuildDecoy(int mhcID);
	void DownregulateMHC();
	bool IsDownregulatingMHC();
	bool IsStealingDecoy();
	bool IsWildType();
	bool IsOnlyAcute()const{return onlyAcute;};

	void SaveBackupVirus(fstream& backupFile);//works
	string RestoreVirus(stringstream& svline);//works

	void SaveParametersVirus(fstream& outfile);
	Gene mhcDecoy;

protected:
	type virusType;
	double mutationRateDownregulation;
	double mutationRateDecoy;
	double viralLoad; //as the increase of the intrinsic death rate
	double lifeTimeVirus; // time that a virus can live in one's organism
	bool onlyAcute;

};



class Host {
public:
	Host(){}; //for creating memory
	Host(int loci_kir, int mhc1, int mhc2, double _mutationRate, bool _tuning,
			int numberOfExtraKirs, Map& kirMap /*MHCGenePool& mhcPool, bool hla, */); //for initialization of the population //works

	//This is where all my troubles come from
	Host(int loci_kir, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist,
			vector<KIRGene>& kirGenesMother, vector<KIRGene>& kirGenesFather, double _mutationRate,
			bool _tuning,int numberOfExtraKirs, Map& kirMap, int mutationType);


	//Host(int loci_kir, Gene* mhcGenesParent, int mhcGenePool, Gene* kirGenesMother, Gene* kirGenesFather, int specificity, double _mutationRate, bool _tuning); //for birth function
	virtual ~Host(){};
	void SetDead(){dead = true;}
	bool IsDead()const{return dead;}
//	bool IsHostToBeTuned()const{return tuning;}

	void InitializeHostParameters(double mutationRate, bool _tuning, int kirloci);
	void SetHostParameters(bool t, double mut, int inftyp, double inftim, double viraldeathm ,double clrtim);
	void EducateKIRs(); //works with activating receptors also!
	//void EducateKIRs(int specificity);//works
	void MutateGenes(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool);


	void ExpressKIRs(int numberOfExpressedKirs);
	double GetAge()const{return age;}
	void SetAge(double number){age = number;}
	void CountFunctionalKIRs();//works
	int CountExpressedKIRs();
	Host& Copy(Host& rightHandSideHost);//works

	double GetIntrinsicDeathRate(vector<double>& rates);
	double GetAgeDependentBirthRate(vector<double>& rates);

	bool IsSusceptible();
	bool IsImmune();
	bool IsChronicInfected();
	bool IsAcuteInfected();
	bool IsIncubating();

	void InfectWith(Virus& nastyVirus, double simulationTime);//works
	void ResetInfection(double simulationTime);
	void ClearInfection(double simulationTime);
	void ClearDecoyWithInhibitoryOnly(int inhibiting_signal, double simulationTime);
	void ClearDecoyWithActivatingOnly(int activating_signal, double simulationTime);
	void ClearDecoyWithActivatingAndInhibitory(int inhibiting_signal, int activating_signal, double simulationTime);

	void SetInfectionType(double simulationTime, double lifeTimeVirus);
	double GetInfectionTime(){return infectionTime;};
	double GetClearanceTime(){return clearanceTime;};
	int GetInfectionType(){return infectionType;};
	double GetMutationRate(){return mutationRateHost;};
	double GetViralDeathRate(){return viralDeathRate;}

	void UpdateParameters(double timeStep, double simulationTime,double lifeTimeVirus);//works
	void SaveGenes(fstream& outfile);//works
	void SaveParameters(fstream& outfile);//works
	void SaveAgeDyingHost(fstream& outfile);

	void SaveBackupHost(fstream& file);//works
	void RestoreHost(const string& sline);//works

	//Gene mhcGenes[LOCI_MHC][TWO];
	//Gene kirGenes[LOCI_KIR][TWO];
	vector<Gene> mhcGenes;
	vector<KIRGene> kirGenes;
	enum state{susceptible, incubating, acute, chronic, immune};

	Virus pathogen;
protected:

	int LOCI_KIR;
	//int LOCI_MHC;
	double age;
	bool dead;
	bool tuning;
//	double tuning_rate;
	double mutationRateHost;

	double intrinsicDeathRate;
	double ageDependentBirthrate;
	double viralDeathRate;
    state infectionType;

    double infectionTime;
    double clearanceTime;
    int immunityTime;
    double ageInfection;
    double ageClearance;

    int inhibitoryKIRs; //number of FUNTIONAL inhibitory / activating receptors
    int activatingKIRs;
    //int MHCsearchinPool;
    //int functionalKIRs;
};

#endif /* HOST_H_ */

/*
 * 23.08.2012: removed int specificity from all function arguments -> all specificity values were taken from genes.GetSpecificity() ....
 *
 * 06.09.2012: changed Host constructors. New argument: kirMap. Changed the mutation-> now it picks a new gene
 *
 * 18.10.2012: changed Host constructors: introduced a new type of mutation-> to PointMutation() and MutateSpecificity()
 */
