/*
 * World.h
 *
 *  Created on: Apr 26, 2011
 *      Author: paola
 */

#ifndef WORLD_H_
#define WORLD_H_
#include "Host.h"//Host.h includes Virus.h and Gene.h which includes Bitstrings.h and Mathfunctions.h

class World {
public:
	World(); //Default constructor
	virtual ~World(){}; //destructor that so far does nothing at all!

	void WriteInfo();
	void LoadParameterFile(const string& fileName);//works
	bool Initialize(); //works
	void CreateBirthAndDeathRates();

	bool Birth(int index, unsigned long int next_id);//, Host& host); 	// Birth() function takes a host at index (index) and picks a random host at a random index and create a "child" that is copied to (host).
										// It returns true or false according if the child should be added to the population or not
	bool Death (int index);
	void Infect(int index);
	void Escape(int index);

	void ShuffleHosts();
	void IntroduceVirus(const string& secondVirus);
	void RemoveDeadHosts();
	void RemoveDeadHosts_HappyNewYear();
	void TrackInfectedIndividuals();
	
	double GetAgeDependentBirthRate(const double age)const;
	double GetIntrinsicDeathRate(const double age, const double viralrate)const;
	
	void Simulate();

	void SaveGenes();//works
	//void SavePopulationSize(int babies,int dead_people);//works
	void SavePopulationSize();//works
	void SaveParameters();//works
	void SaveMap();
	bool SaveAgeDyingHosts(double lastOutfileTime, double lastStopOutfileTime);

	void SaveBackupFile();//works
	void LoadBackupFile(const string& backupName);//works

	vector<Host> hosts; //vector of hosts
	MHCGenePool MHCPool; //vector of integers representing MHC pools
	//GenePool KIRPool;
	Map KIRGenesMap; //map of KIR_id|L|M_id
	Virus nastyVirus;
	Virus downregulatingVirus;
	Virus decoyVirus;
	

	double mutationRate;
	bool education;
	int expressionExtraKIRs;
	bool HLA_C;
	int KIRspecificity;
	int KIRLoci;
	int KIRGeneType;
	int MHCLoci;
	int sizeMHCPool;
	unsigned int initHostPop;
	int mutationType;
	int mutationTypeHost;

protected:
	double birthRate;
	double deathRate;
	double infectionRate;
	double escapeRate;

	double timeStep;
	double timeEnd;
	double simulationTime;
	double outfileRate;
	double backupRate;
	double populationSizeRate;
	int maxHostPop;

	int contactRate;
	double timeIntroducingInfection;
	double timeMHCDownregulation;
	double timeDecoy;
	string secondVirus;

	double deltaVirus;
	double timeInfection;
	double downregulationRate;
	double decoyRate;

	double transmissionRateAcute;
	double transmissionRateChronic;

	vector<int> shuffledHosts;
	vector<int> virtualHosts;

    vector<double> deathRates;
    vector<double> birthRates;

	/*
	int acute_infected;
	int chronic_infected;
	int immune;
	int downregulating;
	int decoy;
	int wildtype;/*/

	int downregulating;
	int downregulating_immune;
	int decoy;
	int decoy_immune;
	int wildtype;
	int wildtype_immune;

	int simpleInfection;
	int doubleInfection;
	int tripleInfection;
	
	int maxNumberOfInfectionsPerHost;
	
	bool isFileOpen;
	ofstream populationSize;
	fstream genesFile;
	fstream parameterFile;
	fstream backupFile;
	fstream notBornChildren;
	fstream dyingHosts;
	fstream mapFile;
};

#endif /* WORLD_H_ */

/*
 * new member variable: sizeMHCpool -> it's initialized in LoadParameterFile()
 */

