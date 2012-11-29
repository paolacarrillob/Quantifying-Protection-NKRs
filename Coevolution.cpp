//============================================================================
// Name        : CoevolutionKIR.cpp
// Author      : P.Carrillo-Bustamante
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "World.h"
using namespace std;

int main(int argc, char*argv[])
{
	/* string parameterFile(argv[1]);
	World oneWorld;
	oneWorld.LoadParameterFile(parameterFile);
	oneWorld.Initialize();

	map< int, pair <int, int> > ::iterator it;
	for(it = oneWorld.KIRGenesMap.GetMap().begin(); it != oneWorld.KIRGenesMap.GetMap().end(); it ++)
	{
		cout << (*it).first << "|" << (*it).second.first <<"|"<< (*it).second.second <<endl;

	}

	int randomHost = RandomNumber(0, oneWorld.hosts.size()-1);
	Host dummy;
	dummy.Copy(oneWorld.hosts.at(0));
	cout << "educating dummy................."<<endl;
	dummy.EducateKIRs();
	dummy.CountFunctionalKIRs();
	int y = dummy.CountExpressedKIRs();
	//cout << "number of FunctionalKIRs: " << x <<endl;
	cout << "number of KIRs: " << y <<endl;
	exit(0);
	dummy.InfectWith(oneWorld.decoyVirus,0.0);
	cout << "virus type: "<<oneWorld.decoyVirus.GetVirusType()<<endl;
	dummy.ClearInfection(1.0);
	cout << "virus type: "<<oneWorld.decoyVirus.GetVirusType()<<endl;
	//oneWorld.Escape(0);
	dummy.ClearInfection(1.0);

	/*
	for(unsigned int i = 0; i< dummy.kirGenes.size(); i++)
	{
		cout << dummy.kirGenes.at(i).GetGeneID()<<"_"<<dummy.kirGenes.at(i).GetGeneSpecificity() << " ";
	}

	cout << endl;

	/*for(unsigned int i=0; i< dummy.kirGenes.size(); i++)
	{
		if(RandomNumberDouble() <0.8)
		{
		dummy.kirGenes.at(i).Mutate();
		}
		cout << dummy.kirGenes.at(i).GetGeneID()<<"_"<<dummy.kirGenes.at(i).GetGeneSpecificity() << " ";

	}
	cout <<endl;

	KIRGene mygene(5);
	cout << mygene.GetGeneID()<< " | " << mygene.GetGeneSpecificity() << " | ";
	mygene.PrintBits();
	mygene.PointMutation();
	cout << mygene.GetGeneID()<< " | " << mygene.GetGeneSpecificity()<< " | ";
	mygene.PrintBits();

	exit(0);

	KIRGene bla(3);
	KIRGene blub(6);

	if(!(bla==blub))
	{
		bla.PrintBits();
		blub.PrintBits();
	}


	KIRGene bla;
	cout << bla.GetGeneID()<< " | " << bla.GetGeneSpecificity()<< " | "<< bla.GetGeneType() << " | " << bla.GetGeneMid() << " | ";
	bla.PrintBits();
	KIRGene blub(5);
	bla.Copy(blub);
	cout << bla.GetGeneID()<< " | " << bla.GetGeneSpecificity()<< " | "<< bla.GetGeneType() << " | "<< bla.GetGeneMid() << " | ";
	bla.PrintBits();
	bla.MutateReceptorType();
	cout << bla.GetGeneID()<< " | " << bla.GetGeneSpecificity()<< " | "<< bla.GetGeneType() << " | "<< bla.GetGeneMid() << " | ";
	bla.PrintBits();



	exit(0);

*/

	//ACTUAL SIMULATION STARTING!
	cerr << "Testing if stderr is redirected to stdoutput" << endl;
	if (argc<2)
	{
		cerr << "Usage: "<<argv[0] << " <Parameter file> <-b Loading Backup file> \n CoevolutionKIR simulates the evolution of the complex KIR system. For correct usage you have to indicate a parameter file and a backup file name. If you are not loading any backup, just give any random name"<< endl;
		exit(-1);
	}
	// delete the existing host_file
	char buffer[512];
	//Ouss: Use boost library here for boost:filesystem to delete files... etc
	cout <<"deleting old files ..."<<endl;
	sprintf(buffer, "rm *.txt");
	system(buffer);

	string parameterFile(argv[1]);
	string backupFile;
	bool loadingBackup=false;
	for(int i=2; i< argc;i++)
	{
		if(strcmp(argv[i],"-b")==0 && argc > i)
		{
			loadingBackup = true;
			backupFile.assign(argv[++i]);
		}
	}

	World theWorld;
	theWorld.LoadParameterFile(parameterFile);
	cout <<"welcome, hello" << endl;
	if(!loadingBackup)
	{
		// initialize host population
		theWorld.Initialize();
	}
	else
	{
		cout << "\n Loading data from backup file: "  << backupFile << endl;
		theWorld.LoadBackupFile(backupFile);
	}
	theWorld.Simulate();
	cout << "bye bye \n";

	return 0;
}


/*
 * 22.08.2012: allow for evolution of specificity of the system:
 *
 * 1. start by creating the KIR pool, with a constraint of a high specificity, e.g., L = 12
 * 2. create "int specificity" as a member variable of the class Gene
 * 3. change everything accordingly: print out functions, copy functions, arguments of the functions
 * 4. Create function Mutate host
 * 4. in every new host, let during mutation, change L to L+1.
 *
 * 23.08.2012: NO, talked to Rob, he had another idea (of course)
 * 1. create RANDOM KIR molecules. Each random KIR has a specificity.
 * 2. mutate by randomly generating a new bitstring
 *
 * -> change initialization!
 * -> Host constructor: upon mutation create a random Bitstring with a random specificity (use BitstringMHC -> othwerwise it will create a KIR that fits one MHC in the pop)
 * -> upon creation do: specificity = randomNumber(1,16)
 *
 * 06.09.2012
 * all done, check other *.h files to see the latest changes!
 *
 * 05.10.2012
 * new changes in GenePool and World:Initialize -> check other *h files
 *
 *12.10.2012
 *new changes in Gene... new functions ->check other *h files
 *
 *18.10.2012
 *changed Host constructors to implement new types of mutations-> check other *h files
 *
 *24.10.2012 fixed bug in LoadBackupFile!
 *
 */
