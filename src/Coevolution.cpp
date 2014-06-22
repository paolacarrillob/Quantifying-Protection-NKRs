//============================================================================
// Name        : CoevolutionKIR.cpp
// Author      : P.Carrillo-Bustamante
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "World.h"
#include <string.h>
#include "Exceptions.h"
using namespace std;

int main(int argc, char*argv[])
{
	
	/*string parameterFile(argv[1]);
	World oneWorld;
	oneWorld.LoadParameterFile(parameterFile);
	oneWorld.Initialize();

	map< int, pair <int, int> > ::iterator it;
	for(it = oneWorld.KIRGenesMap.GetMap().begin(); it != oneWorld.KIRGenesMap.GetMap().end(); it ++)
	{
		cout << (*it).first << "|" << (*it).second.first <<"|"<< (*it).second.second <<endl;

	}
	
	Host dummy;
	dummy.Copy(oneWorld.hosts.at(0));
	dummy.InfectWith(oneWorld.downregulatingVirus, 0.1);	
	//dummy.InfectWith(oneWorld.decoyVirus, 0.1);
	cout << "virus type: "<<oneWorld.nastyVirus.GetVirusType()<<endl;
	cout << "virus type in host: " << dummy.infections.size() <<endl;
	//oneWorld.decoyVirus.PrintParametersVirus();
	//cout <<"printing dummy infection:\n";
	//stupidvirus.PrintParametersVirus();
	
	list<Infection>::iterator it;
	for(it = dummy.infections.begin(); it != dummy.infections.end(); it++)
	{
		cout << "...." <<endl;
		it->PrintParameters();
		dummy.ClearInfection(2.0, (*it));
		it->PrintParameters();
		it->SetInfectionType(2.1);
		it->PrintParameters();
	}
	
*/	
	// To test: copy function of host: is it working??? TransmitInfection of host; Infect event; Escape event.... right now I cannot think of anything else that might be causing the mistake!!!


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
	//sprintf(buffer, "rm *.txt");
	sprintf(buffer, "rm *.log");
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
	try
	{
		theWorld.LoadParameterFile(parameterFile);
	}
	catch (OussException& e ) {
		cout << e.GetErrorData() << endl;
	}
	catch (...)
	{
		cout << "unknown exception thrown" <<endl;
	}
	
	//theWorld.LoadParameterFile(parameterFile);
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
