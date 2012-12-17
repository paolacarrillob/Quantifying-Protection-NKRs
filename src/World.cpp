/*
 * World.cpp
 *
 *  Created on: Apr 26, 2011
 *      Author: paola
 */

#include "World.h"
#include "backup.h"

World::World() {
}

void World :: LoadParameterFile(const string& fileName)
{
	double t_outfile;
	double t_backup;
	double t_popFile;
	bool onlyAcuteInfection;
	kaBackup paramFile(fileName,false);
	paramFile.Load("Maximal population size", &maxHostPop);
	paramFile.Load("Initial population size", &initHostPop);
	paramFile.Load("Time Step", &timeStep);
	paramFile.Load("End of Simulation", &timeEnd);
	paramFile.Load("Time introducing infection", &timeIntroducingInfection);
	paramFile.Load("Time outfile", &t_outfile);
	//paramFile.Load("Time recording", &timeRecording);
	paramFile.Load("Time backup", &t_backup);
	paramFile.Load("Population Time", &t_popFile);
	paramFile.Load("Mutation rate host", &mutationRate);
	paramFile.Load("Mutation type host", &mutationTypeHost);
	paramFile.Load("Contacts per week", &contactRate);
	paramFile.Load("KIR Loci", &KIRLoci);
	paramFile.Load("MHC Loci", &MHCLoci);
	paramFile.Load("HLA-C alleles distribution", &HLA_C);
	paramFile.Load("Number HLA-C alleles",&sizeMHCPool);
	paramFile.Load("MHC-KIR specificity", &KIRspecificity);
	paramFile.Load("Tuning", &education);
	paramFile.Load("Extra Number KIRs to be expressed", &expressionExtraKIRs);
	paramFile.Load("Viral load", &deltaVirus);
	paramFile.Load("Life time of the virus", &timeInfection);
	paramFile.Load("MHC dowregulation rate", &downregulationRate);
	paramFile.Load("Decoy rate", &decoyRate);
	paramFile.Load("Transmission rate acute infection", &transmissionRateAcute);
	paramFile.Load("Transmission rate chronic infection", &transmissionRateChronic);
	paramFile.Load("Type of mutation", &mutationType);
	paramFile.Load("Type of infection", &onlyAcuteInfection);

	const string populationFile("PopulationSize.txt");	//in the future if you would like to...
														//maybe you would think of reading the name of the output file from the parameters file...
														//
														//char populationFileChars[256];
														//paramFile.LoadString("Population Size", populationFileChars, 256);
														//string populationFile(populationFileChars);

	populationSize.open(populationFile.c_str(), ios::out);
	populationSize <<"#time\tpopSize\tbabies\tdead\tacute\tchronic\timmune\tmhc_down\tdecoy\n";


	timeStep= timeStep*WEEK;
	simulationTime = 0.0;
	timeEnd = timeEnd*YEAR;
	//timeRecording = timeRecording * YEAR;
	timeIntroducingInfection = timeIntroducingInfection*YEAR;
	timeMHCDownregulation = timeIntroducingInfection + 5000*YEAR;
/*	if(onlyAcuteInfection)
		timeMHCDownregulation = timeIntroducingInfection + 700*YEAR;
	else
		timeMHCDownregulation = timeIntroducingInfection + 200*YEAR;	 //an acute infection needs longer time
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 //to reach steady state, that is why MHC downregulation
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 //should come even later
*/
	timeDecoy = timeMHCDownregulation + 100000*YEAR;
	outfileRate = 1.0/(t_outfile*YEAR);
	backupRate = 1.0/(t_backup*YEAR);
	populationSizeRate = 1.0/(t_popFile*YEAR);

	birthRate = 0.5*timeStep/YEAR; // every birth event will happen once every four years
	deathRate = timeStep/YEAR; // every death event will happen once every year
	infectionRate = 52.0*timeStep/YEAR; //every infection event will happen every week (pop concerts removed!)
	escapeRate = 0.0001*timeStep/YEAR;

	nastyVirus.SetViralParameters(downregulationRate, decoyRate,deltaVirus,timeInfection, 0, onlyAcuteInfection);
	downregulatingVirus.SetViralParameters(downregulationRate, decoyRate,deltaVirus,timeInfection, 1, onlyAcuteInfection);
	decoyVirus.SetViralParameters(downregulationRate, decoyRate, deltaVirus, timeInfection, 2, onlyAcuteInfection);

	acute_infected = 0;
	chronic_infected = 0;
	immune = 0;

	cout << "Beginning simulation with: \n";
	cout << "End Simulation after " << timeEnd/YEAR << " years \n";
	cout << "Maximal population size " << maxHostPop << " individuals \n";
	cout << "Initial population size " << initHostPop << " individuals \n";
	cout << "Introducing virus after " << timeIntroducingInfection/YEAR << " years \n";
	cout << "Mutation rate host " << mutationRate << "\n";
	cout << "Mutation type host " << mutationTypeHost << "\n";
	cout << "KIR loci " << KIRLoci << "\n";
	cout << "MHC loci " << MHCLoci << "\n";
	cout << "Contacts per week: " << contactRate << "\n";
	cout << "HLA-C alleles distribution " << HLA_C<<"\n";
	cout << "Number of HLA alleles " << sizeMHCPool<<"\n";
	cout << "MHC-KIR specificity " << KIRspecificity<<"\n"; //this is only for the intial pool, basically, the maximum specificity is set with this
	cout << "Tuning "<< education << "\n";
	cout << "Expression extra KIRs "<< expressionExtraKIRs << "\n";
	cout << "Viral load: " << deltaVirus << "\n";
	cout << "Life time virus: " << timeInfection << "\n";
	cout << "MHC dowregulation rate " << downregulationRate <<"\n";
	cout << "Decoy rate "<< decoyRate << "\n";
	cout << "Tranmission rate acute infection "<< transmissionRateAcute << "\n";
}

// initialize host population
bool World::Initialize()
{
	//initialize MHcPool
	MHCPool.FillMHCGenePool(sizeMHCPool);

	//initialize the map of KIR with the first 10 genes to initilaize the population
	for (int i =  0; i <KIRLoci; i++)
	{
		KIRGene dummy_gene(KIRspecificity);
		if (!KIRGenesMap.IsGeneInMap(dummy_gene))
			KIRGenesMap.FillMap(MHCPool, dummy_gene);
	}

	//initialize the population with the MHC genes of the pools
	for (unsigned int i = 0; i< initHostPop; i++)
	{
		Host dummyhost(KIRLoci, MHCLoci, mutationRate, education, expressionExtraKIRs, KIRGenesMap, MHCPool, HLA_C);
		hosts.push_back(dummyhost);
	}

	CreateBirthAndDeathRates();
	return true;
}

/* This function creates a table upon initialization that calculates the birth/death rates for all possible ages*/
void World ::CreateBirthAndDeathRates()
{
	for(int i = 1; i <150; i ++)
	{
		double birth = -1/(1+ exp(i-20)) + 1/(1+exp(i-45));
		birthRates.push_back(birth);

		double death = exp(0.1*i-10.5)+ exp(-0.4*i-8);
		deathRates.push_back(death);
	}

}


/*EVENT functions*/
/*Birth function: creates a child with the haplotype of one parent and another randomly chosen host*/
bool World::Birth(int index, Host& baby_host)
{

	//check if host's age allows him to become a parent
	double ageDependentBirth = hosts.at(index).GetAgeDependentBirthRate(birthRates);
	//cout << "do i get stuck here??? birth event!"<< endl;
	if(RandomNumberDouble()< birthRate*ageDependentBirth*(1-(hosts.size()/(maxHostPop*0.99753))))
	{
		int randomindex = RandomNumber(0,shuffledHosts.size()-1);
		//check if the potential parent is himself
		while(randomindex == index)
		{
			randomindex = RandomNumber(0,shuffledHosts.size()-1);
		}

		//choose randomly which parent is going to "donate" his/her mhc molecule
		int parent;
		if((RandomNumberDouble()<0.5))
			parent = index;
		else
			parent = randomindex;
		Host testHost(KIRLoci,MHCLoci, hosts.at(parent).mhcGenes ,MHCPool, HLA_C,/*KIRPool, */hosts.at(index).kirGenes,hosts.at(randomindex).kirGenes, mutationRate,education,expressionExtraKIRs, KIRGenesMap, mutationTypeHost);
		baby_host.Copy(testHost);
		return true;
	}
	return false;
}

/*Death function: according to an age-dependent rate, a host will be removed from the population*/
bool World::Death(int index)
{
	//cout << "do i get stuck here??? death event!"<< endl;
	double intrinsicDeath;
	intrinsicDeath = deathRate*hosts.at(index).GetIntrinsicDeathRate(deathRates);
	if(RandomNumberDouble()<intrinsicDeath)
	{
		return true;
	}
	return false;
}
/*Infect function: upon contact between two hosts, virus can spread*/

void World::Infect(int index)
{
	//cout << "do i get stuck here??? infect event!"<< endl;
	if(RandomNumberDouble()<infectionRate)
	{
		for(int i=0; i<contactRate; i++)
		{
			//pick random partner that is NOT yourself!
			int randomindex = RandomNumber(0,shuffledHosts.size()-1);
			while (randomindex == index)
				randomindex=RandomNumber(0,shuffledHosts.size()-1);
			//check whether partner is infectious

			//if(hosts.at(randomindex).IsSusceptible() || hosts.at(randomindex).IsImmune())
			if(!hosts.at(randomindex).IsInfected())
				continue;
			//transmit the virus according to the type of infection
			else
			{
				int infectionState = hosts.at(randomindex).GetMainInfectionType();
				switch (infectionState)
				{
				case 1:
				{
					if(RandomNumberDouble()<transmissionRateAcute)
					{
						//get the acute virus
						hosts.at(index).InfectWith(hosts.at(randomindex).GetAcuteInfection(), simulationTime);
					}

				}break;
				case 2:
				{
					if(RandomNumberDouble()<transmissionRateChronic)
					{
						// get the chronic virus
						hosts.at(index).InfectWith(hosts.at(randomindex).GetChronicInfection(), simulationTime);
					}

				}break;
				}
			}
		}
	}
}

void World::Escape(int index)
{
	//cout << "do i get stuck here??? escape event!"<< endl;
	list<Infection>::iterator it;
	for(it = hosts.at(index).infections.begin(); it!= hosts.at(index).infections.end(); it++)
	{
		if(mutationType == 0) //allow for only MHC downregulation to evolve
		{
			if(simulationTime>=timeMHCDownregulation)
			{
				it->pathogen.DownregulateMHC();
			}
		}
		
		if(mutationType == 1) //allow for evolution of directly a decoy virus
		{
			if(simulationTime>=timeMHCDownregulation)
			{
				int hap=RandomNumberDouble()<0.5;
				int mhcID = hosts.at(index).mhcGenes.at(hap).GetGeneID(); //this allows new decoy molecules to arise
				it->pathogen.BuildDecoy(mhcID);
			}
		}
		
		if(mutationType == 2) //allow for first MHC downregulation to evolve and then decoy
		{
			if(simulationTime>=timeMHCDownregulation)
			{
				it->pathogen.DownregulateMHC();
			}
			
			if(simulationTime >= timeDecoy)
			{
				int hap=RandomNumberDouble()<0.5;
				int mhcID = hosts.at(index).mhcGenes.at(hap).GetGeneID(); //this allows new decoy molecules to arise
				it->pathogen.BuildDecoy(mhcID);
			}
		}
		
	}	
}

/*SIMULATION functions*/

void World::ShuffleHosts()
{
	//virtualHosts.clear();
	shuffledHosts.clear();

	for(unsigned int i=0; i<hosts.size(); i++)
	{
		shuffledHosts.push_back(i);
	}

   for(unsigned int i = 0 ; i < hosts.size() ; i ++ )
   {
      int j = RandomNumber(i,hosts.size()-1);
      // swap a[i] and a[j]
      int t = shuffledHosts.at(j);
      shuffledHosts.at(j) = shuffledHosts.at(i);
      shuffledHosts.at(i) = t;
   }

   /*
	while(virtualHosts.size()!=0)
	{
		int randomIndex = RandomNumber(0,virtualHosts.size()-1);
		shuffledHosts.push_back(virtualHosts.at(randomIndex));
		vector<int>:: iterator it = virtualHosts.begin();
		virtualHosts.erase(it + randomIndex);
	}*/
}

void World::Simulate()
{
	double lastPopulationOutfileTime = 0.0;
	double lastOutfileTime = 0.0;
	//double lastStopOutfileTime = timeRecording;
	double lastBackupTime = 0.0;
	double lastAcuteInfectionTime = 0.0;

	isFileOpen = false;
	SaveMap();

	cout << "time Decoy:  "<< timeDecoy << "\n";
	cout << "simulation Time: "<<simulationTime <<"\n";

	while(simulationTime <=timeEnd)
	{
		// printing out the backup files
		if(floor((simulationTime-lastBackupTime)*backupRate)>0)
		{
			cout <<"\tSaving Backup\n";
			SaveBackupFile();
			SaveMap();
			lastBackupTime = simulationTime;
		}


		int number_babies = 0;
		int number_dead_people = 0;
		//cout<< "Time: " << simulationTime/YEAR <<endl;

		ShuffleHosts();
		vector<int>::iterator shuffledHostsit;

		//introduce infection
		//if(floor(timeIntroducingInfection -simulationTime) > 0 && floor(timeIntroducingInfection -simulationTime) < 3.0*WEEK)
			IntroduceVirus();

		for(shuffledHostsit = shuffledHosts.begin(); shuffledHostsit!=shuffledHosts.end(); shuffledHostsit++)
		{
			//let event happen for every random host
			int index = *shuffledHostsit;
			Host babyHost;
			if(Birth(index, babyHost))
			{
				hosts.push_back(babyHost);
				number_babies ++;
			}
			//if(hosts.at(index).IsSusceptible())
			Infect(index);
			//if(hosts.at(index).IsChronicInfected()|| hosts.at(index).IsAcuteInfected())
			Escape(index);
			if(Death(index))
			{
				hosts.at(index).SetDead();
				number_dead_people++;
			}
			
			//clear the infection
			list<Infection>::iterator inf;
			for(inf = hosts.at(index).infections.begin(); inf!=hosts.at(index).infections.end(); inf++) //check for every infection within one host, whether it's time to be cleared
			{
				//Infection nastyInfection = (*inf);
				if(inf->IsAcute())
				{
					if((simulationTime - inf->GetInfectionTime()) == (1.0 + 4.0 *timeInfection)*WEEK)
						hosts.at(index).ClearInfection(simulationTime,(*inf));					
				}

			}
			/*if(hosts.at(index).IsAcuteInfected())
			{

				if((simulationTime -hosts.at(index).GetInfectionTime()) == (1.0 + 4.0*timeInfection)*WEEK)
				{
					hosts.at(index).ClearInfection(simulationTime);
				}
			}
*/
			hosts.at(index).UpdateParameters(timeStep,simulationTime);
			 
		}
		//record those who are dying!
	/*	if(SaveAgeDyingHosts(lastOutfileTime, lastStopOutfileTime))
			lastStopOutfileTime = simulationTime;*/

		RemoveDeadHosts();

		//printing out the population size not every week but with a rate (otherwise the file will turn huge!!!!)
		bool timeToPrintPopulationSize = floor((simulationTime-lastPopulationOutfileTime)*populationSizeRate)>0;
		if(simulationTime == 0.0 || timeToPrintPopulationSize)
		{
			SavePopulationSize(number_babies,number_dead_people);
			lastPopulationOutfileTime = simulationTime;
		}

		//cout<< "Time: " << simulationTime/YEAR <<endl;
		// printing out the gene files
		bool timeToPrintOut = floor((simulationTime-lastOutfileTime)*outfileRate)>0;
		if(simulationTime == 0.0 || timeToPrintOut)
		{
			cout <<"\tPrinting parameters\n";
			cout<< "Time: " << simulationTime/YEAR <<endl;
			SaveGenes();
			SaveParameters();
			populationSize.flush();
			if(populationSize.bad())
				cout << "something bad happened \n";
			lastOutfileTime = simulationTime;
		}

		//re introducing acute infection
		if(simulationTime>timeIntroducingInfection && nastyVirus.IsOnlyAcute())
		{
			if(floor((simulationTime -lastAcuteInfectionTime)*(1.0/(1.0*YEAR)))>0)
			{
				cout << "reintroducing the infection"<<endl;
				IntroduceVirus();
				lastAcuteInfectionTime = simulationTime;
			}
		}

		simulationTime+=timeStep;
	}
	populationSize.close();
	SaveMap();

	cout << "The simulation was carried out with: \n";
	cout << "End Simulation after " << timeEnd/YEAR << " years \n";
	cout << "Maximal population size " << maxHostPop << " individuals \n";
	cout << "Initial population size " << initHostPop << " individuals \n";
	cout << "Introducing virus after " << timeIntroducingInfection/YEAR << " years \n";
	cout << "Mutation rate host " << mutationRate << "\n";
	cout << "Mutation type host " << mutationTypeHost << "\n";
	cout << "Contacts per week: " << contactRate << "\n";
	cout << "KIR Loci: " << KIRLoci << "\n";
	cout << "HLA-C alleles distribution " << HLA_C<<"\n";
	cout << "Number of HLA alleles " << sizeMHCPool<<"\n";
	cout << "Initial MHC-KIR specificity " << KIRspecificity<<"\n"; //this is only for the initial pool, basically, the maximum specificity is set with this
	cout << "Tuning "<< education << "\n";
	cout << "Expression extra KIRs "<< expressionExtraKIRs << "\n";
	cout << "Viral load: " << deltaVirus << "\n";
	cout << "Life time virus: " << timeInfection << "\n";
	cout << "MHC dowregulation rate " << downregulationRate <<"\n";
	cout << "Decoy rate "<< decoyRate << "\n";
	cout << "Tranmission rate acute infection "<< transmissionRateAcute << "\n";
}

void World ::IntroduceVirus()
{
	if(floor(timeIntroducingInfection -simulationTime) > 0 && floor(timeIntroducingInfection -simulationTime) < 3.0*WEEK)
	{
		cout <<"\t Introducing the infection"<<endl;
		for(int i= 0; i<0.05* hosts.size(); i++)
		{
			int randomindex = RandomNumber(0,hosts.size()-1);
			/*while(!hosts.at(randomindex).IsSusceptible())
			{
				randomindex = RandomNumber(0,hosts.size()-1);
			}*/
			hosts.at(randomindex).InfectWith(nastyVirus, simulationTime);
		}
	}



	/*
		//introducing mutant virus: downregulation
	if(floor(timeDecoy -simulationTime) > 0 && floor(timeDecoy -simulationTime) < 3.0*WEEK)
	{
		cout <<"\t Introducing mutant viruses: MHC Downregulation"<<endl;
		for(int i= 0; i<0.5* hosts.size(); i++)
		{
			int randomindex = RandomNumber(0,hosts.size()-1);
			hosts.at(randomindex).InfectWith(downregulatingVirus, simulationTime);
		}
	}


		 * Manually introduce the mutant viruses
		 *

	//introducing mutant virus: downregulation+decoy
	if(floor(timeDecoy -simulationTime) > 0 && floor(timeDecoy -simulationTime) < 3.0*WEEK)
	{
		int randomindex = RandomNumber(0,hosts.size()-1);
		int hap=RandomNumberDouble(0,1)<0.5;
		//int mhcID = MHCPool.GetGenes().at(5);
		int mhcID = hosts.at(randomindex).mhcGenes.at(hap).GetGeneID();
		cout <<"\t Introducing mutant viruses: MHC Downregulation + Decoys"<<endl;
		for(int i= 0; i<0.95*hosts.size(); i++)
		{
			int randomindex2 = RandomNumber(0,hosts.size()-1);
			decoyVirus.mhcDecoy.SetGeneID(mhcID);
			hosts.at(randomindex2).InfectWith(decoyVirus, simulationTime);
		}
	}*/
}

void World::RemoveDeadHosts()
{
	//deleting the hosts
	acute_infected = 0;
	chronic_infected = 0;
	immune = 0;
	downregulating = 0;
	decoy = 0;
	int pos = 0;
	int max = hosts.size();
	while(pos < max)
	{
		if(hosts.at(pos).IsDead())
		{

			//remove the host
			hosts.erase(hosts.begin() + pos);
			max = hosts.size();
		}
		else
		{
			list<Infection>:: iterator it;
			for(it = hosts.at(pos).infections.begin(); it!= hosts.at(pos).infections.end(); it++)
			{
				int inf_type = it->GetInfectionType();
				switch(inf_type)
				{
					case 1:{acute_infected++;}break;
					case 2:{chronic_infected++;}break;
					case 3:{immune++;}break;				
				}
				int virus_type = it->pathogen.GetVirusType();
				switch(virus_type)
				{
					case 1:{downregulating++;}break;
					case 2:{decoy++;}break;
				}
			
			}
			/*
			if(hosts.at(pos).IsAcuteInfected())
				acute_infected++;
			if(hosts.at(pos).IsChronicInfected())
				chronic_infected++;
			if(hosts.at(pos).IsImmune())
				immune ++;
			if(hosts.at(pos).pathogen.IsDownregulatingMHC())
				downregulating ++;
			if(hosts.at(pos).pathogen.IsStealingDecoy())
				decoy ++;*/
		}
		pos ++;
	}
}

/*Functions which do the OUTPUT files*/

/*this function saves the genes of each host -> to keep track of MHC and KIR diversity*/
void World:: SaveGenes()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Genes.txt";
	
	genesFile.open(s.c_str(), ios::out);

	genesFile << "#Host ID\t";
	for(unsigned int i = 0; i< hosts.at(0).mhcGenes.size(); i++)
	{
		genesFile << "Mhc"<< i+1 <<"\t";
	}

	for(unsigned int i = 0; i< hosts.at(0).kirGenes.size(); i++)
	{
		genesFile << "Kir"<< i+1 <<"\t";
	}

	genesFile << "\tfunctional Kirs\texpressed Kirs\n";
	int index = 0;
	vector<Host>::iterator hostIt;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
	{
		genesFile << index << "\t";
		hostIt->SaveGenes(genesFile);
		index ++;
	}

	genesFile.close();
}

/*This function keeps track (and saves) of the population size*/
void World::SavePopulationSize(int babies,int dead_people)
{
	populationSize <<simulationTime/YEAR <<"\t"<< hosts.size() <<"\t"<<babies <<"\t"<<dead_people<<"\t"<<acute_infected << "\t"<< chronic_infected << "\t"<< immune << "\t"<<downregulating<< "\t"<< decoy <<"\n";
}

/*This function keeps track (ans saves) several parameters of the host and virus*/
void World::SaveParameters()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".txt";
	parameterFile.open(s.c_str(), ios::out);

	parameterFile << "#HostID\t" << "Age \t" << "Infection Time \t" << "Infection type \t" << "Virus type\t" << "viral load \t"<< "decoy ID\n";
	vector<Host>::iterator hostIt;
	int index = 0;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
	{
		parameterFile << index <<"\t";
		hostIt->SaveParameters(parameterFile);
		index ++;
	}

	parameterFile.close();
}

void World::SaveMap()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Map.txt";
	mapFile.open(s.c_str(), ios::out);

	map< int, pair <int, int> > ::iterator it;
	for(it = KIRGenesMap.GetMap().begin(); it != KIRGenesMap.GetMap().end(); it ++)
	{
		mapFile << (*it).first << "|" << (*it).second.first <<"|"<< (*it).second.second <<endl;
	}
	mapFile.close();
}

bool World::SaveAgeDyingHosts(double lastOutfileTime, double lastStopOutfileTime)
{

	bool timeToPrintOut = floor((simulationTime-lastOutfileTime)*outfileRate);
	bool timeToStopOut = floor((simulationTime-lastStopOutfileTime)*outfileRate);

	if(timeToPrintOut>0)
	{
		isFileOpen = true;
		stringstream ss;
		ss << simulationTime/YEAR;
		string s(ss.str());
		s+=".Age.txt";
		dyingHosts.open(s.c_str(), ios::out);
		dyingHosts << "#Age \t KIR genes\t virus type\t viral load\t decoy ID\t onlyAcute?\t infection type\t infection time\t clearance time\n";
	}

	if(isFileOpen)
	{
		vector<Host>::iterator hostIt;
		//int index = 0;
		for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
		{
			if (hostIt->IsDead())
			{
				//write in a file their age, gene content, and infection type
				hostIt->SaveAgeDyingHost(dyingHosts);
				dyingHosts << simulationTime/YEAR <<"\n";
			}
			//index++;
		}
		if(timeToStopOut>0)
		{
			isFileOpen = false;
			dyingHosts.close();
			return true;
		}
		else
			return false;
	}
	return false;
}


/*Functions which do the BACKUP files*/

void World::SaveBackupFile()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Backup.data";
	backupFile.open(s.c_str(), ios::out);

	if(!backupFile.is_open())
	{
		throw OussException("Backup: Could not Open File");
	}
	else
	{
		cout <<"Saving backup file\n"<<endl;
	}

	backupFile << simulationTime<< "\t "<< MHCPool.GetPoolSize() << "\t";
	//save mhc gene pool
	for(unsigned int i = 0; i<MHCPool.GetPoolSize(); i++)
	{
		backupFile << MHCPool.GetGenes().at(i) << "\t";
	}
	//save kir gene pool
	/*for(unsigned int i = 0; i<KIRPool.GetPoolSize(); i++)
	{
		backupFile << KIRPool.GetGenes().at(i) << "\t";
	}*/
	backupFile << "\n";

	//save the KIRMap
	KIRGenesMap.SaveBackupMap(backupFile);
	backupFile << "\n";
	vector<Host>::iterator hostIt;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
	{
		hostIt->SaveBackupHost(backupFile);
	}
	backupFile.close();
}

void World ::LoadBackupFile(const string& fileName)
{
	CreateBirthAndDeathRates();
	int l = 0;
	int poolSize = 0;
	int mapSize = 0;
	string sline;
//	string slineIt;
//	string line;
	vector<Host>::iterator hostIt;
	backupFile.open(fileName.c_str(), ios::in);
	if(!backupFile.is_open())
	{
		throw OussException("Backup: Could not Open File");
	}
	while(!backupFile.eof())
	{
		if(l==0)
		{
			getline(backupFile,sline);
			stringstream ssline(sline);
			ssline >>simulationTime;
			ssline >>poolSize;
			//restore gene pools
			for(int i = 0; i<poolSize; i++)
			{
				int mhc_gene;
				ssline>>mhc_gene;
				MHCPool.GetGenes().push_back(mhc_gene);
			}
			/*for(int i = 0; i<kirPoolSize; i++)
			{
				int kir_gene;
				ssline>>kir_gene;
				KIRPool.GetGenes().push_back(kir_gene);
			}*/

			// get the second line with the information of the Map
			string secondline;
			string mline;
			getline(backupFile, secondline);
			stringstream msline (secondline);
			msline >> mapSize;
			for(int i= 0; i<mapSize; i++)
			{
				mline = KIRGenesMap.RestoreMap(msline);
				msline.str() = mline;
			}

		}

		else
		{
			string sline;
			getline(backupFile, sline);
			if(sline.size()!=0)
			{
				//cout <<"sline1:  "<<sline <<endl;
				Host tempHost;
				tempHost.RestoreHost(sline);
				//cout <<"restored:  "<<tempHost.GetAge() <<"\t" << tempHost.GetInfectionType() << "\t" << tempHost.GetInfectionTime() << "\t" <<tempHost.GetClearanceTime() <<endl;
				hosts.push_back(tempHost);

			}
		}
		l++;
	}
	backupFile.close();
}
