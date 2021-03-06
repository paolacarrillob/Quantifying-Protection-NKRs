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
	char buff[256];
	
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
	paramFile.Load("KIR type", &KIRGeneType);
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
	paramFile.Load("Maximal Number of Infections", &maxNumberOfInfectionsPerHost);
	paramFile.LoadString("Second Virus", buff, 256);
	secondVirus.assign(buff);
	
	timeStep= timeStep*WEEK;
	simulationTime = 0.0;
	timeEnd = timeEnd*YEAR;

	timeIntroducingInfection = timeIntroducingInfection*YEAR;
	timeMHCDownregulation = timeIntroducingInfection + 5000*YEAR;

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
	
	//	acute_infected = 0;
	//	chronic_infected = 0;
	//	immune = 0;
	downregulating = 0;
	downregulating_immune = 0;
	decoy = 0;
	decoy_immune = 0;
	wildtype = 0;
	wildtype_immune = 0;

	simpleInfection = 0;
	doubleInfection = 0;
	tripleInfection = 0;
	WriteInfo();
}

void World :: WriteInfo()
{
	cout << "Beginning simulation with: \n";
	cout << "End Simulation after " << timeEnd/YEAR << " years \n";
	cout << "Maximal population size " << maxHostPop << " individuals \n";
	cout << "Initial population size " << initHostPop << " individuals \n";
	cout << "Introducing virus after " << timeIntroducingInfection/YEAR << " years \n";
	cout << "Mutation rate host " << mutationRate << "\n";
	cout << "Mutation type host " << mutationTypeHost << "\n";
	cout << "KIR loci " << KIRLoci << "\n";
	cout << "KIR Gene Type" << KIRGeneType << "\n";
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
		dummy_gene.SetGeneType(KIRGeneType);
		if (!KIRGenesMap.IsGeneInMap(dummy_gene))
			KIRGenesMap.FillMap(MHCPool, dummy_gene);
	}

	//initialize the population with the MHC genes of the pools
	for (unsigned int i = 0; i< initHostPop; i++)
	{
		Host dummyhost(KIRLoci, MHCLoci, mutationRate, education, expressionExtraKIRs, KIRGenesMap, MHCPool, HLA_C, KIRGeneType);
		dummyhost.Set_Host_ID(i);
		hosts.push_back(dummyhost);
	}

	CreateBirthAndDeathRates();
	return true;
}

/* This function creates a table upon initialization that calculates the birth/death rates for all possible ages*/
void World ::CreateBirthAndDeathRates()
{
	try
	{
		for(unsigned long int j=0; j <150*52; j++)
		{
			double i = j/52.0;
			double birth = -1/(1+ exp(i-20)) + 1/(1+exp(i-45));
			birthRates.push_back(birth);

			double death = exp(0.1*i-10.5)+ exp(-0.4*i-8);
			deathRates.push_back(death);
			//cout << i << " "<< birth <<" "<< death << endl;
		}
	}
	catch (...) {
		cout << "caught something" << endl;
	}
}


/*EVENT functions*/
/*Birth function: creates a child with the haplotype of one parent and another randomly chosen host*/
bool World::Birth(int index,unsigned long int next_id)//, Host& baby_host)
{
	if(hosts.size() <=1)
		return false;
	//check if host's age allows him to become a parent
	//cout << "do i get stuck here??? birth event!"<< endl;
	double ageDependentBirth = 0.0;
	try
	{
		//ageDependentBirth = hosts.at(index).GetAgeDependentBirthRate(birthRates);
		ageDependentBirth = GetAgeDependentBirthRate(hosts.at(index).GetAge());
	}
	catch (...)
	{
		cout << "caught something!!!"<< endl;
		cout << "index:" <<index << endl;
		cout << "hosts:" <<hosts.size() << endl;		
	}
	//cout << "do i get stuck here??? birth event!"<< endl;
	if(RandomNumberDouble()< birthRate*ageDependentBirth*(1-(hosts.size()/(maxHostPop*0.99753))))
	{
		int randomindex = RandomNumber(0,shuffledHosts.size()-1);
		//cout << "do i get stuck here??? birth event!"<< endl;
		//check if the potential parent is himself
		while(randomindex == index)
		{
			randomindex = RandomNumber(0,shuffledHosts.size()-1);
		//	cout << "do i get stuck here??? birth event!"<< endl;
		}

		//choose randomly which parent is going to "donate" his/her mhc molecule
		int parent;
		if((RandomNumberDouble()<0.5))
			parent = index;
		else
			parent = randomindex;
		//cout << "do i get stuck here??? birth event!"<< endl;
		Host testHost(KIRLoci,MHCLoci, hosts.at(parent).mhcGenes ,MHCPool, HLA_C,/*KIRPool, */hosts.at(index).kirGenes,hosts.at(randomindex).kirGenes, mutationRate,education,expressionExtraKIRs, KIRGenesMap, mutationTypeHost, KIRGeneType);
		//cout << "do i get stuck here??? birth event!"<< endl;
		testHost.Set_Host_ID(next_id);
		hosts.push_back(testHost);
		//baby_host.Copy(testHost);
		//cout << "do i get stuck here??? birth event!"<< endl;
		return true;
	}
	return false;
}

/*Death function: according to an age-dependent rate, a host will be removed from the population*/
bool World::Death(int index)
{
	//cout << "do i get stuck here??? death event!"<< endl;
	double intrinsicDeath = 0;
	vector<Host>::iterator it = hosts.begin() + index;
	//double intrinsicDeathRate = it->GetIntrinsicDeathRate(deathRates);
	double intrinsicDeathRate = GetIntrinsicDeathRate(it->GetAge(),it->GetViralDeathRate());
	intrinsicDeath = deathRate*intrinsicDeathRate;
	double r = RandomNumberDouble();
	if(r<intrinsicDeath)
	{
		//cout << "r , intrinsicdeathrate, intrinsicdeath:" << r << " , " << intrinsicDeathRate << " , " << intrinsicDeath << (r<intrinsicDeath?" ====> will DIE":"") << endl; 		
		it->SetDead();
		return true;
	}
	return false;
}
/*Infect function: upon contact between two hosts, virus can spread*/

void World::Infect(int index)
{
	//cout << "do i get stuck here??? infect event!"<< endl;
	if(hosts.size() <=1)
		return;
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
				Virus stupidVirus;
				switch (infectionState)
				{
					case 1:
					{
						if(RandomNumberDouble()<transmissionRateAcute)
						{
							//get the acute virus
							stupidVirus.Copy(hosts.at(randomindex).GetAcuteInfection());
							hosts.at(index).InfectWith(stupidVirus, simulationTime, maxNumberOfInfectionsPerHost);
						}

					}break;
					case 2:
					{
						if(RandomNumberDouble()<transmissionRateChronic)
						{
							// get the chronic virus
							stupidVirus.Copy(hosts.at(randomindex).GetChronicInfection());
							hosts.at(index).InfectWith(stupidVirus, simulationTime, maxNumberOfInfectionsPerHost);
						}
					}break;
					//default: cout <<"ERROR!!!!!! shouldn't happen!" <<endl; exit(1);
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
	//populationSize << "blabla\n" <<endl;
	const string populationFile("PopulationSize.log");	//in the future if you would like to...
	//maybe you would think of reading the name of the output file from the parameters file...
	//
	//char populationFileChars[256];
	//paramFile.LoadString("Population Size", populationFileChars, 256);
	//string populationFile(populationFileChars);
	
	populationSize.open(populationFile.c_str());
	populationSize <<"#time\tpopSize\tbabies\tdead\tacute\tchronic\timmune\twt\tmhc_down\tdecoy\n";
	
	double lastPopulationOutfileTime = 0.0;
	double lastOutfileTime = 0.0;
	//double lastStopOutfileTime = timeRecording;
	double lastBackupTime = 0.0;
	double lastAcuteInfectionTime = 0.0;

	isFileOpen = false;
	SaveMap();
	//cout << "time Decoy:  "<< timeDecoy << "\n";
	cout << "simulation Time: "<<simulationTime <<"\n"<<endl;
	unsigned long int id_counter = initHostPop;
	while(simulationTime <=timeEnd)
	{
		// printing out the backup files
		if(floor((simulationTime-lastBackupTime)*backupRate)>0)
		{
			cout <<"\tSaving Backup\n"<<endl;
			SaveBackupFile();
			SaveMap();
			lastBackupTime = simulationTime;
		}
		
//		cout <<"\t1\n"<<endl;


		int number_babies = 0;
		int number_dead_people = 0;
		//cout<< "Time: " << simulationTime/YEAR <<endl;

		ShuffleHosts();
//		cout <<"\t2\n"<<endl;
		vector<int>::iterator shuffledHostsit;

		//introduce infection
		IntroduceVirus(secondVirus);
//		cout <<"\t3\n"<<endl;
		for(shuffledHostsit = shuffledHosts.begin(); shuffledHostsit!=shuffledHosts.end(); shuffledHostsit++)
		{
			//let event happen for every random host
			int index = *shuffledHostsit;
			Host babyHost;
			if(Birth(index,id_counter))
			{
				number_babies++;
				id_counter++;
			}

			Infect(index);

			Escape(index);

			if(Death(index))
			{
				number_dead_people++;
			}
			
			//clear the infection
			list<Infection>::iterator inf;
			for(inf = hosts.at(index).infections.begin(); inf!=hosts.at(index).infections.end(); inf++) //check for every infection within one host, whether it's time to be cleared
			{
				if(inf->IsAcute())
				{
					if((simulationTime - inf->GetInfectionTime()) == (1.0 + 4.0 *timeInfection)*WEEK)
						hosts.at(index).ClearInfection(simulationTime,(*inf));					
				}

			}
			hosts.at(index).UpdateParameters(timeStep,simulationTime);
		}
		RemoveDeadHosts_HappyNewYear();

		//printing out the population size not every week but with a rate (otherwise the file will turn huge!!!!)
		bool timeToPrintPopulationSize = floor((simulationTime-lastPopulationOutfileTime)*populationSizeRate)>0;
		if(simulationTime == 0.0 || timeToPrintPopulationSize)
		{
			//cout << "what's going on?????" <<endl;
			//SavePopulationSize(number_babies,number_dead_people);
			SavePopulationSize();
			lastPopulationOutfileTime = simulationTime;
		}
		if(hosts.size() == 0)
		{
			cout << "Oooops pao is angry and by the way host size is:" << hosts.size() << endl;
			break;
		}
		TrackInfectedIndividuals();

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
				IntroduceVirus(secondVirus);
				lastAcuteInfectionTime = simulationTime;
			}
		}

		simulationTime+=timeStep;
	}
	populationSize.close();
	SaveMap();
	WriteInfo();
}

void World ::IntroduceVirus(const string& secondVirus)
{
	string decoyVirusString("decoy");
	string mhcDownVirusString("mhc down");

	switch(maxNumberOfInfectionsPerHost)
	{
		case 1:
		{
			if(floor(timeIntroducingInfection -simulationTime) > 0 && floor(timeIntroducingInfection -simulationTime) < 3.0*WEEK)
			{
				cout <<"\t Introducing the infection"<<endl;
				for(int i= 0; i<0.05* hosts.size(); i++)
				{
					int randomindex = RandomNumber(0,hosts.size()-1);
					while(hosts.at(randomindex).IsInfected())
					 {
					 randomindex = RandomNumber(0,hosts.size()-1);
					 }
					hosts.at(randomindex).InfectWith(nastyVirus, simulationTime, maxNumberOfInfectionsPerHost);
				}
			}
			
		}break;
		case 2: 
		{
			if(floor(timeDecoy -simulationTime) > 0 && floor(timeDecoy -simulationTime) < 3.0*WEEK)
			{
				if(secondVirus.compare(mhcDownVirusString) == 0)
				{
					cout <<"\t Introducing mutant viruses: MHC Downregulation"<<endl;
					for(int i= 0; i<0.5* hosts.size(); i++)
					{
						int randomindex = RandomNumber(0,hosts.size()-1);
						hosts.at(randomindex).InfectWith(downregulatingVirus, simulationTime, maxNumberOfInfectionsPerHost);
					}
				}
				if(secondVirus.compare(decoyVirusString) == 0)
				{
					int randomindex = RandomNumber(0,hosts.size()-1);
					int hap=RandomNumberDouble(0,1)<0.5;
					int mhcID = hosts.at(randomindex).mhcGenes.at(hap).GetGeneID();
					cout <<"\t Introducing mutant viruses: MHC Downregulation + Decoys"<<endl;
					for(int i= 0; i<0.95*hosts.size(); i++)
					{
						int randomindex2 = RandomNumber(0,hosts.size()-1);
						decoyVirus.mhcDecoy.SetGeneID(mhcID);
						hosts.at(randomindex2).InfectWith(decoyVirus, simulationTime, maxNumberOfInfectionsPerHost);
					}
					
				}
				
			}			
		}break;
		case 3:
		{
			//introducing mutant virus: downregulation
			if(floor(timeDecoy -simulationTime) > 0 && floor(timeDecoy -simulationTime) < 3.0*WEEK)
			{
				cout <<"\t Introducing mutant viruses: MHC Downregulation"<<endl;
				for(int i= 0; i<0.5* hosts.size(); i++)
				{
					int randomindex = RandomNumber(0,hosts.size()-1);
					hosts.at(randomindex).InfectWith(downregulatingVirus, simulationTime, maxNumberOfInfectionsPerHost);
				}
			}
			
			//introducing mutant virus: downregulation+decoy
			if(floor(timeDecoy - simulationTime) > 0 && floor(timeDecoy -simulationTime) < 3.0*WEEK)
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
					hosts.at(randomindex2).InfectWith(decoyVirus, simulationTime, maxNumberOfInfectionsPerHost);
				}
			}			
		}break;
		default: cout <<"ERROR!!!!!! an impossible number if infections! shouldn't happen!" <<endl; exit(1);
	}
}

void World::TrackInfectedIndividuals()
{
	//	acute_infected = 0;
	//	chronic_infected = 0;
	//	immune = 0;
	downregulating = 0;
	downregulating_immune = 0;
	decoy = 0;
	decoy_immune = 0;
	wildtype = 0;
	wildtype_immune = 0;

	simpleInfection = 0;
	doubleInfection = 0;
	tripleInfection = 0;
	
	vector<Host>::iterator it_host = hosts.begin();
	while(it_host != hosts.end())
	{	//check what kind of infections the hosts have
		/*
		list<Infection>:: iterator it;
		for(it = it_host->infections.begin(); it!= it_host->infections.end(); it++)
		{
			int inf_type = it->GetInfectionType();
			switch(inf_type)
			{
				case 0:{acute_infected++;}break;
				case 1:{acute_infected++;}break;
				case 2:{chronic_infected++;}break;
				case 3:{immune++;}break;
					//default: cout <<"ERROR!!!!!! COUNTING OTHER INFECTIONS! shouldn't happen!" <<endl; exit(1);
			}
			int virus_type = it->pathogen.GetVirusType();
			switch(virus_type)
			{
				case 0:
				{
					if(inf_type!=3)
						wildtype++;
				}break;
				case 1:
				{
					if(inf_type != 3) //count them only if they are not immune!
						downregulating++;
				}break;
				case 2:
				{
					if(inf_type !=3)
						decoy++;
				}break;
				default: cout <<"ERROR!!!!!! COUNTING OTHER VIRUSES! shouldn't happen!" <<endl; exit(1);
			}
		}
		it_host++; /*/
		// check how many infections each host has
		int number_infections = it_host->CountInfections();
		switch(number_infections)
		{
		case 0: break;
		case 1: simpleInfection++; break;
		case 2: doubleInfection++; break;
		case 3: tripleInfection++; break;
		default: cout << "ERROR!!!! counting more than three infections per host, this should be IMPOSSIBLE!!!! " <<number_infections <<endl; exit(1);
		}

		//check what kind of infections the hosts have
		list<Infection>:: iterator it;
		for(it = it_host->infections.begin(); it!= it_host->infections.end(); it++)
		{
			int inf_type = it->GetInfectionType();
			int virus_type = it->pathogen.GetVirusType();
			switch(virus_type)
			{
				case 0:
					if(inf_type!=3) //count them only if they are not immune!
						wildtype++;
					else
						wildtype_immune ++;
				break;
				case 1:
					if(inf_type != 3) //count them only if they are not immune!
						downregulating++;
					else
						downregulating_immune++;
				break;
				case 2:
					if(inf_type !=3) //count them only if they are not immune!
						decoy++;
					else
						decoy_immune++;
				break;
				default: cout <<"ERROR!!!!!! COUNTING OTHER VIRUSES! shouldn't happen!" <<endl; exit(1);
			}
		}
		it_host++;
	}
}

void World::RemoveDeadHosts_HappyNewYear()
{	
//	int deathcount=0;
	vector<Host>::iterator it = hosts.begin();
	while(it!=hosts.end())
	{
		if(it->IsDead())
		{
			it = hosts.erase(it);
//			deathcount++;
			continue;
		}
		it++;
	}
//	if(deathcount >0)
//		cout << deathcount << endl;
}

/*Functions which do the OUTPUT files*/

/*this function saves the genes of each host -> to keep track of MHC and KIR diversity*/
void World:: SaveGenes()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Genes.log";
	
	genesFile.open(s.c_str(), ios::out);

	genesFile << "#Host Index\t";
	for(unsigned int i = 0; i< hosts.at(0).mhcGenes.size(); i++)
	{
		genesFile << "Mhc"<< i+1 <<"\t";
	}

	for(unsigned int i = 0; i< hosts.at(0).kirGenes.size(); i++)
	{
		genesFile << "Kir"<< i+1 <<"\t";
	}

	genesFile << "\tfunctional Kirs\texpressed Kirs\t host_id \n";
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
void World::SavePopulationSize()
{
	//populationSize << "hehehehehehe"<<endl;
	//populationSize <<simulationTime/YEAR <<"\t"<< hosts.size() <<"\t"<<acute_infected << "\t"<< chronic_infected << "\t"<< immune << "\t"<<wildtype<<"\t"<<downregulating<< "\t"<< decoy <<"\n";

	populationSize <<simulationTime/YEAR <<"\t"<< hosts.size() <<"\t"<<simpleInfection << "\t"<< doubleInfection<< "\t"<< tripleInfection<< "\t"<<wildtype<<"\t"<<wildtype_immune << "\t"<<downregulating<< "\t"<<downregulating_immune<<"\t"<< decoy <<"\t"<<decoy_immune<<"\n";
}
/*This function keeps track (ans saves) several parameters of the host and virus*/
void World::SaveParameters()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".log";
	parameterFile.open(s.c_str(), ios::out);
	parameterFile << "#HostIndex\t Age\t totalInfections\t infectionTime\t infectionType\t virusType\t viralLoad\t originalViralLoad\t decoyID\t onlyAcute\t Host_ID\n";
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
	s+=".Map.log";
	mapFile.open(s.c_str(), ios::out);

	map< pair <int,int>, pair <int, int> > ::iterator it;
	for(it = KIRGenesMap.GetMap().begin(); it != KIRGenesMap.GetMap().end(); it ++)
	{
		mapFile << (*it).first.first << "|" <<(*it).first.second << "|" << (*it).second.first <<"|"<< (*it).second.second <<endl;
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
		s+=".Age.log";
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

double World :: GetIntrinsicDeathRate(const double age, const double viralDeathRate)const
{
	int a = round(age*52.0);
	if(a >= deathRates.size())
		cout << a << "|" << deathRates.size() << "age" << "ERROR in the deathrate!!!!" <<endl;
	double intrinsicDeathRate = deathRates.at(a) + viralDeathRate;
	// 	cout << "a, age, rate at, intrinsicdeathRates: "<< a << ", " << age << " , " << rates.at(a) << " , "<<intrinsicDeathRate << endl;
	return intrinsicDeathRate;
}

double World :: GetAgeDependentBirthRate(const double age)const
{
	int a = round(age*52.0);
	if(a >= birthRates.size())
		cout << a << "|" << birthRates.size() << "age" << "ERROR in the birth rate!!!!" <<endl;
	double ageDependentBirthrate = birthRates.at(a);
	return ageDependentBirthrate;
}
