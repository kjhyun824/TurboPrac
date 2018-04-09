#include <iostream>
#include<string>

#include "InputCommandLineParser.h"
#include "LevelSub.h"
//#include "TimeUtility.h"
#include "runtimecounter.h"

using namespace std;

std::string dataGraphFileName;
std::string queryGraphFileName;
int k;
int debug;

void help();

int main(int argc, char** argv)
{
	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-q")){
		queryGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-q");
	}

	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-d")){
		dataGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-d");
	}

	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-k")){
		k = stoi(InputCommandLineParser::getCmdOption(argc, argv, "-k"));
	}
	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-p")){
		debug = stoi(InputCommandLineParser::getCmdOption(argc, argv, "-p"));
	}
	if (!InputCommandLineParser::cmdOptionExists(argc, argv, "-q") || !InputCommandLineParser::cmdOptionExists(argc, argv, "-d") || !InputCommandLineParser::cmdOptionExists(argc, argv, "-d")){
		cout << "Wrong Parameters" << endl;
		help();
		exit(1);
	}

//	TimeUtility timer;
//	timer.StartCounterMill();
	Runtimecounter timer;
	timer.start();
	LevelSub levelSub(dataGraphFileName, queryGraphFileName, k, debug);
	timer.stop();
//	std::cout << "Load time(ms):" << timer.GetCounterMill() << std::endl;
	std::cout << "Load time(ms):" << timer.GetRuntime() << std::endl;

//	TimeUtility totaltimer;
//	totaltimer.StartCounterMill();
	Runtimecounter totaltimer;
	totaltimer.start();
	levelSub.subgraphSearch();
//	levelSub.exploreCR();
	totaltimer.stop();
//	std::cout << "Total Cost time(ms): " << totaltimer.GetCounterMill() << std::endl;
	std::cout << "Total Cost time(ms): " << totaltimer.GetRuntime() << std::endl;
	return 0;
}

void help(){
	cout << "Options: " << endl;
	cout << "-d  The datagraph file" << endl;
	cout << "-q  The querygraph file" << endl;
	cout << "-k	 The top k" << endl;
}
