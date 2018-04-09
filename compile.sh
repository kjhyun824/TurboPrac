rm DSQL 
g++ -c InputCommandLineParser.cpp -std=c++11
g++ -c runtimecounter.cpp -std=c++11
g++ -c StringUtility.cpp -std=c++11
g++ -c LevelSub.cpp -std=c++11 -I ../boost_1_58_0/
g++ -o DSQL main.cpp  *.o -std=c++11 -I ../boost_1_58_0/
