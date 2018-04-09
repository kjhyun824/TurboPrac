#include <windows.h>
#include <iostream>


#ifndef TIME_UTILITY
#define TIME_UTILITY

class TimeUtility{
private:
	double PCFreq = 0.0;
	__int64 CounterStart = 0;

public:
	void StartCounterMicro();
	double GetCounterMicro();

	void StartCounterMill();
	double GetCounterMill();
};

#endif