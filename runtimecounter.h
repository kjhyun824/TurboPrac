#ifndef RUNTIMECOUNTER_H_
#define RUNTIMECOUNTER_H_
#include <cstring>
#include <sys/time.h>

class Runtimecounter{
public: //timezone tz;
	timeval t1;
	timeval t2;
	public:
		Runtimecounter();
		void start();
		void stop();
		double GetRuntime();
};

#endif /* RUNTIMECOUNTER_H_ */

