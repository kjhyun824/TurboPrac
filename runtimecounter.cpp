#include "runtimecounter.h"
Runtimecounter::Runtimecounter(){
}

void Runtimecounter::start(){
	gettimeofday(&this->t1, NULL);
}

void Runtimecounter::stop(){
	gettimeofday(&this->t2, NULL);
}

double Runtimecounter::GetRuntime(){
	double t=(double)(t2.tv_sec-t1.tv_sec)*1000.0+(double)(t2.tv_usec-t1.tv_usec)/1000.0;
	return t;
}

