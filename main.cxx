#include<iostream>
#include<fstream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<iomanip>
#include "Vector2d.h"
#include "Particle.h"
#include "constants.h"
#include "maxwellDist.h"
using namespace std;

int main()
{
	maxwellDist *max = new maxwellDist(2, 800.0, 600.0);
	max->setParticleType(4.0, 1.0, 20000);
	max->setTemp(300.0);
	max->initBox();
	
	max->run(1000, 2.0);
	max->write("data.dat");
//	max->print();
}
