#include "Particle.h"
#include<iostream>
using namespace std;

Particle::Particle()
{
	time = 0.0;
	mass = 0.0;
	id = 0;

	vel = Vector2d(0.0, 0.0);
	pos = Vector2d(0.0, 0.0);
}


Particle::Particle(int idd, Vector2d poss, Vector2d vell, double m, double rad, bool coll)
{
	id = idd;
	pos = poss;
	vel = vell;
	mass = m;
	radius = rad;
	collision = coll;
}
