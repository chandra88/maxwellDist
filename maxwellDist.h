//==========================================================
// Author: Chandra S. Sah
// Program: Base class for velocity distribution of idea gas
// Version: 1.0
// Date: June 21, 2015
//==========================================================

#ifndef __maxwellDist__
#define __maxwellDist__

#include<vector>
#include<cstdlib>
#include "Particle.h"
#include "Vector2d.h"
using namespace std;

class maxwellDist
{
	private:
		vector< vector< vector<Particle*> > > partList;
		double temp;					// temperature

		int dim;

		double boxX;					// length of the container
		double boxY;					// width of the container

		int nCellX;
		int nCellY;
		int nCells;

		double cellX;
		double cellY;

		double radius;	
		double mass;
		int nPart;

	public:
		maxwellDist(int d);
		maxwellDist(int d, double x, double y)			{ dim = d; boxX = x; boxY = y; }
		virtual ~maxwellDist();
	
		void initBox();
		void initParticles();
		void setBox(int x, int y)				{ boxX = x; boxY = y; }
		void setParticleType(double r, double m, int n)		{ radius = r; mass = m; nPart = n; } 
		void setTemp(double t)					{ temp = t; }

		double ran();
		bool checkCloseNess(Vector2d p);

		int getDim()		{ return dim; }
		Vector2d getBoxSize()	{ return Vector2d(boxX, boxY); }
		double getRadius()	{ return radius; }
		double getMass()	{ return mass; }
		double getTemp()	{ return temp; }
		Vector2d getCellSize()	{ return Vector2d(nCellX, nCellY); }

		Particle *getClosest(Particle *part1);
		void checkCollision();
		bool checkConservation(Vector2d v1_i, Vector2d v2_i, Vector2d v1_f, Vector2d v2_f, double mass1, double mass2);
		void update(double tStep=1.0);
		void checkBoundary();
		void reset();
		void run(int loop, double tStep=1.0);
		void write(string file);

		void print();
};
#endif
