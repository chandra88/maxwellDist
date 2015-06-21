/*================================================
Author:		Chandra Shekhar Sah
Program:	Particle class to hold a particle
Version:	1.0
Date:		June 14, 2015
=================================================*/

#ifndef __Particle__
#define __Particle__

#include<vector>
#include<cmath>
#include<iostream>
#include "Vector2d.h"
using namespace std;

class Particle 
{
	private:
		int id;					// particle id
		int charge;				// charge
		double mass;				// PDG mass
		Vector2d vel;				// particle current velocity
		Vector2d pos;				// particle current vertex
		double time;				// particle current time
		double radius;				// particle radius
		bool collision;
		int cellX;				// cell number
		int cellY;
	
	public:
		Particle();
		Particle(int id, Vector2d pos, Vector2d vel, double mass, double radius, bool coll);
		Particle(const Particle &part) {}
		virtual ~Particle() { }
	
		void setId(int i)					{ id = i; }
		void setCharge(int i=0)					{ charge = i; }
		void setMass(double mas)				{ mass = mas; }
		void setVel(double vx, double vy)			{ vel = Vector2d(vx, vy); }
		void setVel(Vector2d vect)				{ vel = vect; }
		void setPos(double xx, double yy)			{ pos = Vector2d(xx, yy); }
		void setPos(Vector2d p)					{ pos = p; }
		void setTime(double t)					{ time = t; }
		void setRadius(double r)				{ radius = r; }
		void setCollision(bool col)				{ collision = col; }
		void setCell(int x, int y)				{ cellX = x; cellY = y; }

		int getId()			const { return id; }
		int getCharge()			const { return charge; }
		double getMass()		const { return mass; }
		Vector2d getVel()		const { return vel; }
		Vector2d getPos()		const { return pos; }
		double getTime()		const { return time; }
		double getRadius()		const { return radius; }
		bool isCollided()		{ return collision; }
		Vector2d getCell()		{ return Vector2d(cellX, cellY); }

//		void update(double t);
};
#endif
