/*================================================
Author:		Chandra S. Sah
Program:	2d vector
Verson:		1.0
Date:		June 14, 2015
================================================*/

#ifndef __Vector2d__
#define __Vector2d__
#include<cmath>
using namespace std;

class Vector2d
{
	protected:
		double x, y;

	public:
		Vector2d()				{ x = 0.0; y =0.0; }
		Vector2d(double xx, double yy)		{ x = xx; y = yy; }
		virtual ~Vector2d() { }

		Vector2d operator +(const Vector2d& other);
		Vector2d operator -(const Vector2d& other);
		Vector2d& operator =(const Vector2d& other);
		Vector2d operator *(const double fact);
		bool operator ==(const Vector2d& other);
		bool operator !=(const Vector2d& other);

		double Dot(Vector2d vect)		{ return sqrt(x*vect.x + y*vect.y); }

		void setVector2d(double xx, double yy)	{ x = xx; y = yy; }	
		void setVector2d(Vector2d vect)		{ x = vect.X(); y = vect.Y(); }
		void setX(double xx)			{ x = xx; }
		void setY(double yy)			{ y = yy; }

		double X()				{ return x; }
		double Y()				{ return y; }
		double Theta()				{ return atan2(y, x); }
		double Angle(Vector2d& vec);
		double Mag()				{ return sqrt(x*x + y*y); }
		Vector2d Rotate(double th);
};
#endif
