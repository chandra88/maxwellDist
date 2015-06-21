#include "Vector2d.h"


Vector2d Vector2d::operator +(const Vector2d& other)
{
	Vector2d vect;
	vect.x = this->x + other.x;
	vect.y = this->y + other.y;
	return vect;
}

Vector2d Vector2d::operator -(const Vector2d& other)
{
	Vector2d vect;
	vect.x = this->x - other.x;
	vect.y = this->y - other.y;
	return vect;
}

Vector2d& Vector2d::operator =(const Vector2d& other)
{
	x = other.x;
	y = other.y;
	return *this;
}

bool Vector2d::operator ==(const Vector2d& other)
{
	if(other.x == this->x && other.y == this->y) return true;
	else return false;
}

bool Vector2d::operator !=(const Vector2d& other)
{
	if(other.x != this->x || other.y != this->y) return true;
	else return false;
}

Vector2d Vector2d::operator *(const double fact)
{
	Vector2d vect;
	vect.x = this->x*fact;
	vect.y = this->y*fact;
	return vect;
}

double Vector2d::Angle(Vector2d& other)
{
	return atan2(this->y - other.y, this->x - other.x);
}

Vector2d Vector2d::Rotate(double th)
{
	double xx = this->x*cos(th) - this->y*sin(th);
	double yy = this->x*sin(th) + this->y*cos(th);
	return Vector2d(xx, yy);
}
