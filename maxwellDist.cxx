#include "maxwellDist.h"
#include<iostream>
#include<fstream>
#include<ctime>
#include<cmath>
#include<iomanip>
#include "constants.h"
using namespace std;


maxwellDist::maxwellDist(int d)
{
	time_t seconds;
	seconds = time (NULL);
	unsigned int seed = seconds;
	srand(seed);

	dim = d;
}


maxwellDist::~maxwellDist()
{
	partList.clear();
        vector< vector< vector<Particle*> > >().swap(partList);
}


double maxwellDist::ran()
{
	return double(rand())/double(RAND_MAX);
}


void maxwellDist::initBox()
{
	cellX = 5*radius; 
	cellY = 5*radius;
	nCellX = boxX/cellX + 1;
	nCellY = boxY/cellY + 1;
	nCells = nCellX * nCellY;

	partList.resize(nCellX);
	for(unsigned i=0; i<partList.size(); i++) partList[i].resize(nCellY);
}


void maxwellDist::initParticles()
{
	double therm_vel = sqrt( (dim * temp * 1.38 * 10000)/(mass * 1.67) );
	therm_vel = therm_vel/100;

	int n = 0;
	while(n < nPart) {
		double x = ran() * boxX;
		double y = ran() * boxY;
		Vector2d pos(x, y);
//		if(checkCloseNess(pos)) continue;

		double theta = 2 * pi * ran() - pi;
		double v_x = therm_vel * cos(theta);
		double v_y = therm_vel * sin(theta);
		Vector2d v = Vector2d(v_x, v_y);

		n++;
		Particle *part = new Particle(n, pos, v, 1.0, radius, false);

		int cX = x/cellX;
		int cY = y/cellY;
		part->setCell(cX, cY);

		partList[cX][cY].push_back(part);
	}
	cout << "initialization done" << endl;
}


bool maxwellDist::checkCloseNess(Vector2d pos1)
{
	for(unsigned int i=0; i<partList.size(); i++) {
		for(unsigned int j=0; j<partList[i].size(); j++) {
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				Particle *part = partList[i][j][k];
				Vector2d pos2 = part->getPos();

				double diff = (pos2 - pos1).Mag();
				if(diff <= 2*radius) return true;
			}
		}
	}
	return false;
}


Particle* maxwellDist::getClosest(Particle *part1)
{
	double dist = 99999999;

	int id1 = part1->getId();
	Vector2d pos1 = part1->getPos();
	int n = part1->getCell().X();
	int m = part1->getCell().Y();

	int iMax = max(n-1, 0);
	int iMin = min(n+1, nCellX);

	int jMax = max(m-1, 0);
	int jMin = min(m+1, nCellY);

	int ii = -1;
	int jj = -1;
	int kk = -1;
	for(int i=iMax; i<iMin; i++) {
		for(int j=jMax; j<jMin; j++) {
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				Particle *part2 = partList[i][j][k];
				int id2 = part2->getId();
				Vector2d pos2 = part2->getPos();

				if(id1 == id2 || part2->isCollided()) continue;

				double diff = (pos2 - pos1).Mag();
				if(diff < dist) {
					dist = diff;
					ii = i;
					jj = j;
					kk = k;
				}
			}
		}
	}

	if(ii < 0 || jj < 0 || kk < 0) return NULL;
	else return partList[ii][jj][kk];
}

void maxwellDist::checkCollision()
{
	for(unsigned int i=0; i<partList.size(); i++) {
		for(unsigned int j=0; j<partList[i].size(); j++) {
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				Particle *part1 = partList[i][j][k];
				Particle *part2 = getClosest(part1);
				if(!part2) continue;

				double mass1 = part1->getMass();
				double mass2 = part2->getMass();

				if(part1->isCollided() or part2->isCollided()) continue;
				if(part1->getId() == part2->getId()) continue;

				Vector2d pos1 = part1->getPos();
				Vector2d pos2 = part2->getPos();
				Vector2d rad_vector = pos2 - pos1;

				if(rad_vector.Mag() > 2*radius) continue;

				Vector2d v1 = part1->getVel();
				Vector2d v2 = part2->getVel();

				Vector2d v1_rest = v1 - v2;
				Vector2d v2_rest = v2 - v2;
				Vector2d v1_rest_rot = v1_rest.Rotate(-v1_rest.Theta());
				Vector2d v2_rest_rot = v2_rest;

				rad_vector = rad_vector.Rotate(-v1_rest.Theta());
				double theta2 = pi * ran() - pi/2;

				double v2_rest_rot_mag_f = 2.0 * mass1 * v1_rest_rot.Mag() * cos(theta2) / (mass1 + mass2);
				double v1_rest_rot_mag_f2 = (mass1 * v1_rest_rot.Mag()*v1_rest_rot.Mag() - mass2 * v2_rest_rot_mag_f*v2_rest_rot_mag_f) / mass1;
				double v1_rest_rot_mag_f = sqrt(v1_rest_rot_mag_f2);

				double theta1 = 0.0;
				if(fabs(mass1*v1_rest_rot_mag_f) < 0.00001) theta1 = 0.0;
				else theta1 = asin( -mass2*v2_rest_rot_mag_f*sin(theta2) / (mass1*v1_rest_rot_mag_f) );

				double v1_rest_rot_f_x = v1_rest_rot_mag_f * cos(theta1);
				double v1_rest_rot_f_y = v1_rest_rot_mag_f * sin(theta1);

				double v2_rest_rot_f_x = v2_rest_rot_mag_f * cos(theta2);
				double v2_rest_rot_f_y = v2_rest_rot_mag_f * sin(theta2);

				Vector2d v1_rest_rot_f = Vector2d(v1_rest_rot_f_x, v1_rest_rot_f_y);
				Vector2d v2_rest_rot_f = Vector2d(v2_rest_rot_f_x, v2_rest_rot_f_y);

				if(checkConservation(v1_rest_rot, v2_rest_rot, v1_rest_rot_f, v2_rest_rot_f, mass1, mass2)) {}
				else {
					cout << "conservation in rest frame is not ok" << endl;
					exit(-1);
				}

				// convert to the lab frame -----------------
				Vector2d v1_lab_i = v1_rest_rot.Rotate(v1_rest.Theta()) + v2;
				Vector2d v2_lab_i = v2_rest_rot.Rotate(v1_rest.Theta()) + v2;

				Vector2d v1_lab_f = v1_rest_rot_f.Rotate(v1_rest.Theta()) + v2;
				Vector2d v2_lab_f = v2_rest_rot_f.Rotate(v1_rest.Theta()) + v2;
				//-------------------------------------------

				if(checkConservation(v1_lab_i, v2_lab_i, v1_lab_f, v2_lab_f, mass1, mass2)) {
					part1->setVel(v1_lab_f);
					part2->setVel(v2_lab_f);
					part1->setCollision(true);
					part2->setCollision(true);

					if(isnan(v1_lab_f.Mag()) || isnan(v2_lab_f.Mag())) {
						cout << "nan == " << v1_rest_rot.Mag() << "  " << v1_rest_rot_mag_f << "  " << v2_rest_rot_mag_f << "  " << theta1 << "  " << theta2 << "  " << endl;
				//		exit(-1);
						continue;
					}
				}
				else {
					cout << "conservation in lab frame is not ok" << endl;
					exit(-1);
				}
			}
		}
	}
}


bool maxwellDist::checkConservation(Vector2d v1_i, Vector2d v2_i, Vector2d v1_f, Vector2d v2_f, double mass1, double mass2)
{
	double momX_i = mass1*v1_i.X() + mass2*v2_i.X();
	double momX_f = mass1*v1_f.X() + mass2*v2_f.X();

	double momY_i = mass1*v1_i.Y() + mass2*v2_i.Y();
	double momY_f = mass1*v1_f.Y() + mass2*v2_f.Y();

	double KE_i = 0.5*mass1*v1_i.Mag()*v1_i.Mag() + 0.5*mass2*v2_i.Mag()*v2_i.Mag();
	double KE_f = 0.5*mass1*v1_f.Mag()*v1_f.Mag() + 0.5*mass2*v2_f.Mag()*v2_f.Mag();

	if(fabs(momX_i-momX_f) > 0.001 || fabs(momY_i-momY_f) > 0.001 || fabs(KE_i-KE_f) > 0.001) {
		cout << momX_i-momX_f << "\t" << momY_i-momY_f << "\t" << KE_i-KE_f << endl;
		return false;
	}
	else {
		return true;
	}
}


void maxwellDist::update(double tStep)
{
	for(unsigned int i=0; i<partList.size(); i++) {
		for(unsigned int j=0; j<partList[i].size(); j++) {
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				Particle *part = partList[i][j][k];

				Vector2d vel = part->getVel();
				Vector2d pos = part->getPos();

				double vel_mag = vel.Mag();
				double theta = vel.Theta();
				double new_vel = vel_mag;

				double xx = pos.X() + tStep * new_vel * cos(theta);
				double yy = pos.Y() + tStep * new_vel * sin(theta);
				Vector2d pos_tmp = Vector2d(xx, yy);
				part->setPos(pos_tmp);
			}
		}
	}
}


void maxwellDist::checkBoundary()
{
	for(unsigned int i=0; i<partList.size(); i++) {
		for(unsigned int j=0; j<partList[i].size(); j++) {
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				Particle *part = partList[i][j][k];
				Vector2d pos = part->getPos();
				double x = pos.X();
				double y = pos.Y();
				if(pos.X() >= boxX) x = 0.0;
				else if(pos.X() < 0) x = boxX;
				if(pos.Y() >= boxY) y = 0.0;
				else if(pos.Y() < 0) y = boxY;
				part->setPos(Vector2d(x, y));
			}
		}
	}
}


void maxwellDist::reset()
{
	vector< vector< vector<Particle*> > > temp = partList;
	for(unsigned int i=0; i<temp.size(); i++) {
		for(unsigned int j=0; j<temp[i].size(); j++) {
			for(unsigned int k=0; k<temp[i][j].size(); k++) {
				Particle *part = temp[i][j][k];
				part->setCollision(false);

				Vector2d pos = part->getPos();
				double x = pos.X();
				double y = pos.Y();

				part->setCell(int(x/cellX), int(y/cellY));
				partList[i][j][k] = part;
			}
		}
	}
}


void maxwellDist::run(int nLoops, double tStep)
{
	initParticles();
	for(int i=0; i<nLoops; i++) {
		cout << "loop " << i+1 << "/" << nLoops << endl;
		checkCollision();
		update(tStep);
		checkBoundary();
		reset();
	}
}


void maxwellDist::write(string file)
{
	ofstream out(file.c_str());
	for(unsigned int i=0; i<partList.size(); i++) {
		for(unsigned int j=0; j<partList[i].size(); j++) {
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				out << fixed << setprecision(2) << partList[i][j][k]->getId() << "\t" << partList[i][j][k]->getVel().X() << "\t" << partList[i][j][k]->getVel().Y() << endl;
			}
		}
	}
	out.close();
}


void maxwellDist::print() 
{
	for(unsigned int i=0; i<partList.size(); i++) {
		for(unsigned int j=0; j<partList[i].size(); j++) {
//			cout << partList[i][j].size() << endl;
			for(unsigned int k=0; k<partList[i][j].size(); k++) {
				cout << fixed << setprecision(2) << partList[i][j][k]->getId() << "\t" << partList[i][j][k]->isCollided() << "\t" << partList[i][j][k]->getCell().X() << "\t" << partList[i][j][k]->getCell().Y() << "\t" << partList[i][j][k]->getPos().X() << "\t" << partList[i][j][k]->getPos().Y() << "\t" << partList[i][j][k]->getVel().X() << "\t" << partList[i][j][k]->getVel().Y() << "\t" << partList[i][j][k]->getVel().Mag() << endl;
			}
		}
	}
}
