/*
* Author: ziga
*/

/* This file is part of ultrasound propagation model and is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "source.h"
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <math.h>  

float density = 3.0f;
float fis = PI/5.0f;



Probe::Probe (int pc, int qc, int freq, float amp, float rad, float fic, bool focusedc) {
	//set number of points
	focused = focusedc;
	radius = rad;	
	if(focused) {
		numOfPoints = rad*fis*density;		
	}
	else {
		numOfPoints = radius*sqrt(2.0f - 2*(1.0f - cos(fis)))*density;
	}
	p = pc;
	q = qc;
	minFrequency = 500000; // >= 500 kHz
	maxFrequency = 5000000; // <= 5Mhz
	setFrequency(freq);
	minAmplitude = 10000.0f;
	maxAmplitude = 100000.0f;
	setAmplitude(amp);
	fi = fic;
	//allocate radius array
	setPointMesh(true);
}

void Probe::setFocusedMesh(bool init) {
	float dt = fis/numOfPoints;
	float t = fi;	

	float sint;
	float cost;	
	if(!init) {
		free(pointMesh);
		free(damperMesh);
	}
	pointMesh = (point *)malloc(numOfPoints*sizeof(point));
	damperMesh = (point *)malloc((2*numOfPoints)*sizeof(point));
	
	for(int i = 0; i < numOfPoints; i++) {
		cost = cos(t);
		sint = sin(t);	

		pointMesh[i].x = ((int)(p + (radius*cost))%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		pointMesh[i].y = ((int)(q + (radius*sint))%NUM_FIN_EL +  NUM_FIN_EL)%NUM_FIN_EL;				

		damperMesh[2*i + 0].x = ((int)(p + (radius + 1)*cost)%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		damperMesh[2*i + 0].y = ((int)(q + (radius + 1)*sint)%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;

		damperMesh[2*i + 1].x = ((int)(p + (radius + 2)*cos(t + dt/2))%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		damperMesh[2*i + 1].y = ((int)(q + (radius + 2)*sin(t + dt/2))%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		//update t
		t+=dt;
	}
}

void Probe::setUnfocusedMesh(bool init) {
	float t = 0;
	float vx = radius*((cos(fi + fis) - cos(fi)));
	float vy = radius*((sin(fi + fis) - sin(fi)));	
	float ps = p + radius*(cos(fi));
	float qs  = q + radius*(sin(fi));	
	float psm1 = p + (radius + 1)*(cos(fi));
	float qsm1  = q + (radius + 1)*(sin(fi));
	float psm2 = p + (radius + 2)*(cos(fi));
	float qsm2  = q + (radius + 2)*(sin(fi));		
	float dt = 1.0f/numOfPoints;
	
	if(!init) {
		free(pointMesh);
		free(damperMesh);
	}
	pointMesh = (point *)malloc(numOfPoints*sizeof(point));
	damperMesh = (point *)malloc((2*numOfPoints)*sizeof(point));
	
	for(int i = 0; i < numOfPoints; i++) {
		pointMesh[i].x = ((int)(ps + (t*vx))%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		pointMesh[i].y = ((int)(qs + (t*vy))%NUM_FIN_EL +  NUM_FIN_EL)%NUM_FIN_EL;				

		damperMesh[2*i + 0].x = ((int)(psm1 + (t*vx))%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		damperMesh[2*i + 0].y = ((int)(qsm1 + (t*vy))%NUM_FIN_EL +  NUM_FIN_EL)%NUM_FIN_EL;

		damperMesh[2*i + 1].x = ((int)(psm2 + ((t + dt/2)*vx))%NUM_FIN_EL + NUM_FIN_EL)%NUM_FIN_EL;
		damperMesh[2*i + 1].y = ((int)(qsm2 + ((t + dt/2)*vy))%NUM_FIN_EL +  NUM_FIN_EL)%NUM_FIN_EL;

		t += dt;
	}
}

void Probe::setPointMesh(bool init) {
	if(focused) {
		//set point mesh for focused probe
		setFocusedMesh(init);
	}
	else {
		//set probe mesh for unfocused mesh
		setUnfocusedMesh(init);
	}
}

void Probe::updateValue(float t) {
	value = sin(t*frequency)*amplitude;
}
point * Probe::getPoints() {
	return pointMesh;
}
point * Probe::getDampedPoints() {
	return damperMesh;
}
float Probe::getValue() {
	return value;
}
int Probe::getNumOfPoints() {
	return numOfPoints;
}
void Probe::setFi(int x, int y) {
	x = x - p;
	y = y - q;
	fi = atan2(y,x) - fis/2.0;
}
void Probe::setPQ(int P, int Q) {
	float fih = fi + fis/2.0f;	
	if(focused) {
		p = P - (int)(radius*cos(fih));
		q = Q - (int)(radius*sin(fih));		
	}
	else {
		float vx = radius*((cos(fi + fis) - cos(fi)));
		float vy = radius*((sin(fi + fis) - sin(fi)));	
		p = P - (int)(radius*cos(fi + fis) - 0.5*vx);
		q = Q - (int)(radius*sin(fi + fis) - 0.5*vy);
	}
}
void Probe::setFrequency(int freq) {
	if(freq > maxFrequency) {
		freq = maxFrequency;
	}
	else if(freq < minFrequency) {
		freq = minFrequency;
	}
	frequency = freq;
}
float Probe::getRadius() {
	return radius;
}
point * Probe::getMiddlePoint() {
	return &(pointMesh[numOfPoints/2]);
}
int Probe::getMaxFrequency() {
	return maxFrequency;
}
int Probe::getMinFrequency() {
	return minFrequency;
}
int Probe::getFrequency() {
	return frequency;
}
void Probe::setAmplitude(float amp) {
	if(amp > maxAmplitude) {
		amp = maxAmplitude;
	}
	else if(amp < minAmplitude) {
		amp = minAmplitude;
	}
	amplitude = amp;	
}

float Probe::getAmplitude(){
	return amplitude;
}
float Probe::getMaxAmplitude(){
	return maxAmplitude;
}
float Probe::getMinAmplitude(){
	return minAmplitude;
}
float Probe::getFi() {
	return fi;
}
bool Probe::isFocused() {
	return focused;
}

