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


#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <new>
#include "source.h"

float *k1u;
float *k2u;
float *k3u;
float *k4u;
float *k1v;
float *k2v;
float *k3v;
float *k4v;
float *u;
float *v;
float *c_tissueS; //speed of sound in tissue (values are already squared)
float *tempu;
float *tempv;
float *temp;
float *displayMatrixX;
float *displayMatrixY;
int *damping;
Probe **sources;

int frameCount = 0;
//Number of frames per second
float fps = 0;
int currentTime = 0;
int previousTime = 0; 

bool pause = false;
bool newFocProbe = false;
bool newUnfocProbe = false;
bool removeProbe = false;

float c = 1.0f;
float k2 = 4000000.0f;
float t = 0;
int ms = NUM_FIN_EL*NUM_FIN_EL;
float ts_h = TS/2.0;
float ts_s = TS/6.0;
float hss = H*H;
float ratio = 1.0f/hss;
float dx;	
float dy;

clock_t reset;
GLuint texture;
GLuint probefTexture;
GLuint probeuTexture;


float alfa;
float alfac;
float alfaca;
float alfa_1;
int numSources;

bool holdingSource;
bool holdingFrequencySlider;
float frequencySliderPos;
bool holdingAmpSlider;
float ampSliderPos;
bool roMode;
int sourceId; 

float displayPressure = 0;
int fileSaveIndex = 0;

void addSource(int p, int q, int freq, int amp, float radius, float fi, bool focused) {
	numSources++;
	Probe* src = (Probe *)malloc(sizeof(Probe));
	sources = (Probe ** )realloc(sources, numSources*sizeof(Probe *));	
	//creating new Source
	src = new(src) Probe(p, q, freq, amp, radius, fi, focused);	
	sources[numSources - 1] = src;
}
/**

**/
float getSpeed(unsigned char red, unsigned char green, unsigned char blue) {
	float speed;
	//Brain (Grey Matter) 	1500.0 m/s
	//Brain (White Matter) 	1552.5 m/s
	//Cerebrospinal Fluid (csf)	1504.5 m/s
	//Skin 1624.0
	//Skull Cancellous 	2117.5 m/s
	//Skull Cortical 	2813.7 m/s
	//air 343.0 m/s

	if(red == 0 && green == 255 && blue == 0) {
		//skin
		speed = 1730.0f;
	}
	//185 if mrit2
	else if(red == 229 && green == 229 && blue == 229) {
		//white matter
		speed = 1570.5f;
	}
	else if(red == 51 && green == 51 && blue == 51) {
		//grey matter
		speed = 1570.0f;
	}
	else if(red == 255 && green == 255 && blue == 255) {
		//therapeutic gel
		speed = 1400;
	}
	else if(red == 0 && green == 0 && blue == 255) {
		//cortical skull
		speed = 4080.0;
	}
	else if(red == 0 && green == 255 && blue == 255) {
		//Cerebrospinal Fluid (csf)	
		speed = 1500;
	}
	else {		
		printf("unrecognized class: %d %d %d \n", red, green, blue);
		speed = 0;		
	}
	return speed*speed;
}
int getAttenuationCoefficent(unsigned char red, unsigned char green, unsigned char blue) {
	if(red == 0 && green == 0 && blue == 255) {
		return 160000;
	}
	return 16000;
}

/**
Reads bitmap image and sets sound speed for given tissue
**/
void setTissue(char const* filename) {
    FILE* f = fopen(filename, "rb");
    unsigned char info[54]; //54
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int widthImage = *(int*)&info[18];
    int heightImage = *(int*)&info[22];
	
    int size =  3 * widthImage * heightImage;
    unsigned char* data = (unsigned char*)malloc(size); // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);

	float r1 = (float)widthImage/NUM_FIN_EL;
	float r2 = (float)heightImage/NUM_FIN_EL;
	
	float* speedData = (float *)malloc(widthImage*heightImage*sizeof(float));
	int* dampingData = (int *)malloc(widthImage*heightImage*sizeof(int));
	int j = 0;
	
	for(int i = 0; i < size; i += 3) {
		unsigned char tmp = data[i];
		data[i] = data[i + 2];
		data[i + 2] = tmp;	
		speedData[j] = getSpeed(data[i], data[i+1], data[i + 2]);
		dampingData[j] = getAttenuationCoefficent(data[i], data[i+1], data[i + 2]);
		j++;
	}
	
	//map pixel data to speed data and also damping data
	for(int x= 0; x < NUM_FIN_EL; x++) {
		for(int y = 0; y < NUM_FIN_EL; y++) {
			c_tissueS[y*NUM_FIN_EL + x] = speedData[(int)((int)(y*r2)*widthImage + (int)(x*r1))];
			damping[y*NUM_FIN_EL + x] = dampingData[(int)((int)(y*r2)*widthImage + (int)(x*r1))];
		}
	}
	//free buffer
	free(data);
	
	//free unnecessary speed and damping data
	free(speedData);
	free(dampingData);
}


void setprobefTexture(char const* filename1, char const* filename2) {
    FILE* f1 = fopen(filename1, "rb");
    FILE* f2 = fopen(filename2, "rb");
    unsigned char info1[54]; //54
    unsigned char info2[54]; //54
    fread(info1, sizeof(unsigned char), 54, f1); // read the 54-byte header
    fread(info2, sizeof(unsigned char), 54, f2); // read the 54-byte header

    // extract image height and width from header
    int widthImage1 = *(int*)&info1[18];
    int widthImage2 = *(int*)&info2[18];
    int heightImage1 = *(int*)&info1[22];
    int heightImage2 = *(int*)&info2[22];
	
    int size1 =  3 * widthImage1 * heightImage1;
    int size2 =  3 * widthImage2 * heightImage2;
    unsigned char* data1 = (unsigned char*)malloc(size1); // allocate 3 bytes per pixel
	unsigned char* buffer1 = (unsigned char*)malloc(size1 + widthImage1 * heightImage1);
    unsigned char* data2 = (unsigned char*)malloc(size2); // allocate 3 bytes per pixel
	unsigned char* buffer2 = (unsigned char*)malloc(size2 + widthImage2 * heightImage2);
    fread(data1, sizeof(unsigned char), size1, f1); // read the rest of the data at once
    fread(data2, sizeof(unsigned char), size2, f2); // read the rest of the data at once
    fclose(f1);
    fclose(f2);
	int j = 0;
	for(int i = 0; i < size1; i += 3) {
		buffer1[j] = data1[i+2];
		buffer1[j + 1] = data1[i+1];
		buffer1[j + 2] = data1[i];
		buffer1[j + 3] = 255;
		if(data1[i+2] == 255) {
			buffer1[j + 3] = 0;			
		}
		j += 4;
	}
	j = 0;
	for(int i = 0; i < size2; i += 3) {
		buffer2[j] = data2[i+2];
		buffer2[j + 1] = data2[i+1];
		buffer2[j + 2] = data2[i];
		buffer2[j + 3] = 255;
		if(data2[i + 2] == 255) {
			buffer2[j + 3] = 0;		
		}
		j += 4;
	}	
	
	
	// allocate a texture name
	glGenTextures( 1, &probefTexture );

	// select our current texture
	glBindTexture( GL_TEXTURE_2D, probefTexture );

	// select modulate to mix texture with color for shading
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// when texture area is small, bilinear filter the closest MIP map
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST );
	// when texture area is large, bilinear filter the first MIP map
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

	// if wrap is true, the texture wraps over at the edges (repeat)
	//       ... false, the texture ends at the edges (clamp)
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );

	// build our texture MIP maps
	gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, widthImage1, heightImage1, GL_RGBA, GL_UNSIGNED_BYTE, buffer1 );

	// free buffer
	free(data1);
	free(buffer1);
	
	// allocate a texture name
	glGenTextures( 1, &probeuTexture );

	// select our current texture
	glBindTexture( GL_TEXTURE_2D, probeuTexture );

	// select modulate to mix texture with color for shading
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// when texture area is small, bilinear filter the closest MIP map
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST );
	// when texture area is large, bilinear filter the first MIP map
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

	// if wrap is true, the texture wraps over at the edges (repeat)
	//       ... false, the texture ends at the edges (clamp)
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );

	// build our texture MIP maps
	gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, widthImage2, heightImage2, GL_RGBA, GL_UNSIGNED_BYTE, buffer2 );	
	
	// free buffer
	free(data2);
	free(buffer2);	
}

/**
Reads bitmap image and sets texture
**/
GLuint setTexture(char const* filename) {
	GLuint texture;	
    FILE* f = fopen(filename, "rb");
    unsigned char info[54]; //54
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int widthImage = *(int*)&info[18];
    int heightImage = *(int*)&info[22];
	
    int size =  3 * widthImage * heightImage;
    unsigned char* data = (unsigned char*)malloc(size); // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);
	
	for(int i = 0; i < size; i += 3) {
		unsigned char tmp = data[i];
		data[i] = data[i + 2];
		data[i + 2] = tmp;	
	}
	
	// allocate a texture name
	glGenTextures( 1, &texture );

	// select our current texture
	glBindTexture( GL_TEXTURE_2D, texture );

	// select modulate to mix texture with color for shading
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// when texture area is small, bilinear filter the closest MIP map
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST );
	// when texture area is large, bilinear filter the first MIP map
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

	// if wrap is true, the texture wraps over at the edges (repeat)
	//       ... false, the texture ends at the edges (clamp)
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );

	// build our texture MIP maps
	gluBuild2DMipmaps( GL_TEXTURE_2D, 3, widthImage, heightImage, GL_RGB, GL_UNSIGNED_BYTE, data );

	// free buffer
	free(data);

	return texture;
}

/**
Computes Discrete Laplace. All actual computation is actually done in this method.
**/
void Discrete_Laplace(float *u, float *v, float *du, float *dv) {	
	//du = v
	memcpy(du, v, ms*sizeof(float));
	#pragma omp parallel for
	for(int i = 20; i < NUM_FIN_EL -20; i++) {
		for(int j = 20; j < NUM_FIN_EL - 20; j++) {
			dv[i*NUM_FIN_EL + j] = ratio*c_tissueS[i*NUM_FIN_EL + j]*(u[(i + 1)*NUM_FIN_EL + j] + u[(i - 1)*NUM_FIN_EL + j] + u[i*NUM_FIN_EL + j + 1] + u[i*NUM_FIN_EL + j - 1] - 4*u[i*NUM_FIN_EL + j]) - damping[i*NUM_FIN_EL + j]*v[i*NUM_FIN_EL + j];
			
		}
	}
	
	//padding regions
	#pragma omp parallel for
	for(int i = 1; i < NUM_FIN_EL - 1; i++) {
		for(int j = 1; j < 20; j++) {
			ALPHA(j);
			int je = j + NUM_FIN_EL - 21;
			dv[i*NUM_FIN_EL + j] = ratio*c_tissueS[i*NUM_FIN_EL + j]*(u[(i + 1)*NUM_FIN_EL + j] + u[(i - 1)*NUM_FIN_EL + j] + u[i*NUM_FIN_EL + j + 1] + u[i*NUM_FIN_EL + j - 1] - 4*u[i*NUM_FIN_EL + j]) - ((1.0f - alfaca)*damping[i*NUM_FIN_EL + j] +  alfaca*k2)*v[i*NUM_FIN_EL + j];
			dv[i*NUM_FIN_EL + je] = ratio*c_tissueS[i*NUM_FIN_EL + je]*(u[(i + 1)*NUM_FIN_EL + je] + u[(i - 1)*NUM_FIN_EL + je] + u[i*NUM_FIN_EL + je + 1] + u[i*NUM_FIN_EL + je - 1] - 4*u[i*NUM_FIN_EL + je]) - ((1.0f - alfac)*damping[i*NUM_FIN_EL + je] +  alfac*k2)*v[i*NUM_FIN_EL + je];
		}
		//borders
		dv[i] = ratio*c_tissueS[i]*(u[i + 1] + u[i - 1] + u[i + NUM_FIN_EL]  - 3*u[i]) - k2*v[i];
		dv[(NUM_FIN_EL - 1)*NUM_FIN_EL + i] = ratio*c_tissueS[(NUM_FIN_EL - 1)*NUM_FIN_EL + i]*(u[(NUM_FIN_EL - 1)*NUM_FIN_EL + i + 1] + u[(NUM_FIN_EL - 1)*NUM_FIN_EL + i - 1] + u[(NUM_FIN_EL - 2)*NUM_FIN_EL + i]  - 3*u[(NUM_FIN_EL - 1)*NUM_FIN_EL + i]) - k2*v[(NUM_FIN_EL - 1)*NUM_FIN_EL + i];
		dv[i*NUM_FIN_EL] = ratio*c_tissueS[i*NUM_FIN_EL]*(u[i*NUM_FIN_EL + 1] + u[(i - 1)*NUM_FIN_EL] + u[(i + 1)*NUM_FIN_EL]  - 3*u[i*NUM_FIN_EL]) - k2*v[i*NUM_FIN_EL];
		dv[(i + 1)*NUM_FIN_EL - 1] = ratio*c_tissueS[(i + 1)*NUM_FIN_EL - 1]*(u[(i + 1)*NUM_FIN_EL - 2] + u[(i)*NUM_FIN_EL - 1] + u[(i + 2)*NUM_FIN_EL - 1]  - 3*u[(i + 1)*NUM_FIN_EL - 1]) - k2*v[(i + 1)*NUM_FIN_EL - 1];		
	}
	
	//pading regions
	#pragma omp parallel for
	for(int j = 20; j < NUM_FIN_EL - 20; j++) {
		for(int i = 1; i < 20; i++) {
			ALPHA(i);
			int ie = i + NUM_FIN_EL - 21; 
			dv[i*NUM_FIN_EL + j] = ratio*c_tissueS[i*NUM_FIN_EL + j]*(u[(i + 1)*NUM_FIN_EL + j] + u[(i - 1)*NUM_FIN_EL + j] + u[i*NUM_FIN_EL + j + 1] + u[i*NUM_FIN_EL + j - 1] - 4*u[i*NUM_FIN_EL + j]) - ((1.0f - alfaca)*damping[i*NUM_FIN_EL + j] +  (alfaca)*k2)*v[i*NUM_FIN_EL + j];
			dv[ie*NUM_FIN_EL + j] = ratio*c_tissueS[ie*NUM_FIN_EL + j]*(u[(ie + 1)*NUM_FIN_EL + j] + u[(ie - 1)*NUM_FIN_EL + j] + u[ie*NUM_FIN_EL + j + 1] + u[ie*NUM_FIN_EL + j - 1] - 4*u[ie*NUM_FIN_EL + j]) - ((1.0f - alfac)*damping[ie*NUM_FIN_EL + j] +  alfac*k2)*v[ie*NUM_FIN_EL + j];
		}
	}
	
	//corners
	dv[0] = ratio*(u[1] + u[NUM_FIN_EL] - 2*u[0]) - k2*v[0];
	dv[NUM_FIN_EL - 1] = ratio*(u[NUM_FIN_EL - 2] + u[2*NUM_FIN_EL - 1] - 2*u[NUM_FIN_EL - 1]) - k2*v[NUM_FIN_EL - 1];
	dv[(NUM_FIN_EL - 1)*NUM_FIN_EL] = ratio*(u[(NUM_FIN_EL - 2)*NUM_FIN_EL] + u[(NUM_FIN_EL - 1)*NUM_FIN_EL + 1] - 2*u[(NUM_FIN_EL - 1)*NUM_FIN_EL]) - k2*v[(NUM_FIN_EL - 1)*NUM_FIN_EL];
	dv[ms - 1] = ratio*(u[ms - 2] + u[(NUM_FIN_EL - 1)*NUM_FIN_EL - 1] - 2*u[NUM_FIN_EL*NUM_FIN_EL - 1]) - k2*v[NUM_FIN_EL*NUM_FIN_EL - 1];
}

void calculateFPS() {
    //Increase frame count
    frameCount++;
 
    //Get the number of milliseconds since glutInit called
    //(or first call to glutGet(GLUT ELAPSED TIME)).
    currentTime = glutGet(GLUT_ELAPSED_TIME);
 
    //Calculate time passed
    int timeInterval = currentTime - previousTime;
 
    if(timeInterval > 1000)
    {
        //calculate the number of frames per second
        fps = frameCount / (timeInterval / 1000.0f);
 
        //Set time
        previousTime = currentTime;
 
        //Reset frame count
        frameCount = 0;
    }
}

void updateDampedSource(float * kvalu, float * kvalv) {
	//update source values 
	for(int i = 0; i < numSources; i++) {	
		Probe *s = sources[i];
		point * damped = s->getDampedPoints();
		int numPoints = s->getNumOfPoints();
		int numDamped = 2*numPoints;
		
		for(int j = 0; j < numDamped; j++) {
			kvalu[damped[j].y*NUM_FIN_EL + damped[j].x] = 0.0f;
			kvalv[damped[j].y*NUM_FIN_EL + damped[j].x] = 0.0f;			
		}
	}	
}

void updateSource() {
	//update source values 
	for(int i = 0; i < numSources; i++) {	
		Probe *s = sources[i];
		s->updateValue(t);
		float val = s->getValue();
		point * points = s->getPoints();
		int numPoints = s->getNumOfPoints();

		for(int j = 0; j < numPoints; j++) {
			u[points[j].y*NUM_FIN_EL + points[j].x] = val;
		}
	}
}

/**
Computes Runge-Kutta of 4th order.
**/
void RK4(float *u, float *v, float *tempu, float *tempv) {
	updateSource();
	
	//compute k1
	Discrete_Laplace(u, v, k1u, k1v);
	updateDampedSource(k1u, k1v);
	
	//compute yn + ts/2*k1
	//#pragma omp parallel for
	for(int i = 0; i < ms; i++) {
		tempu[i] = u[i] + ts_h*k1u[i];
		tempv[i] = v[i] + ts_h*k1v[i];
	}
	Discrete_Laplace(tempu, tempv, k2u, k2v);
	updateDampedSource(k2u, k2v);	
	//compute yn + ts/2*k2
	//#pragma omp parallel for
	for(int i = 0; i < ms; i++) {
		tempu[i] = u[i] + ts_h*k2u[i];
		tempv[i] = v[i] + ts_h*k2v[i];
	}
	Discrete_Laplace(tempu, tempv, k3u, k3v);
	updateDampedSource(k3u, k3v);	
	//compute yn + ts*k3
	//#pragma omp parallel for
	for(int i = 0; i < ms; i++) {
		tempu[i] = u[i] + TS*k3u[i];
		tempv[i] = v[i] + TS*k3v[i];
	}
	Discrete_Laplace(tempu, tempv, k4u, k4v);
	updateDampedSource(k4u, k4v);	

	//compute full solution
	#pragma omp parallel for
	for(int i = 0; i < ms; i++) {
		tempu[i] = u[i] + ts_s*(k1u[i] + 2*k2u[i] + 2*k3u[i] + k4u[i]);
		tempv[i] = v[i] + ts_s*(k1v[i] + 2*k2v[i] + 2*k3v[i] + k4v[i]);
	}
}
void drawString(char const *str, float x, float y) {
    glColor3f(1.0f, 1.0f, 1.0f);	
    glRasterPos2f(x, y);		
	
	int len = strlen(str);
	for (int i = 0; i < len; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, (int)str[i]);
	}	
}

void drawMenu() {
			
	//draw menu bar
	glColor3f(0.4, 0.4, 0.4);
	glBegin( GL_QUADS );
		glVertex2d(-1.0, +0.9);
		glVertex2d(+1.0, +0.9);
		glVertex2d(+1.0, +1.0);
		glVertex2d(-1.0, +1.0);
	glEnd();	

	//draw frames per second number
	glColor3f(0.2, 0.2, 0.2);
	glBegin( GL_QUADS );
		glVertex2d(0.72, +0.91);
		glVertex2d(0.99, +0.91);
		glVertex2d(0.99, +0.99);
		glVertex2d(0.72, +0.99);
	glEnd();	
    char ss[15] ;
    sprintf(ss,"%.2f FPS", fps);
	drawString(ss, 0.75f, 0.935f);

	//draw play/pause button
	glColor3f(0.2, 0.2, 0.2);
	glBegin( GL_QUADS );
		glVertex2d(-0.99, +0.91);
		glVertex2d(-0.85, +0.91);
		glVertex2d(-0.85, +0.99);
		glVertex2d(-0.99, +0.99);
	glEnd();
	drawString("|| >>\0", -0.965, 0.935);
	
	//draw new focused probe button
	glColor3f(0.2, 0.2, 0.2);
	glBegin( GL_QUADS );
		glVertex2d(-0.84, +0.91);
		glVertex2d(-0.37, +0.91);
		glVertex2d(-0.37, +0.99);
		glVertex2d(-0.84, +0.99);
	glEnd();
	drawString("New focused probe\0", -0.815, 0.935);
	
	//draw new unfocused probe button
	glColor3f(0.2, 0.2, 0.2);
	glBegin( GL_QUADS );
		glVertex2d(-0.36, +0.91);
		glVertex2d(0.16, +0.91);
		glVertex2d(0.16, +0.99);
		glVertex2d(-0.36, +0.99);
	glEnd();
	drawString("New unfocused probe\0", -0.335, 0.935);	
	
	//draw remove probe button
	glColor3f(0.2, 0.2, 0.2);
	glBegin( GL_QUADS );
		glVertex2d(0.17, +0.91);
		glVertex2d(0.53, +0.91);
		glVertex2d(0.53, +0.99);
		glVertex2d(0.17, +0.99);
	glEnd();
	drawString("Remove probe\0", 0.195, 0.935);	

	//draw clear waves button
	glColor3f(0.2, 0.2, 0.2);
	glBegin( GL_QUADS );
		glVertex2d(0.54, +0.91);
		glVertex2d(0.62, +0.91);
		glVertex2d(0.62, +0.99);
		glVertex2d(0.54, +0.99);
	glEnd();	
	drawString("X\0", 0.565, 0.935);		
	
}

void drawProbes() {

	for(int i = 0; i < numSources; i++) {	
		Probe *s = sources[i];
		point *p = s->getMiddlePoint();
		
		float x = (((float)(p->x))/NUM_FIN_EL)*2 - 1.0f;
		float y = ((((float)(p->y))/NUM_FIN_EL)*2 - 1.0f)*0.95f - 0.05f;
		float fi = PI/2 - s->getFi() - PI/(10);
		float cosfi = cos(fi);
		float sinfi = sin(fi);
		float x1, x2, x3, x4, y1, y2, y3, y4;
		//setup texture mapping
		glEnable( GL_TEXTURE_2D );
		if(s->isFocused()) {
			x1 = - 0.32*cosfi - 0.05*sinfi + x;
			y1 = + 0.32*sinfi - 0.05*cosfi + y;		
			x2 = 0.32*cosfi - 0.05*sinfi + x;
			y2 = -0.32*sinfi - 0.05*cosfi + y;			
			x3 = 0.32*cosfi + 0.2*sinfi + x;
			y3 = -0.32*sinfi + 0.2*cosfi + y;	
			x4 = -0.32*cosfi + 0.2*sinfi + x;
			y4 = 0.32*sinfi + 0.2*cosfi + y;				
			glBindTexture( GL_TEXTURE_2D, probefTexture);
		}
		else {
			x1 = - 0.32*cosfi - 0.025*sinfi + x;
			y1 = + 0.32*sinfi - 0.025*cosfi + y;		
			x2 = 0.32*cosfi - 0.025*sinfi + x;
			y2 = -0.32*sinfi - 0.025*cosfi + y;			
			x3 = 0.32*cosfi + 0.20*sinfi + x;
			y3 = -0.32*sinfi + 0.20*cosfi + y;	
			x4 = -0.32*cosfi + 0.20*sinfi + x;
			y4 = 0.32*sinfi + 0.20*cosfi + y;				
			glBindTexture( GL_TEXTURE_2D, probeuTexture);			
		}
		glColor3f(1.0f, 1.0f, 1.0f);	
		glPushMatrix();
		//glColor4f    (1.0f, 1.0f, 1.0f, 0.5);	
		glBegin( GL_QUADS );
		glTexCoord2d(0.0,0.0); glVertex2d(x1, y1);
		glTexCoord2d(1.0,0.0); glVertex2d(x2, y2);
		glTexCoord2d(1.0,1.0); glVertex2d(x3, y3);
		glTexCoord2d(0.0,1.0); glVertex2d(x4, y4);
		glEnd();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);	//disable texture so it doesn't mess up colors	
		
		glColor3f(0.0f, 0.0f, 0.8f);			
	
		glBegin( GL_QUADS );
			glVertex2d(x - 0.01, y - 0.01);
			glVertex2d(x + 0.01, y - 0.01);
			glVertex2d(x + 0.01, y + 0.01);
			glVertex2d(x - 0.01, y + 0.01);
		glEnd();			
		
	}
}

void drawSliders() {
	if(sourceId >= 0) {
		Probe *s = sources[sourceId];		
		int freq = s->getFrequency();
		int maxFrequency = s->getMaxFrequency();
		int minFrequency = s->getMinFrequency();
		float amplitude = s->getAmplitude();
		float maxAmplitude = s->getMaxAmplitude();
		float minAmplitude = s->getMinAmplitude();
		char ss[20];
		
		glColor3f(1.0, 1.0, 1.0);	
		glRasterPos2f(0.5, 0.835);
		sprintf(ss,"Frequency: %.2f MHz", freq*0.000001);
		int len = strlen(ss);
		for (int i = 0; i < len; i++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (int)ss[i]);
		}		
		
		glRasterPos2f(0.5, 0.735);
		sprintf(ss,"Amplitude: %.2f kPa", amplitude*0.001);
		len = strlen(ss);
		for (int i = 0; i < len; i++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (int)ss[i]);
		}
		
		float pos = (((float)(freq - minFrequency))/(maxFrequency - minFrequency))*0.45 + 0.5;
		frequencySliderPos = pos;
		
		//draw horizontal line
		glColor3f(0.6, 0.6, 0.6);		
		glBegin(GL_QUADS);
			glVertex2d(+0.5, +0.8);
			glVertex2d(+0.95, +0.8);
			glVertex2d(+0.95,+0.81);
			glVertex2d(+0.5,+0.81);			
		glEnd();
			
		
		//draw vertical line
		glColor3f(1.0, 1.0, 1.0);		
		glBegin(GL_QUADS);
			glVertex2d(pos - 0.01, +0.79);
			glVertex2d(pos + 0.01, +0.79);
			glVertex2d(pos + 0.01, +0.82);
			glVertex2d(pos - 0.01, +0.82);			
		glEnd();			
		
		pos = (((float)(amplitude - minAmplitude))/(maxAmplitude - minAmplitude))*0.45 + 0.5;
		ampSliderPos = pos;
		
		//draw horizontal line
		glColor3f(0.6, 0.6, 0.6);		
		glBegin(GL_QUADS);
			glVertex2d(+0.5, +0.7);
			glVertex2d(+0.95, +0.7);
			glVertex2d(+0.95,+0.71);
			glVertex2d(+0.5,+0.71);			
		glEnd();
		
		//draw vertical line
		glColor3f(1.0, 1.0, 1.0);		
		glBegin(GL_QUADS);
			glVertex2d(pos - 0.01, +0.69);
			glVertex2d(pos + 0.01, +0.69);
			glVertex2d(pos + 0.01, +0.72);
			glVertex2d(pos - 0.01, +0.72);			
		glEnd();		
	}
}
void drawTime() {
	char buffer[50];
	int n = sprintf(buffer, "Time: %.3f ms", t*1000);
	glColor3f(1.0, 1.0, 1.0);	
	glRasterPos2f(-0.97, 0.835);	
	for(int i = 0; i < n; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (int)buffer[i]);
	}
}

void drawPressure() {
	if(displayPressure != 0) {
		char buffer[50];
		int n = sprintf(buffer, "Acoustic pressure: %.2f kPa", displayPressure*0.001f);
		glColor3f(1.0, 1.0, 1.0);	
		glRasterPos2f(-0.97, 0.785);	
		for(int i = 0; i < n; i++) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (int)buffer[i]);
		}
	}
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//setup texture mapping
	glEnable( GL_TEXTURE_2D );
	glBindTexture( GL_TEXTURE_2D, texture );

	glPushMatrix();
	glBegin( GL_QUADS );
	glTexCoord2d(0.0,0.0); glVertex2d(-1.0,-1.0);
	glTexCoord2d(1.0,0.0); glVertex2d(+1.0,-1.0);
	glTexCoord2d(1.0,1.0); glVertex2d(+1.0,+0.9);
	glTexCoord2d(0.0,1.0); glVertex2d(-1.0,+0.9);
	glEnd();
	glPopMatrix();
	glDisable(GL_TEXTURE_2D);	//disable texture so it doesn't mess up colors

	int i = 0;
	for(int y = 0; y < NUM_FIN_EL; y++) {
		glBegin(GL_QUAD_STRIP);
		for(int x = 0; x < NUM_FIN_EL; x++) {
			float xpos = displayMatrixX[i];
			float ypos = displayMatrixY[i];
			float ydy = ypos + dy;
			
			glColor4f(GET_COLOR(u[i]));

			glVertex2f(xpos, ydy);		
			glVertex2f(xpos, ypos);			
			i++;
		}
		glEnd();
	}	
	drawProbes();		
	//draw sliders
	if(numSources > 0) {
		drawSliders();		
	}
	//draw menu
	drawMenu();		
	//draw probes
	drawTime();
	drawPressure();
	glutSwapBuffers();
}

void removeSelectedProbe() {
	if(numSources > 0) {
		free((sources[sourceId]));	
		if(sourceId < numSources - 1) {
			sources[sourceId] = sources[numSources - 1];
		}
		else {
			sourceId--;
		}
		numSources--;
	}	
}

void iterate() {
	if(newFocProbe) {
		addSource(NUM_FIN_EL/2, NUM_FIN_EL/2, 60, 1, 235.0f, PI, true);	
		newFocProbe = false;		
	}
	if(newUnfocProbe) {
		addSource(NUM_FIN_EL/2, NUM_FIN_EL/2, 60, 1, 235.0f, PI, false);	
		newUnfocProbe = false;
	}
	if(removeProbe) {
		removeSelectedProbe();
		removeProbe = false;		
	}
	
	
	
	//if pause is false, then play
	if(!pause) {		
		RK4(u, v, tempu, tempv);
		
		//switch u & tempu pointers
		temp = u;
		u = tempu;
		tempu = temp;

		//switch v & tempv pointers
		temp = v;
		v = tempv;
		tempv = temp;
		calculateFPS();
		
		//update t
		t += TS;
	}

	if(clock() - reset > 85) {
		reset = clock();
		glutPostRedisplay();
	}	
}
// set camera position and angle
void setCamera () {
	glMatrixMode(GL_PROJECTION);      // Select the Projection matrix for operation
	glLoadIdentity();                 // Reset Projection matrix
	gluOrtho2D(-1.0f, 1.0f, -1.0f, 1.0f);	
	
}
void initGL() {
	glEnable( GL_BLEND );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	
	glDisable(GL_LIGHTING);
	setCamera();
}

/**
This function is called when window is resized.
**/
void resize(int width, int height) {
    glutReshapeWindow(768, 768);
}

int hoveringOverSource(int x, int y) {
	for(int i = 0; i < numSources; i++) {	
		Probe *s = sources[i];
		point *p = s->getMiddlePoint();
		if(x < p->x + 10 && x > p->x - 10 && y < p->y + 10 && y > p->y - 10) {
			return i;
		}
	}
	return -1;
}

void clearSoundMatrices(bool init) {
	t = 0;
	displayPressure = 0;
	if(!init) {
		free(u);
		free(v);
		free(tempu);
		free(tempv);		
	}
	
	u = (float *)calloc(ms, sizeof(float));
	v = (float *)calloc(ms, sizeof(float));
	tempu = (float *)calloc(ms, sizeof(float));
	tempv = (float *)calloc(ms, sizeof(float));
}

void calculateMaxPressure(int x, int y) {
	float pressure = 0;
	if(!(x < NUM_FIN_EL - 10 && y < NUM_FIN_EL - 10 && x > 10 && y > 10)) {
		return;
	}
	
	for(int i = x - 10; i <  x + 10; i++) {
		for(int j = y - 10; j <  y + 10; j++) {
			if(abs(u[j*NUM_FIN_EL + i]) > pressure) {
				pressure = abs(u[j*NUM_FIN_EL + i]);
			}
		}
	}
	displayPressure = pressure;
	//printf("%f \n", displayPressure);
}
void pressRelease(int button, int state, int x, int y) {
	if(state == 0 && button != 1) {
		
		//check if pause/play button is pressed
		if( x > 0.005*glutGet(GLUT_WINDOW_WIDTH) && x < 0.075*glutGet(GLUT_WINDOW_WIDTH) && y > 0.005*glutGet(GLUT_WINDOW_HEIGHT) && y < 0.045*glutGet(GLUT_WINDOW_HEIGHT)) {
			pause = !pause; //toggle pause
			return;
		}
		//check if new focused probe button is pressed
		if( x > 0.08*glutGet(GLUT_WINDOW_WIDTH) && x < 0.315*glutGet(GLUT_WINDOW_WIDTH) && y > 0.005*glutGet(GLUT_WINDOW_HEIGHT) && y < 0.045*glutGet(GLUT_WINDOW_HEIGHT)) {
			newFocProbe = true;
			return;
		}
		//to do: check if new unfocused probe button is pressed 
		if( x > 0.32*glutGet(GLUT_WINDOW_WIDTH) && x < 0.58*glutGet(GLUT_WINDOW_WIDTH) && y > 0.005*glutGet(GLUT_WINDOW_HEIGHT) && y < 0.045*glutGet(GLUT_WINDOW_HEIGHT)) {
			newUnfocProbe = true;
			return;
		}		
		//check if remove button is pressed
		if( x > 0.585*glutGet(GLUT_WINDOW_WIDTH) && x < 0.765*glutGet(GLUT_WINDOW_WIDTH) && y > 0.005*glutGet(GLUT_WINDOW_HEIGHT) && y < 0.045*glutGet(GLUT_WINDOW_HEIGHT)) {
			removeProbe = true;
			return;
		}
		//check if clear button is pressed
		if( x > 0.77*glutGet(GLUT_WINDOW_WIDTH) && x < 0.81*glutGet(GLUT_WINDOW_WIDTH) && y > 0.005*glutGet(GLUT_WINDOW_HEIGHT) && y < 0.045*glutGet(GLUT_WINDOW_HEIGHT)) {
			clearSoundMatrices(false);
			return;
		}		
		float xGlut = (x*2.0f)/(glutGet(GLUT_WINDOW_WIDTH)) - 1.0f;
		float yGlut = ((y*2.0f)/(glutGet(GLUT_WINDOW_HEIGHT)) - 1.0f)*(-1);
		//check if frequencySlider is clicked
		if(xGlut <= frequencySliderPos + 0.01f && xGlut >= frequencySliderPos - 0.01f && yGlut >= 0.79 && yGlut <=  0.82) {
			holdingFrequencySlider = true;
			return;
		}
		
		//check if amplitude slider is clicked
		if(xGlut <= ampSliderPos + 0.01f && xGlut >= ampSliderPos - 0.01f && yGlut >= 0.69 && yGlut <=  0.72) {
			holdingAmpSlider = true;
			return;
		}
		
		x = ((float)x/(glutGet(GLUT_WINDOW_WIDTH)))*NUM_FIN_EL;
		y = NUM_FIN_EL - (((float)y - glutGet(GLUT_WINDOW_HEIGHT)*0.05)/(glutGet(GLUT_WINDOW_HEIGHT)*0.95))*NUM_FIN_EL;		
		
		//left or right mouse button is pressed
		int id = hoveringOverSource(x,y);
		if(id >= 0) {	
			holdingSource = true;
			roMode = true;
			sourceId = id;		
			if(button == 0) {
				roMode = false;
			}		
			return;
		}
		calculateMaxPressure(x,y);
	}
	else if(state == 1){
		//left mouse button is released
		holdingSource = false;
		holdingFrequencySlider = false;
		holdingAmpSlider = false;
		roMode = false;
		newFocProbe = false;
	}
}

void drag(int x, int y) {
	//drag callback is called
	if(holdingSource) {
		x = ((float)x/(glutGet(GLUT_WINDOW_WIDTH)))*NUM_FIN_EL;
		y = NUM_FIN_EL - (((float)y - glutGet(GLUT_WINDOW_HEIGHT)*0.05)/(glutGet(GLUT_WINDOW_HEIGHT)*0.95))*NUM_FIN_EL;		
		Probe *s = sources[sourceId];
		if(roMode) {
			s->setFi(x, y);
		}
		else {
			s->setPQ(x, y);		
		}
		s->setPointMesh(false);
	}
	else if(holdingAmpSlider) {
		Probe *s = sources[sourceId];
		float xGlut = (x*2.0f)/(glutGet(GLUT_WINDOW_WIDTH)) - 1.0f;
		float minAmp = s->getMinAmplitude();
		float maxAmp = s->getMaxAmplitude();
	
		float amp = ((xGlut - 0.5f)/(0.45f))*(maxAmp - minAmp) + minAmp;
		s->setAmplitude(amp);		
	}
	else if(holdingFrequencySlider) {
		Probe *s = sources[sourceId];
		float xGlut = (x*2.0f)/(glutGet(GLUT_WINDOW_WIDTH)) - 1.0f;
		int minFrequency = s->getMinFrequency();
		int maxFrequency = s->getMaxFrequency();
	
		float freq = ((xGlut - 0.5f)/(0.45f))*(maxFrequency - minFrequency) + minFrequency;
		s->setFrequency((int)freq);
	}
}

void keyboardPress(unsigned char key, int xpos, int ypos) {
	//if right button is pressed, save current values from u matrix horizontal and vertical
	if(key == 's') {		
		FILE * vertical;
		FILE * horizontal;
		
		xpos = ((float)xpos/(glutGet(GLUT_WINDOW_WIDTH)))*NUM_FIN_EL;
		ypos = NUM_FIN_EL - (((float)ypos - glutGet(GLUT_WINDOW_HEIGHT)*0.05)/(glutGet(GLUT_WINDOW_HEIGHT)*0.95))*NUM_FIN_EL;		
		
		printf("Saving vertical and horizontal cut to file %d, %d \n", xpos, ypos);	

		char filename[50];
		sprintf(filename, "./result_data/%d_vertical_acoustic_pressure.txt", fileSaveIndex);	
		vertical = fopen(filename, "w");
		sprintf(filename, "./result_data/%d_horizontal_acoustic_pressure.txt", fileSaveIndex);		
		horizontal = fopen(filename, "w");
		
		if(vertical != NULL) {
			for(int y = 0; y < NUM_FIN_EL; y++) {
				fprintf(vertical, "%f \n", u[y*NUM_FIN_EL + xpos]);
			}		
		}
		if(horizontal != NULL) {
			for(int x = 0; x < NUM_FIN_EL; x++) {
				fprintf(horizontal, "%f \n", u[ypos*NUM_FIN_EL + x]);
			}			
		}
		fclose(horizontal);
		fclose(vertical);
		fileSaveIndex++;
	}
}

int main(int argc, char *argv[]) {
	
	//set num sources to 0
	numSources = 0;

	//add new source
	addSource(NUM_FIN_EL/2, NUM_FIN_EL/2, 60, 1, 235.0f, PI, true);
	
	k1u = (float *)malloc(ms*sizeof(float));
	k2u = (float *)malloc(ms*sizeof(float));
	k3u = (float *)malloc(ms*sizeof(float));
	k4u = (float *)malloc(ms*sizeof(float));
	k1v = (float *)malloc(ms*sizeof(float));
	k2v = (float *)malloc(ms*sizeof(float));
	k3v = (float *)malloc(ms*sizeof(float));
	k4v = (float *)malloc(ms*sizeof(float));
	c_tissueS = (float *)calloc(ms, sizeof(float));
	displayMatrixX = (float *)malloc(ms*sizeof(float));
	displayMatrixY = (float *)malloc(ms*sizeof(float));
	damping = (int *)malloc(ms*sizeof(int));
	
	int i = 0;
	for(int y = 0; y  < NUM_FIN_EL; y++) {
		for(int x = 0; x < NUM_FIN_EL; x++) {
			displayMatrixX[i] = POS(x, (NUM_FIN_EL -1)); 
			displayMatrixY[i] = POS(y, (NUM_FIN_EL - 1))*0.95 - 0.05;
			i++;
		}		
	}
	

	dy = (displayMatrixY[1] - displayMatrixY[NUM_FIN_EL]);
	dx = (displayMatrixX[1] - displayMatrixX[0]);

	omp_set_num_threads(4);	
	
	clearSoundMatrices(true);
	
	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA); //float buffered glut with transparent colors
    glutInitWindowSize(768, 768);
    glutInitWindowPosition(100,100);
    glutCreateWindow("High intensity focused ultrasound");

	
	//load texture
	texture = setTexture("./head.bmp");
	setTissue("./segm_head.bmp");
	setprobefTexture("./probe_focused.bmp", "./probe_unfocused.bmp");
	
	
	initGL();
	
	glutMouseFunc(pressRelease);
	glutKeyboardFunc(keyboardPress);	
	glutMotionFunc(drag);	
    
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutIdleFunc(iterate);
	reset = clock();
    glutMainLoop();
	
	return 0;
}



