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

#define POS(X,Y) ((2.0f*X/Y) - 1.0f)
#define GET_COLOR(V) 0.45, V*V*0.000000001, 0.45, V*V*0.0000000005f
#define ALPHA(I) 	alfa = (I/20.0f); \
					alfa_1 = (1.0f - alfa); \
					alfac = alfa*alfa; \
					alfaca = alfa_1*alfa_1;

#define PI 3.14159265358979323846
#define TS 0.00000005 //time step
#define H 0.0003 //square size
#define NUM_FIN_EL 500 //number of finite elements in one dimension

struct point {
	//x position of point
	int x;
	//y position of point
	int y;
};

class Probe {
	private:
		bool focused;
		int numOfPoints;
		int p;
		int q;
		float amplitude;
		float maxAmplitude;
		float minAmplitude;
		float radius;
		float value;
		float fi;
		int frequency;
		int maxFrequency;
		int minFrequency;
		point * pointMesh;
		point * damperMesh;
		void setFocusedMesh(bool);
		void setUnfocusedMesh(bool);
	
	public:
		Probe(int, int, int, float, float, float, bool);
		void setPointMesh(bool);
		void updateValue(float);
		void setFi(int, int);
		float getFi();
		void setPQ(int, int);
		void setFrequency(int);
		void setAmplitude(float);
		point * getPoints();
		point * getDampedPoints();
		point *getMiddlePoint();
		int getNumOfPoints();
		float getValue();
		float getRadius();
		int getFrequency();
		int getMaxFrequency();
		int getMinFrequency();
		float getAmplitude();
		float getMaxAmplitude();
		float getMinAmplitude();
		bool isFocused();
};