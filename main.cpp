#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
using namespace std;

double myu = 10., lyambda = 1.0, k = 0.5, eps = pow(10., -17.), constH = 0.002, epsHalfH = pow(10., -15.);
double h = constH;
double start = 0., checkpoint;
double x[3], y[3], wholeX[3], wholeY[3], halfX[3], halfY[3];

double precFunc(double x) {
	return sin(k * x);
}
double myFunc(double x, double y) {
	return lyambda * y - lyambda * sin(k * x) + k * cos(k * x);
}
double mySecondFunc(double x, double y) {
	return myu * (y - sin(k * x)) / (x + 1) + k * cos(k * x);
}
ofstream fout1("values1.txt");
ofstream fout2("values2.txt");
ofstream fout3("values3.txt");
ofstream fout4("values4.txt");
ofstream foutH1("step1.txt");
ofstream foutH2("step2.txt");

void solveWithConstantH(ofstream &fout, double (*func)(double, double)) {
	x[0] = start;
	x[1] = x[0] + h;
	x[2] = x[1] + h;
	y[0] = precFunc(x[0]);
	y[1] = precFunc(x[1]);
	y[2] = y[1] + h * func(x[1], y[1]);
	double constPart = y[0] + (h / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
	do {
		checkpoint = y[2];
		y[2] = constPart + (h / 3.) * func(x[2], y[2]);
	} while (abs(1 - y[2] / checkpoint) > eps);
	fout << "{" << setprecision(11) << x[0] << "," << setprecision(11) << y[0] << "}, {" << setprecision(11) << x[1] << "," << setprecision(11) << y[1] << "}, {" << setprecision(11) << x[2] << "," << setprecision(11) << y[2] << "},";
	while (abs(y[2]) < 2) {
		x[0] = x[1];
		x[1] = x[0] + h;
		x[2] = x[1] + h;
		y[0] = y[1]; 
		y[1] = y[2];
		y[2] = y[1] + h * func(x[1], y[1]);
		constPart = y[0] + (h / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
		do {
			checkpoint = y[2];
			y[2] = constPart + (h / 3.) * func(x[2], y[2]);
		} while (abs(1 - y[2] / checkpoint) > eps);
		fout << "{" << setprecision(11) << x[2] << "," << setprecision(11) << y[2] << "},";
	}
}

void solveWithHalfOrDoubleH(ofstream &fout, ofstream &foutH, double (*func)(double, double)) {
	bool checker = false;
	double lastY = 0.;
	double halfH = h / 2;
	x[0] = start;
	x[1] = x[0] + h;
	x[2] = x[1] + h;
	y[0] = precFunc(x[0]);
	y[1] = precFunc(x[1]);
	y[2] = y[1] + h * func(x[1], y[1]);
	double constPart = y[0] + (h / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
	do {
		checkpoint = y[2];
		y[2] = constPart + (h / 3.) * func(x[2], y[2]);
	} while (abs(1 - y[2] / checkpoint) > eps);
	fout << "{" << setprecision(11) << x[0] << "," << setprecision(11) << y[0] << "}, {" << setprecision(11) << x[1] << "," << setprecision(11) << y[1] << "}, ";
	wholeX[0] = x[0];
	wholeX[1] = x[1];
	wholeX[2] = x[2];
	wholeY[0] = y[0];
	wholeY[1] = y[1];
	wholeY[2] = y[2];

	x[0] = (wholeX[1] + wholeX[0]) / 2;
	x[1] = x[0] + halfH;
	x[2] = x[1] + halfH;
	y[0] = precFunc(x[0]);
	y[1] = precFunc(x[1]);
	y[2] = y[1] + halfH * func(x[1], y[1]);
	constPart = y[0] + (halfH / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
	do {
		checkpoint = y[2];
		y[2] = constPart + (halfH / 3.) * func(x[2], y[2]);
	} while (abs(1 - y[2] / checkpoint) > eps);

	halfX[0] = x[0];
	halfX[1] = x[1];
	halfX[2] = x[2];
	halfY[0] = y[0];
	halfY[1] = y[1];
	halfY[2] = y[2];

	x[0] = halfX[1];
	x[1] = halfX[2];
	x[2] = x[1] + halfH;
	y[0] = halfY[1];
	y[1] = halfY[2];
	y[2] = y[1] + halfH * func(x[1], y[1]);
	constPart = y[0] + (halfH / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
	do {
		checkpoint = y[2];
		y[2] = constPart + (halfH / 3.) * func(x[2], y[2]);
	} while (abs(1 - y[2] / checkpoint) > eps);


	double tempHalfY = 0.;
	while (abs(wholeY[2]) < 2) {
		if (abs(16 / 15 * (y[2] - wholeY[2])) < epsHalfH) {
			tempHalfY = halfY[2];
			wholeY[2] += 16 / 15 * (y[2] - wholeY[2]);
			fout << "{" << setprecision(11) << wholeX[2] << "," << setprecision(11) << wholeY[2] << "},";
			foutH << "{" << setprecision(11) << wholeX[2] << "," << setprecision(11) << h << "},";

			h *= 2;
			halfH *= 2;
			halfX[0] = wholeX[1];
			halfX[1] = wholeX[2];
			halfX[2] = halfX[1] + halfH;
			halfY[0] = wholeY[1];
			halfY[1] = wholeY[2];

			halfY[2] = halfY[1] + halfH * func(halfX[1], halfY[1]);
			constPart = halfY[0] + (halfH / 3.) * (func(halfX[0], halfY[0]) + 4 * func(halfX[1], halfY[1]));
			do {
				checkpoint = halfY[2];
				halfY[2] = constPart + (halfH / 3.) * func(halfX[2], halfY[2]);
			} while (abs(1 - halfY[2] / checkpoint) > eps);

			wholeX[1] = wholeX[2];
			wholeX[2] = wholeX[1] + h;
			wholeY[1] = wholeY[2];
			wholeY[2] = wholeY[1] + h * func(wholeX[1], wholeY[1]);
			constPart = wholeY[0] + (h / 3.) * (func(wholeX[0], wholeY[0]) + 4 * func(wholeX[1], wholeY[1]));
			do {
				checkpoint = wholeY[2];
				wholeY[2] = constPart + (h / 3.) * func(wholeX[2], wholeY[2]);
			} while (abs(1 - wholeY[2] / checkpoint) > eps);
			lastY = tempHalfY;
			
			x[0] = halfX[1];
			x[1] = halfX[2];
			x[2] = x[1] + halfH;
			y[0] = halfY[1];
			y[1] = halfY[2];
			y[2] = y[1] + halfH * func(x[1], y[1]);
			constPart = y[0] + (halfH / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
			do {
				checkpoint = y[2];
				y[2] = constPart + (halfH / 3.) * func(x[2], y[2]);
			} while (abs(1 - y[2] / checkpoint) > eps);
			checker = true;
		}
		else {
			double neededX = wholeX[0];
			double neededY = wholeY[0];
			h /= 2;
			halfH /= 2;
			wholeX[0] = halfX[0];
			wholeX[1] = halfX[1];
			wholeX[2] = halfX[2];
			wholeY[0] = halfY[0];
			wholeY[1] = halfY[1];
			wholeY[2] = halfY[2];

			halfX[0] = (wholeX[1] + wholeX[0]) / 2;
			halfX[1] = halfX[0] + halfH;
			halfX[2] = halfX[1] + halfH;
			if (halfX[1] <= constH) {
				halfY[0] = precFunc(halfX[0]);
			}
			else {
				if (checker == true) {
					halfY[0] = lastY;
				}
				else {
					halfY[0] = ((11 * neededY + 72 * wholeY[0] + 45 * wholeY[1]) / 128.) + ((3. / 64.) * halfH * (func(neededX, neededY) + 12 * func(wholeX[0], wholeY[0]) - 3 * func(wholeX[1], wholeY[1])));
				}
			}
			halfY[1] = wholeY[1];
			halfY[2] = halfY[1] + halfH * func(halfX[1], halfY[1]);
			constPart = halfY[0] + (halfH / 3.) * (func(halfX[0], halfY[0]) + 4 * func(halfX[1], halfY[1]));
			do {
				checkpoint = halfY[2];
				halfY[2] = constPart + (halfH / 3.) * func(halfX[2], halfY[2]);
			} while (abs(1 - halfY[2] / checkpoint) > eps);

			x[0] = halfX[1];
			x[1] = halfX[2];
			x[2] = x[1] + halfH;
			y[0] = halfY[1];
			y[1] = halfY[2];
			y[2] = y[1] + halfH * func(x[1], y[1]);
			constPart = y[0] + (halfH / 3.) * (func(x[0], y[0]) + 4 * func(x[1], y[1]));
			do {
				checkpoint = y[2];
				y[2] = constPart + (halfH / 3.) * func(x[2], y[2]);
			} while (abs(1 - y[2] / checkpoint) > eps);
			checker = false;
		}
	}
}

int main() {
	fout1 << fixed << "ListPlot[{";
	solveWithConstantH(fout1,myFunc);
	fout1 << "}, Joined->True, PlotRange->All, PlotStyle->Green]";
	fout1.close();
	fout2 << fixed << "ListPlot[{";
	foutH1 << fixed << "ListPlot[{";
	solveWithHalfOrDoubleH(fout2, foutH1, myFunc);
	fout2 << "}, Joined->True, PlotRange->All, PlotStyle->Purple]";
	foutH1 << "}, Joined->True, PlotRange->All]";
	fout2.close();
	foutH1.close(); 
	fout3 << fixed << "ListPlot[{";
	solveWithConstantH(fout3, mySecondFunc);
	fout3 << "}, Joined->True, PlotRange->All, PlotStyle->Green]";
	fout3.close();
	fout4 << fixed << "ListPlot[{";
	foutH2 << fixed << "ListPlot[{";
	solveWithHalfOrDoubleH(fout4, foutH2, mySecondFunc);
	fout4 << "}, Joined->True, PlotRange->All, PlotStyle->Purple]";
	foutH2 << "}, Joined->True, PlotRange->All]";
	fout4.close();
	foutH2.close();
	return 0;
}