#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "types.hpp"
#include "vec_cal.hpp"

void setParams();
void initAtoms();
void accumProps(int);
void singleStep();
void leapfrogStep(int);
void computeForces();
void evalProps();
void printSummary();

// global variables
real rCut, density, velMag, temperature, deltaT, timeNow;
real uSum, virSum;
vecR cells, initUcell, region, velSum;
int nMol, nDim, stepCount, stepLimit, stepAvg;
Prop kinEnergy, totEnergy, pressure;
Mol *mol;
int *cellList;

int main(int argc, char **argv) {
	//getInputs(argc, argv);

	nDim = 3;
	// input temp. and density from user
	cin >> temperature >> density;
	initUcell = {.x = 20, .y = 20};
	stepLimit = 10000;
	stepAvg = 100;
	deltaT = 0.001;

	setParams();
	mol = new Mol[nMol];
	cellList = new int[vecProd(cells) + nMol];

	initAtoms();
	accumProps(0);

	for (stepCount = 0; stepCount < stepLimit; stepCount++) {
		singleStep();
	}

	delete[] mol;
	delete[] cellList;
	return 0;
}

void setParams() {
	rCut = pow(2, 1/6.0);
	vecScaleCopy(region, 1.0/sqrt(density), initUcell);
	nMol = vecProd(initUcell);
	velMag = sqrt(nDim * (1.0 - 1.0/nMol) * temperature);
}

void initAtoms() {
	vecR c, gap;
	int n = 0;
	vecDiv(gap, region, initUcell);
	for (int ny = 0; ny < initUcell.y; ny++) {
		for (int nx = 0; nx < initUcell.x; nx++) {
			vecSet(c, nx+0.5, ny+0.5);
			vecMul(c, c, gap);
			vecScaleAdd(c, c, -0.5, region);
			mol[n].r = c;
			n++;
		}
	}

	vecSet(velSum, 0, 0);
	for (int i = 0; i < nMol; i++) {
		// radom velocity generator
		real sdf = (rand() % 1000) / 999.9 * 6.283185;
		vecSet(mol[i].vel, cos(sdf), sin(sdf));
		vecScale(mol[i].vel, velMag);
		vecAdd(velSum, velSum, mol[i].vel);
		// assign zero init. acceleration
		vecSet(mol[i].acc, 0, 0);
	}

	// account for COM shift
	for (int i = 0; i < nMol; i++) {
		vecScaleAdd(mol[i].vel, mol[i].vel, -1.0/nMol, velSum);
	}
}

void accumProps(int icode) {
	switch (icode) {
		case 0:
			propZero(kinEnergy);
			propZero(totEnergy);
			propZero(pressure);
			break;
		case 1:
			propAccum(kinEnergy);
			propAccum(totEnergy);
			propAccum(pressure);
			break;
		case 2:
			propAvg(kinEnergy, stepAvg);
			propAvg(totEnergy, stepAvg);
			propAvg(pressure, stepAvg);
			break;
	}
}

void singleStep() {
	timeNow = stepCount * deltaT;

	leapfrogStep(1);
	// apply boundary conditions
	for (int i = 0; i < nMol; i++) {
		vecWrapAll(mol[i].r, region);
	}

	computeForces();
	leapfrogStep(2);
	evalProps();
	accumProps(1);

	if (stepCount % stepAvg == 0) {
		accumProps(2);
		printSummary();
		accumProps(0);
	}
}

void leapfrogStep(int part) {
	if (part == 1) {
		for (int i = 0; i < nMol; i++) {
			vecScaleAdd(mol[i].vel, mol[i].vel, 0.5*deltaT, mol[i].acc);
			vecScaleAdd(mol[i].r, mol[i].r, deltaT, mol[i].vel);
		}
	} else {
		for (int i = 0; i < nMol; i++) {
			vecScaleAdd(mol[i].vel, mol[i].vel, 0.5*deltaT, mol[i].acc);
		}
	}
}

void computeForces() {
	vecR dr;
	real fcVal, rr, rrCut, rri, rri3;

	rrCut = Sqr(rCut);
	for (int i = 0; i < nMol; i++) {
		vecSet(mol[i].acc, 0, 0);
	}

	uSum = 0;
	virSum = 0;
	for (int j1 = 0; j1 < nMol - 1; j1++) {
		for (int j2 = j1+1; j2 < nMol; j2++) {
			vecSub(dr, mol[j1].r, mol[j2].r);
			vecWrapAll(dr, region);
			rr = vecLenSq(dr);
			if (rr < rrCut) {
				rri = 1.0 / rr;
				rri3 = Cub(rri);
				fcVal = 48.0 * rri3 * (rri3 - 0.5) * rri;
				vecScaleAdd(mol[j1].acc, mol[j1].acc, fcVal, dr);
				vecScaleAdd(mol[j2].acc, mol[j2].acc, -fcVal, dr);
				uSum += 4.0 * rri3 * (rri3 - 1.0);
				virSum += fcVal * rr;
			}
		}
	}
}

void evalProps() {
	vecSet(velSum, 0, 0);
	real v2sum = 0;

	for (int i = 0; i < nMol; i++) {
		vecAdd(velSum, velSum, mol[i].vel);
		v2sum += vecLenSq(mol[i].vel);
	}

	kinEnergy.val = 0.5 * v2sum / nMol;
	totEnergy.val = kinEnergy.val + uSum / nMol;
	pressure.val = density * (v2sum + virSum) / (nMol * nDim);
}

void printSummary() {
	std::cout << stepCount << '\t' << timeNow << '\t'
		<< velSum.x/nMol << '\t' << velSum.y/nMol << '\t'
		<< kinEnergy.sum << '\t' << totEnergy.sum << '\t'
		<< pressure.sum << '\n';
}
