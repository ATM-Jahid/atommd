#include <iostream>
#include <math.h>

#include "types.hpp"
#include "vec_cal.hpp"

void setParams(real &, vecR &, int &, real &, real, vecR, int, real);
void initAtoms(Mol *, vecR, vecR, vecR &, int, real);
void accumProps(int, Prop &, Prop &, Prop &, int);

int main(int argc, char **argv) {
	//getInputs(argc, argv);
	//printOutputs(stdout);

	real rCut, density, velMag, temperature;
	vecR initUcell, region, velSum;
	int nMol, nDim, stepCount, stepAvg = 10;

	rCut = 2;
	density = 0.8;
	temperature = 1;
	initUcell = {.x = 20, .y = 20};
	nDim = 2;

	Prop kinEnergy, totEnergy, pressure;
	kinEnergy.val = 0;
	totEnergy.val = 0;
	pressure.val = 0;

	setParams(rCut, region, nMol, velMag,
		density, initUcell, nDim, temperature);

	Mol *mol;
	mol = new Mol[nMol];

	initAtoms(mol, initUcell, region, velSum, nMol, velMag);
	accumProps(0, kinEnergy, totEnergy, pressure, stepAvg);


	delete[] mol;
	return 0;
}

void setParams(real &rCut, vecR &region, int &nMol, real &velMag,
		real density, vecR initUcell, int nDim, real temperature) {

	rCut = pow(2, 1/6.0);
	vecScaleCopy(region, 1/sqrt(density), initUcell);
	nMol = vecProd(initUcell);
	velMag = sqrt(nDim * (1 - 1.0/nMol) * temperature);
}

void initAtoms(Mol *mol, vecR initUcell, vecR region, vecR &velSum, int nMol, real velMag) {
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

void accumProps(int icode, Prop &kinEnergy, Prop &totEnergy, Prop &pressure, int stepAvg) {
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
