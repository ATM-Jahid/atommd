#include <iostream>
#include <cmath>
#include <random>

#include "types.hpp"
#include "vec_cal.hpp"

void setParams();
void initAtoms();
void rescaleVels();
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
int nMol, nDim, stepCount, stepEquil, stepAdjTemp, stepLimit, stepAvg;
Prop kinEnergy, totEnergy, pressure;
Mol *mol;
int *cellList;

int main(int argc, char **argv) {
	//getInputs(argc, argv);

	nDim = 3;
	// input temp. and density from user
	std::cout << "Temperature: ";
	std::cin >> temperature;
	std::cout << "Density: ";
	std::cin >> density;
	rCut = 3;
	initUcell = {5, 5, 5};
	stepLimit = 5000;
	stepEquil = 1000;
	stepAdjTemp = 25;
	stepAvg = 50;
	deltaT = 0.001;

	setParams();
	mol = new Mol[nMol];
	cellList = new int[int(vecProd(cells)) + nMol];

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
	vecScaleCopy(region, 1.0/std::pow(density/4.0, 1/3.0), initUcell);
	nMol = 4 * vecProd(initUcell);
	//velMag = std::sqrt(nDim * (1.0 - 1.0/nMol) * temperature);
}

void initAtoms() {
	vecR c, gap;
	int n = 0;
	vecDiv(gap, region, initUcell);
	for (int nz = 0; nz < initUcell.z; nz++) {
		for (int ny = 0; ny < initUcell.y; ny++) {
			for (int nx = 0; nx < initUcell.x; nx++) {
				vecSet(c, nx+0.25, ny+0.25, nz+0.25);
				vecMul(c, c, gap);
				vecScaleAdd(c, c, -0.5, region);
				for (int j = 0; j < 4; j++) {
					mol[n].r = c;
					switch (j) {
						case 0:
							mol[n].r.x += 0.5 * gap.x;
							mol[n].r.y += 0.5 * gap.y;
							break;
						case 1:
							mol[n].r.y += 0.5 * gap.y;
							mol[n].r.z += 0.5 * gap.z;
							break;
						case 2:
							mol[n].r.z += 0.5 * gap.z;
							mol[n].r.x += 0.5 * gap.x;
							break;
					}
					n++;
				}
			}
		}
	}

	// radom velocity generator
	std::default_random_engine rand_gen;
	std::normal_distribution<real> normal_dist(0.0, 1.0);

	vecSet(velSum, 0, 0, 0);
	for (int i = 0; i < nMol; i++) {
		vecSet(mol[i].vel, normal_dist(rand_gen),
				normal_dist(rand_gen), normal_dist(rand_gen));
		//vecScale(mol[i].vel, velMag);
		vecAdd(velSum, velSum, mol[i].vel);
		// assign zero init. acceleration
		vecSet(mol[i].acc, 0, 0, 0);
	}

	// account for COM shift
	for (int i = 0; i < nMol; i++) {
		vecScaleAdd(mol[i].vel, mol[i].vel, -1.0/nMol, velSum);
	}

	// adjust temperature
	rescaleVels();
}

void rescaleVels() {
	real velSqSum = 0;
	for (int i = 0; i < nMol; i++) {
		velSqSum += vecLenSq(mol[i].vel);
	}

	real lambda = std::sqrt(3 * (nMol - 1) * temperature / velSqSum);
	for (int i = 0; i < nMol; i++) {
		vecScale(mol[i].vel, lambda);
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

	// rescale velocities
	if ((stepCount < stepEquil) && !(stepCount % stepAdjTemp)) {
		rescaleVels();
	}

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
		vecSet(mol[i].acc, 0, 0, 0);
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
	vecSet(velSum, 0, 0, 0);
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
		<< std::sqrt(vecLenSq(velSum))/nMol << '\t'
		<< kinEnergy.sum << '\t' << totEnergy.sum << '\t'
		<< pressure.sum << '\n';
}
