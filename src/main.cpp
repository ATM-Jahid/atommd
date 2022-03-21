#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>

#include "types.hpp"
#include "vec_cal.hpp"

void setParams();
void initAtoms();
void rescaleVels();
void accumProps(int);
void singleStep();
void leapfrogStep(int);
void buildNebrList();
void computeForces();
void evalProps();
void printSummary();

// global variables
real rCut, density, temperature, deltaT, timeNow;
real uSum, virSum;
vecR cells, initUcell, region, velSum;
int nMol, nDim, stepCount, stepEquil, stepAdjTemp, stepLimit, stepAvg;
Prop kinEnergy, totEnergy, pressure;
Mol *mol;
int *cellList;
real dispHi, rNebrShell;
int *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax;
int num_atoms, cell_list = 0, neigh_list = 0;

int main() {
	// program start time
	auto start = std::chrono::system_clock::now();

	nDim = 3;
	rCut = 3;
	stepLimit = 10000;
	stepEquil = 5000;
	stepAdjTemp = 20;
	stepAvg = 50;
	deltaT = 0.001;

	// input from user
	std::cout << "Temperature: ";
	std::cin >> temperature;
	std::cout << "Density: ";
	std::cin >> density;
	std::cout << "Number of atoms: ";
	std::cin >> num_atoms;
	real num_unit_cell = int(std::pow(num_atoms/4, 1/3.0)+0.5);
	initUcell = {num_unit_cell, num_unit_cell, num_unit_cell};

	// which method to use
	std::cout << "Cell subdivision (0 or 1): ";
	std::cin >> cell_list;
	std::cout << "Neighbor list (0 or 1): ";
	std::cin >> neigh_list;

	// for neighbor list
	nebrTabFac = 100;
	rNebrShell = 0.4;
	nebrNow = 1;

	setParams();
	mol = new Mol[nMol];
	cellList = new int[int(vecProd(cells)+0.5) + nMol];
	nebrTab = new int[2*nebrTabMax];

	initAtoms();
	accumProps(0);

	for (stepCount = 0; stepCount < stepLimit; stepCount++) {
		singleStep();
	}

	delete[] mol;
	delete[] cellList;
	delete[] nebrTab;

	// program end time
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
	std::cout << "Wall time: " << elapsed.count() << " seconds\n";

	return 0;
}

void setParams() {
	vecScaleCopy(region, 1.0/std::pow(density/4.0, 1/3.0), initUcell);
	vecScaleCopy(cells, 1.0/rCut, region);
	vecRound(cells);
	nMol = 4 * int(vecProd(initUcell)+0.5);
	nebrTabMax = nebrTabFac * nMol;
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

	// execute this when neigh_list is on
	// and nebrNow is 1
	if (neigh_list && nebrNow) {
		nebrNow = 0;
		dispHi = 0;
		buildNebrList();
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

void buildNebrList() {
	vecR dr, invWid, rs, shift, cc, m1v, m2v;
	vecR vecOffset[] = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1}, {1,0,1},
			{1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, {-1,-1,1}, {0,-1,1}, {1,-1,1}};
	real rrNebr, rr;

	rrNebr = Sqr(rCut + rNebrShell);
	nebrTabLen = 0;

	if (cell_list) {
		/*
		 * CELL SUBDIVISION FOR NEIGHBOR LIST
		 */
		vecDiv(invWid, cells, region);
		// initialize the cell values to -1
		for (int i = nMol; i < nMol + vecProd(cells); i++) {
			cellList[i] = -1;
		}
	
		// make a linked list
		for (int i = 0; i < nMol; i++) {
			vecScaleAdd(rs, mol[i].r, 0.5, region);
			vecMul(cc, rs, invWid);
			vecFloor(cc);
			int c = vecLinear(cc, cells) + nMol;
			cellList[i] = cellList[c];
			cellList[c] = i;
		}
	
		for (int m1z = 0; m1z < cells.z; m1z++) {
			for (int m1y = 0; m1y < cells.y; m1y++) {
				for (int m1x = 0; m1x < cells.x; m1x++) {
					vecSet(m1v, m1x, m1y, m1z);
					int m1 = vecLinear(m1v, cells) + nMol;
					for (int Noff = 0; Noff < 14; Noff++) {
						vecAdd(m2v, m1v, vecOffset[Noff]);
						vecSet(shift, 0, 0, 0);
						cellWrapAll(m2v, shift, cells, region);
						int m2 = vecLinear(m2v, cells) + nMol;
						for (int j1 = cellList[m1]; j1 > -1; j1 = cellList[j1]) {
							for (int j2 = cellList[m2]; j2 > -1; j2 = cellList[j2]) {
								if (m1 != m2 || j1 > j2) {
									vecSub(dr, mol[j1].r, mol[j2].r);
									vecSub(dr, dr, shift);
									rr = vecLenSq(dr);
									if (rr < rrNebr) {
										if (nebrTabLen >= nebrTabMax) {
											std::cout << "too many neighbors!\n";
											exit(0);
										}
										nebrTab[2*nebrTabLen] = j1;
										nebrTab[2*nebrTabLen+1] = j2;
										nebrTabLen++;
									}
								}
							}
						}
					}
				}
			}
		}
	} else {
		/*
		 * ONLY NEIGHBOR LIST
		 */
		for (int j1 = 0; j1 < nMol - 1; j1++) {
			for (int j2 = j1+1; j2 < nMol; j2++) {
				vecSub(dr, mol[j1].r, mol[j2].r);
				vecWrapAll(dr, region);
				rr = vecLenSq(dr);
				if (rr < rrNebr) {
					if (nebrTabLen >= nebrTabMax) {
						std::cout << "too many neighbors!\n";
						exit(0);
					}
					nebrTab[2*nebrTabLen] = j1;
					nebrTab[2*nebrTabLen+1] = j2;
					nebrTabLen++;
				}
			}
		}
	}
}

void computeForces() {
	vecR dr, invWid, rs, shift, cc, m1v, m2v;
	vecR vecOffset[] = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1}, {1,0,1},
			{1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, {-1,-1,1}, {0,-1,1}, {1,-1,1}};
	real fcVal, rr, rrCut, rri, rri3;

	rrCut = Sqr(rCut);
	// resetting the acc. values since they are incremented later on
	for (int i = 0; i < nMol; i++) {
		vecSet(mol[i].acc, 0, 0, 0);
	}
	uSum = 0;
	virSum = 0;

	if (neigh_list) {
		/*
		 * NEIGHBOR LIST
		 */
		for (int i = 0; i < nebrTabLen; i++) {
			int j1 = nebrTab[2*i];
			int j2 = nebrTab[2*i+1];
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
	} else if (cell_list) {
		/*
		 * CELL SUBDIVISION
		 */
		vecDiv(invWid, cells, region);
		// initialize the cell values to -1
		for (int i = nMol; i < nMol + vecProd(cells); i++) {
			cellList[i] = -1;
		}
	
		// make a linked list
		for (int i = 0; i < nMol; i++) {
			vecScaleAdd(rs, mol[i].r, 0.5, region);
			vecMul(cc, rs, invWid);
			vecFloor(cc);
			int c = vecLinear(cc, cells) + nMol;
			cellList[i] = cellList[c];
			cellList[c] = i;
		}
	
		for (int m1z = 0; m1z < cells.z; m1z++) {
			for (int m1y = 0; m1y < cells.y; m1y++) {
				for (int m1x = 0; m1x < cells.x; m1x++) {
					vecSet(m1v, m1x, m1y, m1z);
					int m1 = vecLinear(m1v, cells) + nMol;
					for (int Noff = 0; Noff < 14; Noff++) {
						vecAdd(m2v, m1v, vecOffset[Noff]);
						vecSet(shift, 0, 0, 0);
						cellWrapAll(m2v, shift, cells, region);
						int m2 = vecLinear(m2v, cells) + nMol;
						for (int j1 = cellList[m1]; j1 > -1; j1 = cellList[j1]) {
							for (int j2 = cellList[m2]; j2 > -1; j2 = cellList[j2]) {
								if (m1 != m2 || j1 > j2) {
									vecSub(dr, mol[j1].r, mol[j2].r);
									vecSub(dr, dr, shift);
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
					}
				}
			}
		}
	} else {
		/*
		 * ALL PAIRS
		 */
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
}

void evalProps() {
	vecSet(velSum, 0, 0, 0);
	real v2, v2sum = 0, v2max = 0;

	for (int i = 0; i < nMol; i++) {
		vecAdd(velSum, velSum, mol[i].vel);
		v2 = vecLenSq(mol[i].vel);
		v2sum += v2;
		v2max = std::max(v2max, v2);
	}

	kinEnergy.val = 0.5 * v2sum / nMol;
	totEnergy.val = kinEnergy.val + uSum / nMol;
	pressure.val = density * (v2sum + virSum) / (nMol * nDim);

	dispHi += std::sqrt(v2max) * deltaT;
	if (dispHi > 0.5 * rNebrShell) {
		nebrNow = 1;
	}
}

void printSummary() {
	std::cout << stepCount << '\t' << timeNow << '\t'
		<< std::sqrt(vecLenSq(velSum))/nMol << '\t'
		<< kinEnergy.sum << '\t' << totEnergy.sum << '\t'
		<< pressure.sum << '\n';
}
