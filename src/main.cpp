#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>
#include <string>
#include <fstream>

#include "types.hpp"
#include "vec_cal.hpp"

void setParams();
void initAtoms();
void rescaleVels();
void accumProps(int);
void singleStep(std::string);
void leapfrogStep(int);
void buildNebrList();
void computeForces();
void evalProps();
void printSummary(std::string);
void posDump(std::string);
void evalRdf_AB(std::string);
void printRdf_AB(std::string);
void evalLatticeCorr();
void initDiffusion();
void zeroDiffusion();
void evalDiffusion(std::string);
void accumDiffusion(std::string);
void printMsd(std::string);
void printDiffusion(std::string);
void initVacf();
void zeroVacf();
void evalVacf(std::string);
void accumVacf(std::string);
double integrate(double *, int);
void printVacf(std::string);

// global variables
double rCut, density, temperature, deltaT, timeNow;
double uSum, virSum;
vecR cells, initUcell, region, momSum;
int nDim, nMol, nMolA, nMolB;
int stepCount, stepEquil, stepRun, stepLimit;
int stepAdjTemp, stepAvg, stepDump;
Prop kinEnergy, totEnergy, pressure;
Mol *mol;
int *cellList;
double dispHi, rNebrShell;
int *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax;
int num_atoms, cell_list = 1, neigh_list = 1;
double *histRdfAA, *histRdfBB, *histRdfAB, rangeRdf;
int countRdf, limitRdf, sizeHistRdf, stepRdf;
double latticeCorr;
Tbuff *bufferAA, *bufferBB, *bufferAB;
double *rrDiffAvgAA, *rrDiffAvgBB, *rrDiffAvgAB;
int countDiffAvg, limitDiffAvg, nBuffDiff, nValDiff, stepDiff;
Vbuff *vacBuff;
double *avgAcfVel, intAcfVel;
int countAcfAvg, limitAcfAvg, nBuffAcf, nValAcf, stepAcf;
double nAlpha, nBeta, mass1, mass2, mRatio, Q;
double eps, epsAA = 1.0, epsBB = 0.50, epsAB = 1.5;
double sig, sigAA = 1.0, sigBB = 0.88, sigAB = 0.8;

int main(int argc, char **argv) {
	// program start time
	auto start = std::chrono::system_clock::now();

	// process i/o files
	std::string dot_in(argv[1]);

	nDim = 3;
	rCut = 3;
	stepEquil = 10000;
	stepRun = 10000;
	stepLimit = stepEquil + stepRun;
	stepAdjTemp = 20;
	stepAvg = 50;
	stepDump = 100;

	// rdf parameters
	limitRdf = 200;
	stepRdf = stepRun / limitRdf;
	rangeRdf = 4;
	sizeHistRdf = 200;

	// diffusivity parameters
	stepDiff = 10;
	nValDiff = 500;
	nBuffDiff = 50;
	limitDiffAvg = (stepRun/stepDiff/nValDiff - 1) * nBuffDiff;

	// VACF parameters
	stepAcf = stepDiff;
	nValAcf = nValDiff;
	nBuffAcf = nBuffDiff;
	limitAcfAvg = limitDiffAvg;

	// input from user
	std::ifstream inputFile(dot_in);
	inputFile >> temperature >> density >> num_atoms >> mRatio >> deltaT;
	double num_unit_cell = int(std::pow(num_atoms/4, 1/3.0)+0.5);
	initUcell = {num_unit_cell, num_unit_cell, num_unit_cell};

	// for neighbor list
	nebrTabFac = 100;
	rNebrShell = 0.4;
	nebrNow = 1;

	setParams();
	mol = new Mol[nMol];
	cellList = new int[int(vecProd(cells)+0.5) + nMol];
	nebrTab = new int[2*nebrTabMax];
	histRdfAA = new double[sizeHistRdf];
	histRdfBB = new double[sizeHistRdf];
	histRdfAB = new double[sizeHistRdf];
	rrDiffAvgAA = new double[nValDiff];
	rrDiffAvgBB = new double[nValDiff];
	rrDiffAvgAB = new double[nValDiff];
	bufferAA = new Tbuff[nBuffDiff];
	bufferBB = new Tbuff[nBuffDiff];
	bufferAB = new Tbuff[nBuffDiff];
	for (int nb = 0; nb < nBuffDiff; nb++) {
		bufferAA[nb].orgR = new vecR[nMol];
		bufferAA[nb].rTrue = new vecR[nMol];
		bufferAA[nb].rrDiff = new double[nValDiff];
		bufferBB[nb].orgR = new vecR[nMol];
		bufferBB[nb].rTrue = new vecR[nMol];
		bufferBB[nb].rrDiff = new double[nValDiff];
		bufferAB[nb].orgR = new vecR[nMol];
		bufferAB[nb].rTrue = new vecR[nMol];
		bufferAB[nb].rrDiff = new double[nValDiff];
	}
	avgAcfVel = new double[nValAcf];
	vacBuff = new Vbuff[nBuffAcf];
	for (int nb = 0; nb < nBuffAcf; nb++) {
		vacBuff[nb].acfVel = new double[nValAcf];
		vacBuff[nb].orgVel = new vecR[nMol];
	}

	countRdf = 0;
	initAtoms();
	accumProps(0);
	initDiffusion();
	initVacf();

	for (stepCount = 0; stepCount < stepLimit; stepCount++) {
		singleStep(dot_in);
	}

	delete[] mol;
	delete[] cellList;
	delete[] nebrTab;
	delete[] histRdfAA;
	delete[] histRdfBB;
	delete[] histRdfAB;
	delete[] rrDiffAvgAA;
	delete[] rrDiffAvgBB;
	delete[] rrDiffAvgAB;
	for (int nb = 0; nb < nBuffDiff; nb++) {
		delete[] bufferAA[nb].orgR;
		delete[] bufferAA[nb].rTrue;
		delete[] bufferAA[nb].rrDiff;
		delete[] bufferBB[nb].orgR;
		delete[] bufferBB[nb].rTrue;
		delete[] bufferBB[nb].rrDiff;
		delete[] bufferAB[nb].orgR;
		delete[] bufferAB[nb].rTrue;
		delete[] bufferAB[nb].rrDiff;
	}
	delete[] bufferAA;
	delete[] bufferBB;
	delete[] bufferAB;
	delete[] avgAcfVel;
	for (int nb = 0; nb < nBuffAcf; nb++) {
		delete[] vacBuff[nb].acfVel;
		delete[] vacBuff[nb].orgVel;
	}
	delete[] vacBuff;

	// program end time
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
	std::string dot_out = dot_in.erase(dot_in.length()-2).append("out");
	std::ofstream outputFile;
	outputFile.open(dot_out, std::ofstream::app);
	outputFile << "Wall time: " << elapsed.count() << " seconds\n";
	outputFile.close();

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
	mass2 = 1.0;
	mass1 = mass2 * mRatio;

	nMolA = 0;
	nMolB = 0;
	for (int n = 0; n < nMol; n++) {
		if (n % 5 == 0) {
			mol[n].type = 2;
			mol[n].mass = mass2;
			nMolB++;
		} else {
			mol[n].type = 1;
			mol[n].mass = mass1;
			nMolA++;
		}
	}
	nAlpha = nMolA / double(nMol);
	nBeta = nMolB / double(nMol);
	Q = 1.0/(nAlpha*nBeta)*Sqr((mRatio*nAlpha)+nBeta);

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
	std::normal_distribution<double> normal_dist(0.0, 1.0);

	vecSet(momSum, 0, 0, 0);
	for (int i = 0; i < nMol; i++) {
		vecSet(mol[i].vel, normal_dist(rand_gen),
				normal_dist(rand_gen), normal_dist(rand_gen));
		vecScaleAdd(momSum, momSum, mol[i].mass, mol[i].vel);
		// assign zero init. acceleration
		vecSet(mol[i].acc, 0, 0, 0);
	}

	// account for COM shift
	for (int i = 0; i < nMol; i++) {
		vecScaleAdd(mol[i].vel, mol[i].vel, -1.0/mol[i].mass/nMol, momSum);
	}

	// adjust temperature
	rescaleVels();
}

void rescaleVels() {
	double mv2sum = 0;
	for (int i = 0; i < nMol; i++) {
		mv2sum += mol[i].mass * vecLenSq(mol[i].vel);
	}

	double lambda = std::sqrt(3 * (nMol - 1) * temperature / mv2sum);
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

void singleStep(std::string dot_in) {
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
		evalLatticeCorr();
		printSummary(dot_in);
		accumProps(0);
	}

	if (stepCount % stepDump == 0) {
		posDump(dot_in);
	}

	if (stepCount >= stepEquil && (stepCount - stepEquil) % stepRdf == 0) {
		evalRdf_AB(dot_in);
	}

	if (stepCount >= stepEquil && (stepCount - stepEquil) % stepDiff == 0) {
		evalDiffusion(dot_in);
	}

	if (stepCount >= stepEquil && (stepCount - stepEquil) % stepAcf == 0) {
		evalVacf(dot_in);
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
	double rrNebr, rr;

	rrNebr = Sqr(rCut + rNebrShell);
	nebrTabLen = 0;

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
}

void computeForces() {
	vecR dr;
	double fcVal, rr, rrCut;

	rrCut = Sqr(rCut);
	// resetting the acc. values since they are incremented later on
	for (int i = 0; i < nMol; i++) {
		vecSet(mol[i].acc, 0, 0, 0);
	}
	uSum = 0;
	virSum = 0;

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
			if (mol[j1].type == 1 && mol[j2].type == 1) {
				eps = epsAA;
				sig = sigAA;
			} else if (mol[j1].type == 2 && mol[j2].type == 2) {
				eps = epsBB;
				sig = sigBB;
			} else if ((mol[j1].type == 1 && mol[j2].type == 2)
				|| (mol[j1].type == 2 && mol[j2].type == 1)) {
				eps = epsAA;
				sig = sigBB;
			}
			fcVal = 48.0 * eps * std::pow(sig, 12) / std::pow(rr, 7)
				- 24.0 * eps * std::pow(sig, 6) / std::pow(rr, 4);
			vecScaleAdd(mol[j1].acc, mol[j1].acc, fcVal/mol[j1].mass, dr);
			vecScaleAdd(mol[j2].acc, mol[j2].acc, -fcVal/mol[j2].mass, dr);
			uSum += 4.0 * eps * std::pow(sig, 12) / std::pow(rr, 6)
				- 4.0 * eps * std::pow(sig, 6) / std::pow(rr, 3);
			virSum += fcVal * rr;
		}
	}
}

void evalProps() {
	vecSet(momSum, 0, 0, 0);
	double v2, v2sum = 0, v2max = 0;

	for (int i = 0; i < nMol; i++) {
		vecScaleAdd(momSum, momSum, mol[i].mass, mol[i].vel);
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

void printSummary(std::string dot_in) {
	std::string dot_out = dot_in.erase(dot_in.length()-2).append("out");
	std::ofstream outputFile;
	outputFile.open(dot_out, std::ofstream::app);
	outputFile << stepCount << '\t' << timeNow << '\t'
		<< std::sqrt(vecLenSq(momSum))/nMol << '\t'
		<< kinEnergy.sum << '\t' << totEnergy.sum << '\t'
		<< pressure.sum << '\t' << latticeCorr << '\n';
	outputFile.close();
}

void posDump(std::string dot_in) {
	std::string dot_dump = dot_in.erase(dot_in.length()-2).append("dump");
	std::ofstream dumpFile;
	dumpFile.open(dot_dump, std::ofstream::app);
	dumpFile << "ITEM: TIMESTEP\n" << timeNow << '\n'
		<< "ITEM: NUMBER OF ATOMS\n" << nMol << '\n'
		<< "ITEM: BOX BOUNDS pp pp pp\n"
		<< -0.5*region.x << ' ' << 0.5*region.x << '\n'
		<< -0.5*region.y << ' ' << 0.5*region.y << '\n'
		<< -0.5*region.z << ' ' << 0.5*region.z << '\n'
		<< "ITEM: ATOMS id type x y z\n";
	for (int i = 0; i < nMol; i++) {
		dumpFile << i+1 << ' ' << mol[i].type << ' '
			<< mol[i].r.x << ' ' << mol[i].r.y << ' ' << mol[i].r.z << '\n';
	}
	dumpFile.close();
}

void evalRdf_AB(std::string dot_in) {
	vecR dr;
	double deltaR, normFacAA, normFacBB, normFacAB, rr;

	if (countRdf == 0) {
		for (int n = 0; n < sizeHistRdf; n++) {
			histRdfAA[n] = 0;
			histRdfAB[n] = 0;
			histRdfBB[n] = 0;
		}
	}
	deltaR = rangeRdf / sizeHistRdf;

	for (int j1 = 0; j1 < nMol - 1; j1++) {
		for (int j2 = j1 + 1; j2 < nMol; j2++) {
			vecSub(dr, mol[j1].r, mol[j2].r);
			vecWrapAll(dr, region);
			rr = vecLenSq(dr);
			if (rr < Sqr(rangeRdf)) {
				int n = std::sqrt(rr) / deltaR;
				if (mol[j1].type == 1 && mol[j2].type == 1) {
					histRdfAA[n]++;
				} else if (mol[j1].type == 2 && mol[j2].type == 2) {
					histRdfBB[n]++;
				} else if ((mol[j1].type == 1 && mol[j2].type == 2)
					|| (mol[j1].type == 2 && mol[j2].type == 1)) {
					histRdfAB[n]++;
				}
			}
		}
	}

	countRdf++;
	if (countRdf == limitRdf) {
		normFacAA = vecProd(region)
			/ (2.0 * 3.141592654 * Cub(deltaR) * Sqr(nMolA) * countRdf);
		normFacBB = vecProd(region)
			/ (2.0 * 3.141592654 * Cub(deltaR) * Sqr(nMolB) * countRdf);
		normFacAB = vecProd(region)
			/ (4.0 * 3.141592654 * Cub(deltaR) * (nMolA*nMolB) * countRdf);
		for (int n = 0; n < sizeHistRdf; n++) {
			histRdfAA[n] *= normFacAA / Sqr(n + 0.5);
			histRdfBB[n] *= normFacBB / Sqr(n + 0.5);
			histRdfAB[n] *= normFacAB / Sqr(n + 0.5);
		}
		printRdf_AB(dot_in);
		countRdf = 0;
	}
}

void printRdf_AB(std::string dot_in) {
	std::string dot_rdf = dot_in.erase(dot_in.length()-2).append("rdf");
	std::ofstream rdfFile;
	rdfFile.open(dot_rdf, std::ofstream::app);

	rdfFile << "RDF AA BB AB\n";
	for (int n = 0; n < sizeHistRdf; n++) {
		double rb = (n + 0.5) * rangeRdf / sizeHistRdf;
		rdfFile << rb << '\t'
			<< histRdfAA[n] << '\t'
			<< histRdfBB[n] << '\t'
			<< histRdfAB[n] << '\n';
	}

	rdfFile.close();
}

void evalLatticeCorr() {
	vecR kVec;
	double si = 0, sr = 0, t;

	kVec.x = 2.0 * 3.141592654 * initUcell.x / region.x;
	kVec.y = - kVec.x;
	kVec.z = kVec.x;

	for (int n = 0; n < nMol; n++) {
		t = vecDot(kVec, mol[n].r);
		sr += std::cos(t);
		si += std::sin(t);
	}

	latticeCorr = std::sqrt(Sqr(sr) + Sqr(si)) / nMol;
}

void initDiffusion() {
	for (int nb = 0; nb < nBuffDiff; nb++) {
		bufferAA[nb].count = -nb * nValDiff / nBuffDiff;
		bufferBB[nb].count = -nb * nValDiff / nBuffDiff;
		bufferAB[nb].count = -nb * nValDiff / nBuffDiff;
	}
	zeroDiffusion();
}

void zeroDiffusion() {
	countDiffAvg = 0;
	for (int j = 0; j < nValDiff; j++) {
		rrDiffAvgAA[j] = 0;
		rrDiffAvgBB[j] = 0;
		rrDiffAvgAB[j] = 0;
	}
}

void evalDiffusion(std::string dot_in) {
	vecR dr, rSum;
	for (int nb = 0; nb < nBuffDiff; nb++) {
		if (bufferAA[nb].count == 0) {
			for (int n = 0; n < nMol; n++) {
				if (mol[n].type == 1) {
					bufferAA[nb].orgR[n] = mol[n].r;
					bufferAA[nb].rTrue[n] = mol[n].r;
				} else if (mol[n].type == 2) {
					bufferBB[nb].orgR[n] = mol[n].r;
					bufferBB[nb].rTrue[n] = mol[n].r;
				}
			}
		}
		if (bufferAA[nb].count >= 0) {
			vecSet(rSum, 0, 0, 0);
			int ni = bufferAA[nb].count;
			bufferAA[nb].rrDiff[ni] = 0;
			bufferBB[nb].rrDiff[ni] = 0;
			bufferAB[nb].rrDiff[ni] = 0;
			for (int n = 0; n < nMol; n++) {
				if (mol[n].type == 1) {
					vecSub(dr, bufferAA[nb].rTrue[n], mol[n].r);
					vecDiv(dr, dr, region);
					vecRound(dr);
					vecMul(dr, dr, region);
					vecAdd(bufferAA[nb].rTrue[n], mol[n].r, dr);
					vecSub(dr, bufferAA[nb].rTrue[n], bufferAA[nb].orgR[n]);
					bufferAA[nb].rrDiff[ni] += vecLenSq(dr);
					vecAdd(rSum, rSum, dr);
				} else if (mol[n].type == 2) {
					vecSub(dr, bufferBB[nb].rTrue[n], mol[n].r);
					vecDiv(dr, dr, region);
					vecRound(dr);
					vecMul(dr, dr, region);
					vecAdd(bufferBB[nb].rTrue[n], mol[n].r, dr);
					vecSub(dr, bufferBB[nb].rTrue[n], bufferBB[nb].orgR[n]);
					bufferBB[nb].rrDiff[ni] += vecLenSq(dr);
				}
			}
			bufferAB[nb].rrDiff[ni] = vecLenSq(rSum);
		}
		bufferAA[nb].count++;
	}

	accumDiffusion(dot_in);
}

void accumDiffusion(std::string dot_in) {
	double facAA, facBB, facAB;
	for (int nb = 0; nb < nBuffDiff; nb++) {
		if (bufferAA[nb].count == nValDiff) {
			for (int j = 0; j < nValDiff; j++) {
				rrDiffAvgAA[j] += bufferAA[nb].rrDiff[j];
				rrDiffAvgBB[j] += bufferBB[nb].rrDiff[j];
				rrDiffAvgAB[j] += bufferAB[nb].rrDiff[j];
			}
			bufferAA[nb].count = 0;
			countDiffAvg++;
			if (countDiffAvg == limitDiffAvg) {
				printMsd(dot_in);
				facAA = 1.0 / (nDim * 2 * nMolA * stepDiff * deltaT * limitDiffAvg);
				facBB = 1.0 / (nDim * 2 * nMolB * stepDiff * deltaT * limitDiffAvg);
				facAB = Q / (nDim * 2 * nMol * stepDiff * deltaT * limitDiffAvg);
				for (int k = 1; k < nValDiff; k++) {
					rrDiffAvgAA[k] *= facAA / k;
					rrDiffAvgBB[k] *= facBB / k;
					rrDiffAvgAB[k] *= facAB / k;
				}
				printDiffusion(dot_in);
				zeroDiffusion();
			}
		}
	}
}

void printMsd(std::string dot_in) {
	std::string dot_msd = dot_in.erase(dot_in.length()-2).append("msd");
	std::ofstream msdFile;
	msdFile.open(dot_msd, std::ofstream::app);

	double tVal;
	msdFile << "MSD AA BB AB\n";
	for (int j = 0; j < nValDiff; j++) {
		tVal = j * stepDiff * deltaT;
		msdFile << tVal << '\t'
			<< rrDiffAvgAA[j] / limitDiffAvg / nMolA << '\t'
			<< rrDiffAvgBB[j] / limitDiffAvg / nMolB << '\t'
			<< rrDiffAvgAB[j] / limitDiffAvg / nMol << '\n';
	}

	msdFile.close();
}

void printDiffusion(std::string dot_in) {
	std::string dot_dfs = dot_in.erase(dot_in.length()-2).append("dfs");
	std::ofstream dfsFile;
	dfsFile.open(dot_dfs, std::ofstream::app);

	double tVal;
	dfsFile << "Diffusion AA BB AB\n";
	for (int j = 0; j < nValDiff; j++) {
		tVal = j * stepDiff * deltaT;
		dfsFile << tVal << '\t'
			<< rrDiffAvgAA[j] << '\t'
			<< rrDiffAvgBB[j] << '\t'
			<< rrDiffAvgAB[j] << '\n';
	}

	dfsFile.close();
}

void initVacf() {
	for (int nb = 0; nb < nBuffAcf; nb++) {
		vacBuff[nb].count = -nb * nValAcf / nBuffAcf;
	}
	zeroVacf();
}

void zeroVacf() {
	countAcfAvg = 0;
	for (int j = 0; j < nValAcf; j++) {
		avgAcfVel[j] = 0;
	}
}

void evalVacf(std::string dot_in) {
	for (int nb = 0; nb < nBuffAcf; nb++) {
		if (vacBuff[nb].count == 0) {
			for (int n = 0; n < nMol; n++) {
				vacBuff[nb].orgVel[n] = mol[n].vel;
			}
		}
		if (vacBuff[nb].count >= 0) {
			int ni = vacBuff[nb].count;
			vacBuff[nb].acfVel[ni] = 0;
			for (int n = 0; n < nMol; n++) {
				vacBuff[nb].acfVel[ni] += vecDot(vacBuff[nb].orgVel[n], mol[n].vel);
			}
		}
		vacBuff[nb].count++;
	}

	accumVacf(dot_in);
}

void accumVacf(std::string dot_in) {
	double fac;
	for (int nb = 0; nb < nBuffAcf; nb++) {
		if (vacBuff[nb].count == nValAcf) {
			for (int j = 0; j < nValAcf; j++) {
				avgAcfVel[j] += vacBuff[nb].acfVel[j];
			}
			vacBuff[nb].count = 0;
			countAcfAvg++;
			if (countAcfAvg == limitAcfAvg) {
				fac = stepAcf * deltaT / (nDim * nMol * limitAcfAvg);
				intAcfVel = fac * integrate(avgAcfVel, nValAcf);
				for (int k = 1; k < nValAcf; k++) {
					avgAcfVel[k] /= avgAcfVel[0];
				}
				avgAcfVel[0] = 1;
				printVacf(dot_in);
				zeroVacf();
			}
		}
	}
}

double integrate(double *f, int nf) {
	double s = 0.5 * (f[0] + f[nf-1]);
	for (int i = 1; i < nf-1; i++) {
		s += f[i];
	}
	return s;
}

void printVacf(std::string dot_in) {
	std::string dot_acf = dot_in.erase(dot_in.length()-2).append("acf");
	std::ofstream acfFile;
	acfFile.open(dot_acf, std::ofstream::app);

	double tVal;
	acfFile << "VACF\n";
	for (int j = 0; j < nValAcf; j++) {
		tVal = j * stepAcf * deltaT;
		acfFile << tVal << '\t' << avgAcfVel[j] << '\n';
	}
	acfFile << "VACF integral: " << intAcfVel << '\n';

	acfFile.close();
}
