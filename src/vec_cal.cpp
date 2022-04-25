#include <cmath>
#include <algorithm>
#include "types.hpp"

real Sqr(real x) {
	return x * x;
}

real Cub(real x) {
	return x * x * x;
}

void vecSet(vecR &vec, real cx, real cy, real cz) {
	vec.x = cx;
	vec.y = cy;
	vec.z = cz;
}

void vecRound(vecR &vec) {
	vec.x = int(std::round(vec.x));
	vec.y = int(std::round(vec.y));
	vec.z = int(std::round(vec.z));
}

void vecFloor(vecR &vec) {
	vec.x = int(vec.x);
	vec.y = int(vec.y);
	vec.z = int(vec.z);
}

void vecAdd(vecR &sum, vecR u, vecR v) {
	sum.x = u.x + v.x;
	sum.y = u.y + v.y;
	sum.z = u.z + v.z;
}

void vecSub(vecR &diff, vecR u, vecR v) {
	diff.x = u.x - v.x;
	diff.y = u.y - v.y;
	diff.z = u.z - v.z;
}

void vecMul(vecR &mul, vecR u, vecR v) {
	mul.x = u.x * v.x;
	mul.y = u.y * v.y;
	mul.z = u.z * v.z;
}

void vecDiv(vecR &div, vecR u, vecR v) {
	div.x = u.x / v.x;
	div.y = u.y / v.y;
	div.z = u.z / v.z;
}

real vecProd(vecR u) {
	return u.x * u.y * u.z;
}

real vecDot(vecR u, vecR v) {
	return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

real vecLenSq(vecR u) {
	return vecDot(u, u);
}

void vecScale(vecR &vec, real s) {
	vec.x *= s;
	vec.y *= s;
	vec.z *= s;
}

void vecScaleCopy(vecR &vec, real s, vecR u) {
	vec.x = s * u.x;
	vec.y = s * u.y;
	vec.z = s * u.z;
}

void vecScaleAdd(vecR &vec, vecR u, real s, vecR v) {
	vec.x = u.x + s * v.x;
	vec.y = u.y + s * v.y;
	vec.z = u.z + s * v.z;
}

int vecLinear(vecR u, vecR v) {
	return (u.z * v.y * v.x) + (u.y * v.x) + u.x + 0.5;
}

void vecWrapAll(vecR &vec, vecR region) {
	if (vec.x >= 0.5 * region.x) {
		vec.x -= region.x;
	} else if (vec.x < -0.5 * region.x) {
		vec.x += region.x;
	}
	if (vec.y >= 0.5 * region.y) {
		vec.y -= region.y;
	} else if (vec.y < -0.5 * region.y) {
		vec.y += region.y;
	}
	if (vec.z >= 0.5 * region.z) {
		vec.z -= region.z;
	} else if (vec.z < -0.5 * region.z) {
		vec.z += region.z;
	}
}

void cellWrapAll(vecR &vec, vecR &shift, vecR cells, vecR region) {
	if (vec.x >= cells.x) {
		vec.x = 0;
		shift.x = region.x;
	} else if (vec.x < 0) {
		vec.x = cells.x - 1;
		shift.x = - region.x;
	}
	if (vec.y >= cells.y) {
		vec.y = 0;
		shift.y = region.y;
	} else if (vec.y < 0) {
		vec.y = cells.y - 1;
		shift.y = - region.y;
	}
	if (vec.z >= cells.z) {
		vec.z = 0;
		shift.z = region.z;
	} else if (vec.z < 0) {
		vec.z = cells.z - 1;
		shift.z = - region.z;
	}
}

void propZero(Prop &prop) {
	prop.sum = 0;
	prop.sum2 = 0;
}

void propAccum(Prop &prop) {
	prop.sum += prop.val;
	prop.sum2 += Sqr(prop.val);
}

void propAvg(Prop &prop, int n) {
	prop.sum /= n;
	prop.sum2 = std::sqrt(std::max(prop.sum2/n - Sqr(prop.sum), 0.0));
}
