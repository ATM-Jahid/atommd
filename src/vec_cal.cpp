#include <math.h>
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
	prop.sum2 = sqrt(std::max(prop.sum2/n - Sqr(prop.sum), 0.0));
}
