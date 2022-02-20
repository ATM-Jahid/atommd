#include "types.hpp"

real Sqr(real x) {
	return x * x;
}

real Cub(real x) {
	return x * x * x;
}

void vecSet(vecR &vec, real cx, real cy) {
	vec.x = cx;
	vec.y = cy;
}

void vecSetAll(vecR &vec, real c) {
	vecSet(vec, c, c);
}

void vecSetZero(vecR &vec) {
	vecSetAll(vec, 0);
}

void vecAdd(vecR &sum, vecR u, vecR v) {
	sum.x = u.x + v.x;
	sum.y = u.y + v.y;
}

void vecSub(vecR &diff, vecR u, vecR v) {
	diff.x = u.x - v.x;
	diff.y = u.y - v.y;
}

void vecMul(vecR &mul, vecR u, vecR v) {
	mul.x = u.x * v.x;
	mul.y = u.y * v.y;
}

void vecDiv(vecR &div, vecR u, vecR v) {
	div.x = u.x / v.x;
	div.y = u.y / v.y;
}

real vecProd(vecR u) {
	return u.x * u.y;
}

real vecDot(vecR u, vecR v) {
	return (u.x * v.x) + (u.y * v.y);
}

real vecLenSq(vecR u) {
	return vecDot(u, u);
}

void vecScale(vecR &vec, real s) {
	vec.x *= s;
	vec.y *= s;
}

void vecScaleCopy(vecR &vec, real s, vecR u) {
	vec.x = s * u.x;
	vec.y = s * u.y;
}

void vecScaleAdd(vecR &vec, vecR u, real s, vecR v) {
	vec.x = u.x + s * v.x;
	vec.y = u.y + s * v.y;
}

void vecSelfScaleAdd(vecR &u, real s, vecR v) {
	vecScaleAdd(u, u, s, v);
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
}
