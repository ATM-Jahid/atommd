real Sqr(real);
real Cub(real);

void vecSet(vecR &, real, real);
void vecAdd(vecR &, vecR, vecR);
void vecSub(vecR &, vecR, vecR);
void vecMul(vecR &, vecR, vecR);
void vecDiv(vecR &, vecR, vecR);

real vecProd(vecR);
real vecDot(vecR, vecR);
real vecLenSq(vecR);

void vecScale(vecR &, real);
void vecScaleCopy(vecR &, real, vecR);
void vecScaleAdd(vecR &, vecR, real, vecR);

void vecWrapAll(vecR &, vecR);

void propZero(Prop &);
void propAccum(Prop &);
void propAvg(Prop &, int);
