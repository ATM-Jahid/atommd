double Sqr(double);
double Cub(double);

void vecSet(vecR &, double, double, double);
void vecRound(vecR &);
void vecFloor(vecR &);
void vecAdd(vecR &, vecR, vecR);
void vecSub(vecR &, vecR, vecR);
void vecMul(vecR &, vecR, vecR);
void vecDiv(vecR &, vecR, vecR);

double vecProd(vecR);
double vecDot(vecR, vecR);
double vecLenSq(vecR);

void vecScale(vecR &, double);
void vecScaleCopy(vecR &, double, vecR);
void vecScaleAdd(vecR &, vecR, double, vecR);
int vecLinear(vecR, vecR);

void vecWrapAll(vecR &, vecR);
void cellWrapAll(vecR &, vecR &, vecR, vecR);

void propZero(Prop &);
void propAccum(Prop &);
void propAvg(Prop &, int);
