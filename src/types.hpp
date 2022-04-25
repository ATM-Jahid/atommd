typedef double real;

typedef struct {
	real x, y, z;
} vecR;

typedef struct {
	vecR r, vel, acc;
} Mol;

typedef struct {
	real val, sum, sum2;
} Prop;

typedef struct {
	vecR *orgR, *rTrue;
	real *rrDiff;
	int count;
} Tbuff;
