typedef struct {
	double x, y, z;
} vecR;

typedef struct {
	vecR r, vel, acc;
	double mass;
	int type;
} Mol;

typedef struct {
	double val, sum, sum2;
} Prop;

typedef struct {
	vecR *orgR, *rTrue;
	double *rrDiff;
	int count;
} Tbuff;

typedef struct {
	vecR *orgVel;
	double *acfVel;
	int count;
} Vbuff;
