double** Create2DMatrix(int, int);
int** Create2DMatrix_INT(int, int);
double*** Create3DMatrix(int, int, int);
int*** Create3DMatrix_INT(int, int, int);
double**** Create4DMatrix(int, int, int, int);
void Delete4DMatrix(double****&, int, int, int, int);
void Delete3DMatrix(double***&, int, int, int);
void Delete3DMatrix_INT(int***&, int, int, int);
void Delete2DMatrix(double**&, int, int);

static void swap4(void *);
static void swap8(void *);

