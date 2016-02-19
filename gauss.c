/* 
 * Original author:  UNKNOWN
 *
 * Modified:         Kai Shen (January 2010)
 */
#ifndef _REENTRANT
#define _REENTRANT		/* basic 3-lines for threads */
#endif

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>


/* #define DEBUG */

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}

/* Solve the equation:
 *   matrix * X = R
 */

double **matrix, *X, *R;

/* Pre-set solution. */

double *X__;

/* Number of threads */

int task_num = 1;

/* size of the matrix */

int nsize = 0;

/* pivot for all threads - set by thread 1 */

int t = 0;

double pivotval_global;

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
void
barrier (int expect)
{
    static int arrived = 0;
     
    pthread_mutex_lock (&mut);  //lock
    arrived++;
    
    if (arrived < expect)
        pthread_cond_wait (&cond, &mut);
    else {
        arrived = 0;		// reset the barrier before broadcast is important
        pthread_cond_broadcast (&cond);
    }
    
    pthread_mutex_unlock (&mut);	//unlock
}

/* Initialize the matirx. */

int initMatrix(const char *fname)
{
    FILE *file;
    int l1, l2, l3;
    double d;
    int nsize;
    int i, j;
    double *tmp;
    char buffer[1024];

    if ((file = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "The matrix file open error\n");
        exit(-1);
    }
    
    /* Parse the first line to get the matrix size. */
    if(fgets(buffer, 1024, file))
	// maybe add check 
    sscanf(buffer, "%d %d %d", &l1, &l2, &l3);
    nsize = l1;
#ifdef DEBUG
    fprintf(stdout, "matrix size is %d\n", nsize);
#endif

    /* Initialize the space and set all elements to zero. */
    matrix = (double**)malloc(nsize*sizeof(double*));
    assert(matrix != NULL);
    tmp = (double*)malloc(nsize*nsize*sizeof(double));
    assert(tmp != NULL);    
    for (i = 0; i < nsize; i++) {
        matrix[i] = tmp;
        tmp = tmp + nsize;
    }
    for (i = 0; i < nsize; i++) {
        for (j = 0; j < nsize; j++) {
            matrix[i][j] = 0.0;
        }
    }

    /* Parse the rest of the input file to fill the matrix. */
    for (;;) {
	if(fgets(buffer, 1024, file))
		// maybe add check
	sscanf(buffer, "%d %d %lf", &l1, &l2, &d);
	if (l1 == 0) break;

	matrix[l1-1][l2-1] = d;
#ifdef DEBUG
	fprintf(stdout, "row %d column %d of matrix is %e\n", l1-1, l2-1, matrix[l1-1][l2-1]);
#endif
    }

    fclose(file);
    return nsize;
}

/* Initialize the right-hand-side following the pre-set solution. */

void initRHS(int nsize)
{
    int i, j;

    X__ = (double*)malloc(nsize * sizeof(double));
    assert(X__ != NULL);
    for (i = 0; i < nsize; i++) {
	X__[i] = i+1;
    }

    R = (double*)malloc(nsize * sizeof(double));
    assert(R != NULL);
    for (i = 0; i < nsize; i++) {
	R[i] = 0.0;
	for (j = 0; j < nsize; j++) {
	    R[i] += matrix[i][j] * X__[j];
	}
    }
}

/* Initialize the results. */

void initResult(int nsize)
{
    int i;

    X = (double*)malloc(nsize * sizeof(double));
    assert(X != NULL);
    for (i = 0; i < nsize; i++) {
	X[i] = 0.0;
    }
}

/* Get the pivot - the element on column with largest absolute value. */

void getPivot(int nsize, int currow)
{
    int i, pivotrow;

    pivotrow = currow;
    for (i = currow+1; i < nsize; i++) {
	if (fabs(matrix[i][currow]) > fabs(matrix[pivotrow][currow])) {
	    pivotrow = i;
	}
    }
    if (fabs(matrix[pivotrow][currow]) == 0.0) {
        fprintf(stderr, "The matrix is singular\n");
        exit(-1);
    }
    
    if (pivotrow != currow) {
#ifdef DEBUG
	//printf("pivot row at step %d is %d\n", currow, pivotrow);
#endif
        for (i = currow; i < nsize; i++) {
            SWAP(matrix[pivotrow][i],matrix[currow][i]);
        }
        SWAP(R[pivotrow],R[currow]);

    }
}

void
errexit (const char *err_str)
{
    fprintf (stderr, "%s", err_str);
    exit (1);
}


/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */
void computeGauss(int task_id)
{
    int i, j, k;
    double pivotval;
	
    for (i = 0; i < nsize; i++) { 
	if(task_id == 0) {
		getPivot(nsize,i);
		pivotval_global = matrix[i][i];
	
	    /* Scale the main row. */	
	    if (pivotval_global != 1.0) {
	        matrix[i][i] = 1.0;
	        for (j = i + 1; j < nsize; j++) {
		    matrix[i][j] /= pivotval_global;
	        }  
	       R[i] /= pivotval_global; 
	    }	
	} 
	barrier(task_num);
	/* Factorize the rest of the matrix. */
	for (j = i +1 + task_id; j < nsize; j = j + task_num ) {
	pivotval = matrix[j][i];
        matrix[j][i] = 0.0;  	
	    for (k = i + 1; k < nsize; k++) {
                matrix[j][k] -= pivotval * matrix[i][k];
            }          
            R[j] -= pivotval * R[i];
	}
        barrier(task_num);
   }
}

/* Solve the equation. */

void solveGauss(int nsize)
{
    int i, j;

    X[nsize-1] = R[nsize-1];
    for (i = nsize - 2; i >= 0; i --) {
        X[i] = R[i];
        for (j = nsize - 1; j > i; j--) {
            X[i] -= matrix[i][j] * X[j];
        }
    }

#ifdef DEBUG
    fprintf(stdout, "X = [");
    for (i = 0; i < nsize; i++) {
        fprintf(stdout, "%.6f ", X[i]);
    }
    fprintf(stdout, "];\n");
#endif
}

extern char *optarg;


void *
work_thread (void *lp)
{
    int task_id = *((int *) lp);
    struct timeval start, finish;
	
    int begin, end;    
    begin = (nsize * task_id) / task_num + 1;
    end = (nsize * (task_id + 1)) / task_num;
  
    barrier (task_num);
    
    // start timer if first thread
    if(task_id==0)
        gettimeofday (&start, NULL);
    
    // function for each thread
    computeGauss(task_id);

    //computeGauss(task_id);
    barrier(task_num);
    
    // once all threads have completed, check the time 
    gettimeofday (&finish, NULL);
    
    // print out the time for only one thread
    if(task_id==0)
	fprintf(stdout, "Time:  %f seconds\n", (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec)*0.000001);
}


int main(int argc, char *argv[])
{
    int i;
    struct timeval start, finish;
    double error;
    
    pthread_attr_t attr;
    pthread_t *tid;
    int *id;
    
    if (argc < 2) {
	fprintf(stderr, "usage: %s <matrixfile>\n", argv[0]);
	exit(-1);
    }
    
    // for getting the threads
    if(argc == 3) {
        task_num = strtol(argv[2], NULL, 10);
    }

    nsize = initMatrix(argv[1]);
    initRHS(nsize);
    initResult(nsize);
   
    // create threads
    id = (int *) malloc (sizeof (int) * task_num);
    tid = (pthread_t *) malloc (sizeof (pthread_t) * task_num);
    if (!id || !tid)
        errexit ("out of shared memory");
    pthread_attr_init (&attr);
    pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);
    for (i = 1; i < task_num; i++) {
        id[i] = i;
        pthread_create (&tid[i], &attr, work_thread, &id[i]);
    }

    id[0]=0;
    work_thread(&id[0]);
    // wait for all threads to finish
    for (i = 1; i < task_num; i++)
        pthread_join (tid[i], NULL);
	    
    solveGauss(nsize);
    
    error = 0.0;
    for (i = 0; i < nsize; i++) {
	double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i])/X__[i]);
	if (error < error__) {
		error = error__;
	}
    }
    fprintf(stdout, "Error: %e\n", error);
}

 
