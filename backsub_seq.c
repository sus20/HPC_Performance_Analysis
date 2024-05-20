/* 
 * Backward substitution
 * sequential version, row oriented
 * MPI time measurement
 *
 * H.Moritsch 2021-05-05
 */

#include <stdio.h>
#include <stdlib.h>
// #include <mpi.h>

#define MAXN 32 
#define MAXK 10

/*-----------------------------------------------------------*/
int main(argc, argv)
/*-----------------------------------------------------------*/
int argc;
char **argv;

{
double U[MAXN][MAXN]; 
double b[MAXN];
double x[MAXN];

int n = MAXN;
int i, j, k;

double sum, totalsum;

double time_init, time_backsub;
double t_start, t_stop;

// initialization

// t_start = MPI_Wtime();

for(i=0;i<n;i++) {
    for(j=0;j<i;j++)
        U[i][j] = 0.0;
    for(j=i;j<n;j++)
        U[i][j]= ((double)(i+1)/(double)(j+1));
    }

// t_stop = MPI_Wtime();
// time_init = t_stop - t_start;

// backward substituion

totalsum = 0.0;

// t_start = MPI_Wtime();

for (k=0;k<MAXK;k++) {   // do several times to get longer execution times

    // generate (new) right hand side

    for(i=0;i<n;i++)
        b[i] = ( (double)i + (double)k/(double)MAXK ) * (double)n;

    // perform backward substitution (column oriented)

    for (i = n - 1; i >= 0; i--) {
        x[i] = b[i];

        for (j = i + 1; j < n; j++) 
            // update x[i]
            x[i] -= U[i][j] * x[j];

        x[i] /= U[i][i];
        }

    // calculate sum of elements of x

    sum = 0.0;
    for (i=0; i<n; i++) {
        sum += x[i];
        }
    printf("sum:      %12.3f\n", sum); // comment this out for time measurement

    totalsum += sum;
        
    } // k-loop

// t_stop = MPI_Wtime();
// time_backsub = t_stop - t_start;

printf("totalsum: %12.3f\n", totalsum);
// printf("init:     %12.3f msec\n", time_init*1000.0);
// printf("backward: %12.3f msec\n", time_backsub*1000.0);

return 0;
}

