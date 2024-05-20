/* 
 * Backward substitution
 * MPI version, row oriented 
 * point to point communication vectorized
 * augmented multiplication
 *
 * H. Moritsch 2022-05-05
 */

/* 
 mpicc -DMAXN=2048 -DMAXK=100 -DNP=32 -DAUGN=100 -DAUGM=1 -o backsub_mpi_row -O3 backsub_mpi_row.c multadd.c
 
 sbatch submit.sh:
 #!/bin/bash
 #SBATCH -n 32
 mpirun -np 32 ./backsub_mpi_row
 
 results on alma:
 p00 final: MAXN=2048 MAXK=100 NP=32 AUGN=100 AUGM=1 totalsum:  417767412.125 t_total_p0:  2830.132 msec

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "multadd.h"
#define min(a, b) (a<b? a : b)
#define max(a, b) (a>b? a : b)

/* MAXN must be a multiple of NP! */

/* #define LOG 1 */

double  U[MAXN/NP][MAXN];   // (block, *) distributed
double  b[MAXN/NP];         // (block) distributed
double  x[MAXN];            // replicated

int     rank, size;
/*-----------------------------------------------------------*/
int main(argc, argv)
/*-----------------------------------------------------------*/
int argc;
char **argv;

{
int n = MAXN;
int lb, ub, len, last;
char processor[MPI_MAX_PROCESSOR_NAME];

int i, j, k;
int p, m;
int count;

double sum, totalsum;
double t_total, t_total_p0;

MPI_Status status;

/* MPI initialization */

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

if (size != NP) {
    printf("*** np must be %d vs. %d ***\n", NP, size);
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
    }

MPI_Get_processor_name(processor, &len);
printf("p%03d/%d on %s\n", rank, size, processor);

/* local segment */

len = n/size;
lb = rank*len;
ub = (rank+1)*len-1;

printf("p%02d: [%d .. %d] %d\n", rank, lb, ub, len);

last = size-1;

/* initialization */

for(i=lb; i<=ub; i++) {
    for(j=0; j<i; j++)
        U[i-lb][j] = 0.0;
    for(j=i; j<n; j++)
        U[i-lb][j]= ((double)(i+1)/(double)(j+1));
    }

MPI_Barrier(MPI_COMM_WORLD);
t_total = MPI_Wtime();

if (rank == 0) {
    t_total_p0 = t_total;
    totalsum = 0.0;
    }

for (k=0; k<MAXK; k++) {   // do it several times to get longer execution times

    MPI_Barrier(MPI_COMM_WORLD);

    /* generate right hand side */

    for(i=lb; i<=ub; i++) 
        // b[i-lb] = ( (double)i + (double)k/(double)MAXK ) * (double)n; // different values in each k-iteration
        b[i-lb] = ( (double)i + 1.0/(double)MAXK ) * (double)n; // same values in each k-iteration

    /* backward substitution */

    /* initialize x */

    for(i=lb; i<=ub; i++) 
        x[i] = b[i-lb];

    count = 0;

    /* first solution */

    if (rank == last) {
        x[n-1] /= U[n-1-lb][n-1];
        #ifdef LOG
            printf("p%02d: x[%d] = x[%d] / U[%d, %d] %d!\n", rank, n-1, n-1, n-1, n-1, ++count);
        #endif
        }

    /*
    j-loop with a check for ((j+1) mod len == 0)
    -----------------------------------------------------------------------------------
    for (j = n-1; j > lb; j--) { 

        if ((j+1) % len == 0) {     // if we need a(n additional) nonlocal segment of x
            q = (j+1)/len - 1;      // from higher-rank process q:
            if (q > rank) {     
                // receive x[q*len:(q+1)*len-1] from process q 
                }
            }
        ...
        }
    -----------------------------------------------------------------------------------
    is transformed into outer m-loop and j-loop with adjusted bounds:
    */
        
    for (m = size-1; m >= rank; m--) { 

        if (rank < min(m, last)) {   
            // we need a(n additional) nonlocal segment of x from higher-rank process m

            #ifdef LOG
                printf("p%02d: receive x[%d:%d] from process %d\n", rank, m*len, (m+1)*len-1, m); 
            #endif
            MPI_Recv(x+m*len, len, MPI_DOUBLE, m, 0, MPI_COMM_WORLD, &status);
            }

/* 
 *      if (rank == last) 
 *          seq_start = MPI_Wtime(); // measurement of pure sequential part
 */

        for (j = (m+1)*len-1; j > max(m*len-1, lb); j--) { 
            #ifdef LOG
                printf("p%02d: j=%d\n", rank, j);
            #endif
    
            for (i = lb; i < min(j, ub+1); i++) {    
    
                /* update x[i] */

                if (AUGN == 0)
                    x[i] -= U[i-lb][j] * x[j];
                else 
                    x[i] -= multadd(U[i-lb][j], x[j], AUGN, AUGM);
                #ifdef LOG
                    printf("p%02d:\ti=%d x[%d] = x[%d] - U[%d, %d] * x[%d] %d.\n", rank, i, i, i, i, j, j, ++count);
                #endif
                }
    
            if ((j-1)/len == rank) {  // j-1 div len
                x[j-1] /= U[j-1-lb][j-1];
                #ifdef LOG
                    printf("p%02d:\tx[%d] = x[%d] / U[%d, %d] done %d.\n", rank, j-1, j-1, j-1, j-1, ++count);
                #endif
                }

            } // j-loop
        
/* 
 *      if (rank == last) 
 *          seq_stop = MPI_Wtime(); 
 */

        } // m-loop

    /* send x[lb:ub] to all lower-rank processes */

    for (p=rank; p>0; p--) {
        #ifdef LOG
            printf("p%02d: send x[%d:%d] to process %d\n", rank, lb, ub, p-1);
        #endif
        MPI_Send(x+lb, len, MPI_DOUBLE, p-1, 0, MPI_COMM_WORLD);
        }

    /* calculate sum(x) */

    if (rank == 0) {
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum += x[i];
            // printf("p%02d: sum:      %12.3f\n", rank, sum); // comment this out for time measurement
            }
        totalsum += sum;
        }

    } // k-loop

t_total = MPI_Wtime() - t_total;

// printf("p%02d: t_total %12.3f msec\n", rank, t_total * 1000.0);

MPI_Barrier(MPI_COMM_WORLD);
if (rank == 0) {
    t_total_p0 = MPI_Wtime() - t_total_p0;
    printf("p%02d final: MAXN=%d MAXK=%d NP=%d AUGN=%d AUGM=%d totalsum:%15.3f t_total_p0:%10.3f msec\n", 
        rank, MAXN, MAXK, NP, AUGN, AUGM, totalsum, t_total_p0 * 1000.0);
    }

MPI_Finalize();

return 0;
}

