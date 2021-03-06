#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "bhtree.h"
#define density 0.0005
//
//  benchmarking program
//
int main(int argc, char **argv)
{
    //int navg, nabsavg = 0;
    //double dmin, absmin = 1.0, davg, absavg = 0.0;
    //double rdavg, rdmin;
    //int rnavg;

    //
    //  process command line parameters
    //
    if (find_option(argc, argv, "-h") >= 0)
    {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n");
        printf("-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("test1 current rank %d\n", rank);
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen(sumname, "a") : NULL;

    particle_t *particles = (particle_t *)malloc(n * sizeof(particle_t));

    //define a new type PARTICLE which contains 6 MPI_DOUBLE data
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);

    //
    //  set up the data partitioning across processors
    //  n is the number of particle, n_proc is the number of thread (why plus n_proc here?)
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int *)malloc((n_proc + 1) * sizeof(int));
    for (int i = 0; i < n_proc + 1; i++)
        partition_offsets[i] = min(i * particle_per_proc, n);

    int *partition_sizes = (int *)malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; i++)
        partition_sizes[i] = partition_offsets[i + 1] - partition_offsets[i];

    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t *)malloc(nlocal * sizeof(particle_t));

    //init the array to gather all the root value
    particle_t *allroot = (particle_t *)calloc(n_proc, sizeof(particle_t));

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //

    //every thread init their own part
    //make sure the points fall into low left to upper right
    Quad *initq = new (Quad);
    double localsize;
    //get the total size
    localsize = sqrt(density * particle_per_proc);

    //printf("localsize %f\n", localsize);
    initq->cl = localsize;

    //partition

    if (n_proc == 1)
    {
        initq->llx = 0;
        initq->lly = 0;
    }
    else if (n_proc == 2)
    {
        initq->llx = 0;
        initq->lly = 0 + rank * localsize;
    }
    else if (n_proc == 4)
    {
        initq->llx = 0 + (rank % 2) * localsize;
        initq->lly = 0 + (rank / 2) * localsize;
    }
    else if (n_proc == 8)
    {
        initq->llx = 0 + (rank % 4) * localsize;
        initq->lly = 0 + (rank / 4) * localsize;
    }
    else if (n_proc == 16)
    {
        initq->llx = 0 + (rank % 4) * localsize;
        initq->lly = 0 + (rank / 4) * localsize;
    }

    init_particles_inthread(nlocal, initq, local);

    //printf("rank value %d, x  %f, y  %f\n",rank,initq->llx,initq->lly);

    //printf("test current rank %d\n", rank);
    //
    //if (rank == 0)
    //{
    //    init_particles(n, particles);
    //}

    //MPI_Scatterv(particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD);

    //
    //  simulate a number of time steps
    //
    //  Output for visualization
#ifdef DEBUG_MPI  
    char outname [100];
#endif
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++)
    {
        //navg = 0;
        //dmin = 1.0;
        //davg = 0.0;

        //1 get speecific local elements
        //MPI_Allgatherv(local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD);


        //2 build tree (tree should be buit for every iteration)

        //init position should relate to the rank

        BHTreeNode *bhtree = BHTree(initq);
        for (int i = 0; i < nlocal; i++)
        {
            if (Quadcontains(initq, (local + i)->x, (local + i)->y) == true)
            {
                BHTinsert(bhtree, local + i);
            }
        }


        //printf("rank %d step %d insert ok\n",rank, step);

        //printf("debug rank %d step %d, insert ok\n",rank,step);

        //get root node from other partition
        //printf("threadid %d %f %f\n", rank, bhtree->particle->x, bhtree->particle->y);
        //attention, recieve number of recieve value is 1
        if(bhtree->particle==NULL){
            //insert a psudo node
            particle_t*p=new(particle_t);
            p->x=initq->llx;
            p->y=initq->lly;
            bhtree->particle=p;
        }
        
        //printf("rank %d step %d parx %f pary %f alltogether start\n",rank, step, bhtree->particle->x,bhtree->particle->y);
        MPI_Allgather(bhtree->particle, 1, PARTICLE, allroot, 1, PARTICLE, MPI_COMM_WORLD);

        
        //MPI_Barrier(MPI_COMM_WORLD);

        //printf("rank %d step %d alltogether ok\n",rank, step);
       
        //printf("thread %d step %d all together ok\n",rank,step);
        //3 compute force (for the particle in other mpi thread, use root value to take place )
        //an array is needed to store all the root value from other mpi thread
        //compute the force in local tree

        for (int i = 0; i < nlocal; i++)
        {
            applyForceTree(local + i, bhtree);
            for (int j = 0; j < n_proc; j++)
            {
                if (j != rank)
                {
                    applyForceTwoParticle((local + i), allroot+j);
                }
            }
        }

        //printf("rank %d step %d force caculation ok\n",rank, step);

        // printf("debug rank %d step %d, apply force ok\n",rank, step);

        //printf("rank %d applyforce ok\n", rank);
        //compute the force from the other mpi thread
        //collect the root value of other mpi thread

        //MPI_Allgather(void* send_data,int send_count,MPI_Datatype send_datatype,void* recv_data,int recv_count,MPI_Datatype recv_datatype,MPI_Comm communicator)

        //MPI_Barrier(MPI_COMM_WORLD);

        //if (bhtree->particle != NULL)
        //{

        // }

        /*
        if (rank == 0)
        {
            for (int i = 0; i < n_proc; i++)
            {
                printf("index %d x %f y %f\n", i, (allroot + i)->x, (allroot + i)->y);
            }
        }
*/
        //
        //  collect all global data locally (not good idea to do)
        //

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //

        //
        //  compute all forces
        //
        /*
        for (int i = 0; i < nlocal; i++)
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < n; j++)
                apply_force(local[i], particles[j], &dmin, &davg, &navg);
        }

        if (find_option(argc, argv, "-no") == -1)
        {

            MPI_Reduce(&davg, &rdavg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&navg, &rnavg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&dmin, &rdmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                //
                // Computing statistical data
                //
                if (rnavg)
                {
                    absavg += rdavg / rnavg;
                    nabsavg++;
                }
                if (rdmin < absmin)
                    absmin = rdmin;
            }
        }
        */

        //
        //  move particles, only consider the particles in current thread
        //

        //printf("apply force ok\n");
        //4 move
        for (int i = 0; i < nlocal; i++)
        {
            //for n=2 and n=8 the region is not squre, consider this
            //printf("index %d x %f y %f ax %f ay %f\n", i, local[i].x, local[i].y, local[i].ax, local[i].ay);
            move_mpi(step, local[i], localsize * n_proc);
        }

        //printf("debug rank %d step %d, move mpi ok\n",rank, step);

        //if (rank == 0)
        //{
        //    printf("index %d\n", step);
        //}
        if (find_option(argc, argv, "-no") == -1)
        {
            if (fsave && (step % SAVEFREQ) == 0)
            {
                save(fsave, nlocal, local);
            }
        }

#ifdef DEBUG_MPI
    // For visualization
    sprintf(outname, "out/fout-rank%d-%05d.txt", rank, step );
    FILE *fout = fopen( outname, "w" );
    save( fout, nlocal, local);
    fclose( fout );
#endif
        //printf("rank %d step %d finish\n",rank, step);
    }

    simulation_time = read_timer() - simulation_time;

    if (rank == 0)
    {
        printf("rank id (%d) n = %d, simulation time = %g seconds", rank, n, simulation_time);

        if (find_option(argc, argv, "-no") == -1)
        {
            //if (nabsavg)
            //    absavg /= nabsavg;
            //
            //  -the minimum distance absmin between 2 particles during the run of the simulation
            //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
            //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
            //
            //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
            //
            //printf(", absmin = %lf, absavg = %lf", absmin, absavg);
            //if (absmin < 0.4)
            //    printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
            //if (absavg < 0.8)
            //    printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
        }
        printf("\n");

        //
        // Printing summary data
        //
        if (fsum)
            fprintf(fsum, "%d %d %g\n", n, n_proc, simulation_time);
    }

    //
    //  release resources
    //
    if (fsum)
        fclose(fsum);
    free(partition_offsets);
    free(partition_sizes);
    free(local);
    free(particles);
    if (fsave)
        fclose(fsave);

    MPI_Finalize();

    return 0;
}
