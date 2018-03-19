#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,num_threads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

// #define DEBUG_LA
#ifdef DEBUG_LA
    // For visualization
    char *outname = (char*) malloc(32 * sizeof(char));
#endif

    // Size of one side of a 2D bin square
    int bin_count = get_bin_count();
    // bin_t * binArr = make_bin(n);
    make_bin(n);
    clear_bins(0, bin_count);

    // #pragma omp parallel private(dmin) 
    // {
    // num_threads = omp_get_num_threads();
    // int thread_id = omp_get_thread_num();
    // Index start of bins for this thread
    // int thread_start = (bin_count+1)/num_threads*thread_id;
    // Index end of bins for this thread
    // Max function added to account for odd # of bins
    // int thread_end = min((bin_count+1)/num_threads*(thread_id+1), bin_count);
    // printf("This thread num = %i, start = %i, end = %i\r\n", thread_id, thread_start, thread_end);
    for( int step = 0; step < 1000; step++ )
    {
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute all forces
        //
        // Clear bins out to redo in case of move
        // clear_bins(thread_start, thread_end);
        #pragma omp for
        for( int i = 0; i < bin_count; i++ ) 
            clear_bin(i);

        // #pragma omp barrier


        // Make bins and set particles into bins
        // #pragma omp master
        #pragma omp for
        for(int pCount = 0; pCount < n; pCount++ )
        {
            set_bin(particles[pCount]);
        }

        // #pragma omp barrier
        // printf("Computing forces... step: %i\r\n", step);
        //
        //  compute forces
        //
        // #pragma omp for reduction (+:navg) reduction(+:davg)
        #pragma omp for
        for (int i=0; i < bin_count; i++)
        {
            apply_force_bin(i,&dmin,&davg,&navg);
        }

        /*
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }
        */
        
		
        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	      if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }

#ifdef DEBUG_LA
          // Output for visualization
          #pragma omp master
          {
          sprintf( outname, "out/fout-%05d.txt", step );
          FILE *fout = fopen( outname, "w" );
          save( fout, n, particles);
          fclose( fout );
          }
#endif
    // }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,num_threads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,num_threads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
