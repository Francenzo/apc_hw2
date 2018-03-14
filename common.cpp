#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "common.h"

using namespace std;

double size;
// Amount of bins in one row
int block_row_count;
// Vector of bins
vector< vector<particle_t*> > binVec;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
    // Divide length of side by cutoff length
    block_row_count = (size/cutoff);
    printf("size = %f, b = %i\r\n", size, block_row_count);
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
// Get total number of bins
//
int get_bin_count() 
{
    return (block_row_count*block_row_count);
}

//
// Create vector of bins
//
void make_bin(int n) 
{
    int b_square = block_row_count*block_row_count;
    binVec = vector< vector<particle_t *> >(b_square);
}


//
// Puts particles into their bins
//
void set_bin(particle_t & particle) 
{
    int frac_x = particle.x/size;
    int frac_y = particle.y/size;
    int binNum = (frac_x * block_row_count ) + (frac_y * block_row_count * block_row_count );
    // printf("binNum = %i, max = %i, fx = %f, fy = %f, cx = %i, cy = %i\r\n", binNum, block_row_count*block_row_count, particle.x/size, particle.y/size, (int)( (particle.x/size) * block_row_count ), (int) ( (particle.y/size) * block_row_count *block_row_count ));
    particle.binNum = binNum;
    binVec.at(binNum).push_back(&particle);
}

//
// Move particle into new bin
//
void move_bin(particle_t & particle) 
{
    vector<particle_t*> bin = binVec.at(particle.binNum);

    // Look for particle in bin and remove it before transfer
    for(int iCount = 0; iCount < bin.size(); iCount++)
    {
        if ((&particle) == bin.at(iCount))
        {
            bin.erase(bin.begin()+iCount);
            set_bin(particle);
            return;
        }
    }

    printf("Error: Could not find particle in bin.\r\n");
    exit(-1);

}

//
// Apply force to all particles in current bin
// Include surrounding bins to address cutoff
//
void apply_force_bin(int binNum, double *dmin, double *davg, int *navg)
{
    vector<particle_t*> bin = binVec.at(binNum);

    // All surrounding bins in a 3x3 square
    int binsToCheck[] = {   binNum - block_row_count - 1,
                            binNum - block_row_count,
                            binNum - block_row_count + 1,
                            binNum - 1,
                            binNum,
                            binNum + 1,
                            binNum + block_row_count - 1,
                            binNum + block_row_count,
                            binNum + block_row_count + 1
                        };

    // Set acceleration of all particles in the bin to 0
    for (int i = 0 ; i < bin.size(); i++)
    {
        bin.at(i)->ax = bin.at(i)->ay = 0;
    }

    // Address each bin fully before moving to next
    for (int i = 0; i < 9; i++)
    {
        int compareBinNum = binsToCheck[i];
        if (compareBinNum >= 0 && compareBinNum < block_row_count*block_row_count)
        {
            vector<particle_t*> compare_bin = binVec.at(compareBinNum);
            for (int j = 0; j < bin.size(); j++)
            {
                for (int k = 0; k < compare_bin.size(); k++)
                {
                    apply_force(*bin.at(j), *compare_bin.at(k),dmin,davg,navg);
                }
            }
        }
    }
}

//
// Print all sizes of bins
// Debug function
//
void print_bins()
{
    for( int i = 0; i < binVec.size(); i++ )
    {
        printf("bin[%i].size = %d,\r\n", i, (int) binVec.at(i).size());
    }
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
	if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }

    // Check if particle moved to new bin
    int frac_x = p.x/size;
    int frac_y = p.y/size;
    int binNum = (frac_x * block_row_count ) + (frac_y * block_row_count * block_row_count );
    if (p.binNum != binNum)
    {
        // printf("wtf binnum = %i, p.binnum = %i\r\n", binNum, p.binNum);
        move_bin(p);
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
