#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include "stdio.h"
#include "quad.h"

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

#include <vector>

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  int binNum;
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;


// quard using
// Returns true if the body is in quadrant q.
bool particlein(particle_t* p,Quad *q);

// return a new particle_t that represents the center-of-mass of the two bodies a and b. 
particle_t* aggregate(particle_t* a, particle_t* b);
particle_t *newparticle(Quad *q);

//
//  timing routines
//
double read_timer( );

//
//  bin block routines
//
int get_bin_count();
void make_bin(int n);
void clear_bins(int begin, int end);
void clear_bin(int index);
void set_bin(particle_t & particle);
void move_bin(particle_t & particle);
void apply_force_bin(int binNum, double *dmin, double *davg, int *navg);
void print_bins();

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void init_particles_inthread(int num, Quad *initq,particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move(particle_t &p);
void move_mpi(particle_t &p, double size);


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
