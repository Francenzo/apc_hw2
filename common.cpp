#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include "quad.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass 0.01
#define cutoff 0.01
#define min_r (cutoff / 100)
#define dt 0.0005

//
//  timer
//
double read_timer()
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if (!initialized)
    {
        gettimeofday(&start, NULL);
        initialized = true;
    }
    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size(int n)
{
    size = sqrt(density * n);
}

void init_particles_inthread(int num, Quad *initq, particle_t *p)
{
    //printf("size %f\n",initq->cl);

    double LOx = initq->llx;
    double LOy = initq->lly;
    double HIx = LOx + initq->cl;
    double HIy = LOy + initq->cl;
    srand(time(0));
    for (int i = 0; i < num; i++)
    {

        double fxfinal = LOx + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (HIx - LOx)));

        double fyfinal = LOy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (HIy - LOy)));
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = fxfinal;
        p[i].y = fyfinal;

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48() * 2 - 1;
        p[i].vy = drand48() * 2 - 1;

        p[i].ax = 0;
        p[i].ay = 0;

        //printf("init %f %f\n", p[i].x, p[i].y);
    }
}

//
//  Initialize the particle positions and velocities
//
void init_particles(int n, particle_t *p)
{
    srand48(time(NULL));

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n + sx - 1) / sx;

    int *shuffle = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++)
        shuffle[i] = i;

    for (int i = 0; i < n; i++)
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48() % (n - i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n - i - 1];

        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size * (1. + (k % sx)) / (1 + sx);
        p[i].y = size * (1. + (k / sx)) / (1 + sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48() * 2 - 1;
        p[i].vy = drand48() * 2 - 1;

        printf("init %f %f\n", p[i].x, p[i].y);
    }
    free(shuffle);
}

//
//  interact two particles
//
void apply_force(particle_t &particle, particle_t &neighbor, double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if (r2 > cutoff * cutoff)
        return;
    if (r2 != 0)
    {
        if (r2 / (cutoff * cutoff) < *dmin * (*dmin))
        {
            *dmin = sqrt(r2) / cutoff;
        }
        (*davg) += sqrt(r2) / cutoff;
        (*navg)++;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    //
    //  very simple short-range repulsive force
    //
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move(particle_t &p, double size)
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //

    
    if (p.ax == 0 && p.ay == 0)
    {
        return;
    }
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    //
    //  bounce from walls
    //
    while (p.x < 0 || p.x > size)
    {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
        //printf("while 1 x (%f) y (%f) vx (%f) vy (%f)\n",p.x,p.y, p.vx,p.vy);
    }
    while (p.y < 0 || p.y > size)
    {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
        // printf("while 2 x (%f) y (%f) vx (%f) vy (%f)\n",p.x,p.y, p.vx,p.vy);
    }
}

//
//  I/O routines
//
void save(FILE *f, int n, particle_t *p)
{
    static bool first = true;
    if (first)
    {
        fprintf(f, "%d %g\n", n, size);
        first = false;
    }
    for (int i = 0; i < n; i++)
        fprintf(f, "%g %g\n", p[i].x, p[i].y);
}

//
//  command line option processing
//
int find_option(int argc, char **argv, const char *option)
{
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], option) == 0)
            return i;
    return -1;
}

int read_int(int argc, char **argv, const char *option, int default_value)
{
    int iplace = find_option(argc, argv, option);
    if (iplace >= 0 && iplace < argc - 1)
        return atoi(argv[iplace + 1]);
    return default_value;
}

char *read_string(int argc, char **argv, const char *option, char *default_value)
{
    int iplace = find_option(argc, argv, option);
    if (iplace >= 0 && iplace < argc - 1)
        return argv[iplace + 1];
    return default_value;
}

bool particlein(particle_t *p, Quad *q)
{
    double cl = q->cl;
    if ((q->llx + cl > p->x) && (q->lly + cl > p->y) && (p->x > q->llx) && (p->y > q->lly))
    {
        return true;
    }

    return false;
}

particle_t *aggregate(particle_t *a, particle_t *b)
{
    particle_t *newp = new (particle_t);
    //assume mass is 0.01
    newp->x = 0.5 * (a->x + b->x);
    newp->y = 0.5 * (a->y + b->y);

    return newp;
}
