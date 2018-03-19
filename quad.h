#ifndef __CS267_QUAD_H__
#define __CS267_QUAD_H__

typedef struct Quad{
    //low left x, y
    //length of the cell
    double llx;
    double  lly;
    double cl;
}Quad;


//Returns true if (x, y) is in the region, and false otherwise.
bool Quadcontains(Quad *q, double x, double y);


//These four methods create and return a new Quad representing a sub-quadrant of the invoking quadrant.
Quad* NW(Quad *q);
Quad* NE(Quad *q);
Quad* SW(Quad *q);
Quad* SE(Quad *q);

#endif