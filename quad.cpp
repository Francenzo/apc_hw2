#include "quad.h"
#include "stdio.h"



bool Quadcontains(Quad *q, double x, double y)
{
    if (q == NULL)
    {
        printf("the quad is null\n");
        return false;
    }
    double cl = q->cl;
    //up right x, y
    double urx = q->llx + cl;
    double ury = q->lly + cl;
    if (urx >= x && ury >= y && x>=q->llx && y>=q->lly)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Quad *NW(Quad *q)
{
    if (q == NULL)
    {
        //printf("the quad is null in NW\n");
        return NULL;
    }
    //printf("funcNW %f %f %f\n",q->llx,q->lly,q->cl);
    Quad *newq = new (Quad);
    
    newq->llx = q->llx;
    newq->lly = q->lly + (q->cl / 2.0);
    newq->cl = 0.5 * q->cl;

    //printf("funcNW newq %f %f %f\n",newq->llx,newq->lly,q->cl);
    return newq;
}

Quad *NE(Quad *q)
{
    if (q == NULL)
    {
        printf("the quad is null in NE\n");
        return NULL;
    }
    Quad *newq = new (Quad);
    
    newq->llx = q->llx + (q->cl / 2.0);
    newq->lly = q->lly + (q->cl / 2.0);
    newq->cl = 0.5 * q->cl;
    return newq;
}

Quad *SW(Quad *q)
{
    if (q == NULL)
    {
        printf("the quad is null in SW\n");
        return NULL;
    }
    Quad *newq = new (Quad);
    
    newq->llx = q->llx;
    newq->lly = q->lly;
    newq->cl = 0.5 * q->cl;
    return newq;
}

Quad *SE(Quad *q)
{
    if (q == NULL)
    {
        printf("the quad is null in SE\n");
        return NULL;
    }
    Quad *newq = new (Quad);
    
    newq->llx = q->llx + (q->cl/ 2.0);
    newq->lly = q->lly;
    newq->cl = 0.5 * q->cl;
    return newq;
}