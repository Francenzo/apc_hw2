
#include "bhtree.h"
#include "stdio.h"
#include "iostream"
#include <queue>
#include <math.h>
#define cutoff 0.01
#define min_r (cutoff / 100)
#define mass 0.01

using namespace std;

BHTreeNode *BHTree(Quad *q)
{

    BHTreeNode *bhnode = new (BHTreeNode);
    bhnode->quad = q;
    bhnode->particle = NULL;
    bhnode->NE = NULL;
    bhnode->NW = NULL;
    bhnode->SW = NULL;
    bhnode->SE = NULL;
    return bhnode;
}

bool EdgeNode(BHTreeNode *bht)
{
    if (bht->particle == NULL)
    {
        return true;
    }

    if (bht->NE == NULL && bht->NW == NULL && bht->SE == NULL && bht->SW == NULL)
    {
        return true;
    }
    return false;
}

void BHTinsert(BHTreeNode *bht, particle_t *p)
{
    double px = p->x;
    double py = p->y;

    //printf("insert p %f %f\n", px, py);

    if (Quadcontains(bht->quad, px, py) == false)

    {
        //printf("step0\n");
        return;
    }
    else
    {
        //if external (accurate nodes)
        if (EdgeNode(bht) == true)
        {
            if (bht->particle == NULL)
            {
                //printf("add particle directly\n");
                bht->particle = p;
            }
            else
            {

                //not internal, divided further
                //printf("bht quad %f %f\n", bht->quad->llx, bht->quad->lly);
                Quad *nwq = NW(bht->quad);
                Quad *neq = NE(bht->quad);
                Quad *swq = SW(bht->quad);
                Quad *seq = SE(bht->quad);
                //printf("first pa nwq %f %f %f\n", nwq->llx, nwq->lly, nwq->cl);
                bht->NW = BHTree(nwq);
                bht->NE = BHTree(neq);
                bht->SW = BHTree(swq);
                bht->SE = BHTree(seq);
                //aggregate pa pb
                particle_t *pa = bht->particle;
                particle_t *pb = p;

                //put aggregated data into current bht
                particle_t *integratep = aggregate(pa, pb);
                bht->particle = integratep;
                //insert pa and pb into specific for cell
                //printf("pa nwq %f %f %f particle %f %f\n", nwq->llx, nwq->lly, nwq->cl, pa->x, pa->y);
                if (Quadcontains(nwq, pa->x, pa->y) == true)
                {
                    //printf("nwq\n");
                    //printf("pa nwq %f %f %f particle %f %f\n", nwq->llx, nwq->lly, nwq->cl, pa->x, pa->y);
                    BHTinsert(bht->NW, pa);
                }
                if (Quadcontains(neq, pa->x, pa->y) == true)
                {
                    //printf("neq\n");
                    //printf("pa neq %f %f %f particle %f %f\n", neq->llx, neq->lly, neq->cl, pa->x, pa->y);
                    BHTinsert(bht->NE, pa);
                }
                if (Quadcontains(swq, pa->x, pa->y) == true)
                {
                    //printf("swq\n");
                    //printf("pa swq %f %f %f particle %f %f\n", swq->llx, swq->lly, swq->cl, pa->x, pa->y);
                    BHTinsert(bht->SW, pa);
                }
                if (Quadcontains(seq, pa->x, pa->y) == true)
                {
                    //printf("seq\n");
                    //printf("pa seq %f %f %f particle %f %f\n", seq->llx, seq->lly, seq->cl, pa->x, pa->y);
                    BHTinsert(bht->SE, pa);
                }

                //insert pb
                //printf("pb nwq %f %f particle %f %f\n", nwq->llx, nwq->lly, pb->x, pb->y);
                if (Quadcontains(nwq, pb->x, pb->y) == true)
                {
                    BHTinsert(bht->NW, pb);
                }
                if (Quadcontains(neq, pb->x, pb->y) == true)
                {
                    BHTinsert(bht->NE, pb);
                }
                if (Quadcontains(swq, pb->x, pb->y) == true)
                {
                    BHTinsert(bht->SW, pb);
                }
                if (Quadcontains(seq, pb->x, pb->y) == true)
                {
                    BHTinsert(bht->SE, pb);
                }

                return;
            }
        }
        else
        {
            //internal particles
            //printf("update particle\n");
            bht->particle->x = 0.5 * (bht->particle->x + p->x);
            bht->particle->y = 0.5 * (bht->particle->y + p->y);

            //find more accurate one
            if (bht->NW != NULL)
            {
                //printf("step1\n");
                BHTinsert(bht->NW, p);
            }
            if (bht->NE != NULL)
            {
                //printf("step2\n");
                BHTinsert(bht->NE, p);
            }
            if (bht->SE != NULL)
            {
                //printf("step3\n");
                BHTinsert(bht->SE, p);
            }
            if (bht->SW != NULL)
            {
                //printf("step4\n");
                BHTinsert(bht->SW, p);
            }
        }
    }
}

void printBT(BHTreeNode *bnode, queue<BHTreeNode *> &q)
{

    //printf("r %f %f %f\n", bnode->quad->llx, bnode->quad->lly, bnode->quad->cl);
    if (bnode->NW != NULL)
    {
        q.push(bnode->NW);
    }
    if (bnode->NE != NULL)
    {
        q.push(bnode->NE);
    }
    if (bnode->SW != NULL)
    {
        q.push(bnode->SW);
    }
    if (bnode->SE != NULL)
    {
        q.push(bnode->SE);
    }

    BHTreeNode *qb = q.front();

    if (qb->particle != NULL)
    {
        printf("loc %f %f %f particle %f %f\n", qb->quad->llx, qb->quad->lly, qb->quad->cl, qb->particle->x, qb->particle->y);
    }
    else
    {
        printf("null node\n");
    }

    q.pop();
    if (q.size() == 0)
    {
        return;
    }
    else
    {
        printBT(q.front(), q);
    }
}

void applyForceTree(particle_t *particle, BHTreeNode *bht)
{
    BHTreeNode *curr = bht;
    double dx = curr->quad->llx - particle->x;
    double dy = curr->quad->lly - particle->y;
    double r2 = dx * dx + dy * dy;
    if (r2 > cutoff * cutoff)
    {
        //do not compute the under branch
        return;
    }
    else
    {

        r2 = fmax(r2, min_r * min_r);
        double r = sqrt(r2);

        //
        //  very simple short-range repulsive force
        //
        double coef = (1 - cutoff / r) / r2 / mass;
        particle->ax += coef * dx;
        particle->ay += coef * dy;
        //range the child node
        if (curr->NW != NULL)
        {
            applyForceTree(particle, curr->NW);
        }
        if (curr->NE != NULL)
        {
            applyForceTree(particle, curr->NE);
        }
        if (curr->SE != NULL)
        {
            applyForceTree(particle, curr->SE);
        }
        if (curr->SW != NULL)
        {
            applyForceTree(particle, curr->SW);
        }
    }
}

void BHTfree(BHTreeNode *bht)
{

    if (bht->NE != NULL)
    {
        BHTfree(bht->NE);
    }
    if (bht->NW != NULL)
    {
        BHTfree(bht->NW);
    }
    if (bht->SE != NULL)
    {
        BHTfree(bht->SE);
    }
    if (bht->SW != NULL)
    {
        BHTfree(bht->SW);
    }

    delete (bht);
}

//test using
/*
int main()
{
    Quad *q = new (Quad);
    q->cl = 1;
    q->llx = 0;
    q->lly = 0;

    cout << Quadcontains(q, 5.0, 5.0) << endl;
    cout << Quadcontains(q, 0.25, 0.25) << endl;

    Quad *nwq = NW(q);
    Quad *neq = NE(q);
    Quad *swq = SW(q);
    Quad *seq = SE(q);
    printf("nwq %f %f %f\n", nwq->llx, nwq->lly, nwq->cl);
    printf("neq %f %f %f\n", neq->llx, neq->lly, neq->cl);
    printf("swq %f %f %f\n", swq->llx, swq->lly, swq->cl);
    printf("seq %f %f %f\n", seq->llx, seq->lly, seq->cl);

    //test bh tree

    BHTreeNode *bnode = BHTree(q);

    printf("node quadra %f %f %f\n", bnode->quad->llx, bnode->quad->lly, bnode->quad->cl);

    //test tree insert

    particle_t *parray = (particle_t *)malloc(sizeof(particle_t) * 10);
    parray[0].x = 0.2;
    parray[0].y = 0.2;

    parray[1].x = 0.2;
    parray[1].y = 0.6;

    parray[2].x = 0.6;
    parray[2].y = 0.6;

    parray[3].x = 0.9;
    parray[3].y = 0.9;

    parray[4].x = 0.9;
    parray[4].y = 0.2;

    int i = 0;
    queue<BHTreeNode *> que;

    for (i = 0; i < 4; i++)
    {
        BHTinsert(bnode, parray + i);

        printf("\n");

        if (que.size() == 0)
        {
            que.push(bnode);

            printBT(bnode, que);
        }

        printf("\n");
    }

    BHTfree(bnode);
}
*/