#ifndef __CS267_BHTREE_H__
#define __CS267_BHTREE_H__

#include "common.h"
#include "quad.h"


typedef struct BHTreeNode {
    particle_t * particle;     // body or aggregate body stored in this node
    Quad *quad;     // square region that the tree represents
    BHTreeNode* NW;     // tree representing northwest quadrant
    BHTreeNode* NE;     // tree representing northeast quadrant
    BHTreeNode* SW;     // tree representing southwest quadrant
    BHTreeNode* SE;     // tree representing southeast quadrant
}BHTreeNode;


//create a Barnes-Hut tree with no bodies, representing the given quadrant.
BHTreeNode* BHTree(Quad *q);

//add the particle p to the invoking Barnes-Hut tree.
void BHTinsert(BHTreeNode*bht,particle_t *p);

//approximate the net force acting on Body b from all the bodies in the invoking Barnes-Hut tree, and update b's force accordingly.
void BHTupdateForce(BHTreeNode*bht,particle_t *p);


#endif