#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H

typedef struct {
    int i;
    int j;
} Pair;


typedef struct {
    Pair *pairs;
    int count;
} NeighborList;


NeighborList compute_neighbor_list(const double *positions, int N, double cutoff);


void free_neighbor_list(NeighborList *nl);


#endif

