#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "neighbor_list.h"

NeighborList compute_neighbor_list(const double *position, int N, double cutoff) {
    NeighborList nl;
    nl.pairs = NULL;
    nl.count = 0;

    int count = 0;

    for (int i = 0; i < N; i++) {

        double xi = position[3*i];
        double yi = position[3*i+1];
        double zi = position[3*i+2];

        for (int j = i+1; j < N; j++){
            double dx = xi - position[3*j];
            double dy = yi - position[3*j+1];
            double dz = zi - position[3*j+2];

            double dist = dx*dx + dy*dy + dz*dz;

            if (dist <= cutoff*cutoff) {
                count++;
            }
        }
    }

    if (count == 0) {
        return nl;
    }


    nl.pairs = (Pair*)malloc(count * sizeof(Pair));
    if (nl.pairs == NULL) {
        return nl;
    }
    nl.count = count;


    int idx = 0;
    for (int i = 0; i < N; i++) {

        double xi = position[3*i];
        double yi = position[3*i+1];
        double zi = position[3*i+2];

        for (int j = i+1; j < N; j++) {
            double dx = xi - position[3*j];
            double dy = yi - position[3*j+1];
            double dz = zi - position[3*j+2];

            double dist = dx*dx + dy*dy + dz*dz;
            
            if (dist <= cutoff*cutoff) {
                nl.pairs[idx].i = i;
                nl.pairs[idx].j = j;
                idx++;
            }
        }
    }


    return nl;
}


void free_neighbor_list(NeighborList *nl) {
    if (nl->pairs != NULL) {
        free(nl->pairs);
        nl->pairs = NULL;
    }
    nl->count = 0;
}



