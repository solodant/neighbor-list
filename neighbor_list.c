#include <stdio.h>
#include <stdlib.h>
#include "neighbor_list.h"


CutoffSpec cutoff_global(double value) {
    CutoffSpec spec;
    spec.type = CUTOFF_GLOBAL;
    spec.data.global = value;
    return spec;
}


CutoffSpec cutoff_per_atom(const double *values) {
    CutoffSpec spec;
    spec.type = CUTOFF_PER_ATOM;
    spec.data.per_atom = values;  
    return spec;
}


static double get_cutoff_sum(const CutoffSpec *spec, int i, int j) {
    switch (spec->type) {
        case CUTOFF_GLOBAL:
            return spec->data.global; 
        
        case CUTOFF_PER_ATOM:
            return spec->data.per_atom[i] + spec->data.per_atom[j];
    }
    return 0.0;
}


NeighborList primitive_neighbor_list(const double *position, int N, const NeighborListConfig *config) {
    NeighborList nl;
    nl.pairs = NULL;
    nl.count = 0;

    const CutoffSpec *cutoff_spec = config->cutoff_spec;
    int self_interaction = config->self_interaction;

    int count = 0;

    if (self_interaction) {
        count += N;  
    }

    for (int i = 0; i < N; i++) {

        double xi = position[3*i];
        double yi = position[3*i+1];
        double zi = position[3*i+2];

        for (int j = i+1; j < N; j++){
            double dx = xi - position[3*j];
            double dy = yi - position[3*j+1];
            double dz = zi - position[3*j+2];

            double dist = dx*dx + dy*dy + dz*dz;
            double cutoff_sum = get_cutoff_sum(cutoff_spec, i, j);

            if (dist <= cutoff_sum*cutoff_sum) {
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

    if (self_interaction) {
        for (int i = 0; i < N; i++) {
            nl.pairs[idx].i = i;
            nl.pairs[idx].j = i;
            idx++;
        }
    }

    for (int i = 0; i < N; i++) {

        double xi = position[3*i];
        double yi = position[3*i+1];
        double zi = position[3*i+2];

        for (int j = i+1; j < N; j++) {
            double dx = xi - position[3*j];
            double dy = yi - position[3*j+1];
            double dz = zi - position[3*j+2];

            double dist = dx*dx + dy*dy + dz*dz;
            double cutoff_sum = get_cutoff_sum(cutoff_spec, i, j);
            
            if (dist <= cutoff_sum*cutoff_sum) {
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



