#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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


static double distance_sq_pbc(int i, int j, const double *position, const int *pbc, const double *cell) {
        double dx = position[3*i] - position[3*j];
        double dy = position[3*i+1] - position[3*j+1];
        double dz = position[3*i+2] - position[3*j+2]; 

        if (pbc[0]) {
            double Lx = cell[0];
            dx -= round(dx / Lx) * Lx;
        }

        if (pbc[1]) {
        double Ly = cell[4];  
        dy -= round(dy / Ly) * Ly;
        }   

        if (pbc[2]) {
        double Lz = cell[8];  
        dz -= round(dz / Lz) * Lz;
        }
    
    return dx*dx + dy*dy + dz*dz;
}


NeighborList primitive_neighbor_list(const double *position, int N, const NeighborListConfig *config) {
    NeighborList nl;
    nl.pairs = NULL;
    nl.count = 0;

    const CutoffSpec *cutoff_spec = config->cutoff_spec;
    int self_interaction = config->self_interaction;
    int bothways = config->bothways;
    const int *pbc = config->pbc;
    const double *cell = config->cell;


    int count = 0;

    if (self_interaction) {
        count += N;  
    }

    for (int i = 0; i < N; i++) {

        for (int j = i+1; j < N; j++){

            double dist = distance_sq_pbc(i, j, position, config->pbc, config->cell);

            double cutoff_sum = get_cutoff_sum(cutoff_spec, i, j);

            if (dist <= cutoff_sum*cutoff_sum) {
                if (bothways) {
                    count += 2;
                } else {
                    count += 1;
                }
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

        for (int j = i+1; j < N; j++) {

            double dist = distance_sq_pbc(i, j, position, config->pbc, config->cell);

            double cutoff_sum = get_cutoff_sum(cutoff_spec, i, j);
            
            if (dist <= cutoff_sum*cutoff_sum) {
                nl.pairs[idx].i = i;
                nl.pairs[idx].j = j;
                idx++;

                if (bothways) {
                    nl.pairs[idx].i = j;
                    nl.pairs[idx].j = i;
                    idx++;
                }
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



