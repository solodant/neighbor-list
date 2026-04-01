#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "neighbor_list.h"


static int compare_pairs(const void *a, const void *b) {
    const Pair *pa = (const Pair*)a;
    const Pair *pb = (const Pair*)b;
    if (pa->i != pb->i) return pa->i - pb->i;
    return pa->j - pb->j;
}


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


NeighborListObject* neighborlist_create(const double *cutoffs, int natoms, 
                                        int self_interaction, int bothways,
                                        double skin, int sorted) {

    NeighborListObject *nl = (NeighborListObject*)malloc(sizeof(NeighborListObject));
    if (nl == NULL) {
        return NULL;
    }

    nl->cutoffs = (double*)malloc(natoms * sizeof(double));
    if (nl->cutoffs == NULL) {
        free(nl);
        return NULL;
    }
    for (int i = 0; i < natoms; i++) {
        nl->cutoffs[i] = cutoffs[i];
    }


    nl->natoms = natoms;
    nl->self_interaction = self_interaction;
    nl->bothways = bothways;
    nl->skin = skin;
    nl->sorted = sorted;
    nl->cached_nl = NULL;
    nl->nupdates = 0;
    nl->last_positions = NULL;


    for (int i = 0; i < 3; i++) {
        nl->last_pbc[i] = 0;
    }
    for (int i = 0; i < 9; i++) {
        nl->last_cell[i] = (i % 4 == 0) ? 1.0 : 0.0;  
    }
    

    return nl;
}


void neighborlist_free(NeighborListObject *nl) {
    if (nl == NULL) return;
    
    if (nl->cutoffs != NULL) {
        free(nl->cutoffs);
    }
    
    if (nl->cached_nl != NULL) {
        free_neighbor_list(nl->cached_nl);
        free(nl->cached_nl);
    }

    if (nl->last_positions != NULL) {
        free(nl->last_positions);
    }
    
    free(nl);
}


void neighborlist_update(NeighborListObject *nl, const double *positions, 
                         const int *pbc, const double *cell) {
    if (nl == NULL) return;


    int need_rebuild = 0;

    if (nl->nupdates == 0) {
        need_rebuild = 1;
    } else {
        if (pbc[0] != nl->last_pbc[0] || pbc[1] != nl->last_pbc[1] || pbc[2] != nl->last_pbc[2]) {
            need_rebuild = 1;
        }

        if (!need_rebuild) {
            for (int k = 0; k < 9; k++) {
                if (cell[k] != nl->last_cell[k]) {
                    need_rebuild = 1;
                    break;
                }
            }
        }

        if (!need_rebuild && nl->last_positions != NULL) {
            double max_move_sq = 0.0;
            for (int i = 0; i < nl->natoms; i++) {
                double dx = positions[3*i] - nl->last_positions[3*i];
                double dy = positions[3*i+1] - nl->last_positions[3*i+1];
                double dz = positions[3*i+2] - nl->last_positions[3*i+2];
                double move_sq = dx*dx + dy*dy + dz*dz;
                if (move_sq > max_move_sq) max_move_sq = move_sq;
            }
            if (sqrt(max_move_sq) > nl->skin) {
                need_rebuild = 1;
            }
        }
    }


    if (!need_rebuild) {
        return;
    }


    CutoffSpec cutoff_spec = cutoff_per_atom(nl->cutoffs);


    NeighborListConfig config = {
        .cutoff_spec = &cutoff_spec,
        .self_interaction = nl->self_interaction,
        .bothways = nl->bothways,
        .pbc = {pbc[0], pbc[1], pbc[2]},
        .cell = {cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]}
    };


    NeighborList new_nl = primitive_neighbor_list(positions, nl->natoms, &config);


    if (nl->sorted && new_nl.count > 0) {
        qsort(new_nl.pairs, new_nl.count, sizeof(Pair), compare_pairs);
    }


    if (nl->cached_nl != NULL) {
        free_neighbor_list(nl->cached_nl);
        free(nl->cached_nl);
    }


    nl->cached_nl = (NeighborList*)malloc(sizeof(NeighborList));
    if (nl->cached_nl == NULL) {
        free_neighbor_list(&new_nl);
        return;
    }
    *nl->cached_nl = new_nl;


    if (nl->last_positions != NULL) {
        free(nl->last_positions);
    }
    nl->last_positions = (double*)malloc(3 * nl->natoms * sizeof(double));
    if (nl->last_positions != NULL) {
        for (int i = 0; i < 3 * nl->natoms; i++) {
            nl->last_positions[i] = positions[i];
        }
    }


    for (int i = 0; i < 3; i++) {
        nl->last_pbc[i] = pbc[i];
    }
    for (int i = 0; i < 9; i++) {
        nl->last_cell[i] = cell[i];
    }


    nl->nupdates++;
}


AtomNeighbors neighborlist_get_neighbors(NeighborListObject *nl, int atom) {
    AtomNeighbors result;
    result.indices = NULL;
    result.count = 0;


    if (nl == NULL || nl->cached_nl == NULL || nl->cached_nl->pairs == NULL) {
        return result;
    }

    NeighborList *nlist = nl->cached_nl;


    int count = 0;
    for (int k = 0; k < nlist->count; k++) {
        if (nlist->pairs[k].i == atom) {
            count++;
        }
    }

    if (count == 0) {
        return result;
    }


    result.indices = (int*)malloc(count * sizeof(int));
    if (result.indices == NULL) {
        return result;
    }
    result.count = count;

    int idx = 0;
    for (int k = 0; k < nlist->count; k++) {
        if (nlist->pairs[k].i == atom) {
            result.indices[idx++] = nlist->pairs[k].j;
        }
    }

    return result;
}


void atom_neighbors_free(AtomNeighbors *neighbors) {
    if (neighbors->indices != NULL) {
        free(neighbors->indices);
        neighbors->indices = NULL;
    }
    neighbors->count = 0;
}


int neighborlist_get_nupdates(NeighborListObject *nl) {
    return nl->nupdates;
}


void free_neighbor_list(NeighborList *nl) {
    if (nl->pairs != NULL) {
        free(nl->pairs);
        nl->pairs = NULL;
    }
    nl->count = 0;
}
