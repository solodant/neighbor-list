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


typedef enum {
    CUTOFF_GLOBAL,
    CUTOFF_PER_ATOM,
} CutoffType;


typedef struct {
    CutoffType type;
    union {
        double global;
        const double *per_atom;
    } data;
} CutoffSpec;


CutoffSpec cutoff_global(double value);
CutoffSpec cutoff_per_atom(const double *values);


typedef struct {
    const CutoffSpec *cutoff_spec;
    int self_interaction;
} NeighborListConfig;



NeighborList primitive_neighbor_list(
    const double *positions,
    int N,
    const NeighborListConfig *config
);


void free_neighbor_list(NeighborList *nl);


#endif

