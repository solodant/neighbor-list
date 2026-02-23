#include <stdio.h>
#include "neighbor_list.h"

int main() {
    double positions[] = {
        0,0,0, 
        1,0,0, 
        2.5,1,0, 
    };
    int N = 3;
    double cutoff = 1.5;

    NeighborList nl = compute_neighbor_list(positions, N, cutoff);

    printf("Found %d pairs:\n", nl.count);
    for (int k = 0; k < nl.count; k++) {
        printf("  %d -- %d\n", nl.pairs[k].i, nl.pairs[k].j);
    }

    free_neighbor_list(&nl);
    return 0;
}