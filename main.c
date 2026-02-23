#include <stdio.h>
#include "neighbor_list.h"

int main() {
    // Тестовые данные: 3 атома
    double positions[] = {
        0.0, 0.0, 0.0,   // атом 0
        1.0, 0.0, 0.0,   // атом 1
        2.5, 0.0, 0.0    // атом 2
    };
    int N = 3;

    // ТЕСТ 1: Глобальный cutoff
    printf("=== Тест 1: Глобальный cutoff ===\n");
    CutoffSpec global = cutoff_global(1.5);  // все атомы с радиусом 1.5
    
    NeighborList nl1 = compute_neighbor_list(positions, N, &global);
    printf("Найдено %d пар:\n", nl1.count);
    for (int k = 0; k < nl1.count; k++) {
        printf("  %d -- %d\n", nl1.pairs[k].i, nl1.pairs[k].j);
    }
    free_neighbor_list(&nl1);

    // ТЕСТ 2: Индивидуальные радиусы
    printf("\n=== Тест 2: Индивидуальные радиусы ===\n");
    double radii[] = {0.9, 0.9, 1.6};  // у атома 2 радиус больше
    CutoffSpec per_atom = cutoff_per_atom(radii);
    
    NeighborList nl2 = compute_neighbor_list(positions, N, &per_atom);
    printf("Найдено %d пар:\n", nl2.count);
    for (int k = 0; k < nl2.count; k++) {
        printf("  %d -- %d\n", nl2.pairs[k].i, nl2.pairs[k].j);
    }
    free_neighbor_list(&nl2);

    printf("\n=== Тест 3: Глобальный cutoff 1.5 (исправленный) ===\n");
    CutoffSpec global_fixed = cutoff_global(1.5);
    NeighborList nl3 = compute_neighbor_list(positions, N, &global_fixed);
    printf("Пар: %d\n", nl3.count);
    for (int k = 0; k < nl3.count; k++)
    printf("  %d-%d\n", nl3.pairs[k].i, nl3.pairs[k].j);
    free_neighbor_list(&nl3);

    return 0;
}