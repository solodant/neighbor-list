/**
 * @file neighbor_list.h
 * @brief Neighbor list library for atomic systems
 * 
 * Implements neighbor search algorithms similar to ASE neighborlist.
 * Supports global and per-atom cutoffs, periodic boundary conditions,
 * lazy updates (skin), and sorted neighbor lists.
 * 
 * @author Anton Solodukhin
 * @date 2026
 */


#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H


/**
 * @brief A pair of atom indices representing a neighbor bond
 */
typedef struct {
    int i; /**< First atom index */
    int j; /**< Second atom index */
} Pair;


/**
 * @brief Result of neighbor search containing all pairs
 */
typedef struct {
    Pair *pairs; /**< Array of neighbor pairs */
    int count;   /**< Number of pairs */
} NeighborList;


/**
 * @brief Type of cutoff specification
 */
typedef enum {
    CUTOFF_GLOBAL,   /**< Single cutoff value for all atoms */
    CUTOFF_PER_ATOM, /**< Individual cutoff radius per atom */
} CutoffType;


/**
 * @brief Flexible cutoff specification (global or per-atom)
 */
typedef struct {
    CutoffType type; /**< Type of this specification */
    union {
        double global;          /**< Value for CUTOFF_GLOBAL */
        const double *per_atom; /**< Array for CUTOFF_PER_ATOM (size natoms) */
    } data;
} CutoffSpec;


/**
 * @brief Create a global cutoff specification
 * @param value Single cutoff value applied to all atoms
 * @return CutoffSpec configured for global cutoff
 */
CutoffSpec cutoff_global(double value);


/**
 * @brief Create a per-atom cutoff specification
 * @param values Array of cutoff radii (one per atom)
 * @return CutoffSpec configured for per-atom cutoff
 * @note The array is not copied; caller must ensure it remains valid
 */
CutoffSpec cutoff_per_atom(const double *values);


/**
 * @brief Configuration for primitive_neighbor_list
 */
typedef struct {
    const CutoffSpec *cutoff_spec; /**< Cutoff specification */
    int self_interaction;          /**< 1 to include (i,i) pairs, 0 otherwise */
    int bothways;                  /**< 1 to return both (i,j) and (j,i), 0 for i<j only */
    int pbc[3];                    /**< Periodic boundary flags (1=periodic, 0=non-periodic) */
    double cell[9];                /**< Unit cell vectors (3x3 matrix, row-major) */
} NeighborListConfig;


/**
 * @brief Core neighbor search function (O(N²) implementation)
 * 
 * Computes all neighbor pairs for the given atomic configuration.
 * For periodic systems, applies minimum image convention (orthogonal cells only).
 * 
 * @param positions Array of coordinates [x0,y0,z0, x1,y1,z1, ...]
 * @param natoms Number of atoms
 * @param config Configuration parameters (cutoff, PBC, etc.)
 * @return NeighborList containing all pairs (memory allocated internally)
 * 
 * @note Call free_neighbor_list() to release memory
 * @warning Current implementation is O(N²); use with large systems cautiously
 */
NeighborList primitive_neighbor_list(const double *positions, int natoms,
                                     const NeighborListConfig *config);


/**
 * @brief Free memory allocated by primitive_neighbor_list
 * @param nl NeighborList to free (can be NULL)
 */                                     
void free_neighbor_list(NeighborList *nl);


/**
 * @brief Result of neighborlist_get_neighbors
 */
typedef struct {
    int *indices; /**< Array of neighbor indices */
    int count;    /**< Number of neighbors */
} AtomNeighbors;


/**
 * @brief High-level neighbor list object with state caching
 * 
 * Analogous to ASE's NeighborList class. Stores parameters and caches results
 * for lazy updates using the skin parameter.
 */
typedef struct {
    /* Parameters (set at creation) */
    double *cutoffs;         /**< Copy of cutoff radii (owned by this object) */
    int natoms;              /**< Number of atoms */
    int self_interaction;    /**< Self-interaction flag */
    int bothways;            /**< Bothways flag */
    double skin;             /**< Skin distance for lazy updates */
    int sorted;              /**< Sort neighbor pairs flag */

    /* Cached state */
    NeighborList *cached_nl; /**< Cached neighbor list result */
    int nupdates;            /**< Number of rebuilds performed */

    /* Saved state for lazy update checks */
    int last_pbc[3];         /**< Last PBC flags used */
    double last_cell[9];     /**< Last cell matrix used */
    double *last_positions;  /**< Last positions used (for skin check) */
} NeighborListObject;


/**
 * @brief Create a new NeighborListObject
 * 
 * @param cutoffs Array of cutoff radii (size natoms)
 * @param natoms Number of atoms
 * @param self_interaction 1 to include (i,i) pairs, 0 otherwise
 * @param bothways 1 to return both orientations, 0 for i<j only
 * @param skin Skin distance for lazy updates (0 = always rebuild)
 * @param sorted 1 to sort neighbor pairs, 0 otherwise
 * @return New NeighborListObject, or NULL on allocation failure
 * 
 * @note The cutoffs array is copied; caller can free it after creation
 */
NeighborListObject* neighborlist_create(const double *cutoffs, int natoms, 
                                        int self_interaction, int bothways,
                                        double skin, int sorted);

            
/**
 * @brief Free a NeighborListObject and all associated memory
 * @param nl Object to free (can be NULL)
 */
void neighborlist_free(NeighborListObject *nl);


/**
 * @brief Update the neighbor list with new atomic positions
 * 
 * Rebuilds the neighbor list only if:
 *   - First update
 *   - PBC flags changed
 *   - Cell matrix changed
 *   - Any atom moved more than skin distance
 * 
 * @param nl NeighborListObject
 * @param positions New atomic positions [x0,y0,z0, ...]
 * @param pbc New PBC flags [x,y,z] (1=periodic)
 * @param cell New cell matrix (3x3, row-major)
 */
void neighborlist_update(NeighborListObject *nl, const double *positions,
                         const int *pbc, const double *cell);


/**
 * @brief Get neighbors for a specific atom
 * 
 * @param nl NeighborListObject (must have been updated)
 * @param atom Atom index
 * @return AtomNeighbors containing neighbor indices (memory allocated)
 * 
 * @note Call atom_neighbors_free() to release memory
 */
AtomNeighbors neighborlist_get_neighbors(NeighborListObject *nl, int atom);


/**
 * @brief Free memory allocated by neighborlist_get_neighbors
 * @param neighbors AtomNeighbors to free (can be NULL)
 */
void atom_neighbors_free(AtomNeighbors *neighbors);


/**
 * @brief Get the number of updates performed
 * @param nl NeighborListObject
 * @return Number of times the neighbor list was rebuilt
 */
int neighborlist_get_nupdates(NeighborListObject *nl);


#endif /* NEIGHBOR_LIST_H */
