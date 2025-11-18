#include "partitioner.h"
#include "nodelist.h"
#include "geometry.h"

#include <stdio.h>  // For fprintf, stderr
#include <stdlib.h> // For qsort
#include <string.h> // For memset (if needed)

/**
 * @brief Main partitioning function.
 */
void partition_data(NodeList* all_nodes, int* secure_indices, int num_secure_indices,
                    NodeList* healthy_points, NodeList* infected_points) {
    
    // --- 1. Handle Edge Cases ---
    if (num_secure_indices == 0) {
        fprintf(stderr, "Warning: No secure indices provided. All nodes classified as infected.\n");
        for (size_t i = 0; i < all_nodes->size; ++i) {
            add_point(infected_points, all_nodes->points[i]);
        }
        return;
    }

    // --- 2. Extract unique secure nodes ---
    NodeList* secure_nodes_sorted = init_list(num_secure_indices);
    for (int i = 0; i < num_secure_indices; ++i) {
        add_point(secure_nodes_sorted, all_nodes->points[secure_indices[i]]);
    }
    qsort(secure_nodes_sorted->points, secure_nodes_sorted->size, sizeof(Point), compare_points);
    
    NodeList* secure_nodes_unique = get_unique_points(secure_nodes_sorted);
    free_list(secure_nodes_sorted); // Don't need this anymore

    if (secure_nodes_unique->size < 3) {
        fprintf(stderr, "Warning: Only %zu unique secure nodes. "
                        "Cannot form a 2D secure region. All nodes classified as infected.\n", 
                        secure_nodes_unique->size);
        for (size_t i = 0; i < all_nodes->size; ++i) {
            add_point(infected_points, all_nodes->points[i]);
        }
        free_list(secure_nodes_unique);
        return;
    }

    // --- 3. Compute the Convex Hull ---
    NodeList* hull = compute_convex_hull(secure_nodes_unique);
    
    if (hull == NULL) {
         // This can happen if all unique points (>=3) are collinear
        fprintf(stderr, "Warning: All secure nodes are collinear. "
                        "Cannot form a 2D secure region. All nodes classified as infected.\n");
        for (size_t i = 0; i < all_nodes->size; ++i) {
            add_point(infected_points, all_nodes->points[i]);
        }
        free_list(secure_nodes_unique);
        return;
    }

    // --- 4. Classify All Points ---
    for (size_t i = 0; i < all_nodes->size; ++i) {
        if (is_inside(hull, all_nodes->points[i])) {
            add_point(healthy_points, all_nodes->points[i]);
        } else {
            add_point(infected_points, all_nodes->points[i]);
        }
    }

    // --- 5. Cleanup ---
    free_list(secure_nodes_unique);
    free_list(hull);
}