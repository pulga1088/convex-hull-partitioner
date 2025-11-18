#include <stdio.h>
#include "nodelist.h"
#include "partitioner.h"

// --- --- --- --- --- --- --- --- ---
// --- Example Usage ---
// --- --- --- --- --- --- --- --- ---
int main() {
    
    // 1. Define the data
    Point node_data[] = {
        {1, 1}, {1, 6}, {3, 4}, {4, 2}, {4, 7},
        {6, 5}, {7, 1}, {8, 8}, {10, 3}, {11, 6},
        {2, 8}, {5, 10}, {9, 9}
    };
    int num_nodes = sizeof(node_data) / sizeof(Point);
    
    NodeList* all_nodes = init_list(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        add_point(all_nodes, node_data[i]);
    }
    
    // Indices of nodes known to be secure
    int secure_indices[] = {0, 1, 4, 6, 7, 11}; 
    int num_secure = sizeof(secure_indices) / sizeof(int);
    
    printf("--- Running Node Partitioner (C Version) ---\n");

    // 2. Prepare lists and run the function
    NodeList* healthy = init_list(num_nodes);
    NodeList* infected = init_list(num_nodes);
    
    partition_data(all_nodes, secure_indices, num_secure, healthy, infected);

    // 3. Print results
    printf("\nTotal nodes: %d\n", num_nodes);
    
    printf("\nHealthy nodes (%zu):\n", healthy->size);
    for (size_t i = 0; i < healthy->size; ++i) {
        printf("  (%.1f, %.1f)\n", healthy->points[i].x, healthy->points[i].y);
    }
        
    printf("\nInfected nodes (%zu):\n", infected->size);
    for (size_t i = 0; i < infected->size; ++i) {
        printf("  (%.1f, %.1f)\n", infected->points[i].x, infected->points[i].y);
    }

    // 4. Free all allocated memory
    free_list(all_nodes);
    free_list(healthy);
    free_list(infected);

    return 0;
}