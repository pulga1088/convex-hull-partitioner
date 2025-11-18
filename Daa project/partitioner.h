#ifndef PARTITIONER_H
#define PARTITIONER_H

#include "nodelist.h" // Needs NodeList definition

/**
 * @brief Main partitioning function.
 * Takes all nodes and secure indices, and populates the healthy/infected lists.
 */
void partition_data(NodeList* all_nodes, int* secure_indices, int num_secure_indices,
                    NodeList* healthy_points, NodeList* infected_points);

#endif // PARTITIONER_H