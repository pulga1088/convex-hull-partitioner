#ifndef NODELIST_H
#define NODELIST_H

#include <stdlib.h> // For size_t
#include <math.h>   // For fabs

// Epsilon for floating point comparisons
#define EPSILON 1e-9

// --- Data Structures ---

typedef struct {
    double x, y;
} Point;

typedef struct {
    Point* points;
    size_t size;
    size_t capacity;
} NodeList;

// --- NodeList (Dynamic Array) Helper Function Declarations ---

/**
 * @brief Initializes a new NodeList with a given initial capacity.
 */
NodeList* init_list(size_t initial_capacity);

/**
 * @brief Adds a point to the NodeList, resizing if necessary.
 */
void add_point(NodeList* list, Point p);

/**
 * @brief Frees all memory associated with the NodeList.
 */
void free_list(NodeList* list);

/**
 * @brief Comparison function for qsort, sorts by x, then y.
 */
int compare_points(const void* a, const void* b);

/**
 * @brief Gets a new list containing only the unique points from the input.
 * Assumes input list is sorted.
 */
NodeList* get_unique_points(NodeList* sorted_list);

#endif // NODELIST_H