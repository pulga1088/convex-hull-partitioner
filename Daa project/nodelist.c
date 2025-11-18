#include "nodelist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * @brief Initializes a new NodeList with a given initial capacity.
 */
NodeList* init_list(size_t initial_capacity) {
    NodeList* list = (NodeList*)malloc(sizeof(NodeList));
    if (!list) {
        perror("Failed to allocate NodeList");
        exit(EXIT_FAILURE);
    }
    // Use 1 as min capacity to avoid 0-size malloc
    size_t real_capacity = initial_capacity > 0 ? initial_capacity : 10;
    list->points = (Point*)malloc(real_capacity * sizeof(Point));
    if (!list->points) {
        perror("Failed to allocate points array");
        free(list);
        exit(EXIT_FAILURE);
    }
    list->size = 0;
    list->capacity = real_capacity;
    return list;
}

/**
 * @brief Adds a point to the NodeList, resizing if necessary.
 */
void add_point(NodeList* list, Point p) {
    if (list->size >= list->capacity) {
        size_t new_capacity = list->capacity * 2;
        Point* new_points = (Point*)realloc(list->points, new_capacity * sizeof(Point));
        if (!new_points) {
            perror("Failed to reallocate points array");
            // Original list->points is still valid, but we exit
            exit(EXIT_FAILURE);
        }
        list->points = new_points;
        list->capacity = new_capacity;
    }
    list->points[list->size++] = p;
}

/**
 * @brief Frees all memory associated with the NodeList.
 */
void free_list(NodeList* list) {
    if (list) {
        free(list->points);
        free(list);
    }
}

/**
 * @brief Comparison function for qsort, sorts by x, then y.
 */
int compare_points(const void* a, const void* b) {
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    if (fabs(p1->x - p2->x) > EPSILON) {
        return (p1->x > p2->x) - (p1->x < p2->x);
    }
    if (fabs(p1->y - p2->y) > EPSILON) {
        return (p1->y > p2->y) - (p1->y < p2->y);
    }
    return 0;
}

/**
 * @brief Gets a new list containing only the unique points from the input.
 * Assumes input list is sorted.
 */
NodeList* get_unique_points(NodeList* sorted_list) {
    if (sorted_list->size == 0) {
        return init_list(0);
    }
    NodeList* unique_list = init_list(sorted_list->size);
    add_point(unique_list, sorted_list->points[0]);

    for (size_t i = 1; i < sorted_list->size; ++i) {
        if (compare_points(&sorted_list->points[i-1], &sorted_list->points[i]) != 0) {
            add_point(unique_list, sorted_list->points[i]);
        }
    }
    return unique_list;
}