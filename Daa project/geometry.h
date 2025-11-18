#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "nodelist.h" // Needs Point and NodeList definitions

// --- Geometric Helper Function Declarations ---

/**
 * @brief Calculates the 2D cross product of vectors (o->p1) and (o->p2).
 * [Image of 2D cross product formula and vector diagram]
 * > 0: p2 is left of vector o->p1 (counter-clockwise)
 * < 0: p2 is right of vector o->p1 (clockwise)
 * = 0: p1, p2, and o are collinear
 */
double cross_product(Point o, Point p1, Point p2);

/**
 * @brief Computes the convex hull of a set of points using the Monotone Chain algorithm.
 * Returns a new NodeList containing the hull vertices in counter-clockwise order.
 */
NodeList* compute_convex_hull(NodeList* nodes);

/**
 * @brief Checks if a point is inside or on the boundary of a convex hull.
 * 
 */
int is_inside(NodeList* hull, Point p);

#endif // GEOMETRY_H