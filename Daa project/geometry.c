    #include "geometry.h"
    #include <stdlib.h> // For qsort
    #include <math.h>   // For fabs

    /**
    * @brief Calculates the 2D cross product of vectors (o->p1) and (o->p2).
    */
    double cross_product(Point o, Point p1, Point p2) {
        return (p1.x - o.x) * (p2.y - o.y) - (p1.y - o.y) * (p2.x - o.x);
    }

    /**
    * @brief Computes the convex hull of a set of points using the Monotone Chain algorithm.
    */
    NodeList* compute_convex_hull(NodeList* nodes) {
        if (nodes->size < 3) {
            return NULL; // Not enough points
        }

        // Sort points lexicographically
        qsort(nodes->points, nodes->size, sizeof(Point), compare_points);

        NodeList* hull = init_list(nodes->size * 2);
        
        // Build lower hull
        for (size_t i = 0; i < nodes->size; ++i) {
            while (hull->size >= 2 && 
                cross_product(hull->points[hull->size - 2], hull->points[hull->size - 1], nodes->points[i]) <= -EPSILON) {
                hull->size--; // Pop
            }
            add_point(hull, nodes->points[i]);
        }

        // Build upper hull
        size_t lower_hull_size = hull->size;
        for (size_t i = nodes->size - 1; i > 0; --i) {
            while (hull->size > lower_hull_size && 
                cross_product(hull->points[hull->size - 2], hull->points[hull->size - 1], nodes->points[i-1]) <= -EPSILON) {
                hull->size--; // Pop
            }
            add_point(hull, nodes->points[i-1]);
        }

        // Remove the last point, which is a duplicate of the first
        if (hull->size > 1) {
            hull->size--;
        }

        // Check for collinear case (hull will have < 3 points)
        if (hull->size < 3) {
            free_list(hull);
            return NULL;
        }

        return hull;
    }

    /**
    * @brief Checks if a point is inside or on the boundary of a convex hull.
    */
    int is_inside(NodeList* hull, Point p) {
        if (hull == NULL || hull->size < 3) {
            return 0; // Cannot be inside a line or a point
        }

        for (size_t i = 0; i < hull->size; ++i) {
            Point p1 = hull->points[i];
            Point p2 = hull->points[(i + 1) % hull->size]; // Wrap around
            
            double cp = cross_product(p1, p2, p);
            
            // If point is to the "right" of the edge (with tolerance)
            if (cp > EPSILON) {
                return 0; // Outside
            }
        }
        return 1; // Inside or on boundary
    }