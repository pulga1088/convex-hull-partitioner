// Required includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#define COLOR_RESET "\x1b[0m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_RED "\x1b[31m"
#define COLOR_CYAN "\x1b[36m"
#define COLOR_BLUE "\x1b[34m"

// --- Data Structures ---

// Structure to represent a 2D point
typedef struct
{
    double x;
    double y;
} Point;

// Structure to represent a line equation (Ax + By + C = 0)
typedef struct
{
    double A;
    double B;
    double C;
} LineEquation;

// Structure to hold partition results
typedef struct
{
    Point *healthy_points;
    int healthy_count;
    Point *infected_points;
    int infected_count;
} PartitionResult;

// Structure to hold convex hull result
typedef struct
{
    Point *vertices;
    int vertex_count;
    LineEquation *equations;
    int equation_count;
    bool is_valid;
} ConvexHull;

typedef enum
{
    CELL_EMPTY = 0,
    CELL_HULL,
    CELL_SECURE,
    CELL_HEALTHY,
    CELL_INFECTED
} CellType;

typedef struct
{
    CellType type;
} PlotCell;

// --- Utility Functions ---

/**
 * Compares two points for equality with tolerance
 */
bool points_equal(Point p1, Point p2, double tolerance)
{
    return fabs(p1.x - p2.x) < tolerance && fabs(p1.y - p2.y) < tolerance;
}

bool point_in_list(Point p, Point *list, int count, double tolerance)
{
    for (int i = 0; i < count; i++)
    {
        if (points_equal(p, list[i], tolerance))
        {
            return true;
        }
    }
    return false;
}

bool is_secure_point(Point p, Point *all_nodes, int *secure_indices, int secure_count, double tolerance)
{
    if (!secure_indices || secure_count == 0)
    {
        return false;
    }
    for (int i = 0; i < secure_count; i++)
    {
        if (points_equal(p, all_nodes[secure_indices[i]], tolerance))
        {
            return true;
        }
    }
    return false;
}

void draw_line_on_grid(PlotCell *grid, int width, int height, int x0, int y0, int x1, int y1)
{
    int dx = abs(x1 - x0);
    int sx = (x0 < x1) ? 1 : -1;
    int dy = -abs(y1 - y0);
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx + dy;

    while (true)
    {
        if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height)
        {
            PlotCell *cell = &grid[y0 * width + x0];
            if (cell->type == CELL_EMPTY)
            {
                cell->type = CELL_HULL;
            }
        }
        if (x0 == x1 && y0 == y1)
        {
            break;
        }
        int e2 = 2 * err;
        if (e2 >= dy)
        {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx)
        {
            err += dx;
            y0 += sy;
        }
    }
}

void render_ascii_plot(Point *all_nodes, int node_count,
                       int *secure_indices, int secure_count,
                       PartitionResult *partition, ConvexHull *hull)
{
    if (!hull || !hull->is_valid)
    {
        printf("\nConvex hull plot unavailable (invalid hull).\n");
        return;
    }

    double min_x = DBL_MAX;
    double min_y = DBL_MAX;
    double max_x = -DBL_MAX;
    double max_y = -DBL_MAX;

    for (int i = 0; i < node_count; i++)
    {
        if (all_nodes[i].x < min_x)
            min_x = all_nodes[i].x;
        if (all_nodes[i].x > max_x)
            max_x = all_nodes[i].x;
        if (all_nodes[i].y < min_y)
            min_y = all_nodes[i].y;
        if (all_nodes[i].y > max_y)
            max_y = all_nodes[i].y;
    }

    for (int i = 0; i < hull->vertex_count; i++)
    {
        Point v = hull->vertices[i];
        if (v.x < min_x)
            min_x = v.x;
        if (v.x > max_x)
            max_x = v.x;
        if (v.y < min_y)
            min_y = v.y;
        if (v.y > max_y)
            max_y = v.y;
    }

    int margin = 1;
    int width = (int)ceil(max_x - min_x) + 1 + (margin * 2);
    int height = (int)ceil(max_y - min_y) + 1 + (margin * 2);

    if (width < 3)
    {
        width = 3;
    }
    if (height < 3)
    {
        height = 3;
    }

    PlotCell *grid = (PlotCell *)calloc(width * height, sizeof(PlotCell));
    if (!grid)
    {
        fprintf(stderr, "Error: Unable to allocate plot grid.\n");
        return;
    }

    for (int i = 0; i < hull->vertex_count; i++)
    {
        Point a = hull->vertices[i];
        Point b = hull->vertices[(i + 1) % hull->vertex_count];
        int x0 = (int)llround(a.x - min_x) + margin;
        int y0 = (int)llround(a.y - min_y) + margin;
        int x1 = (int)llround(b.x - min_x) + margin;
        int y1 = (int)llround(b.y - min_y) + margin;
        draw_line_on_grid(grid, width, height, x0, y0, x1, y1);
    }

    double tolerance = 1e-9;
    for (int i = 0; i < node_count; i++)
    {
        int gx = (int)llround(all_nodes[i].x - min_x) + margin;
        int gy = (int)llround(all_nodes[i].y - min_y) + margin;
        if (gx < 0 || gx >= width || gy < 0 || gy >= height)
        {
            continue;
        }

        bool infected = point_in_list(all_nodes[i], partition->infected_points,
                                      partition->infected_count, tolerance);
        bool secure = is_secure_point(all_nodes[i], all_nodes, secure_indices,
                                      secure_count, tolerance);

        PlotCell *cell = &grid[gy * width + gx];
        if (infected)
        {
            cell->type = CELL_INFECTED;
        }
        else if (secure)
        {
            cell->type = CELL_SECURE;
        }
        else
        {
            cell->type = CELL_HEALTHY;
        }
    }

    printf("\nASCII plot (integer grid, origin bottom-left):\n\n");

    for (int y = height - 1; y >= 0; y--)
    {
        if (y >= margin && y < height - margin)
        {
            double actual_y = min_y + (y - margin);
            printf("%6.1f | ", actual_y);
        }
        else
        {
            printf("       | ");
        }

        for (int x = 0; x < width; x++)
        {
            PlotCell cell = grid[y * width + x];
            const char *color = COLOR_RESET;
            char symbol = ' ';
            switch (cell.type)
            {
            case CELL_HULL:
                color = COLOR_BLUE;
                symbol = '.';
                break;
            case CELL_SECURE:
                color = COLOR_GREEN;
                symbol = 'S';
                break;
            case CELL_HEALTHY:
                color = COLOR_CYAN;
                symbol = 'H';
                break;
            case CELL_INFECTED:
                color = COLOR_RED;
                symbol = 'X';
                break;
            default:
                symbol = ' ';
                break;
            }
            printf("%s%c%s", color, symbol, COLOR_RESET);
        }
        printf("\n");
    }

    printf("       + ");
    for (int x = 0; x < width; x++)
    {
        printf("-");
    }
    printf("\n");

    printf("         X range: %.1f to %.1f\n", min_x, max_x);
    printf("         Y range: %.1f to %.1f\n", min_y, max_y);

    printf("\nLegend:\n");
    printf("  %sS%s Secure node\n", COLOR_GREEN, COLOR_RESET);
    printf("  %sH%s Healthy node\n", COLOR_CYAN, COLOR_RESET);
    printf("  %sX%s Infected node\n", COLOR_RED, COLOR_RESET);
    printf("  %s.%s Convex hull edge\n", COLOR_BLUE, COLOR_RESET);

    free(grid);
}
/**
 * Removes duplicate points from an array
 * Returns the new count of unique points
 */
int remove_duplicates(Point *points, int count, Point *unique_points)
{
    if (count == 0)
        return 0;

    int unique_count = 0;
    double tolerance = 1e-12;

    for (int i = 0; i < count; i++)
    {
        bool is_duplicate = false;
        for (int j = 0; j < unique_count; j++)
        {
            if (points_equal(points[i], unique_points[j], tolerance))
            {
                is_duplicate = true;
                break;
            }
        }
        if (!is_duplicate)
        {
            unique_points[unique_count++] = points[i];
        }
    }

    return unique_count;
}

/**
 * Computes the cross product of vectors OA and OB where O = p0, A = p1, B = p2
 * Returns > 0 for counter-clockwise turn, < 0 for clockwise turn, 0 for collinear
 */
double cross_product(Point p0, Point p1, Point p2)
{
    return (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
}

/**
 * Comparator function for qsort to sort points by polar angle
 * (used in Graham's scan algorithm)
 */
Point pivot_point;

int compare_polar_angle(const void *a, const void *b)
{
    Point *p1 = (Point *)a;
    Point *p2 = (Point *)b;

    double cross = cross_product(pivot_point, *p1, *p2);

    if (fabs(cross) < 1e-12)
    {
        // Collinear points - sort by distance from pivot
        double dist1 = (p1->x - pivot_point.x) * (p1->x - pivot_point.x) +
                       (p1->y - pivot_point.y) * (p1->y - pivot_point.y);
        double dist2 = (p2->x - pivot_point.x) * (p2->x - pivot_point.x) +
                       (p2->y - pivot_point.y) * (p2->y - pivot_point.y);
        return (dist1 < dist2) ? -1 : (dist1 > dist2) ? 1
                                                      : 0;
    }

    return (cross > 0) ? -1 : 1;
}

/**
 * Computes the convex hull using Graham's scan algorithm
 */
ConvexHull compute_convex_hull(Point *points, int count)
{
    ConvexHull hull;
    hull.vertices = NULL;
    hull.vertex_count = 0;
    hull.equations = NULL;
    hull.equation_count = 0;
    hull.is_valid = false;

    if (count < 3)
    {
        return hull;
    }

    // Find the point with the lowest y-coordinate (and leftmost if tied)
    int min_idx = 0;
    for (int i = 1; i < count; i++)
    {
        if (points[i].y < points[min_idx].y ||
            (fabs(points[i].y - points[min_idx].y) < 1e-12 && points[i].x < points[min_idx].x))
        {
            min_idx = i;
        }
    }

    // Swap the pivot point to the first position
    Point temp = points[0];
    points[0] = points[min_idx];
    points[min_idx] = temp;

    // Set the pivot for sorting
    pivot_point = points[0];

    // Sort the remaining points by polar angle
    qsort(points + 1, count - 1, sizeof(Point), compare_polar_angle);

    // Check if all points are collinear
    bool all_collinear = true;
    for (int i = 2; i < count; i++)
    {
        if (fabs(cross_product(points[0], points[1], points[i])) > 1e-12)
        {
            all_collinear = false;
            break;
        }
    }

    if (all_collinear)
    {
        return hull;
    }

    // Perform Graham's scan
    Point *stack = (Point *)malloc(count * sizeof(Point));
    int stack_size = 0;

    stack[stack_size++] = points[0];
    stack[stack_size++] = points[1];

    for (int i = 2; i < count; i++)
    {
        while (stack_size > 1 &&
               cross_product(stack[stack_size - 2], stack[stack_size - 1], points[i]) <= 0)
        {
            stack_size--;
        }
        stack[stack_size++] = points[i];
    }

    // Store hull vertices
    hull.vertices = (Point *)malloc(stack_size * sizeof(Point));
    hull.vertex_count = stack_size;
    for (int i = 0; i < stack_size; i++)
    {
        hull.vertices[i] = stack[i];
    }

    // Compute equations for each edge of the hull
    hull.equations = (LineEquation *)malloc(stack_size * sizeof(LineEquation));
    hull.equation_count = stack_size;

    for (int i = 0; i < stack_size; i++)
    {
        Point p1 = hull.vertices[i];
        Point p2 = hull.vertices[(i + 1) % stack_size];

        // For edge from p1 to p2, the outward normal is (-(p2.y - p1.y), p2.x - p1.x)
        // Equation: A*x + B*y + C = 0
        // We want points on the left (inside) to satisfy A*x + B*y + C <= 0
        hull.equations[i].A = -(p2.y - p1.y);
        hull.equations[i].B = p2.x - p1.x;
        hull.equations[i].C = -hull.equations[i].A * p1.x - hull.equations[i].B * p1.y;

        // Normalize so that the equation points outward
        // Check with the centroid of the hull
        double cx = 0, cy = 0;
        for (int j = 0; j < stack_size; j++)
        {
            cx += hull.vertices[j].x;
            cy += hull.vertices[j].y;
        }
        cx /= stack_size;
        cy /= stack_size;

        double test = hull.equations[i].A * cx + hull.equations[i].B * cy + hull.equations[i].C;
        if (test > 0)
        {
            // Flip the equation so it points outward
            hull.equations[i].A = -hull.equations[i].A;
            hull.equations[i].B = -hull.equations[i].B;
            hull.equations[i].C = -hull.equations[i].C;
        }
    }

    hull.is_valid = true;

    free(stack);
    return hull;
}

/**
 * Tests if a point is inside or on the boundary of the convex hull
 */
bool is_point_inside_hull(Point p, ConvexHull *hull, double tolerance)
{
    if (!hull->is_valid)
    {
        return false;
    }

    for (int i = 0; i < hull->equation_count; i++)
    {
        double value = hull->equations[i].A * p.x +
                       hull->equations[i].B * p.y +
                       hull->equations[i].C;
        if (value > tolerance)
        {
            return false;
        }
    }

    return true;
}

/**
 * Frees memory allocated for a ConvexHull structure
 */
void free_convex_hull(ConvexHull *hull)
{
    if (hull->vertices)
    {
        free(hull->vertices);
        hull->vertices = NULL;
    }
    if (hull->equations)
    {
        free(hull->equations);
        hull->equations = NULL;
    }
}

/**
 * Frees memory allocated for a PartitionResult structure
 */
void free_partition_result(PartitionResult *result)
{
    if (result->healthy_points)
    {
        free(result->healthy_points);
        result->healthy_points = NULL;
    }
    if (result->infected_points)
    {
        free(result->infected_points);
        result->infected_points = NULL;
    }
}

/**
 * Partitions a list of 2D nodes into "healthy" and "infected" lists
 * based on a convex hull computed from a subset of secure nodes.
 *
 * The "secure region" is the convex hull of the secure nodes.
 *
 * "Healthy" points are inside or on the boundary of this hull.
 * "Infected" points are outside the hull.
 *
 * Args:
 *     all_nodes: Array of Point structures representing all points.
 *     node_count: Number of nodes in all_nodes.
 *     secure_indices: Array of indices corresponding to the secure nodes in all_nodes.
 *     secure_count: Number of indices in secure_indices.
 *     out_hull: Optional pointer to store the computed convex hull (may be NULL).
 *
 * Returns:
 *     PartitionResult structure containing healthy_points and infected_points arrays
 */
PartitionResult partition_data(Point *all_nodes, int node_count,
                               int *secure_indices, int secure_count,
                               ConvexHull *out_hull)
{
    PartitionResult result;
    result.healthy_points = NULL;
    result.healthy_count = 0;
    result.infected_points = NULL;
    result.infected_count = 0;

    if (out_hull)
    {
        out_hull->vertices = NULL;
        out_hull->vertex_count = 0;
        out_hull->equations = NULL;
        out_hull->equation_count = 0;
        out_hull->is_valid = false;
    }

    // --- 1. Handle Edge Cases ---

    // Case 1: No secure indices provided
    if (secure_count == 0 || secure_indices == NULL)
    {
        fprintf(stderr, "Warning: No secure indices provided. All nodes classified as infected.\n");
        result.infected_points = (Point *)malloc(node_count * sizeof(Point));
        result.infected_count = node_count;
        for (int i = 0; i < node_count; i++)
        {
            result.infected_points[i] = all_nodes[i];
        }
        if (out_hull)
        {
            out_hull->is_valid = false;
        }
        return result;
    }

    // Extract secure nodes
    Point *secure_nodes = (Point *)malloc(secure_count * sizeof(Point));
    for (int i = 0; i < secure_count; i++)
    {
        secure_nodes[i] = all_nodes[secure_indices[i]];
    }

    // Use only unique secure nodes to build the hull
    Point *secure_nodes_unique = (Point *)malloc(secure_count * sizeof(Point));
    int unique_count = remove_duplicates(secure_nodes, secure_count, secure_nodes_unique);

    // Case 2: Not enough unique points to form a 2D area (hull)
    if (unique_count < 3)
    {
        fprintf(stderr, "Warning: Only %d unique secure nodes. "
                        "Cannot form a 2D secure region. Classifying all nodes as infected.\n",
                unique_count);
        result.infected_points = (Point *)malloc(node_count * sizeof(Point));
        result.infected_count = node_count;
        for (int i = 0; i < node_count; i++)
        {
            result.infected_points[i] = all_nodes[i];
        }
        if (out_hull)
        {
            out_hull->is_valid = false;
        }
        free(secure_nodes);
        free(secure_nodes_unique);
        return result;
    }

    // --- 2. Compute the Convex Hull of Secure Nodes ---

    // Create a working copy of unique secure nodes for hull computation
    Point *hull_points = (Point *)malloc(unique_count * sizeof(Point));
    for (int i = 0; i < unique_count; i++)
    {
        hull_points[i] = secure_nodes_unique[i];
    }

    ConvexHull hull = compute_convex_hull(hull_points, unique_count);

    // Case 3: All secure nodes (>=3) are collinear
    if (!hull.is_valid)
    {
        fprintf(stderr, "Warning: All secure nodes are collinear. "
                        "Cannot form a 2D secure region. Classifying all nodes as infected.\n");
        result.infected_points = (Point *)malloc(node_count * sizeof(Point));
        result.infected_count = node_count;
        for (int i = 0; i < node_count; i++)
        {
            result.infected_points[i] = all_nodes[i];
        }
        if (out_hull)
        {
            out_hull->is_valid = false;
        }
        free(secure_nodes);
        free(secure_nodes_unique);
        free(hull_points);
        free_convex_hull(&hull);
        return result;
    }

    // --- 3. Classify All Points ---

    // Add a small tolerance for floating point errors
    // This ensures points *on* the boundary are included.
    double tolerance = 1e-12;

    // Allocate temporary arrays for classification
    Point *temp_healthy = (Point *)malloc(node_count * sizeof(Point));
    Point *temp_infected = (Point *)malloc(node_count * sizeof(Point));
    int healthy_count = 0;
    int infected_count = 0;

    // Test each point
    for (int i = 0; i < node_count; i++)
    {
        if (is_point_inside_hull(all_nodes[i], &hull, tolerance))
        {
            temp_healthy[healthy_count++] = all_nodes[i];
        }
        else
        {
            temp_infected[infected_count++] = all_nodes[i];
        }
    }

    // --- 4. Return the Partitioned Lists ---

    // Allocate exact-sized arrays for results
    if (healthy_count > 0)
    {
        result.healthy_points = (Point *)malloc(healthy_count * sizeof(Point));
        for (int i = 0; i < healthy_count; i++)
        {
            result.healthy_points[i] = temp_healthy[i];
        }
    }
    result.healthy_count = healthy_count;

    if (infected_count > 0)
    {
        result.infected_points = (Point *)malloc(infected_count * sizeof(Point));
        for (int i = 0; i < infected_count; i++)
        {
            result.infected_points[i] = temp_infected[i];
        }
    }
    result.infected_count = infected_count;

    // Cleanup
    free(secure_nodes);
    free(secure_nodes_unique);
    free(hull_points);
    free(temp_healthy);
    free(temp_infected);

    if (out_hull)
    {
        *out_hull = hull;
    }
    else
    {
        free_convex_hull(&hull);
    }

    return result;
}

// --- --- --- --- --- --- --- --- ---
// --- Example Usage ---
// --- --- --- --- --- --- --- --- ---
int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    Point all_nodes[] = {
        {1, 1}, {1, 6}, {3, 4}, {4, 2}, {4, 7}, {6, 5}, {7, 1}, {8, 8}, {10, 3}, {11, 6}, {2, 8}, {5, 10}, {9, 9}};
    int node_count = sizeof(all_nodes) / sizeof(all_nodes[0]);

    int secure_indices[] = {0, 1, 4, 6, 7, 11};
    int secure_count = sizeof(secure_indices) / sizeof(secure_indices[0]);

    printf("--- Running Node Partitioner ---\n");

    ConvexHull secure_hull = {0};
    PartitionResult result = partition_data(all_nodes, node_count,
                                            secure_indices, secure_count,
                                            &secure_hull);

    printf("\nTotal nodes: %d\n", node_count);

    printf("\n%sSecure nodes (%d):%s\n", COLOR_GREEN, secure_count, COLOR_RESET);
    for (int i = 0; i < secure_count; i++)
    {
        Point p = all_nodes[secure_indices[i]];
        printf("  %sS%s (%.1f, %.1f)\n", COLOR_GREEN, COLOR_RESET, p.x, p.y);
    }

    printf("\nHealthy nodes (%d):\n", result.healthy_count);
    double tolerance = 1e-9;
    for (int i = 0; i < result.healthy_count; i++)
    {
        Point p = result.healthy_points[i];
        bool secure_member = is_secure_point(p, all_nodes, secure_indices, secure_count, tolerance);
        const char *color = secure_member ? COLOR_GREEN : COLOR_CYAN;
        const char *label = secure_member ? "S" : "H";
        printf("  %s%s%s (%.1f, %.1f)\n", color, label, COLOR_RESET, p.x, p.y);
    }

    printf("\n%sInfected nodes (%d):%s\n", COLOR_RED, result.infected_count, COLOR_RESET);
    for (int i = 0; i < result.infected_count; i++)
    {
        Point p = result.infected_points[i];
        printf("  %sX%s (%.1f, %.1f)\n", COLOR_RED, COLOR_RESET, p.x, p.y);
    }

    render_ascii_plot(all_nodes, node_count,
                      secure_indices, secure_count,
                      &result, &secure_hull);

    free_convex_hull(&secure_hull);
    free_partition_result(&result);

    return 0;
}
