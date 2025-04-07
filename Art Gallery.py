

# #NEW DRAWING POLYGON###

# import matplotlib.pyplot as plt
# import networkx as nx
# import numpy as np
# import math

# # Global variables
# p_vertices = []
# graph = nx.Graph()

# def on_click(event):
#     """Capture user clicks and draw straight lines to form the polygon."""
#     if event.xdata is None or event.ydata is None:
#         return  # Ignore clicks outside the plot

#     new_point = (float(round(event.xdata, 2)), float(round(event.ydata, 2)))  # Round for consistency
#     print(f"Clicked: {new_point}")  # Debugging print

#     # Close the polygon if clicking near the first point
#     if len(p_vertices) > 2:
#         first_point = p_vertices[0]
#         dist=math.sqrt((new_point[0]-first_point[0])**2+(new_point[1]-first_point[1])**2)
#         # if np.linalg.norm(np.array(first_point) - np.array(new_point)) < 0.2:
#         if dist<0.3:
#             p_vertices.append(first_point)  # Snap to start
#             draw_polygon()
#             return

#     # Check if the new point is close to an existing point (to merge nodes)
#     for point in p_vertices:
#         # if np.linalg.norm(np.array(existing_point) - np.array(new_point)) < 0.2:
#         dist=math.sqrt((new_point[0]-point[0])**2+(new_point[1]-point[1])**2)
#         # if np.linalg.norm(np.array(first_point) - np.array(new_point)) < 0.2:
#         if dist<0.3:
#             new_point = point  # Snap to existing node
#             break

#     p_vertices.append(new_point)

#     # Redraw dynamically
#     ax.cla()  # Clears only the contents, not the entire axis
#     ax.set_xlim(0, 15)
#     ax.set_ylim(0, 15)
#     ax.set_title("Click to draw straight-line polygon")

#     # ax.scatter(x_vals, y_vals, 'm')  # Mark points
#     x_vals=[]
#     y_vals=[]

#     for i in p_vertices:
#         x_vals.append(i[0])
#         y_vals.append(i[1])

#     if len(p_vertices)>1:
#         ax.plot(x_vals, y_vals, 'k-')  # Blue straight lines
       
#     ax.scatter(x_vals, y_vals, c='m', s=50)  # Mark points

#     # fig.canvas.draw_idle()  # Forces a refresh
#     # plt.pause(0.01)  # Helps Matplotlib update properly
#     fig.canvas.draw()
#     fig.canvas.flush_events()

# def draw_polygon():
#     """Convert the drawn polygon into a NetworkX graph and visualize it."""
#     graph.clear()
    
#     # Remove duplicates while keeping order
#     unique_nodes = []
#     for i in p_vertices:
#         if i not in unique_nodes:
#             unique_nodes.append(i)
#     # # Add nodes
#     for i in range(len(unique_nodes)):
#         x, y = unique_nodes[i]  # Get coordinates at index i
#         graph.add_node(i, pos=(x, y))  

#     # Add edges
#     for i in range(len(unique_nodes) - 1):
#         graph.add_edge(i, i + 1)
#     graph.add_edge(len(unique_nodes) - 1, 0)  # Close the polygon

#     pos = nx.get_node_attributes(graph, 'pos')

#     plt.clf()
#     nx.draw(graph, pos, with_labels=True, node_color='pink', edge_color='gray', node_size=500)
#     plt.title("Final Drawn Polygon")
#     plt.show()

# # Create Matplotlib figure
# fig, ax = plt.subplots()
# ax.set_xlim(0, 15)
# ax.set_ylim(0, 15)
# ax.set_title("Click to draw straight-line polygon")

# # Connect click event
# fig.canvas.mpl_connect("button_press_event", on_click)

# plt.show()

######################################################################################################
def main_polygon_generation():
    import matplotlib.pyplot as plt
    import networkx as nx
    import numpy as np
    import math
    from itertools import combinations

    # Global variables
    p_vertices = []
    graph = nx.Graph()

    def do_intersect(p1, p2, p3, p4):
        """Helper function to check if line segments p1p2 and p3p4 intersect."""
        def ccw(a, b, c):
            """Check if the points a, b, c are in counter-clockwise order."""
            return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

        # Two line segments (p1p2) and (p3p4) intersect if and only if the segments are "straddling" each other
        return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)

    def check_self_intersection(points):
        """Check if any pair of non-adjacent edges of the polygon intersect."""
        # Generate the list of edges (pairs of consecutive points)
        edges = [(points[i], points[i + 1]) for i in range(len(points) - 1)]
        edges.append((points[-1], points[0]))  # Close the polygon by adding the last edge

        # Check all pairs of non-adjacent edges
        for (p1, p2), (p3, p4) in combinations(edges, 2):
            # Skip adjacent edges (including the edge formed between first and last point)
            if p2 == p3 or p1 == p4:
                continue
            if do_intersect(p1, p2, p3, p4):
                return True  # Intersection found, meaning it's a self-intersecting polygon
        return False  # No intersections

    def on_click(event):
        """Capture user clicks and draw straight lines to form the polygon."""
        if event.xdata is None or event.ydata is None:
            return  # Ignore clicks outside the plot

        new_point = (float(round(event.xdata, 2)), float(round(event.ydata, 2)))  # Round for consistency
        print(f"Clicked: {new_point}")  # Debugging print


        # Close the polygon if clicking near the first point
        if len(p_vertices) > 2:
            first_point = p_vertices[0]
            dist=math.sqrt((new_point[0]-first_point[0])**2+(new_point[1]-first_point[1])**2)
            # if np.linalg.norm(np.array(first_point) - np.array(new_point)) < 0.2:
            if dist<0.3:
                p_vertices.append(first_point)  # Snap to start
                draw_polygon()
                return p_vertices

        # Check if the new point is close to an existing point (to merge nodes)
        for point in p_vertices:
            # if np.linalg.norm(np.array(existing_point) - np.array(new_point)) < 0.2:
            dist=math.sqrt((new_point[0]-point[0])**2+(new_point[1]-point[1])**2)
            # if np.linalg.norm(np.array(first_point) - np.array(new_point)) < 0.2:
            if dist<0.3:
                new_point = point  # Snap to existing node
                break

        p_vertices.append(new_point)

        # Check if the polygon remains valid (no self-intersections)
        if len(p_vertices) > 2 and check_self_intersection(p_vertices):
            print("Invalid polygon: self-intersection detected.")
            p_vertices.pop()  # Remove the last point to reject the click
            return

        # Redraw dynamically
        ax.cla()  # Clears only the contents, not the entire axis
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15)
        ax.set_title("Click to draw straight-line polygon")

        x_vals = [pt[0] for pt in p_vertices]
        y_vals = [pt[1] for pt in p_vertices]

        if len(p_vertices) > 1:
            ax.plot(x_vals, y_vals, 'k-')  # Blue straight lines

        ax.scatter(x_vals, y_vals, c='m', s=50)  # Mark points

        fig.canvas.draw()
        fig.canvas.flush_events()

    def draw_polygon():
        """Convert the drawn polygon into a NetworkX graph and visualize it."""
        graph.clear()

        # Remove duplicates while keeping order
        unique_nodes = []
        for i in p_vertices:
            if i not in unique_nodes:
                unique_nodes.append(i)

        # Add nodes
        for i in range(len(unique_nodes)):
            x, y = unique_nodes[i]  # Get coordinates at index i
            graph.add_node(i, pos=(x, y))

        # Add edges
        for i in range(len(unique_nodes) - 1):
            graph.add_edge(i, i + 1)
        graph.add_edge(len(unique_nodes) - 1, 0)  # Close the polygon

        pos = nx.get_node_attributes(graph, 'pos')

        plt.clf()
        nx.draw(graph, pos, with_labels=True, node_color='pink', edge_color='gray', node_size=500)
        plt.title("Final Drawn Polygon")
        plt.show()

    # Create Matplotlib figure
    fig, ax = plt.subplots()
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)
    ax.set_title("Click to draw straight-line polygon")

    # Connect click event
    fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()
    return p_vertices



#################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


# Function to check if a point is inside a triangle
def is_point_in_triangle(pt, v1, v2, v3):
    # Barycentric coordinates method
    def sign(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    d1 = sign(pt, v1, v2)
    d2 = sign(pt, v2, v3)
    d3 = sign(pt, v3, v1)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_neg and has_pos)


# Function to check if a triangle is convex
def is_convex(a, b, c):
    return (c[0] - b[0]) * (a[1] - b[1]) - (c[1] - b[1]) * (a[0] - b[0]) > 0


# Ear clipping triangulation function
def ear_clip_triangulate(polygon):
    vertices = polygon.copy()
    triangles = []

    while len(vertices) > 3:
        for i in range(len(vertices)):
            prev_idx = (i - 1) % len(vertices)
            next_idx = (i + 1) % len(vertices)

            prev_vertex = vertices[prev_idx]
            current_vertex = vertices[i]
            next_vertex = vertices[next_idx]

            if is_convex(prev_vertex, current_vertex, next_vertex):
                ear = True
                for j in range(len(vertices)):
                    if j != prev_idx and j != i and j != next_idx:
                        if is_point_in_triangle(vertices[j], prev_vertex, current_vertex, next_vertex):
                            ear = False
                            break

                if ear:
                    triangles.append([prev_vertex, current_vertex, next_vertex])
                    vertices.pop(i)
                    break

    # The remaining polygon should be a triangle
    triangles.append(vertices)

    return triangles


# Plotting function
def plot_polygon_and_triangulation(polygon, triangles):
    fig, ax = plt.subplots()

    # Plot the polygon
    polygon_patch = Polygon(polygon, closed=True, fill=None, edgecolor='black')
    ax.add_patch(polygon_patch)

    # Plot the triangles
    for triangle in triangles:
        triangle_patch = Polygon(triangle, closed=True, fill=True, edgecolor='red', facecolor='lightblue', alpha=0.5)
        ax.add_patch(triangle_patch)

    ax.set_xlim(min([x[0] for x in polygon]) - 1, max([x[0] for x in polygon]) + 1)
    ax.set_ylim(min([x[1] for x in polygon]) - 1, max([x[1] for x in polygon]) + 1)
    ax.set_aspect('equal', 'box')
    plt.show()

# Example polygon (a simple concave polygon)
# polygon = [(1, 1), (4, 1), (4, 4), (2, 3), (1, 4)]

# # Apply ear clipping triangulation
# triangles = ear_clip_triangulate(polygon)

# # Plot the result
# plot_polygon_and_triangulation(polygon, triangles)

vertices = main_polygon_generation()
print(vertices)
# triangles = ear_clip_triangulate(vertices)

# # Plot the result
# plot_polygon_and_triangulation(vertices, triangles)

