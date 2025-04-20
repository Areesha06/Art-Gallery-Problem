import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
import numpy as np
import math
from itertools import combinations

# --- Polygon Generation Code ---
def main_polygon_generation():
    p_vertices = []
    graph = nx.Graph()

    def do_intersect(p1, p2, p3, p4):
        def ccw(a, b, c):
            return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

        return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)

    def check_self_intersection(points):
        if len(points) < 2:
            return False

        new_edge = (points[-2], points[-1])
        for i in range(len(points) - 3):
            edge = (points[i], points[i + 1])
            if do_intersect(edge[0], edge[1], new_edge[0], new_edge[1]):
                return True
        return False

    def on_click(event):
        if event.xdata is None or event.ydata is None:
            return

        new_point = (float(round(event.xdata, 2)), float(round(event.ydata, 2)))

        if len(p_vertices) > 2:
            first_point = p_vertices[0]
            dist = math.sqrt((new_point[0] - first_point[0])**2 + (new_point[1] - first_point[1])**2)
            if dist < 0.3:
                p_vertices.append(first_point)
                draw_polygon()
                return p_vertices

        for point in p_vertices:
            dist = math.sqrt((new_point[0] - point[0])**2 + (new_point[1] - point[1])**2)
            if dist < 0.3:
                new_point = point
                break

        p_vertices.append(new_point)

        if len(p_vertices) > 2 and check_self_intersection(p_vertices):
            p_vertices.pop()
            return

        ax.cla()
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15)
        ax.set_title("Click to draw straight-line polygon")

        x_vals = [pt[0] for pt in p_vertices]
        y_vals = [pt[1] for pt in p_vertices]

        if len(p_vertices) > 1:
            ax.plot(x_vals, y_vals, 'k-')

        ax.scatter(x_vals, y_vals, c='m', s=50)

        fig.canvas.draw()
        fig.canvas.flush_events()

    def draw_polygon():
        graph.clear()
        unique_nodes = []
        for i in p_vertices:
            if i not in unique_nodes:
                unique_nodes.append(i)

        for i in range(len(unique_nodes)):
            x, y = unique_nodes[i]
            graph.add_node(i, pos=(x, y))

        for i in range(len(unique_nodes) - 1):
            graph.add_edge(i, i + 1)
        graph.add_edge(len(unique_nodes) - 1, 0)

        pos = nx.get_node_attributes(graph, 'pos')

        plt.clf()
        nx.draw(graph, pos, with_labels=True, node_color='pink', edge_color='gray', node_size=500)
        plt.title("Final Drawn Polygon")
        plt.show()

    fig, ax = plt.subplots()
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)
    ax.set_title("Click to draw straight-line polygon")

    fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()
    return p_vertices[:-1]

# --- Triangulation and Coloring Code ---
def is_ccw(polygon):
    area = 0.0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        area += (x2 - x1) * (y2 + y1)
    return area < 0

def is_convex(a, b, c):
    return ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) > 0

def is_point_in_triangle(p, a, b, c):
    def sign(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    b1 = sign(p, a, b) < 0.0
    b2 = sign(p, b, c) < 0.0
    b3 = sign(p, c, a) < 0.0
    return b1 == b2 == b3

def ear_clipping_triangulation(polygon_coords):
    triangles = []
    vertices = polygon_coords[:]

    if not is_ccw(vertices):
        vertices.reverse()

    while len(vertices) > 3:
        ear_found = False
        for i in range(len(vertices)):
            prev = vertices[i - 1]
            curr = vertices[i]
            nxt = vertices[(i + 1) % len(vertices)]

            if is_convex(prev, curr, nxt):
                triangle = [prev, curr, nxt]
                contains_point = any(
                    p not in triangle and is_point_in_triangle(p, *triangle)
                    for p in vertices
                )
                if not contains_point:
                    triangles.append(triangle)
                    del vertices[i]
                    ear_found = True
                    break
        if not ear_found:
            raise ValueError("No ear found — the polygon might be invalid or not simple.")

    triangles.append(vertices)
    return triangles

def build_graph_from_triangles(triangles):
    G = nx.Graph()
    point_to_id = {}
    id_to_point = {}
    counter = 0

    for tri in triangles:
        for point in tri:
            if point not in point_to_id:
                point_to_id[point] = counter
                id_to_point[counter] = point
                counter += 1

    for pid, coord in id_to_point.items():
        G.add_node(pid, pos=coord)

    for tri in triangles:
        ids = [point_to_id[p] for p in tri]
        G.add_edge(ids[0], ids[1])
        G.add_edge(ids[1], ids[2])
        G.add_edge(ids[2], ids[0])

    return G, id_to_point

def chvatal_3_coloring(graph):
    coloring = {}
    for node in sorted(graph.nodes):
        used = {coloring.get(n) for n in graph.neighbors(node) if n in coloring}
        for color in range(3):
            if color not in used:
                coloring[node] = color
                break
    return coloring

def draw_triangulated_graph(G):
    pos = nx.get_node_attributes(G, 'pos')
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500)
    plt.title("Triangulated Polygon as Graph")
    plt.show()

def draw_colored_graph(G, coloring):
    color_map = {0: 'Red', 1: 'Green', 2: 'Blue'}

    # Fallback coloring for uncolored nodes
    for node in G.nodes:
        if node not in coloring:
            coloring[node] = 0

    node_colors = [mcolors.to_rgba(color_map[coloring[n]], alpha=0.6) for n in G.nodes()]
    pos = nx.get_node_attributes(G, 'pos')
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, node_color=node_colors, with_labels=True, edge_color='gray', node_size=500)
    plt.title("Coloring Based on Chvátal's Art Gallery Theorem")
    plt.show()

def color_triangulated_graph(triangles):
    G, _ = build_graph_from_triangles(triangles)
    draw_triangulated_graph(G)
    coloring = chvatal_3_coloring(G)
    draw_colored_graph(G, coloring)
    return G, coloring

def show_min_guard_positions(G, coloring):
    import matplotlib.patches as mpatches

    color_map = {0: 'Red', 1: 'Green', 2: 'Blue'}

    # Fallback coloring
    for node in G.nodes:
        if node not in coloring:
            coloring[node] = 0

    # Determine best colors
    color_counts = {c: list(coloring.values()).count(c) for c in range(3)}
    min_count = min(color_counts.values())
    best_colors = [c for c, count in color_counts.items() if count == min_count]

    best_color_names = [color_map[c] for c in best_colors]
    result_color_str = '/'.join(best_color_names)

    pos = nx.get_node_attributes(G, 'pos')

    # Assign colors for visualization
    node_colors = []
    for n in G.nodes():
        c = coloring[n]
        if c in best_colors:
            node_colors.append(mcolors.to_rgba(color_map[c], alpha=0.8))
        else:
            node_colors.append(mcolors.to_rgba('gray', alpha=0.4))

    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, node_color=node_colors, with_labels=True, edge_color='gray', node_size=500)

    # Text box info
    text = (
        f"Minimum guards needed by Chvátal's theorem: {math.ceil(len(coloring)/3)}\n"
        f"Actual minimum guards needed: {min_count}\n"
        f"Guards should be placed at color(s): {result_color_str}"
    )

    plt.gcf().text(0.05, 0.01, text, fontsize=10, va='bottom', ha='left', bbox=dict(facecolor='white', alpha=0.7))
    plt.title("Minimum Guard Positions Highlighted")
    plt.show()


if __name__ == "__main__":
    polygon_coords = main_polygon_generation()
    triangles = ear_clipping_triangulation(polygon_coords)
    G, coloring = color_triangulated_graph(triangles)
    show_min_guard_positions(G, coloring)


