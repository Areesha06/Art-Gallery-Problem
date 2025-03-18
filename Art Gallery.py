# import networkx as nx
# import matplotlib.pyplot as plt
# import itertools 
# import random

# # Function to create a simple polygon graph
# def create_polygon_graph(vertices):
#     G = nx.Graph()
#     n = len(vertices)
    
#     for i in range(n):
#         G.add_edge(i, (i + 1) % n)  # Connecting each vertex to the next (cyclic)
    
#     return G

# # Function to perform graph-based triangulation
# def triangulate_polygon(G):
#     n = len(G.nodes)
#     triangulated_graph = G.copy()
    
#     # Simple ear clipping method for triangulation
#     for i in range(n - 2):
#         triangulated_graph.add_edge(i, i + 2)
    
#     return triangulated_graph

# # Function to perform 3-coloring using a greedy algorithm
# def three_coloring(G):
#     colors = {}
#     available_colors = [0, 1, 2]
    
#     for node in G.nodes:
#         neighbor_colors = {colors[neighbor] for neighbor in G.neighbors(node) if neighbor in colors}
        
#         for color in available_colors:
#             if color not in neighbor_colors:
#                 colors[node] = color
#                 break
    
#     return colors

# # Function to find guard placement based on Chvátal's Theorem
# def place_guards(G, colors):
#     guard_color = min(set(colors.values()), key=lambda c: list(colors.values()).count(c))
#     guards = [node for node in colors if colors[node] == guard_color]
#     return guards

# # Example: Define a polygon (hexagon for simplicity)
# n = 9  # Number of vertices
# vertices = [(random.randint(0, 10), random.randint(0, 10)) for _ in range(n)]
# polygon_graph = create_polygon_graph(vertices)

# # Triangulate the polygon
# t_graph = triangulate_polygon(polygon_graph)

# # Perform 3-coloring
# colors = three_coloring(t_graph)

# # Determine guards' placement
# guards = place_guards(t_graph, colors)

# # Visualize the polygon and guards
# plt.figure(figsize=(8, 6))
# pos = {i: vertices[i] for i in range(n)}
# nx.draw(t_graph, pos, with_labels=True, node_color=[colors[i] for i in t_graph.nodes], cmap=plt.cm.Set1, edge_color='black')
# nx.draw_networkx_nodes(t_graph, pos, nodelist=guards, node_color='red', label='Guards', node_size=300)
# plt.legend()
# plt.title("Chvátal's Art Gallery Theorem Simulation")
# plt.show()

# print(f"Number of guards needed: {len(guards)}")
# print(f"Guards placed at vertices: {guards}")


# import matplotlib.pyplot as plt
# import networkx as nx
# import numpy as np
# from scipy.spatial import Delaunay

# graph = nx.Graph()
# points = []

# def onclick(event):
#     if event.xdata is None or event.ydata is None:
#         return
#     x, y = event.xdata, event.ydata
#     points.append((x, y))
#     graph.add_node(len(points) - 1, pos=(x, y))
    
#     plt.clf()
#     draw_graph()

# def draw_graph():
#     pos = nx.get_node_attributes(graph, 'pos')
#     nx.draw(graph, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500)
#     plt.scatter(*zip(*points), color='red')
#     plt.draw()

# def triangulate():
#     if len(points) < 3:
#         print("Need at least 3 points for triangulation.")
#         return
    
#     tri = Delaunay(points)
#     for simplex in tri.simplices:
#         for i in range(3):
#             graph.add_edge(simplex[i], simplex[(i + 1) % 3])
    
#     plt.clf()
#     draw_graph()
#     color_triangulation(tri)

# def color_triangulation(tri):
#     colors = ["red", "blue", "green"]
#     color_map = {}
    
#     for simplex in tri.simplices:
#         used_colors = set()
#         for i in simplex:
#             if i in color_map:
#                 used_colors.add(color_map[i])
        
#         available_colors = [c for c in colors if c not in used_colors]
#         for i in simplex:
#             if i not in color_map:
#                 color_map[i] = available_colors.pop()
    
#     pos = nx.get_node_attributes(graph, 'pos')
#     node_colors = [color_map[i] for i in graph.nodes()]
#     nx.draw(graph, pos, with_labels=True, node_color=node_colors, edge_color='gray', node_size=500)
#     plt.draw()

# def main():
#     fig, ax = plt.subplots()
#     ax.set_xlim(0, 10)
#     ax.set_ylim(0, 10)
#     ax.set_title("Click to create nodes")
#     fig.canvas.mpl_connect('button_press_event', onclick)
#     plt.show()
    
#     triangulate()
#     plt.show()

# if __name__ == "__main__":
#     main()


import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.patches as patches

graph = nx.Graph()
guards = []
colors = ['blue', 'orange', 'purple', 'cyan', 'magenta']  # Different colors for different guards

def onclick(event):
    if event.xdata is None or event.ydata is None:
        return
    x, y = event.xdata, event.ydata
    
    guards.append((x, y))
    graph.add_node(len(guards) - 1, pos=(x, y))
    
    plt.clf()
    draw_graph()

def draw_graph():
    pos = nx.get_node_attributes(graph, 'pos')
    nx.draw(graph, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500)
    plt.scatter(*zip(*guards), color='red')
    
    for i, guard in enumerate(guards):
        plt.scatter(*guard, color='green', marker='X', s=200, label=f'Guard {i+1}')
        draw_guard_view(guard, colors[i % len(colors)])
    
    plt.legend()
    plt.draw()

def draw_guard_view(guard, color):
    guard_x, guard_y = guard
    radius = 5  # Set guard view radius
    angle = np.linspace(-np.pi / 2, np.pi / 2, 30)  # 180-degree view
    
    arc_x = guard_x + radius * np.cos(angle)
    arc_y = guard_y + radius * np.sin(angle)
    
    polygon = np.vstack(([guard_x, guard_y], np.column_stack((arc_x, arc_y))))
    plt.fill(polygon[:, 0], polygon[:, 1], color=color, alpha=0.3, label=f'View {color}')

def triangulate():
    if len(guards) < 3:
        print("Need at least 3 points for triangulation.")
        return
    
    tri = Delaunay(guards)
    for simplex in tri.simplices:
        for i in range(3):
            graph.add_edge(simplex[i], simplex[(i + 1) % 3])
    
    plt.clf()
    draw_graph()
    color_triangulation(tri)

def color_triangulation(tri):
    tri_colors = ["red", "blue", "green"]
    color_map = {}
    
    for simplex in tri.simplices:
        used_colors = set()
        for i in simplex:
            if i in color_map:
                used_colors.add(color_map[i])
        
        available_colors = [c for c in tri_colors if c not in used_colors]
        for i in simplex:
            if i not in color_map:
                color_map[i] = available_colors.pop()
    
    pos = nx.get_node_attributes(graph, 'pos')
    node_colors = [color_map[i] for i in graph.nodes()]
    nx.draw(graph, pos, with_labels=True, node_color=node_colors, edge_color='gray', node_size=500)
    plt.draw()

def main():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("Click to place guards")
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    triangulate()
    plt.show()

if __name__ == "__main__":
    main()
