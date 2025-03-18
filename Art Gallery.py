###DRAWING POLYGON-USER INPUT###

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# Global variables
polygon_vertices = []
graph = nx.Graph()

def on_click(event):
    """Capture user clicks and draw straight lines to form the polygon."""
    if event.xdata is None or event.ydata is None:
        return  # Ignore clicks outside the plot

    new_point = (round(event.xdata, 2), round(event.ydata, 2))  # Round for consistency

    # Close the polygon if clicking near the first point
    if len(polygon_vertices) > 2:
        first_point = polygon_vertices[0]
        if np.linalg.norm(np.array(first_point) - np.array(new_point)) < 0.2:
            polygon_vertices.append(first_point)  # Snap to start
            draw_polygon()
            return

    # Check if the new point is close to an existing point (to merge nodes)
    for existing_point in polygon_vertices:
        if np.linalg.norm(np.array(existing_point) - np.array(new_point)) < 0.2:
            new_point = existing_point  # Snap to existing node
            break

    polygon_vertices.append(new_point)

    # Draw the lines dynamically
    ax.clear()
    # ax.set_xlim(0, 10)
    # ax.set_ylim(0, 10)
    ax.set_title("Click to draw straight-line polygon, click near start to close")

    if len(polygon_vertices) > 1:
        x_vals, y_vals = zip(*polygon_vertices)
        ax.plot(x_vals, y_vals, 'b-')  # Blue straight lines

    ax.scatter(*zip(*polygon_vertices), color='red')  # Mark points
    plt.draw()

def draw_polygon():
    """Convert the drawn polygon into a NetworkX graph and visualize it."""
    graph.clear()
    
    # Use dictionary to remove duplicates while maintaining order
    unique_nodes = list(dict.fromkeys(polygon_vertices))

    # Add nodes
    for i, (x, y) in enumerate(unique_nodes):
        graph.add_node(i, pos=(x, y))

    # Add edges
    for i in range(len(unique_nodes) - 1):
        graph.add_edge(i, i + 1)
    graph.add_edge(len(unique_nodes) - 1, 0)  # Close the polygon

    pos = nx.get_node_attributes(graph, 'pos')

    plt.clf()
    nx.draw(graph, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500)
    plt.title("Final Drawn Polygon")
    plt.show()

# Create Matplotlib figure
fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_title("Click to draw straight-line polygon")

# Connect click event
fig.canvas.mpl_connect("button_press_event", on_click)

plt.show()

