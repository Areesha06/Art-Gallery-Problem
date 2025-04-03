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


#NEW DRAWING POLYGON###
mport matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math

# Global variables
p_vertices = []
graph = nx.Graph()

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
            return

    # Check if the new point is close to an existing point (to merge nodes)
    for point in p_vertices:
        # if np.linalg.norm(np.array(existing_point) - np.array(new_point)) < 0.2:
        dist=math.sqrt((new_point[0]-point[0])**2+(new_point[1]-point[1])**2)
        # if np.linalg.norm(np.array(first_point) - np.array(new_point)) < 0.2:
        if dist<0.3:
            new_point = point  # Snap to existing node
            break

    p_vertices.append(new_point)

    # Redraw dynamically
    ax.cla()  # Clears only the contents, not the entire axis
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)
    ax.set_title("Click to draw straight-line polygon")

    # ax.scatter(x_vals, y_vals, 'm')  # Mark points
    x_vals=[]
    y_vals=[]

    for i in p_vertices:
        x_vals.append(i[0])
        y_vals.append(i[1])

    if len(p_vertices)>1:
        ax.plot(x_vals, y_vals, 'k-')  # Blue straight lines
       
    ax.scatter(x_vals, y_vals, c='m', s=50)  # Mark points

    # fig.canvas.draw_idle()  # Forces a refresh
    # plt.pause(0.01)  # Helps Matplotlib update properly
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
    # # Add nodes
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

