import matplotlib.pyplot as plt     # module for drawing the graphs
import matplotlib.colors as mcolors # for using the colours in our graph 
import networkx as nx # the interface for plotting and visualizing the graph
import numpy as np 
import math #for operations like ciel and floor

def ccw(a, b, c): # checks if our 3 points are in counterclockwise direction. This will help us check convexity later and allow us to check if our lines intersect or not
    return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])


def intersect(p1, p2, p3, p4):  # uses the counterclockwaise function to determine if the lines intersect or not (formed by vectors) using cross product

    return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)

def edge_valid(new_edge_start, new_edge_end, existing_vertices): # allows us to see if we can add a new edge in our polygon witout intersecting the other edges
    for i in range(len(existing_vertices) - 1):
        p1, p2 = existing_vertices[i], existing_vertices[i + 1]
        if (p1 == new_edge_start or p2 == new_edge_start or p1 == new_edge_end or p2 == new_edge_end):
            continue
        if intersect(new_edge_start, new_edge_end, p1, p2):
            return False    # we cannot add that edge
    return True

def draw_polygon(graph, vertices): # this draws oir graph on networkX
    graph.clear() # makes sure that screen is complelty new
    nodes = [] # we remove any repeated vertices
    for i in vertices:
        if i not in nodes:
            nodes.append(i)

    for i in range(len(nodes)):
        x, y = nodes[i] # we get x and y values of our point 
        graph.add_node(i, pos=(x, y)) # adds our node to the graph at position (x, y)

    for i in range(len(nodes) - 1):  # we add edges in our  graph (-1 because our graph is closed and will always have 1 less edge than the number of vertices)
        graph.add_edge(i, i + 1)  # add edge between current and next point
    graph.add_edge(len(nodes) - 1, 0) # this adds the very last edge to close shape from the last vertex to very first vertex

    pos = nx.get_node_attributes(graph, 'pos') 
    # netwrokX syntax to plot graph
    plt.clf()
    nx.draw(graph, pos, with_labels=True, node_color='pink', edge_color='gray', node_size=500)
    plt.title("Your Art Gallery!")
    plt.show()

def main_polygon_generation(): # our main function which uses all the helper functions
    vertices = []
    graph = nx.Graph()  # initialize the vertex list and graph

    def on_click(event): # mouse clicking syntax, each time i click on the screen this fucntion runs

        if event.xdata is None or event.ydata is None: # if i have not cliked inside the graph screen, there is no data of the click
            return
        
        x = event.xdata
        y = event.ydata

        new = (float(round(x, 2)), float(round(y, 2))) # this is our current point (rounded to 2 decimal places to use less memory)


        if len(vertices) > 2:  # a polygon cannot be made with only 2 points
            # the very first point is required to check if we are closing our polygon or not
            first = vertices[0] 
            x0 = first[0]
            y0 = first[1]
            x1 = new[0]
            y1 = new[1]

            dist = math.sqrt((x1 - x0)**2 + (y1 - y0)**2) # distance formula for finding the distance between the our current point and first point 

            if dist < 0.3: # if we are very close to the first point, we are trying to close the polygon

                if edge_valid(vertices[-1], first, vertices) == False: # if we at a point where closing the polygon will cause edges to intersect then give error message

                    print("Kindly make a simple polygon! Art gallery theorem works on Simple Polygons.")
                    return
                vertices.append(first) # if there are no intersections, then append the first point again so we have a loop 
                draw_polygon(graph, vertices) # no we can simply draw our graph
                return vertices # we need this data for triangulation and coloring
    
        for point in vertices: # a simple check if we are trying to add a point near an already cliked point
            dist = math.sqrt((new[0] - point[0])**2 + (new[1] - point[1])**2)
            if dist < 0.3:
                new = point
                break

        if len(vertices) > 0: # this checks if our currently clicked point does not cause any overlapping edges

            if edge_valid(vertices[-1], new, vertices) == False:
                print("Kindly make a simple polygon! Art gallery theorem works on Simple Polygons.")
                return

        vertices.append(new) # finally add the new point clicked into the set of vertices

        ax.cla()
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15) # clears the screen and sets the axes on the screen
        ax.set_title("Generate your Art Gallery!")

        x_coords = []
        y_coords = []
        for node in vertices:
            x_coords.append(node[0])
            y_coords.append(node[1])


        if len(vertices) > 1: # there needs to be at least 2 points to draw an edge
            ax.plot(x_coords, y_coords, 'k-') # draws black lines

        ax.scatter(x_coords, y_coords, c='m', s=50) # draws magenta points with size 50

        fig.canvas.draw() # updates screen to draw graph
        fig.canvas.flush_events() # proper screen updating


    # draws our final graph on fig window of area ax
    fig, ax = plt.subplots()
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)
    ax.set_title("Make your own Art Gallery!")

    fig.canvas.mpl_connect("button_press_event", on_click) # run on_click when i click mouse on screen
    plt.show()
    return vertices[:-1] # we dont need last vertex beacuse its a cycle and its the same as first

def counterclockwhole(polygon): # this checks if our entire shape has points in counterclockwise direction
    area = 0.0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i] # x and y values of our current point
        x2, y2 = polygon[(i + 1) % n] # x and y values of next point. % allows us  to loop back so that we can close the shape
        area += (x2 - x1) * (y2 + y1)
    nature = area < 0
    return nature

def convexity(a, b, c): # checks if our points have a convex angle or not. This is important for the ear clipping algorithm, as we need to find points that stick out (convex angles) to cut off as triangles.
    ear = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) > 0  # cross product
    return ear

def sign(p1, p2, p3): # tells us which side of the line our point p1 is on of the line from p2 to p3 by cross product
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def point_pos(p, a, b, c): # tells us if we have apoint which lies inside the traingle made by a b c

    # if points are on left side of the lines  ab, bc, ca
    b1 = sign(p, a, b) < 0.0 
    b2 = sign(p, b, c) < 0.0
    b3 = sign(p, c, a) < 0.0
    return b1 == b2 == b3 # point is on the same side (left or right) of all three lines iff all three have neagtive cross product

def earclip(polygon_coords): # breaks our polygon into ears
    triangles = []
    vertices = polygon_coords[:] # so that original vertex arent affected

    if counterclockwhole(vertices) == False:
        vertices.reverse() # make them in counterclockwise direction

    while len(vertices) > 3: # stop when we have three vertices left, as a traingle has to have 3 vertices
        ear_found = False
        for i in range(len(vertices)): # keep ietrating over our vertices and take theree together
            prev = vertices[i - 1]
            curr = vertices[i]
            nxt = vertices[(i + 1) % len(vertices)] # % so that we can go back to first vertex again

            if convexity(prev, curr, nxt) == True: # too remove ear, it has to be a convex ear
                triangle = [prev, curr, nxt]
                contains_point = False 
                for p in vertices: # checks if we have point inside ear
                    if p not in triangle:
                        if point_pos(p, triangle[0], triangle[1], triangle[2]):
                            contains_point = True
                            break


                if contains_point == False: # if no point inside ear then we can clip it
                    triangles.append(triangle)
                    del vertices[i]
                    ear_found = True
                    break
        if not ear_found:
            raise ValueError("No ear found — the polygon might be invalid or not simple.")

    triangles.append(vertices) # when finally 3 points are left, directly add to the list of traingles
    return triangles

def build_graph_from_triangles(triangles): # builds the graph in a way that i can colour it later
    G = nx.Graph() # initialize empty graph
    point_to_id = {} # dict to store each number a value 
    id_to_point = {} # opposite dict so we can display it later
    counter = 0

    for tri in triangles: # loop through each triangle
        for point in tri: # loop through the 3 points
            if point not in point_to_id: # we add it to our dct to make graph
                point_to_id[point] = counter
                id_to_point[counter] = point # keys to retrive the nodes later
                counter += 1 # increase counter for next point

    for pid, coord in id_to_point.items(): # loops through all the points and their numbers
        G.add_node(pid, pos=coord) # add node in the graph at position = cordinate

    for tri in triangles: # loop through traingles
        ids = [] # get their ids in the dct
        for p in tri:
            ids.append(point_to_id[p])
        G.add_edge(ids[0], ids[1])
        G.add_edge(ids[1], ids[2]) # now add lines between all three points
        G.add_edge(ids[2], ids[0])

    # print(point_to_id)
    # print(id_to_point)

    return G, id_to_point

def colorgraph_fisk(triangles, id_to_point):# colours the graph using Fisk's algorithm
    dual_G = nx.Graph() # initialize dual graph
    diagonals = {}  # which traingles share a diagonal
    for i in range(len(triangles)): # loop through all triangles
        dual_G.add_node(i) #first add node
        for j in range(i + 1, len(triangles)): # loop through all traingles after current triangle so we can compare the coloring
            shared_vertices = [] # get the common vertices
            for p1 in triangles[i]:
                for p2 in triangles[j]:
                    if p1 == p2: # Triangles share a diagonal 
                        shared_vertices.append(p1) 
            if len(shared_vertices) == 2: # they share a diagonal if they exactly have 2 common points
                dual_G.add_edge(i, j) # so we add an edge over the shared diagonal, i.e. from triangle i to j
                u, v = tuple(sorted(shared_vertices, key=lambda p: (p[0], p[1]))) # sorting because we want our vertices in the exact order they were clicked
                diagonals[(u, v)] = (i, j) # add in our dict what diagonal belongs to what traingle 
    print("Dual graph:", list(dual_G.edges()))
    print("dias:", list(diagonals.items())) # debugging check

    G = nx.Graph() # initalize our new graph
    point_to_id = {} # dict to go from point to its id
    for i, p in id_to_point.items():
        point_to_id[p] = i 

    for pid, coord in id_to_point.items(): # loops through the ids to get points(coordinates)
        G.add_node(pid, pos=coord)

    for tri in triangles: #add edges between these points/nodes
        ids = [point_to_id[p] for p in tri]
        G.add_edge(ids[0], ids[1])
        G.add_edge(ids[1], ids[2])
        G.add_edge(ids[2], ids[0])


    vertex_colors = {}
    # Start with the first triangle (index 0)
    tri = triangles[0]
    v0, v1, v2 = [point_to_id[p] for p in tri]
    vertex_colors[v0] = 0  # Red
    vertex_colors[v1] = 1  # Green
    vertex_colors[v2] = 2  # Blue
    print(f"Colored triangle 0: {v0}=0, {v1}=1, {v2}=2")

    # Traverse the dual graph using BFS to process triangles
    visited_triangles = {0}
    queue = [(0, tri)]
    while queue:
        tri_idx, current_tri = queue.pop(0)
        # Find adjacent triangles in the dual graph
        for next_tri_idx in dual_G.neighbors(tri_idx):
            if next_tri_idx in visited_triangles:
                continue
            next_tri = triangles[next_tri_idx]
            visited_triangles.add(next_tri_idx)

            # Identify the shared diagonal
            shared_vertices = set(current_tri) & set(next_tri)
            u, v = tuple(sorted(shared_vertices, key=lambda p: (p[0], p[1])))
            u_id, v_id = point_to_id[u], point_to_id[v]

            # Get the third vertex in the new triangle
            third_vertex = (set(next_tri) - shared_vertices).pop()
            third_id = point_to_id[third_vertex]

            # Ensure u and v have different colors (they should already be colored)
            u_color = vertex_colors[u_id]
            v_color = vertex_colors[v_id]
            if u_color == v_color:
                # This shouldn’t happen if coloring is correct; force a fix
                v_neighbor_colors = {vertex_colors.get(n) for n in G.neighbors(v_id) if n in vertex_colors}
                v_allowed = {0, 1, 2} - {u_color} - v_neighbor_colors
                vertex_colors[v_id] = min(v_allowed) if v_allowed else (u_color + 1) % 3
                v_color = vertex_colors[v_id]
                print(f"Forced different color on diagonal ({u_id}, {v_id}): {u_id}={u_color}, {v_id}={v_color}")

            # Color the third vertex to complete {0, 1, 2} in the new triangle
            tri_colors = {u_color, v_color}
            missing_color = (set([0, 1, 2]) - tri_colors).pop()
            third_neighbor_colors = {vertex_colors.get(n) for n in G.neighbors(third_id) if n in vertex_colors}
            if missing_color not in third_neighbor_colors:
                vertex_colors[third_id] = missing_color
            else:
                # Conflict with neighbors; try another color and flag for fixing
                allowed = {0, 1, 2} - third_neighbor_colors
                vertex_colors[third_id] = min(allowed) if allowed else missing_color
                print(f"Conflict at node {third_id}: assigned {vertex_colors[third_id]} (may need fixing)")

            print(f"Colored triangle {next_tri_idx}: {u_id}={u_color}, {v_id}={v_color}, {third_id}={vertex_colors[third_id]}")
            queue.append((next_tri_idx, next_tri))

    # Step 4: Validate
    print("Initial coloring:", vertex_colors)
    for tri_idx, tri in enumerate(triangles):
        tri_vertices = [point_to_id[p] for p in tri]
        colors = {vertex_colors.get(v) for v in tri_vertices}
        if len(colors) != 3:
            print(f"Error: Triangle {tri_idx} does not have three distinct colors: {colors}")

    for u, v in G.edges():
        if vertex_colors[u] == vertex_colors[v]:
            print(f"Error: Adjacent vertices {u} and {v} have the same color: {vertex_colors[u]}")

    # Step 5: Fix violations with a better strategy
    iteration = 0
    while True:
        violations = [(u, v) for u, v in G.edges() if vertex_colors[u] == vertex_colors[v]]
        if not violations:
            break
        u, v = violations[0]
        print(f"Iteration {iteration}: Fixing violation at edge ({u}, {v})")

        # Recolor v, considering all triangles it’s in
        neighbor_colors = {vertex_colors.get(n) for n in G.neighbors(v) if n in vertex_colors}
        v_point = id_to_point[v]
        triangles_with_v = [i for i, tri in enumerate(triangles) if v_point in tri]

        # Try each color for v
        for c in range(3):
            if c in neighbor_colors:
                continue
            vertex_colors[v] = c
            # Check all triangles containing v
            valid = True
            for tri_idx in triangles_with_v:
                tri = triangles[tri_idx]
                tri_vertices = [point_to_id[p] for p in tri]
                colors = {vertex_colors.get(w) for w in tri_vertices}
                if len(colors) != 3:
                    valid = False
                    break
            if valid:
                break
        else:
            # No color works perfectly; minimize damage
            allowed = {0, 1, 2} - neighbor_colors
            vertex_colors[v] = min(allowed) if allowed else (vertex_colors[u] + 1) % 3
            print(f"Couldn’t fully satisfy constraints for {v}; assigned {vertex_colors[v]}")

        iteration += 1
        if iteration > len(G.nodes()):
            print("Warning: Could not resolve all violations")
            break

    # Step 6: Final validation
    print("Final coloring:", vertex_colors)
    for tri_idx, tri in enumerate(triangles):
        tri_vertices = [point_to_id[p] for p in tri]
        colors = {vertex_colors.get(v) for v in tri_vertices}
        if len(colors) != 3:
            print(f"Error after fix: Triangle {tri_idx} does not have three distinct colors: {colors}")
    for u, v in G.edges():
        if vertex_colors[u] == vertex_colors[v]:
            print(f"Error after fix: Adjacent vertices {u} and {v} have the same color: {vertex_colors[u]}")

    return vertex_colors


def drawtrigraph(G):
    pos = nx.get_node_attributes(G, 'pos')
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500)
    plt.title("Triangulated Polygon as Graph")
    plt.show()

def drawCgraph(G, coloring):
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
    G, id_to_point = build_graph_from_triangles(triangles)
    drawtrigraph(G)
    coloring = colorgraph_fisk(triangles, id_to_point)
    drawCgraph(G, coloring)
    return G, coloring



def show_min_guard_positions(G, coloring):

    color_map = {0: 'Red', 1: 'Green', 2: 'Blue'}

    for node in G.nodes:
        if node not in coloring:
            coloring[node] = 0

    color_counts = {c: list(coloring.values()).count(c) for c in range(3)}
    min_count = min(color_counts.values())
    best_colors = [c for c, count in color_counts.items() if count == min_count]

    best_color_names = [color_map[c] for c in best_colors]
    result_color_str = '/'.join(best_color_names)

    pos = nx.get_node_attributes(G, 'pos')

    node_colors = []
    for n in G.nodes():
        c = coloring[n]
        if c in best_colors:
            node_colors.append(mcolors.to_rgba(color_map[c], alpha=0.8))
        else:
            node_colors.append(mcolors.to_rgba('gray', alpha=0.4))

    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, node_color=node_colors, with_labels=True, edge_color='gray', node_size=500)

    text = (
        f"Minimum guards needed by Chvátal's theorem: {math.floor(len(coloring)/3)}\n"
        f"Actual minimum guards needed: {min_count}\n"
        f"Guards should be placed at color(s): {result_color_str}"
    )
    plt.gcf().text(0.05, 0.01, text, fontsize=10, va='bottom', ha='left', bbox=dict(facecolor='white', alpha=0.7))
    plt.title("Minimum Guard Positions Highlighted")
    plt.show()


if __name__ == "__main__":
    polygon_coords = main_polygon_generation()
    triangles = earclip(polygon_coords)
    G, coloring = color_triangulated_graph(triangles)
    show_min_guard_positions(G, coloring)



