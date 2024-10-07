from graphviz import Digraph

# Create a new directed graph
dot = Digraph()

# Add nodes for level 1 (root)
dot.node('A', '')

# Add nodes for level 2
dot.node('B', '')
dot.node('C', '')

# Add nodes for level 3
dot.node('D', '')
dot.node('E', '')

# Add edges from level 1 to level 2
dot.edge('A', 'B')
dot.edge('A', 'C')

# Add edges from level 2 to level 3
dot.edge('B', 'D')
dot.edge('B', 'E')
dot.edge('C', 'E')

# Save the source code to a file and render it
dot.attr(rankdir='LR')
#dot.render('simple_tree_no_labels', format='png', view=True)

import matplotlib.pyplot as plt
import numpy as np

# Define the range of x values
x = np.linspace(-10, 10, 400)

# Compute the y values for the parabola y = x^2
y = x**2

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(x, y, label='$y = x^2$', color='b')

# Add title and labels
plt.title('Simple Concave Upward Parabola')
plt.xlabel('x')
plt.ylabel('y')

# Add a legend
plt.legend()

# Show the plot
plt.show()
