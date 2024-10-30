class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []

def parse_newick(newick):
    root = TreeNode('')
    current_node = root
    stack = []

    for char in newick:
        if char == '(':
            current_node.children.append(TreeNode(''))
            stack.append(current_node)
            current_node = current_node.children[-1]
        elif char == ',':
            stack[-1].children.append(TreeNode(''))
            current_node = stack[-1].children[-1]
        elif char == ')':
            current_node = stack.pop()
        elif char != ';':
            current_node.name += char

    return root.children[0]

def get_leaf_names(node, names):
    if not node:
        return

    if not node.children:
        names.append(node.name)
    else:
        for child in node.children:
            get_leaf_names(child, names)

def split_tree_into_clusters(newick_str, num_clusters):
    root = parse_newick(newick_str)
    leaf_names = []
    get_leaf_names(root, leaf_names)
    num_leaves = len(leaf_names)

    # Calculate the number of leaves per cluster
    leaves_per_cluster = num_leaves // num_clusters
    remainder = num_leaves % num_clusters

    # Assign leaves to clusters
    clusters = []
    cluster_index = 0
    leaves_assigned = 0
    cluster_size = 0
    for leaf_name in leaf_names:
        if cluster_size == 0:
            cluster_size = leaves_per_cluster + (1 if remainder > 0 else 0)
            remainder -= 1

            # Start a new cluster
            clusters.append([])
            cluster_index += 1

        clusters[cluster_index - 1].append(leaf_name)
        leaves_assigned += 1
        cluster_size -= 1

    return clusters