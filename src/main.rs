/*
SCRIPT extract_extant_script
--------------------------------------------------------------------------------
This script takes:
--------------------------------------------------------------------------------
Output:
*/

#[macro_use]
extern crate pest_derive;
#[allow(dead_code)]
use pest::Parser;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::env;
use std::fs::File;
extern crate regex;

#[derive(Clone, Debug)]


struct Node {
    name: String,
    left_child: Option<Box<Node>>,
    right_child: Option<Box<Node>>,
    parent: Option<usize>, // Using index for parent
    depth: Option<f64>,
    length: f64,
}
#[derive(Clone, Debug)]
struct FlatNode {
    name: String,
    left_child: Option<usize>,
    right_child: Option<usize>,
    parent: Option<usize>,
    depth: Option<f64>,
    length: f64,
}
#[derive(Parser)]
#[grammar = "newick.pest"]
struct NewickParser;

// --------------------------------

fn give_depth(node: &mut Node, depth: f64) {
    /*
    Recursively gives the depth of nodes in a "Node" arborescent tree.
    ----------------------------------
    Input:
    - A node, and its depth in the tree
    ----------------------------------
    Output:
    - None. The tree is modified in place via a mutable reference (&mut)
    */
    node.depth = Some(depth);
    if let Some(left_child) = &mut node.left_child {
        give_depth(left_child, depth + left_child.length);
    }
    if let Some(right_child) = &mut node.right_child {
        give_depth(right_child, depth + right_child.length);
    }
}
// The following two functions are used to convert a newick string into an arborescent Tree structure, where each "Node" object owns its two children in a Box object.
fn newick_to_tree(pair: pest::iterators::Pair<Rule>) -> Vec<Node> {
    /*
    The input .nwk file of the script is a concatenation of newick trees, separated by ";".
    This function returns a vector containing each of these trees.
    --------------------------------
    Input:
    - A "pair" object of type "newick" (see newick.pest grammar) which is a newick string representing one single tree.
    ----------------------------------
    Output:

    */
    let mut vec_trees:Vec<Node> = Vec::new();
    for inner in pair.into_inner() {
        let tree = handle_pair(inner);
        vec_trees.push(tree.unwrap());
    }
    return vec_trees
}
fn handle_pair(pair: pest::iterators::Pair<Rule>) -> Option<Node> {
    /* 
    Recursively parses a newick string representing a single newick tree, according to the grammar contained in newick.pest.
    ----------------------------------
    Input:
    - A "pair" object representing a part of this string.
    ----------------------------------
    Output:
    - Either a node if the input string represents a node, or nothing in cases where it does not (e.g. if the string is the length of a node).
    */
    match pair.as_rule() {
        Rule::newick => {
            // newick = { subtree ~ ";" }
            // For now, we only do the implementation for one single tree.
            // This case is trivial : we just call handle_pair on the only tree we found.
            for inner in pair.into_inner() {
                if inner.as_rule() == Rule::subtree {
                    return handle_pair(inner);
                }
            }
            None
        },
        Rule::subtree => {
            // subtree = { leaf | internal }
            // Subtree is like the choice between an inner node and a leaf. Either way, we need to pass it to handle_pair.
            handle_pair(pair.into_inner().next().unwrap())
        },
        Rule::leaf => {
            // Choose default values for the name and length of the leaf.
            // The defaults are an empty string and a 0.0 length.
            let mut name = String::new();
            let mut length = 0.0;

            // leaf = { NAME? ~ ":"? ~ LENGTH? }
            // Therefore, we use a match statement to handle the cases NAME and LENGTH because ":" is unimportant.
            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    }
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    },
                    
                    _ => {} // Ignore other rules
                }
            }
            let node = Node {
                name: name,
                left_child: None,
                right_child: None,
                parent: None,
                depth: None,
                length: length,
            };
            Some(node)
        },
        Rule::internal => {
            // internal = { "(" ~ subtree ~ "," ~ subtree ~ ")" ~ NAME? ~ ":"? ~ LENGTH? }
            
            // Initialize default values for the name and length of the internal node.
            // The defaults are an empty string and a 0.0 length.
            let mut name = String::new();
            let mut length = 0.0;
        
            let mut first_subtree = None;
            let mut second_subtree = None;
        
            // Iterate through the inner rules without assuming their order
            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::subtree => {
                        let subtree = handle_pair(inner_pair).unwrap();
                        if first_subtree.is_none() {
                            first_subtree = Some(subtree);
                        } else {
                            second_subtree = Some(subtree);
                        }
                    },
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    },
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    },
                    
                    _ => {} // Ignore other rules
                }
            }
        
            let node = Node {
                name,
                left_child: first_subtree.map(Box::new),
                right_child: second_subtree.map(Box::new),
                parent: None,
                depth: None,
                length,
            };
            Some(node)
        }
        Rule::NAME | Rule::LENGTH => {
            // We should never directly handle these outside their containing rules.
            None
        },
    }
}

// Convert a Tree in "Node" form to a Newick string.
fn node_to_newick(node: &Node) -> String {
    /* Takes a node and returns the corresponding subtree in Newick format.
        --------------------------------
        INPUT:
            - node: the node to convert to Newick format.
        OUTPUT:
            - the Newick representation of the subtree rooted at node.
        Warning: rounds the lengths to 6 decimal places.
    */
    if let (Some(left_child), Some(right_child)) = (&node.left_child, &node.right_child) {
        // This is an internal node with both left and right children.
        format!(
            "({},{}){}:{:.6}",
            node_to_newick(left_child),
            node_to_newick(right_child),
            node.name,
            node.length
        )
    } else {
        // This is a leaf node.
        format!("{}:{:.6}", node.name, node.length)
    }
}
// Convert from FlatNode to Node
fn flat_to_node(flat_tree: &[FlatNode], index: usize, parent_index: Option<usize>) -> Option<Node> {
    /*
    This function converts a flat tree into an arborescent tree recursively.
    To use it, give the flat_tree, as well as the index of the root. The vector will be traversed recursively, 
    following the descendants of each node being examined.
    ----------------------------------
    Input: 
    - A flat tree, the index of a node in the flat tree, and the index of its parent.
    ----------------------------------
    Output:
    - The corresponding node in the arborescent tree.
    ----------------------------------
    Warning: can bug if applied directly to a non-root node, which will be mistaken for a root node.
    */
    let flat_node = &flat_tree[index];
    let left_child = flat_node.left_child.and_then(|i| {
        flat_to_node(flat_tree, i, Some(index)).map(Box::new)
    });
    let right_child = flat_node.right_child.and_then(|i| {
        flat_to_node(flat_tree, i, Some(index)).map(Box::new)});
    
    Some(Node {
        name: flat_node.name.clone(),
        left_child: left_child,
        right_child: right_child,
        parent: parent_index,
        depth: flat_node.depth,
        length: flat_node.length,
    })
}
// Convert from Node to FlatNode
fn node_to_flat(node: &Node, flat_tree: &mut Vec<FlatNode>, parent: Option<usize>) -> usize {
    /* Transforms the arborescent tree into a "linear tree", which is just a vector of nodes with parent and descendants.
    ----------------------------------
    Input:
    - The root node, which contains the whole arborescent tree
    - The flat_tree vector, which is to be filled by the function
    ----------------------------------
    Output:
    - The index of the node currently being added to flat_tree (usize because the indexing of Rust vectors can only be usize).
    */
    let index = flat_tree.len();
    flat_tree.push(FlatNode {
        name: node.name.clone(),
        left_child: None,  // Will fill this in a moment
        right_child: None, // Will fill this in a moment
        parent: parent,
        depth: node.depth,
        length: node.length,
    });

    if let Some(left) = &node.left_child {
        let left_index = node_to_flat(left, flat_tree, Some(index));
        flat_tree[index].left_child = Some(left_index);
    }

    if let Some(right) = &node.right_child {
        let right_index = node_to_flat(right, flat_tree, Some(index));
        flat_tree[index].right_child = Some(right_index);
    }

    index
}
// Remove any given leaf from the tree
fn change_tree(flat_tree: &mut Vec<FlatNode>, index: usize){
    /* 
    Takes a flat tree and an index of leaf to remove from the flat tree.
    Removes the corresponding leaf as well as its parent.
    Each application of this function should take a correct phylogenetic tree,
    and output an object representing a correct phylogenetic tree along with isolated nodes.
    ----------------------------------
    Input:
    - A flat tree, and the index of a leaf to remove from the tree.
    ----------------------------------
    Output:
    - The flat tree, with the leaf removed.
    */

    // --------------------------------
    // Node indexes
    let parent_index = flat_tree[index].parent.expect("The root is apparently a leaf.");
    let sister_index;
    if flat_tree[parent_index].left_child.unwrap() == index {
        sister_index = flat_tree[parent_index].right_child.unwrap()} else 
        {sister_index = flat_tree[parent_index].left_child.unwrap()};

    let grandparent_index_opt = flat_tree[parent_index].parent;




        // The leave and its parent are removed from the tree.
        flat_tree[parent_index].parent = None;

        // The sister of the leaf becomes the child of the grandparent, because its parent has lost one of its child and hence disappears.
        flat_tree[sister_index].parent = grandparent_index_opt;


        // Changing attributes of the parent of the removed leaf is unimportant, since it will never be referenced again, and hence
        // will not be seen in a tree traversal from the root.

        // Changing the child of the grandparent which is the parent to the sister.
        if let Some(grandparent_index) = grandparent_index_opt {
            if flat_tree[grandparent_index].left_child == Some(parent_index) {
                flat_tree[grandparent_index].left_child = Some(sister_index)} else
                {flat_tree[grandparent_index].right_child = Some(sister_index)
                };
            }
        // None of the depths have to be updated. The lengths are not updated either, since they will
        // be reconstructed from the depths at the end.
    }
// Apply change_tree to remove all unsampled leaves
fn remove_all_unsampled(flat_tree: &mut Vec<FlatNode>, list_indexes: &Vec<usize>) {
    /*
    Takes a flat tree and a vector of indexes of leaves to remove from the tree.
    Removes all the corresponding leaves as well as their parents.
    Acts in place and does not return anything.
    ----------------------------------
    Input:
    - A flat tree, and a vector of indexes of leaves to remove from the tree.
    ----------------------------------
    Output:
    - The flat tree, with the leaves removed.
    */
    for index in list_indexes {
        change_tree(flat_tree, *index);
    }
}
// Extract leaves from the tree
fn find_all_leaves(flat_tree: &Vec<FlatNode>) -> Vec<usize> {
    // Leaves are exactly nodes without children, so they are easy to find.
    let mut leaves: Vec<usize> = Vec::new();
    for (i, node) in flat_tree.iter().enumerate() {
        if node.left_child.is_none() && node.right_child.is_none() {
            leaves.push(i);
        }
    }
    leaves
}

fn find_deepest_nodes(flat_tree: &Vec<FlatNode>, nb_leaves: &usize) -> Vec<usize> {
    // Sort the nodes by depth
    let new_flat_tree = flat_tree.clone();

    let mut temp_iter = new_flat_tree
        .iter()
        .enumerate()
        .filter(|(_, node)| node.left_child.is_none() && node.right_child.is_none())
        .map(|(i, node)| (i, node.depth.unwrap()))
        .collect::<Vec<(usize, f64)>>();
    temp_iter.sort_by(|(_, depth1), (_, depth2)| depth1.partial_cmp(depth2).unwrap());
    // Now take the indices of the deepest nodes
    let mut deepest_nodes: Vec<usize> = Vec::new();
    for (i, _) in temp_iter.iter().rev().take(*nb_leaves) {
        deepest_nodes.push(*i);
    }
    deepest_nodes
}

// Find the leaves to be removed (Complement of leaves to be kept)
fn leaves_to_be_removed(leaves: &Vec<usize>, sampled_leaves: &Vec<usize>) -> Vec<usize> {
    // This function takes a vector of all the leaves in the tree, and a vector of the sampled leaves.
    // It returns a vector of the indexes of the leaves to be removed.
    let mut leaves_to_be_removed: Vec<usize> = Vec::new();
    for leaf in leaves {
        if !sampled_leaves.contains(leaf) {
            leaves_to_be_removed.push(*leaf);
        }
    }
    leaves_to_be_removed
}
// Give the index of the root of the tree
fn find_root(flat_tree: &Vec<FlatNode>, true_leaf: usize) -> usize {
    /*
    Uses an index of a leaf to find the root of the tree.
     */
    let mut current_node = true_leaf;
    let mut current_parent = flat_tree[current_node].parent;
    while current_parent.is_some() {
        current_node = current_parent.unwrap();
        current_parent = flat_tree[current_node].parent;
    }
    current_node
}
// Find the leaves in the gene tree corresponding to the leaves in the species tree
// Corrects lengths of the tree after removing some nodes
fn depths_to_lengths(node: &mut Node, parent_depth: f64) {
    /*
    Computes the lengths of nodes from the depths.
    Necessary because the gene transfer functions only modify the depth of the
    transfered node and not lengths for convenience.
    ----------------------------------
    Input:
    - A "Node" object representing an arborescent tree.
    - The depth of its parent
    ----------------------------------
    Output:
    - None. The tree is modified in place via a mutable reference.
    */
    let depth = node.depth.unwrap();
    if node.parent.is_some() {
        node.length = depth - parent_depth;
    }
    if let Some(left_child) = &mut node.left_child {
        depths_to_lengths(left_child, depth);
    }
    if let Some(right_child) = &mut node.right_child {
        depths_to_lengths(right_child, depth);
    }
}
// Returns and saves a gene tree with sampled leaves

fn species_tree_sample_to_string(species_tree_path: &str, output_dir: &str, nb_leaves: &usize) -> Result<(String, Vec<String>, Vec<String>), io::Error> {
    /*
    Samples the given number of nodes in the species tree, and returns the list of names of sampled nodes, as well as unsampled nodes.
    1. Convert the species tree to a flat tree.
    2. Sample the leaves.
    3. Reconstruct the sampled species tree, by SPR moves, disconnecting unsampled leaves.
        3.1. Apply the function remove_all_unsampled to the flat tree.
        3.2. Use any of the leaves in sampled_leaves to find the root of the new tree.
        3.3. Convert the flat tree to a classical tree by using the root index. This will effectively ignore all the nodes which are now unsampled.
    4. Convert the species tree to a Newick string.
    ----------------------------------
    Input:
    - The path to the species tree file.
    - The number of sampled nodes.
    - The path to the output directory.
    ----------------------------------
    Output:
    - The Newick string of the sampled species tree.
    - The list of names of sampled nodes.
    - The list of names of unsampled nodes.
    */
    // -1. Ensure the output directory exists
    let output_path = Path::new(output_dir);
    if !output_path.exists() {
        fs::create_dir_all(output_path)?;
    }

    // 0. Open the species tree and convert it to a flat tree.
    let species_tree = fs::read_to_string(species_tree_path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let species_tree = species_tree.trim();
    let mut species_tree = NewickParser::parse(Rule::newick, species_tree)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut species_tree = newick_to_tree(species_tree.next().unwrap());

    let root_length = species_tree[0].length.clone();
    let _ = give_depth(&mut species_tree[0], root_length);
    let mut flat_tree: Vec<FlatNode> = Vec::new();
    // Populate the flat tree with the species tree.
    let _ = node_to_flat(&species_tree[0], &mut flat_tree, None);

    // 1. Sample the leaves.
    let sampled_leaves =find_deepest_nodes(&flat_tree, nb_leaves);
    // 2. Construct the species tree by removing the unsampled leaves from the species tree.
    // 2.1. Apply the function remove_all_unsampled to the flat tree.
    let leaves_to_be_removed = leaves_to_be_removed(&find_all_leaves(&flat_tree), &sampled_leaves);
    remove_all_unsampled(&mut flat_tree, &leaves_to_be_removed);
    // 2.2. Use any of the leaves in sampled_leaves to find the root of the new tree.
    let root_of_species_tree = find_root(&flat_tree, sampled_leaves[0]);
    // 2.3. Convert the flat tree to a classical tree by using the root index. This will effectively ignore all the nodes which are now unsampled.
    let mut reconstructed_tree = flat_to_node(&flat_tree, root_of_species_tree, None)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to convert flat tree to node tree"))?;
    
    // The lengths of the nodes are not correct, so we need to update them.
    let root_depth = reconstructed_tree.depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    depths_to_lengths(&mut reconstructed_tree, root_depth);

    // 3. Convert the species tree to a Newick string.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the species tree as a Newick string in the output directory
    let species_filename = format!("extant_species_tree.nwk");
    let species_filename = Path::new(output_dir).join(species_filename);
    let mut species_file = File::create(species_filename)?;
    species_file.write_all(reconstructed_newick.as_bytes())?;

    // 4. Return the Newick string of the sampled species tree, as well as the list of names of sampled nodes, and the list of names of unsampled nodes.
    let sampled_leaves_names: Vec<String> = sampled_leaves.iter().map(|i| flat_tree[*i].name.clone()).collect();
    let leaves_to_be_removed_names: Vec<String> = leaves_to_be_removed.iter().map(|i| flat_tree[*i].name.clone()).collect();
    Ok((reconstructed_newick, sampled_leaves_names, leaves_to_be_removed_names))
}


fn main() {
    // Read the arguments
    let args: Vec<String> = env::args().collect();

    // Ensure the correct number of arguments are provided
    if args.len() != 4 {
        eprintln!("Usage: {} <species_tree_path> <n_extant_nodes> <output_dir>", args[0]);
        eprintln!("Received arguments: {:?}", args);
        panic!("Error with the input arguments! See error above.");
    }

    let species_tree_path = &args[1];
    let n_extant = match args[2].parse::<usize>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Error: n_extant_nodes must be an integer. Received: {}", args[4]);
            eprintln!("All arguments: {:?}", args);
            return;
        }
    };
    let output_dir = &args[3];



    // Sample the species tree
    let result = species_tree_sample_to_string(species_tree_path, output_dir, &n_extant);
    match result {
        Ok((_, sampled_names, removed_names)) => (sampled_names, removed_names),
        Err(e) => {
            eprintln!("Error during species tree sampling: {}", e);
            eprintln!("Species Tree Path: {}", species_tree_path);
            eprintln!("Number of Sampled Nodes: {}", n_extant);
            eprintln!("Output Directory: {}", output_dir);
            return;
        }
    };
}