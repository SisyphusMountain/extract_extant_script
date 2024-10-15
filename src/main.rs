use std::env;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::Path;

extern crate regex;

use pest::Parser;
use newick_parser::node::{FlatTree, TraversalOrder};
use newick_parser::newick::{newick_to_tree, node_to_newick, NewickParser, Rule};

/// Removes a given leaf from the flat tree.
///
/// This function takes a flat tree and an index of a leaf to remove. It removes the corresponding
/// leaf as well as its parent. Each application of this function should maintain a correct
/// phylogenetic tree, outputting an object representing a correct phylogenetic tree along with
/// isolated nodes.
///
/// # Arguments
///
/// * `flat_tree` - A mutable reference to the flat tree.
/// * `index` - The index of the leaf to remove from the tree.
fn change_tree(flat_tree: &mut FlatTree, index: usize) {
    // Node indexes
    let parent_index = flat_tree[index]
        .parent
        .expect("The root is apparently a leaf.");
    let sister_index = if flat_tree[parent_index].left_child.unwrap() == index {
        flat_tree[parent_index].right_child.unwrap()
    } else {
        flat_tree[parent_index].left_child.unwrap()
    };

    let grandparent_index_opt = flat_tree[parent_index].parent;

    // The leaf and its parent are removed from the tree.
    flat_tree[parent_index].parent = None;

    // The sister of the leaf becomes the child of the grandparent.
    flat_tree[sister_index].parent = grandparent_index_opt;

    // Change the child of the grandparent from parent to sister.
    if let Some(grandparent_index) = grandparent_index_opt {
        if flat_tree[grandparent_index].left_child == Some(parent_index) {
            flat_tree[grandparent_index].left_child = Some(sister_index);
        } else {
            flat_tree[grandparent_index].right_child = Some(sister_index);
        }
    }
    // Depths and lengths will be recalculated later.
}

/// Removes all unsampled leaves from the flat tree.
///
/// This function takes a flat tree and a vector of indexes of leaves to remove from the tree.
/// It removes all the corresponding leaves as well as their parents. The function acts in place
/// and does not return anything.
///
/// # Arguments
///
/// * `flat_tree` - A mutable reference to the flat tree.
/// * `list_indexes` - A vector of indexes of leaves to remove from the tree.
fn remove_all_unsampled(flat_tree: &mut FlatTree, list_indexes: &Vec<usize>) {
    for index in list_indexes {
        change_tree(flat_tree, *index);
    }
}

/// Extracts all leaves from the flat tree.
///
/// Leaves are nodes without children.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the flat tree.
///
/// # Returns
///
/// A vector containing the indexes of all leaves in the tree.
fn find_all_leaves(flat_tree: &FlatTree) -> Vec<usize> {
    flat_tree
        .iter(TraversalOrder::PreOrder)
        .enumerate()
        .filter(|(_, node)| node.left_child.is_none() && node.right_child.is_none())
        .map(|(i, _)| i)
        .collect()
}

/// Finds the deepest nodes (leaves) in the flat tree.
///
/// This function sorts the leaves by their depth in descending order and returns the top `nb_leaves` deepest nodes.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the flat tree.
/// * `nb_leaves` - The number of deepest leaves to find.
///
/// # Returns
///
/// A vector containing the indexes of the deepest leaves.
fn find_deepest_nodes(flat_tree: &FlatTree, nb_leaves: usize) -> Vec<usize> {
    // Find the deepest leaves.
    let mut leaves_with_depths: Vec<(usize, f64)> = flat_tree
        .iter(TraversalOrder::PreOrder)
        .enumerate()
        .filter(|(_, node)| node.left_child.is_none() && node.right_child.is_none())
        .map(|(i, node)| (i, node.depth.unwrap()))
        .collect();

    // Sort by depth in descending order.
    leaves_with_depths.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Take the top `nb_leaves` deepest nodes.
    leaves_with_depths.iter().take(nb_leaves).map(|(i, _)| *i).collect()
}

/// Finds the leaves to be removed from the tree.
///
/// This function computes the complement of the sampled leaves, i.e., the leaves that are not in the sampled list.
///
/// # Arguments
///
/// * `leaves` - A vector of all the leaves in the tree.
/// * `sampled_leaves` - A vector of the sampled leaves.
///
/// # Returns
///
/// A vector containing the indexes of the leaves to be removed.
fn leaves_to_be_removed(leaves: &Vec<usize>, sampled_leaves: &Vec<usize>) -> Vec<usize> {
    leaves
        .iter()
        .filter(|leaf| !sampled_leaves.contains(leaf))
        .cloned()
        .collect()
}

/// Determines the root index of the tree using a given leaf.
///
/// This function traverses from a leaf up to the root by following parent pointers.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the flat tree.
/// * `true_leaf` - The index of a leaf node in the tree.
///
/// # Returns
///
/// The index of the root node in the flat tree.
fn find_root(flat_tree: &FlatTree, true_leaf: usize) -> usize {
    let mut current_node = true_leaf;
    let mut current_parent = flat_tree[current_node].parent;
    while let Some(parent) = current_parent {
        current_node = parent;
        current_parent = flat_tree[current_node].parent;
    }
    current_node
}

/// Samples the species tree and returns a Newick string along with sampled and removed leaf names.
///
/// This function performs the following steps:
/// 1. Converts the species tree to a flat tree.
/// 2. Samples the deepest leaves.
/// 3. Removes unsampled leaves from the tree.
/// 4. Reconstructs the tree and updates node lengths based on depths.
/// 5. Converts the reconstructed tree to a Newick string and saves it.
///
/// # Arguments
///
/// * `species_tree_path` - The path to the species tree file in Newick format.
/// * `output_dir` - The output directory where the sampled tree will be saved.
/// * `nb_leaves` - The number of leaves to sample.
///
/// # Returns
///
/// A `Result` containing:
/// - The Newick string of the sampled species tree.
/// - A vector of names of sampled leaves.
/// - A vector of names of unsampled (removed) leaves.
///
/// If an error occurs, an `io::Error` is returned.
fn species_tree_sample_to_string(
    species_tree_path: &str,
    output_dir: &str,
    nb_leaves: usize,
) -> Result<(String, Vec<String>, Vec<String>), io::Error> {
    // Ensure the output directory exists
    let output_path = Path::new(output_dir);
    if !output_path.exists() {
        fs::create_dir_all(output_path)?;
    }

    // Open the species tree and convert it to a flat tree.
    let species_tree_str = fs::read_to_string(species_tree_path)?;
    let species_tree_str = species_tree_str.trim();

    // Parse the Newick string, mapping the pest error to io::Error
    let mut pairs = NewickParser::parse(Rule::newick, species_tree_str)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

    let mut node_tree = newick_to_tree(
        pairs
            .next()
            .expect("Error converting the Newick file")
    )
    .pop()
    .expect("Error: no tree found");

    // Assign depths
    node_tree.zero_root_length();
    node_tree.assign_depths(0.0);

    // Convert to FlatTree
    let mut flat_tree = node_tree.to_flat_tree();

    // Sample the leaves.
    let sampled_leaves = find_deepest_nodes(&flat_tree, nb_leaves);

    // Construct the species tree by removing the unsampled leaves.
    let leaves = find_all_leaves(&flat_tree);
    let leaves_to_be_removed = leaves_to_be_removed(&leaves, &sampled_leaves);
    remove_all_unsampled(&mut flat_tree, &leaves_to_be_removed);

    // Update the root of the flat tree.
    let root_of_species_tree = find_root(&flat_tree, sampled_leaves[0]);
    flat_tree.root = root_of_species_tree;

    // Convert the flat tree back to a Node tree.
    let mut reconstructed_tree = flat_tree.to_node();

    // Update lengths based on depths.
    let root_depth = reconstructed_tree
        .depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    reconstructed_tree.depths_to_lengths(root_depth);

    // Convert the tree to Newick format.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the species tree to the output directory.
    let species_filename = Path::new(output_dir).join("extant_species_tree.nwk");
    let mut species_file = File::create(species_filename)?;
    species_file.write_all(reconstructed_newick.as_bytes())?;

    // Return the Newick string and the lists of sampled and removed leaf names.
    let sampled_leaves_names: Vec<String> = sampled_leaves
        .iter()
        .map(|i| flat_tree[*i].name.clone())
        .collect();
    let leaves_to_be_removed_names: Vec<String> = leaves_to_be_removed
        .iter()
        .map(|i| flat_tree[*i].name.clone())
        .collect();

    Ok((reconstructed_newick, sampled_leaves_names, leaves_to_be_removed_names))
}

fn main() {
    // Read the arguments
    let args: Vec<String> = env::args().collect();
    // This script takes the n most recent nodes, samples them from a tree, and returns the sampled tree.
    // If we know the species tree has n extant nodes, we can sample the n most recent nodes to get the extant species tree.
    // Ensure the correct number of arguments are provided
    if args.len() != 4 {
        eprintln!(
            "Usage: {} <species_tree_path> <n_extant_nodes> <output_dir>",
            args[0]
        );
        eprintln!("Received arguments: {:?}", args);
        panic!("Error with the input arguments! See error above.");
    }

    let species_tree_path = &args[1];
    let n_extant = match args[2].parse::<usize>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!(
                "Error: n_extant_nodes must be an integer. Received: {}",
                args[2]
            );
            eprintln!("All arguments: {:?}", args);
            return;
        }
    };
    let output_dir = &args[3];

    // Sample the species tree
    let result = species_tree_sample_to_string(species_tree_path, output_dir, n_extant);
    match result {
        Ok((_, sampled_names, removed_names)) => {
            // You can use sampled_names and removed_names if needed
            println!("Sampled Leaves: {:?}", sampled_names);
            println!("Removed Leaves: {:?}", removed_names);
        }
        Err(e) => {
            eprintln!("Error during species tree sampling: {}", e);
            eprintln!("Species Tree Path: {}", species_tree_path);
            eprintln!("Number of Sampled Nodes: {}", n_extant);
            eprintln!("Output Directory: {}", output_dir);
        }
    };
}
