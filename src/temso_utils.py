import numpy as np
import pandas as pd
import scanpy as sc
import torch
import logging
import numpy as np
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import pairwise_distances
import networkx as nx
import re

from scipy.sparse.csc import csc_matrix
from scipy.sparse.csr import csr_matrix
from scipy.stats import pearsonr

from . import Mapper as mo
from . import utils as ut

from Bio import Phylo
from Bio.Phylo.Newick import Clade

logging.getLogger().setLevel(logging.INFO)


def pp_adatas(adata_sc, adata_sp, genes=None, gene_to_lowercase = True):
    """
    Pre-process AnnDatas so that they can be mapped. Specifically:
    - Remove genes that all entries are zero
    - Find the intersection between adata_sc, adata_sp and given marker gene list, save the intersected markers in two adatas
    - Calculate density priors and save it with adata_sp

    Args:
        adata_sc (AnnData): single cell data
        adata_sp (AnnData): spatial expression data
        genes (List): Optional. List of genes to use. If `None`, all genes are used.
    
    Returns:
        update adata_sc by creating `uns` `training_genes` `overlap_genes` fields 
        update adata_sp by creating `uns` `training_genes` `overlap_genes` fields and creating `obs` `rna_count_based_density` & `uniform_density` field
    """

    # remove all-zero-valued genes
    sc.pp.filter_genes(adata_sc, min_cells=1)
    sc.pp.filter_genes(adata_sp, min_cells=1)

    if genes is None:
        # Use all genes
        genes = adata_sc.var.index
               
    # put all var index to lower case to align
    if gene_to_lowercase:
        adata_sc.var.index = [g.lower() for g in adata_sc.var.index]
        adata_sp.var.index = [g.lower() for g in adata_sp.var.index]
        genes = list(g.lower() for g in genes)

    adata_sc.var_names_make_unique()
    adata_sp.var_names_make_unique()
    

    # Refine `marker_genes` so that they are shared by both adatas
    genes = list(set(genes) & set(adata_sc.var.index) & set(adata_sp.var.index))
    # logging.info(f"{len(genes)} shared marker genes.")

    adata_sc.uns["training_genes"] = genes
    adata_sp.uns["training_genes"] = genes
    logging.info(
        "{} training genes are saved in `uns``training_genes` of both single cell and spatial Anndatas.".format(
            len(genes)
        )
    )

    # Find overlap genes between two AnnDatas
    overlap_genes = list(set(adata_sc.var.index) & set(adata_sp.var.index))
    # logging.info(f"{len(overlap_genes)} shared genes.")

    adata_sc.uns["overlap_genes"] = overlap_genes
    adata_sp.uns["overlap_genes"] = overlap_genes
    logging.info(
        "{} overlapped genes are saved in `uns``overlap_genes` of both single cell and spatial Anndatas.".format(
            len(overlap_genes)
        )
    )

    # Calculate uniform density prior as 1/number_of_spots
    adata_sp.obs["uniform_density"] = np.ones(adata_sp.X.shape[0]) / adata_sp.X.shape[0]
    logging.info(
        f"uniform based density prior is calculated and saved in `obs``uniform_density` of the spatial Anndata."
    )

    # Calculate rna_count_based density prior as % of rna molecule count
    rna_count_per_spot = np.array(adata_sp.X.sum(axis=1)).squeeze()
    adata_sp.obs["rna_count_based_density"] = rna_count_per_spot / np.sum(rna_count_per_spot)
    logging.info(
        f"rna count based density prior is calculated and saved in `obs``rna_count_based_density` of the spatial Anndata."
    )
        

def adata_to_cluster_expression(adata, cluster_label, scale=True, add_density=True):
    """
    Convert an AnnData to a new AnnData with cluster expressions. Clusters are based on `cluster_label` in `adata.obs`.  The returned AnnData has an observation for each cluster, with the cluster-level expression equals to the average expression for that cluster.
    All annotations in `adata.obs` except `cluster_label` are discarded in the returned AnnData.
    
    Args:
        adata (AnnData): single cell data
        cluster_label (String): field in `adata.obs` used for aggregating values
        scale (bool): Optional. Whether weight input single cell by # of cells in cluster. Default is True.
        add_density (bool): Optional. If True, the normalized number of cells in each cluster is added to the returned AnnData as obs.cluster_density. Default is True.

    Returns:
        AnnData: aggregated single cell data

    """
    try:
        value_counts = adata.obs[cluster_label].value_counts(normalize=True)
    except KeyError as e:
        raise ValueError("Provided label must belong to adata.obs.")
    unique_labels = value_counts.index
    new_obs = pd.DataFrame({cluster_label: unique_labels})
    adata_ret = sc.AnnData(obs=new_obs, var=adata.var, uns=adata.uns)

    X_new = np.empty((len(unique_labels), adata.shape[1]))
    for index, l in enumerate(unique_labels):
        if not scale:
            X_new[index] = adata[adata.obs[cluster_label] == l].X.mean(axis=0)
        else:
            X_new[index] = adata[adata.obs[cluster_label] == l].X.sum(axis=0)
    adata_ret.X = X_new

    if add_density:
        adata_ret.obs["cluster_density"] = adata_ret.obs[cluster_label].map(
            lambda i: value_counts[i]
        )

    return adata_ret


def map_cells_to_space(
    adata_sc,
    adata_sp,
    cv_train_genes=None,
    cluster_label=None,
    mode="cells",
    device="cpu",
    learning_rate=0.1,
    num_epochs=1000,
    scale=True,
    lambda_d=0,
    lambda_g1=1,
    lambda_g2=0,
    lambda_r=0,
    lambda_count=1,
    lambda_f_reg=1,
    target_count=None,
    random_state=None,
    verbose=True,
    density_prior='rna_count_based',
):
    """
    Map single cell data (`adata_sc`) on spatial data (`adata_sp`).
    
    Args:
        adata_sc (AnnData): single cell data
        adata_sp (AnnData): gene spatial data
        cv_train_genes (list): Optional. Training gene list. Default is None.
        cluster_label (str): Optional. Field in `adata_sc.obs` used for aggregating single cell data. Only valid for `mode=clusters`.
        mode (str): Optional. Tangram mapping mode. Currently supported: 'cell', 'clusters', 'constrained'. Default is 'cell'.
        device (string or torch.device): Optional. Default is 'cpu'.
        learning_rate (float): Optional. Learning rate for the optimizer. Default is 0.1.
        num_epochs (int): Optional. Number of epochs. Default is 1000.
        scale (bool): Optional. Whether weight input single cell data by the number of cells in each cluster, only valid when cluster_label is not None. Default is True.
        lambda_d (float): Optional. Hyperparameter for the density term of the optimizer. Default is 0.
        lambda_g1 (float): Optional. Hyperparameter for the gene-voxel similarity term of the optimizer. Default is 1.
        lambda_g2 (float): Optional. Hyperparameter for the voxel-gene similarity term of the optimizer. Default is 0.
        lambda_r (float): Optional. Strength of entropy regularizer. An higher entropy promotes probabilities of each cell peaked over a narrow portion of space. lambda_r = 0 corresponds to no entropy regularizer. Default is 0.
        lambda_count (float): Optional. Regularizer for the count term. Default is 1. Only valid when mode == 'constrained'
        lambda_f_reg (float): Optional. Regularizer for the filter, which promotes Boolean values (0s and 1s) in the filter. Only valid when mode == 'constrained'. Default is 1.
        target_count (int): Optional. The number of cells to be filtered. Default is None.
        random_state (int): Optional. pass an int to reproduce training. Default is None.
        verbose (bool): Optional. If print training details. Default is True.
        density_prior (str, ndarray or None): Spatial density of spots, when is a string, value can be 'rna_count_based' or 'uniform', when is a ndarray, shape = (number_spots,). This array should satisfy the constraints sum() == 1. If None, the density term is ignored. Default value is 'rna_count_based'.

    Returns:
        a cell-by-spot AnnData containing the probability of mapping cell i on spot j.
        The `uns` field of the returned AnnData contains the training genes.
    """

    # check invalid values for arguments
    if lambda_g1 == 0:
        raise ValueError("lambda_g1 cannot be 0.")

    if (type(density_prior) is str) and (
        density_prior not in ["rna_count_based", "uniform", None]
    ):
        raise ValueError("Invalid input for density_prior.")

    if density_prior is not None and (lambda_d == 0 or lambda_d is None):
        lambda_d = 1

    if lambda_d > 0 and density_prior is None:
        raise ValueError("When lambda_d is set, please define the density_prior.")

    if mode not in ["clusters", "cells", "constrained"]:
        raise ValueError('Argument "mode" must be "cells", "clusters" or "constrained')

    if mode == "clusters" and cluster_label is None:
        raise ValueError("A cluster_label must be specified if mode is 'clusters'.")

    if mode == "constrained" and not all([target_count, lambda_f_reg, lambda_count]):
        raise ValueError(
            "target_count, lambda_f_reg and lambda_count must be specified if mode is 'constrained'."
        )

    if mode == "clusters":
        adata_sc = adata_to_cluster_expression(
            adata_sc, cluster_label, scale, add_density=True
        )

    # Check if training_genes key exist/is valid in adatas.uns
    if not set(["training_genes", "overlap_genes"]).issubset(set(adata_sc.uns.keys())):
        raise ValueError("Missing tangram parameters. Run `pp_adatas()`.")

    if not set(["training_genes", "overlap_genes"]).issubset(set(adata_sp.uns.keys())):
        raise ValueError("Missing tangram parameters. Run `pp_adatas()`.")

    assert list(adata_sp.uns["training_genes"]) == list(adata_sc.uns["training_genes"])

    # get training_genes
    if cv_train_genes is None:
        training_genes = adata_sc.uns["training_genes"]
    elif cv_train_genes is not None:
        if set(cv_train_genes).issubset(set(adata_sc.uns["training_genes"])):
            training_genes = cv_train_genes
        else:
            raise ValueError(
                "Given training genes list should be subset of two AnnDatas."
            )

    logging.info("Allocate tensors for mapping.")
    # Allocate tensors (AnnData matrix can be sparse or not)

    if isinstance(adata_sc.X, csc_matrix) or isinstance(adata_sc.X, csr_matrix):
        S = np.array(adata_sc[:, training_genes].X.toarray(), dtype="float32",)
    elif isinstance(adata_sc.X, np.ndarray):
        S = np.array(adata_sc[:, training_genes].X.toarray(), dtype="float32",)
    else:
        X_type = type(adata_sc.X)
        logging.error("AnnData X has unrecognized type: {}".format(X_type))
        raise NotImplementedError

    if isinstance(adata_sp.X, csc_matrix) or isinstance(adata_sp.X, csr_matrix):
        G = np.array(adata_sp[:, training_genes].X.toarray(), dtype="float32")
    elif isinstance(adata_sp.X, np.ndarray):
        G = np.array(adata_sp[:, training_genes].X, dtype="float32")
    else:
        X_type = type(adata_sp.X)
        logging.error("AnnData X has unrecognized type: {}".format(X_type))
        raise NotImplementedError

    if not S.any(axis=0).all() or not G.any(axis=0).all():
        raise ValueError("Genes with all zero values detected. Run `pp_adatas()`.")

    d_source = None

    # define density_prior if 'rna_count_based' is passed to the density_prior argument:
    d_str = density_prior
    if type(density_prior) is np.ndarray:
        d_str = "customized"

    if density_prior == "rna_count_based":
        density_prior = adata_sp.obs["rna_count_based_density"]

    # define density_prior if 'uniform' is passed to the density_prior argument:
    elif density_prior == "uniform":
        density_prior = adata_sp.obs["uniform_density"]

    if mode == "cells":
        d = density_prior

    if mode == "clusters":
        d_source = np.array(adata_sc.obs["cluster_density"])

    if mode in ["clusters", "constrained"]:
        if density_prior is None:
            d = adata_sp.obs["uniform_density"]
            d_str = "uniform"
        else:
            d = density_prior
        if lambda_d is None or lambda_d == 0:
            lambda_d = 1

    # Choose device
    device = torch.device(device)  # for gpu

    if verbose:
        print_each = 100
    else:
        print_each = None

    if mode in ["cells", "clusters"]:
        hyperparameters = {
            "lambda_d": lambda_d,  # KL (ie density) term
            "lambda_g1": lambda_g1,  # gene-voxel cos sim
            "lambda_g2": lambda_g2,  # voxel-gene cos sim
            "lambda_r": lambda_r,  # regularizer: penalize entropy
            "d_source": d_source,
        }

        logging.info(
            "Begin training with {} genes and {} density_prior in {} mode...".format(
                len(training_genes), d_str, mode
            )
        )
        mapper = mo.Mapper(
            S=S, G=G, d=d, device=device, random_state=random_state, **hyperparameters,
        )

        # TODO `train` should return the loss function

        mapping_matrix, training_history = mapper.train(
            learning_rate=learning_rate, num_epochs=num_epochs, print_each=print_each,
        )

    # constrained mode
    elif mode == "constrained":
        hyperparameters = {
            "lambda_d": lambda_d,  # KL (ie density) term
            "lambda_g1": lambda_g1,  # gene-voxel cos sim
            "lambda_g2": lambda_g2,  # voxel-gene cos sim
            "lambda_r": lambda_r,  # regularizer: penalize entropy
            "lambda_count": lambda_count,
            "lambda_f_reg": lambda_f_reg,
            "target_count": target_count,
        }

        logging.info(
            "Begin training with {} genes and {} density_prior in {} mode...".format(
                len(training_genes), d_str, mode
            )
        )

        mapper = mo.MapperConstrained(
            S=S, G=G, d=d, device=device, random_state=random_state, **hyperparameters,
        )

        mapping_matrix, F_out, training_history = mapper.train(
            learning_rate=learning_rate, num_epochs=num_epochs, print_each=print_each,
        )

    logging.info("Saving results..")
    adata_map = sc.AnnData(
        X=mapping_matrix,
        obs=adata_sc[:, training_genes].obs.copy(),
        var=adata_sp[:, training_genes].obs.copy(),
    )

    if mode == "constrained":
        adata_map.obs["F_out"] = F_out

    # Annotate cosine similarity of each training gene
    G_predicted = adata_map.X.T @ S
    cos_sims = []
    for v1, v2 in zip(G.T, G_predicted.T):
        norm_sq = np.linalg.norm(v1) * np.linalg.norm(v2)
        cos_sims.append((v1 @ v2) / norm_sq)

    df_cs = pd.DataFrame(cos_sims, training_genes, columns=["train_score"])
    df_cs = df_cs.sort_values(by="train_score", ascending=False)
    adata_map.uns["train_genes_df"] = df_cs

    # Annotate sparsity of each training genes
    ut.annotate_gene_sparsity(adata_sc)
    ut.annotate_gene_sparsity(adata_sp)
    adata_map.uns["train_genes_df"]["sparsity_sc"] = adata_sc[
        :, training_genes
    ].var.sparsity
    adata_map.uns["train_genes_df"]["sparsity_sp"] = adata_sp[
        :, training_genes
    ].var.sparsity
    adata_map.uns["train_genes_df"]["sparsity_diff"] = (
        adata_sp[:, training_genes].var.sparsity
        - adata_sc[:, training_genes].var.sparsity
    )

    adata_map.uns["training_history"] = training_history

    return adata_map

def generate_spatial_map(
    X,
    mapping_matrix,
    v,
    method = "MAP"
):
    """
    Map single cell data (`adata_sc`) on spatial data (`adata_sp`).
    
    Args:
        X (ndarray): single cell data
        mapping_matrix (ndarray): cell to location mapping matrix
        v (ndarray): spatial coords of the given spots
        method (str): Optional. The method to map the cells, default is max a posteriori estimation.
        
    Returns:
        a cell by coords array indicating the inferred spatial location of cells
    """

    if method == "MAP":
        max_col_indices = np.argmax(mapping_matrix, axis=1)
        v_pred = v[max_col_indices]
        M_pred = np.matmul(mapping_matrix.transpose(), X)
    else:
        #mapping_matrix[mapping_matrix < 0.1] = 0
        v_pred = np.matmul(mapping_matrix,v)
        M_pred = np.matmul(mapping_matrix.transpose(), X)

    
    return v_pred,M_pred

def Benchmark_temso(
    X,
    Y,
    mapping_matrix,
    v_gt,
    v,
    d_gt,
    method = "MAP"
):
    """
    Map single cell data (`adata_sc`) on spatial data (`adata_sp`).
    
    Args:
        X (ndarray): single cell data
        Y (ndarray): spatial transcriptomic data
        mapping_matrix (ndarray): cell to location mapping matrix
        v_gt (ndarray): spatial coords of the cells
        v (ndarray): spatial coords of the given spots
        method (str): Optional. The method to map the cells, default is max a posteriori estimation.
        
    Returns:
        MSE_st (int): MSE between the inferred spatial transcriptomics and the ground truth
        MSE_coords (int): MSE between the inferred spatial coordinates and the ground truth
        N_mismatches (int): The number of inconsistent spatial coordinates between the inferred coordinates and the ground truth
        j_index (int): The jaccard index between the knn graphs of cells' true coordinates and inferred coordinates(k = 20)
    """

    v_pred,M_pred = generate_spatial_map(X,mapping_matrix,v,method)

    MSE_st = ((M_pred - Y)**2).mean()
    MSE_coords = ((v_pred - v_gt)**2).mean()
    N_mismatches = np.sum(v_pred != v_gt)
    L1_coords = abs(v_pred - v_gt).mean()

    pairwise_gt = pairwise_distances(v_gt, v_gt, metric='euclidean')
    pairwise_pred = pairwise_distances(v_pred,v_pred, metric="euclidean")
 
    knn_infer = knn_graph(pairwise_distances(v_pred, v_pred, metric='euclidean'),20)
    knn_true = knn_graph(pairwise_distances(v_gt,v_gt, metric="euclidean"),20)
 
    j_sum = 0
    for node in knn_true.nodes():
        j_sum += jaccard_similarity(knn_true.edges(node),knn_infer.edges(node))
 
    j_index_individual = j_sum/knn_true.number_of_nodes()
    j_index_whole = jaccard_similarity(knn_true.edges(),knn_infer.edges())

    pearson = pearsonr(pairwise_gt.flatten(),pairwise_pred.flatten()).statistic
    if d_gt is not None:
        a_squared = np.sum(v_pred**2, axis=1,keepdims= True)  # N x 1
        b_squared = np.sum(v**2, axis=1,keepdims= True).T  # 1 x M
        distances = a_squared + b_squared - 2 * np.dot(v_pred, v.T)  # N x M
   
        nearest_neighbors = np.argmin(distances, axis=1)
        density_dist = np.bincount(nearest_neighbors, minlength=v.shape[0]).astype(float)
        MSE_density = ((density_dist - d_gt)**2).mean()
    else:
        MSE_density = 0  

    
   
    return MSE_st,MSE_coords,N_mismatches,j_index_individual,j_index_whole,L1_coords,pearson,MSE_density


def define_clones(spatial_map,nclones = 4):
    ncells = spatial_map.shape[0]
    clone_size = ncells/nclones

    spatial_map['clonalid'] = spatial_map['cellid']
    for i in range(nclones):
        max_id = (i+1) * clone_size
        min_id = i * clone_size + 1
        spatial_map.loc[(spatial_map.index >= min_id) & (spatial_map['cellid'] <= max_id), 'clonalid'] = i

    return spatial_map

def knn_graph(pairwise_distances, k):
    # Construct KNN graph
    knn_graph = kneighbors_graph(pairwise_distances, n_neighbors=k, mode='connectivity', include_self=False)
    
    # Convert sparse matrix to networkx graph
    knn_graph_networkx = nx.from_scipy_sparse_array(knn_graph)

    return knn_graph_networkx

def jaccard_similarity(g, h):
    i = set(g).intersection(h)
    return round(len(i) / (len(g) + len(h) - len(i)),3)


def kernel_product(x,y, mode = "gaussian", s = 1.):
    x_i = x.unsqueeze(1)
    y_j = y.unsqueeze(0)
    xmy = ((x_i-y_j)**2).sum(2)

    if   mode == "gaussian" : K = torch.exp( - xmy/s**2)
    elif mode == "laplace"  : K = torch.exp( - torch.sqrt(xmy + (s**2)))
    elif mode == "energy"   : K = torch.pow(   xmy + (s**2), -.25 )

    return torch.t(K)


def reroot_balance(tree):
    ncells = tree.root.count_terminals()
    node = tree.root
    while node.count_terminals() > ncells/2:
        ncells_bound = 0
        for children_node in node.clades:
            ncells_child = children_node.count_terminals()
            if ncells_child >= ncells_bound:
                node = children_node
                ncells_bound = ncells_child
        
    
    tree.root_with_outgroup(node)
    
    new_clades = tree.root.clades
    reg = 0
    for i in range(len(new_clades)):
    #     # Take the last two children and create a new internal node
        children_node = new_clades[i]
        num_leaves = children_node.count_terminals()
        if num_leaves > reg:
            reg = num_leaves
            biggest_clade_index = i
            biggest_clade_length = children_node.branch_length

    new_clade_children_index = [i for i, x in enumerate(new_clades) if i != biggest_clade_index]
    new_clade_children = [x for i, x in enumerate(new_clades) if i != biggest_clade_index]
    for clade_delete in reversed(new_clade_children_index):
        new_clades.pop(clade_delete)
    new_clade = Clade(clades=new_clade_children,branch_length = biggest_clade_length)
        
    # Add this new internal node as a child of the original clade
    new_clades.append(new_clade)

    return tree


def Mask_clone(cloneid):
    n = cloneid.shape[0]
    Z_clone = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if cloneid[i] == cloneid[j]:
                Z_clone[i, j] = 1
    return Z_clone

def Define_clones(tree,num_clones = 4):
    if not tree.is_bifurcating():
        raise ValueError("Input tree is not bifurcating.")
    
    if len(tree.root.clades)!=2:
        print("Input tree root has more than 2 children. Will run reroot and binarize the tree.")
        tree = reroot_balance(tree)
    tree = reroot_balance(tree)

    if num_clones%2 != 0:
        raise ValueError("Input number of clones should be a power of 2.")
    

    Num_generation = np.log2(num_clones)
    node_list = [tree.root]
    for i in range(int(Num_generation)):
        old_node_list = node_list
        node_list = []
        for node in old_node_list:
            children = node.clades
            node_list.extend(children)
    
    cloneid = np.array([0] * tree.root.count_terminals())
    for i in range(num_clones):
        clone_leaf_list = [int(re.findall(r'\d+',x.name)[0])-1 for x in node_list[i].get_terminals()]
        cloneid[clone_leaf_list] = i
    
    Z_clone = Mask_clone(cloneid)
    
    return cloneid,Z_clone