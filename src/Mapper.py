import numpy as np
import logging
import torch
from torch.nn.functional import softmax, cosine_similarity, hardshrink,kl_div
from torch.nn import PairwiseDistance

def compute_nearest_neighbors(source, target,k=5):
    """
    Compute nearest neighbors between source and target datasets.
    """
    # Compute pairwise distances between source and target points
    distances = torch.cdist(source, target)  # Pairwise L2 distance
    # Get the indices of the k-nearest neighbors for each source point in the target
    source_to_target_knn = distances.topk(k, dim=1, largest=False).indices
    # Get the indices of the k-nearest neighbors for each target point in the source
    target_to_source_knn = distances.topk(k, dim=0, largest=False).indices.t()  # Transpose for consistency
    
    
    return source_to_target_knn, target_to_source_knn

def mutual_nearest_neighbors(source, target,k = 5):
    """
    Identify mutual nearest neighbors between source and target.
    """
    # Get k-nearest neighbors for source -> target and target -> source
    source_to_target_knn, target_to_source_knn = compute_nearest_neighbors(source, target, k)

    # Mutual k-nearest neighbors: points where source and target belong to each other's k-nearest neighbors
    mutual_knn_indices = []
    for i in range(source.shape[0]):
        for j in source_to_target_knn[i]:  # Iterate over the k nearest neighbors of source[i]
            if i in target_to_source_knn[j]:  # Check if i is also in the k-nearest neighbors of target[j]
                mutual_knn_indices.append((i, j.item()))  # Collect mutual knn pair
    
    return mutual_knn_indices

def initialize_mapping_matrix(source, target, mutual_nn_indices):
    """
    Initialize mapping matrix using mutual nearest neighbors.
    """
    if len(mutual_nn_indices) == 0:
        raise ValueError("No mutual nearest neighbors found.")

    n_samples_source = source.shape[0]
    n_samples_target = target.shape[0]

    # Initialize the mapping matrix with zeros
    W = torch.zeros(n_samples_source, n_samples_target)

    # Set W[i, j] = 1 for each mutual nearest neighbor pair (i, j)
    for i, j in mutual_nn_indices:
        W[i, j] = 1

    return W

def bin_tensor(input_tensor, num_bins, bin_size):
    # Ensure the input tensor is a 1D tensor
    if input_tensor.ndim != 1:
        raise ValueError("input_tensor must be a 1D tensor")

    # Determine the range of values in the tensor
    min_val = input_tensor.min().item()
    max_val = input_tensor.max().item()

    # Create bin edges based on the bin size
    bin_edges = torch.arange(min_val, max_val + bin_size, bin_size)
    
    # Compute the bin indices for each value in the tensor
    bin_indices = torch.bucketize(input_tensor, bin_edges, right=False)
    
    # Initialize a tensor to hold the binned counts
    binned_tensor = torch.zeros(num_bins, dtype=torch.float32)

    # Count occurrences in each bin
    for i in range(num_bins):
        binned_tensor[i] = (bin_indices == i).sum().float()
    
    return binned_tensor

def compute_binned_distribution(M,v,num_bins,device="cpu"):
    x = v[:,0]
    y = v[:,1]

    min_x = x.min().item()
    max_x = x.max().item()
    min_y = y.min().item()
    max_y = y.max().item()

    bin_size_x = (max_x - min_x)/num_bins
    bin_size_y = (max_y - min_y)/num_bins

    M_binned_x = torch.zeros(M.shape[0], num_bins+1,device=device)
    M_binned_y = torch.zeros(M.shape[0], num_bins+1,device=device)

    bin_edges_x = torch.arange(min_x, max_x + bin_size_x, bin_size_x,device=device)
    bin_edges_y = torch.arange(min_y, max_y + bin_size_y, bin_size_y,device=device)

    bin_indices_x = torch.bucketize(v[:,0], bin_edges_x, right=False)
    bin_indices_y = torch.bucketize(v[:,1], bin_edges_y, right=False)

    for j in range(M.shape[1]):
        prob = M[:,j]
        M_binned_x[:,bin_indices_x[j]] = M_binned_x[:,bin_indices_x[j]] + prob
        M_binned_y[:,bin_indices_y[j]] = M_binned_y[:,bin_indices_y[j]] + prob
    
    return M_binned_x,M_binned_y



def compute_pairwise_distances(x, y):
    x_norm = (x**2).sum(1).view(-1, 1)
    y_t = torch.transpose(y, 0, 1)
    y_norm = (y**2).sum(1).view(1, -1)
    
    dist = x_norm + y_norm - 2.0 * torch.mm(x, y_t)
    return torch.sqrt(torch.clamp(dist, 0.0, np.inf))


def sample_above_threshold(tensor, threshold, v,device="cpu"):
    # Check that the tensor is 2D
    if tensor.dim() != 2:
        raise ValueError("Input tensor must be 2D.")
    
    sampled_values_x = []
    sampled_values_y = []
    
    for row in tensor:
        # Get values above the threshold
        indices_above_threshold = torch.nonzero(row > threshold).flatten()
        
        if indices_above_threshold.numel() > 0:  # Ensure there are values above the threshold
            # Randomly sample one value from the filtered values
            sampled_value_x = v[torch.randint(0, indices_above_threshold.numel(), (1,)).item(),0]
            sampled_value_y = v[torch.randint(0, indices_above_threshold.numel(), (1,)).item(),1]
            sampled_values_x.append(sampled_value_x)
            sampled_values_y.append(sampled_value_y)
        else:
            # If no values are above the threshold, you can choose how to handle this
            sampled_values_x.append(v[torch.randint(0,v.shape[0],(1,)),0])  # or you could use torch.tensor(float('nan')) or a specific fallback value
            sampled_values_y.append(v[torch.randint(0,v.shape[0],(1,)),1])
    
    return torch.tensor(sampled_values_x,device = device),torch.tensor(sampled_values_y,device = device)

class UnimodalLoss(torch.nn.Module):
    def __init__(self):
        super(UnimodalLoss, self).__init__()
    
    def forward(self, distribution):
        # distribution: Tensor of shape (batch_size, num_classes)
        
        # Compute the differences between consecutive elements
        non_empty_mask = distribution.abs().sum(dim=0).bool()
        diffs = torch.diff(distribution[:,non_empty_mask], dim=1)
        
        # Identify where the distribution starts decreasing
        decreasing = diffs < 0
        
        # Identify violations: an increase after a decrease
        violation_mask = decreasing[:, :-1] & (diffs[:, 1:] > 0)

        # Sum the penalty and average over the batch
        unimodal_penalty = violation_mask.float().sum(dim=1).mean()
        
        return unimodal_penalty

class GaussianDensityLoss(torch.nn.Module):
    def __init__(self, target_variance=None, device = "cpu"):
        super(GaussianDensityLoss, self).__init__()
        self.target_variance = target_variance if target_variance is not None else 0.1
        self.device = device
    
    def forward(self, distribution, weights=None,target_means = None):
        # distribution: Tensor of shape (number_of_locs, 2), each row is (x, y)
        # weights: Tensor of shape (sample_size,number_of_locs), representing the density weights (optional)

        if weights is None:
            weights = torch.ones(distribution.size(0))

        # Separate x and y coordinates
        x = distribution[:, 0]
        y = distribution[:, 1]

        # Compute weighted means
        mean_x = torch.sum(weights * x,dim=1)
        mean_y = torch.sum(weights * y,dim=1)

        # Compute target means
        target_x = x[torch.argmax(weights,dim=1)]
        target_y = y[torch.argmax(weights,dim=1)]

        #calculate weighted mean of the top 5 
        weights_sorted, sorted_indices  = torch.sort(weights, dim=1, descending=True)
        top5_weights = weights_sorted[:, :10]
        top_indices = sorted_indices[:,:10]
        indices_choice = top_indices[torch.arange(top_indices.shape[0]),torch.randint(0,10,(top_indices.shape[0],))]
        target_x = x[indices_choice]
        target_y = y[indices_choice]
        #target_x = torch.sum(top5_weights * x[sorted_indices[:,:10]],dim=1)
        #target_y = torch.sum(top5_weights * y[sorted_indices[:,:10]],dim=1)
        #target_x,target_y = sample_above_threshold(weights,0.1,distribution,device=self.device)

        #target_x = target_means[:,0]
        #target_y = target_means[:,1]

        # Compute weighted variances
        var_x = torch.sum(weights * (x[:,None] - mean_x.repeat(x.shape[0],1)).t()**2,dim=1)
        var_y = torch.sum(weights * (y[:,None] - mean_y.repeat(y.shape[0],1)).t()**2,dim=1)

        # Compute weighted covariance (if needed for a more general Gaussian)
        # cov_xy = torch.sum(weights * (x - mean_x) * (y - mean_y))

        # Gaussian likelihood loss
        gaussian_likelihood_loss = (
            ((mean_x - target_x)**2 + (mean_y - target_y)**2) / (2 * self.target_variance)
        ).mean()

        # Regularize variance to encourage small variance
        variance_loss = ((var_x + var_y)/2).mean()
        
        # Combine losses
        loss = gaussian_likelihood_loss + variance_loss

        return loss
    

class PeakLoss(torch.nn.Module):
    def __init__(self, top_n = 5,device = "cpu"):
        super(PeakLoss, self).__init__()
        self.device = device
        self.top_n = top_n
    
    def forward(self,distribution, weights=None,spot_dist = None):
        # distribution: Tensor of shape (number_of_locs, 2), each row is (x, y)
        # weights: Tensor of shape (sample_size,number_of_locs), representing the density weights (optional)

        # Separate x and y coordinates
        x = distribution[:, 0]
        y = distribution[:, 1]

        # Compute weighted means
        mean_x = torch.sum(weights * x,dim=1)
        mean_y = torch.sum(weights * y,dim=1)

        # Compute weighted variances
        var_x = torch.sum(weights * (x[:,None] - mean_x.repeat(x.shape[0],1)).t()**2,dim=1)
        var_y = torch.sum(weights * (y[:,None] - mean_y.repeat(y.shape[0],1)).t()**2,dim=1)

        # Regularize variance to encourage small variance
        variance_loss = ((var_x + var_y)/2).mean()
        

        #calculate weighted mean of the top n
        weights_sorted, sorted_indices  = torch.sort(weights, dim=1, descending=True)
        top5_weights_relative = torch.div(weights_sorted[:, :self.top_n].t(),weights_sorted[:,0])
        top5_weights_relative = top5_weights_relative.t()
        #sorted_indices[top5_weights_relative < 0.1] = -1
        top_indices = sorted_indices[:,:self.top_n]

        loss = 0

        for i in range(self.top_n):
            for j in range(i+1,self.top_n):
                top_indices_i = top_indices[:,i]
                top_indices_j = top_indices[:,j]

                if_include = torch.logical_and(top5_weights_relative[:,i] >= 0.1, top5_weights_relative[:,j] >= 0.1)
                loss += torch.sum(spot_dist[top_indices_i,top_indices_j] * if_include)

        return variance_loss

# def kl_divergence(p, q):
#     """Calculate the Kullback-Leibler divergence D_KL(P || Q)"""
#     # Avoid log(0) by adding a small epsilon

#     return torch.sum(p * torch.log(p / q), dim=-1)

# def jensen_shannon_divergence(p, q):
#     """Calculate the Jensen-Shannon divergence between two probability distributions"""
#     # Ensure p and q are probability distributions
#     m = 0.5 * (p + q)
#     return 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)

def pairwise_jsd(prob_matrix):
    """Calculate the pairwise Jensen-Shannon divergence distance matrix"""
    n = prob_matrix.shape[0]
    # Add a small epsilon to avoid log(0)
    epsilon = 1e-12
    
    # Broadcast the probability matrix to compute pairwise JSD in a vectorized manner
    prob_matrix_expanded = prob_matrix.unsqueeze(0)  # Shape: (1, N, D)
    prob_matrix_tiled = prob_matrix.unsqueeze(1)     # Shape: (N, 1, D)
    
    # Pairwise average distribution M = 0.5 * (P + Q)
    m = 0.5 * (prob_matrix_expanded + prob_matrix_tiled)  # Shape: (N, N, D)
    
    # Compute KL divergences in a batch
    p_kl_m = torch.sum(prob_matrix_expanded * torch.log((prob_matrix_expanded + epsilon) / (m + epsilon)), dim=-1)
    q_kl_m = torch.sum(prob_matrix_tiled * torch.log((prob_matrix_tiled + epsilon) / (m + epsilon)), dim=-1)
    
    # Jensen-Shannon divergence
    jsd_matrix = 0.5 * (p_kl_m + q_kl_m)
    
    return jsd_matrix

def cell_uniformity_loss(coords, bandwidth=0.5, epsilon=1e-8):
    """
    Loss function based on KL divergence to encourage a uniform distribution of 2D coordinates.
    
    Args:
        coords (torch.Tensor): A tensor of shape (N, 2) where N is the number of 2D points.
        bandwidth (int): The bandwidth for gaussian kernel estimation.
        epsilon (float): A small value to prevent division by zero or log(0).
        
    Returns:
        torch.Tensor: A scalar loss value encouraging a uniform distribution.
    """
    pairwise_distances = torch.cdist(coords,coords) 
    N = pairwise_distances.size(0)
    
    # Step 2: Apply Gaussian kernel to estimate KDE
    kde = torch.exp(-pairwise_distances**2 / (2 * bandwidth**2))
    
    # Step 3: Normalize KDE to form an empirical probability distribution
    empirical_dist = kde.sum(dim=1)  # Sum over neighbors to get density at each point
    empirical_dist = empirical_dist / empirical_dist.sum()  # Normalize to make it a distribution
    
    # Step 4: Target uniform distribution (1/N probability for each point)
    uniform_dist = torch.ones_like(empirical_dist) / N
    
    # Step 5: Add epsilon to avoid log(0) in KL divergence
    empirical_dist = empirical_dist + epsilon
    uniform_dist = uniform_dist + epsilon
    
    # Step 6: Compute KL divergence between empirical and uniform distribution
    kl_loss = kl_div(empirical_dist.log(), uniform_dist, reduction='batchmean')

    #penalizes small distances
    mask = ~torch.eye(N, dtype=bool, device=coords.device)
    pairwise_distances = pairwise_distances[mask]
    
    # Add a small epsilon to avoid division by zero
    epsilon = 1e-6
    inv_distances = 1.0 / (pairwise_distances + epsilon)

    # The loss increases when points are too close together
    dist_loss = torch.mean(inv_distances)
    
    return kl_loss + dist_loss



def Density_loss_new(coords,coords_spots,d,epsilon=1e-8):
    

    #New density loss
    # Step 1: Compute the pairwise squared Euclidean distances between a and b
    a_squared = torch.sum(coords**2, axis =1).unsqueeze(1)  # N x 1
    b_squared = torch.sum(coords_spots**2, axis =1).unsqueeze(0)  # 1 x M
    distances = a_squared + b_squared - 2 * torch.matmul(coords, coords_spots.T)  # N x M
    
    # Step 2: For each point in a, find the nearest point in b (minimum distance)
    nearest_neighbors = torch.argmin(distances, dim=1)  # N x 1, index of nearest point in b for each point in a
    
    # Step 3: Compute density distribution by counting the nearest neighbors for each point in b
    density_dist = torch.bincount(nearest_neighbors, minlength=coords_spots.size(0)).float()  # M x 1
    density_dist = density_dist + epsilon
    density_loss = kl_div(density_dist.log(), d, reduction='batchmean')

    return density_loss

class Mapper(torch.nn.Module):
    """
    The optimizer for TemSOMap.
    After initialization, the optimizer is run with the 'train' method.
    """

    def __init__(
        self,
        X,
        Y,
        Z,
        Z_clone,
        v,
        v_gt=None,
        d=None,
        d_source=None,
        lambda_g1=1.0,
        lambda_d=0,
        lambda_g2=0,
        lambda_r=0,
        lambda_l=0,
        lambda_c=0,
        lambda_u=0,
        lambda_cu=0,
        device="cpu",
        adata_map=None,
        random_state=None,
    ):
        """
        Instantiate the TemSOMap optimizer.

        Args:
            X (ndarray): Single nuclei matrix, shape = (number_cell, number_genes).
            Y (ndarray): Spatial transcriptomics matrix, shape = (number_spots, number_genes).
                Spots can be single cells or they can contain multiple cells.
            Z (ndarray): Pairwise lineage barcode distances of cells, shape = (number_cell, number_cell)
            Z_clone (ndarray): Clonal information of cells, shape = (number_cell, number_cell)
            v (ndarray):coordinate info of voxels, shape = (number_spots, 2)
            d (ndarray): Spatial density of cells, shape = (number_spots,). If not provided, the density term is ignored.
                This array should satisfy the constraints d.sum() == 1.
            d_source (ndarray): Density of single cells in single cell clusters. To be used when S corresponds to cluster-level expression.
                This array should satisfy the constraint d_source.sum() == 1.
            lambda_g1 (float): Optional. Strength of Tangram loss function. Default is 1.
            lambda_d (float): Optional. Strength of density regularizer. Default is 0.
            lambda_g2 (float): Optional. Strength of voxel-gene regularizer. Default is 0.
            lambda_r (float): Optional. Strength of entropy regularizer. An higher entropy promotes
                              probabilities of each cell peaked over a narrow portion of space.
                              lambda_r = 0 corresponds to no entropy regularizer. Default is 0.
            device (str or torch.device): Optional. Device is 'cpu'.
            adata_map (scanpy.AnnData): Optional. Mapping initial condition (for resuming previous mappings). Default is None.
            random_state (int): Optional. pass an int to reproduce training. Default is None.
        """
        super().__init__()
        self.X = torch.tensor(X, device=device, dtype=torch.float32)
        self.Y = torch.tensor(Y, device=device, dtype=torch.float32)
        self.Z = torch.tensor(Z, device=device, dtype=torch.float32)

        #self.Z = torch.exp(-self.Z/1**2)

        self.clonal_loss_enabled = Z_clone is not None
        if self.clonal_loss_enabled:
            self.Z_clone = torch.tensor(Z_clone, device=device, dtype=torch.float32)
        self.v = torch.tensor(v, device=device, dtype=torch.float32)

        #self.v_gt = torch.tensor(v_gt,device=device, dtype=torch.float32)

        self.target_density_enabled = d is not None
        if self.target_density_enabled:
            self.d = torch.tensor(d, device=device, dtype=torch.float32)

        self.source_density_enabled = d_source is not None
        if self.source_density_enabled:
            self.d_source = torch.tensor(d_source, device=device, dtype=torch.float32)

        self.spot_dist = torch.cdist(self.v,self.v)

        self.lambda_d = lambda_d
        self.lambda_g1 = lambda_g1
        self.lambda_g2 = lambda_g2
        self.lambda_r = lambda_r
        self.lambda_l = lambda_l
        self.lambda_c = lambda_c
        self.lambda_u = lambda_u
        self.lambda_cu = lambda_cu
        self._density_criterion = torch.nn.KLDivLoss(reduction="batchmean")
        self._density_criterion_lineage = torch.nn.KLDivLoss(reduction="batchmean")
        #self.Unimodal_criterion = GaussianDensityLoss(device=device)
        self.Unimodal_criterion = PeakLoss(top_n = 20)

        self.random_state = random_state

        if adata_map is None:
            if self.random_state:
                np.random.seed(seed=self.random_state)
            self.M = np.random.normal(0, 1, (X.shape[0], Y.shape[0]))
        else:
            raise NotImplemented
            self.M = adata_map.X  # doesn't work. maybe apply inverse softmax

        # self.M = torch.tensor(
        #     self.M, device=device, requires_grad=True, dtype=torch.float32
        # )

        self.M = torch.nn.Parameter(torch.tensor(np.random.normal(0, 1, (X.shape[0], Y.shape[0])), device=device , dtype=torch.float32))
        # Find mutual nearest neighbors
        #mutual_nn_indices = mutual_nearest_neighbors(self.X, self.Y,k=10)

        # Initialize mapping matrix
        #self.M = torch.nn.Parameter(torch.tensor(initialize_mapping_matrix(self.X, self.Y, mutual_nn_indices), device=device , dtype=torch.float32))

    def _loss_fn(self,device="cpu", verbose=True):
        """
        Evaluates the loss function.

        Args:
            verbose (bool): Optional. Whether to print the loss results. If True, the loss for each individual term is printed as:
                density_term, gene-voxel similarity term, voxel-gene similarity term. Default is True.

        Returns:
            Tuple of 6 Floats: Total loss, gv_loss, vg_loss, kl_reg, entropy_reg, lineage_reg
        """
        M_probs = softmax(self.M, dim=1)
        #print(f"printing M_probs {M_probs}.")

        # 4. density term
        if self.target_density_enabled and self.source_density_enabled:
            d_pred = torch.log(
                self.d_source @ M_probs
            )  # KL wants the log in first argument
            density_term = self.lambda_d * self._density_criterion(d_pred, self.d)
        elif self.target_density_enabled and not self.source_density_enabled:
            d_pred = torch.log(
                M_probs.sum(axis=0) / self.M.shape[0]
            )  # KL wants the log in first argument
            density_term = self.lambda_d * self._density_criterion(d_pred, self.d)
        else:
            density_term = 0

        # 1. expression term
        G_pred = torch.matmul(M_probs.t(), self.X)
        gv_term = cosine_similarity(G_pred, self.Y, dim=0).mean()
        vg_term = cosine_similarity(G_pred, self.Y, dim=1).mean()
        expression_term = self.lambda_g1 * gv_term + self.lambda_g2 * vg_term

        # 2. entropy regularization term
        # perform Gaussian smoothing to get smoothed M_probs
        sigma = 0.1
        gaussian_kernel = torch.exp(-self.spot_dist ** 2 / (2 * sigma ** 2))
        gaussian_kernel /= gaussian_kernel.sum(dim=1, keepdim=True)
        M_probs_smoothed = torch.matmul(M_probs, gaussian_kernel)
        M_probs_smoothed /= M_probs_smoothed.sum(dim=1, keepdim=True)
        regularizer_term = self.lambda_r * (torch.log(M_probs_smoothed) * M_probs_smoothed).mean()

        #new implementation of the regularizer_term
        #M_binned_x,M_binned_y = compute_binned_distribution(M_probs,self.v,100,device=device)
        #regularizer_bin_term =  self.lambda_u * self.Unimodal_criterion(M_binned_x) + self.lambda_u * self.Unimodal_criterion(M_binned_y)
        if self.lambda_u != 0:
            regularizer_bin_term = self.lambda_u * self.Unimodal_criterion(self.v,M_probs,self.spot_dist)
        else:
            regularizer_bin_term = 0

        # 3. new lineage term
        self.v_pred = M_probs_smoothed @ self.v # torch.matmul(M_probs, self.v)
        #euclidean distance
        self.Z_pred = torch.cdist(self.v_pred,self.v_pred) 
        #Alternative ways to calculate pairwise distance
        #self.Z_pred = torch.cdist(M_binned_x,M_binned_x) + torch.cdist(M_binned_y,M_binned_y)
        # Pairwise JSD
        #self.Z_pred = pairwise_jsd(M_probs_smoothed)

        # #RBF kernel, Gaussian
        s = 1
        #self.Z_pred = torch.exp(-self.Z_pred/s**2)

        # #RBF kernel, Laplacian
        # s = 1
        # self.Z_pred = torch.exp(-torch.sqrt(self.Z_pred + s**2))

        # #RBF kernel, energy
        # s = 1
        # self.Z_pred = torch.pow(self.Z_pred + s**2,-.25)

        # lineage_term = self.lambda_l * self._density_criterion_lineage(self.Z_pred / (torch.sum(self.Z_pred, dim=1, keepdim = True) + 1e-10), 
        #                                       self.Z / (torch.sum(self.Z, dim=1, keepdim = True) + 1e-10))

        #MSE instead of KL
        P_z = self.Z / (torch.sum(self.Z, dim=1, keepdim = True) + 1e-6)
        P_zpred = self.Z_pred / (torch.sum(self.Z_pred, dim=1, keepdim = True) + 1e-6)        
        lineage_term = self.lambda_l *((P_z - P_zpred)**2).sum()

        # 4. lineage clustering (clone) regularization term
        if self.clonal_loss_enabled:
            clone_term = self.lambda_c * (torch.sum(self.Z_clone * P_zpred) - torch.sum((1 - self.Z_clone) * P_zpred))/(torch.sum(P_zpred))
        else:
            clone_term = 0

        # 5. cell uniformity regularization term
        if self.lambda_cu != 0:
            uniformity_term = self.lambda_cu * cell_uniformity_loss(self.v_pred,bandwidth = 0.2)
        else:
            uniformity_term = 0

        # 6. new density loss
        if self.lambda_d != 0:
            density_term = self.lambda_d * Density_loss_new(self.v_pred,self.v,self.d)
        else:
            density_term = 0

        main_loss = (gv_term / self.lambda_g1).tolist()
        kl_reg = (
            (density_term / self.lambda_d).tolist()
            if density_term != 0
            else np.nan
        )
        vg_reg = (vg_term / self.lambda_g2).tolist()

        entropy_reg = (regularizer_term / self.lambda_r).tolist()

        lineage_reg = (lineage_term / self.lambda_l).tolist()

        clone_reg = (
            (clone_term / self.lambda_c).tolist()
            if clone_term != 0
            else np.nan)
        
        bin_reg = (
            (regularizer_bin_term / self.lambda_u).tolist()
            if regularizer_bin_term != 0
            else np.nan)
        
        uniformity_reg = (
            (uniformity_term / self.lambda_cu).tolist()
            if uniformity_term != 0
            else np.nan)
        
        if verbose:

            term_numbers = [main_loss, vg_reg, kl_reg, entropy_reg, lineage_reg, clone_reg,bin_reg,uniformity_reg]
            term_names = ["Score", "VG reg", "KL reg", "Entropy reg", "Lineage reg", "Clone reg", "Unimodal reg", "Uniformity reg"]

            d = dict(zip(term_names, term_numbers))
            clean_dict = {k: d[k] for k in d if not np.isnan(d[k])}
            msg = []
            for k in clean_dict:
                m = "{}: {:.3f}".format(k, clean_dict[k])
                msg.append(m)

            print(str(msg).replace("[", "").replace("]", "").replace("'", ""))
        
        total_loss = -expression_term - regularizer_term + lineage_term + density_term + clone_term + regularizer_bin_term + uniformity_term

        return total_loss, gv_term, vg_term, density_term, regularizer_term, lineage_term, clone_term, regularizer_bin_term, uniformity_term

    def train(self,device="cpu", num_epochs=1000, learning_rate=0.1, print_each=100):
        """
        Run the optimizer and returns the mapping outcome.

        Args:
            num_epochs (int): Number of epochs.
            learning_rate (float): Optional. Learning rate for the optimizer. Default is 0.1.
            print_each (int): Optional. Prints the loss each print_each epochs. If None, the loss is never printed. Default is 100.

        Returns:
            output (ndarray): The optimized mapping matrix M (ndarray), with shape (number_cells, number_spots).
            training_history (dict): loss for each epoch
        """
        if self.random_state:
            torch.manual_seed(seed=self.random_state)
        self.optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)
        
        if print_each:
            logging.info(f"Printing scores every {print_each} epochs.")

        keys = ["total_loss", "gv_loss", "vg_reg", "kl_reg", "entropy_reg", "lineage_reg", "clone_reg","unimodal_reg","uniformity_reg"]
        values = [[] for i in range(len(keys))]
        training_history = {key: value for key, value in zip(keys, values)}
        for t in range(num_epochs):
            if print_each is None or t % print_each != 0:
                run_loss = self._loss_fn(device=device,verbose=False)
            else:
                run_loss = self._loss_fn(device=device,verbose=True)

            loss = run_loss[0]

            for i in range(len(keys)):
                training_history[keys[i]].append(str(run_loss[i]))

            
            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.parameters(), 1)
            self.optimizer.step()

            
        # take final softmax w/o computing gradients
        with torch.no_grad():
            output = softmax(self.M, dim=1).cpu().numpy()
            return output, training_history
