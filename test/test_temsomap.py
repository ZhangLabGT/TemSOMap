from src import Mapper
from src import temso_utils
import torch
import logging
#import squidpy as sq
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import os
from scipy.spatial import distance
from biotite.sequence.phylo import neighbor_joining
from Bio import Phylo
from io import StringIO

from multiprocessing import Process
from multiprocessing import Semaphore

def define_clones(spatial_map,nclones = 4):
    ncells = spatial_map.shape[0]
    clone_size = ncells/nclones
 
    spatial_map['clonalid'] = spatial_map['cellid']
    for i in range(nclones):
        max_id = (i+1) * clone_size
        min_id = i * clone_size + 1
        spatial_map.loc[(spatial_map['cellid'] >= min_id) & (spatial_map['cellid'] <= max_id), 'clonalid'] = i
 
    return spatial_map

def load_spatedsim(header,ncells,mr,rd,rm,run):
    content = ["spatial_meta","spatial_counts","spatial_map","meta","expr","cm"]
    params = "_ncells_" + str(ncells) + "_mr_" + str(mr) + "_rd_" + str(rd) + "_rm_" + str(rm) + "_run_" + str(run)

    items = []
    for prefix in content:
        if prefix in ["cm"]:
            suffix = ".txt"
        else:
            suffix = ".csv"

        path = os.path.join(header,prefix+params+suffix)
        if prefix in ["spatial_meta","meta"]:
            sep = ","
        elif prefix in ["spatial_map"]:
            sep = "\t"
        else:
            sep = ' '

        items.append(pd.read_csv(path,sep = sep))
    
    items.append(params)
    
    return items


def Temso_run(header,ncells,mr,rd,rm,run,lambda_l = 0.01,training_genes = range(500),lambda_d = 0.5,lambda_g1 = 1,lambda_g2 = 1,lambda_r = 0.1, lambda_c = 0.01, lambda_u = 0.001, lambda_cu = 0.001,
             d_source = None,learning_rate = 0.001,random_state = np.random.seed(seed=1),num_epochs = 10000,device = "cuda:3",print_each = 100,out_dir = "/project/xpan78/results_temso_ablation/"):
    hyperparameters = {
    "ncells" : ncells,
    "mr" : mr,
    "rd" : rd,
    "rm" : rm,  
    "run" : run
    }

    spatial_meta,spatial_count,spatial_map,meta,expr,cm,params_run = load_spatedsim(header = header, **hyperparameters)

    spatial_map.set_index('cellid')
    spatial_map = spatial_map.reindex(spatial_map['cellid'].map(dict(zip(meta.cellID.values, range(len(meta.cellID.values))))).sort_values().index)
    spatial_map = spatial_map.reset_index() 
    spatial_map = define_clones(spatial_map)

    spatial_map = temso_utils.define_clones(spatial_map,nclones = 4)

    d_hamming = distance.squareform(distance.pdist(cm,metric = "hamming"))
    tree = neighbor_joining(d_hamming)
    tree_newick = tree.to_newick()
    tree_phylo = Phylo.read(StringIO(tree_newick),"newick")

    Z = [[0 for _ in range(ncells)] for _ in range(ncells)]
    for leaf1 in tree.leaves:
        for leaf2 in tree.leaves:
            Z[leaf1.index][leaf2.index] = tree.get_distance(leaf1.index,leaf2.index)

    Z = np.array(Z)

    X = expr.to_numpy()
    Y = spatial_count.to_numpy()

    v = spatial_meta[["x", "y"]]
    v = v.to_numpy()

    Z_clone = temso_utils.Define_clones(tree=tree_phylo)

    d = spatial_meta.density.values

    device = torch.device(device)  # for gpu


    print_each = 100

    hyperparameters = {
        "lambda_d": lambda_d,  # KL (ie density) term
        "lambda_g1": lambda_g1,  # gene-voxel cos sim
        "lambda_g2": lambda_g2,  # voxel-gene cos sim
        "lambda_r": lambda_r,  # regularizer: penalize entropy
        "lambda_l": lambda_l,  # regularizer: penalize lineage dissimilarity
        "lambda_c": lambda_c,  # regularizer: penalize clone dissimilarity
        "lambda_u": lambda_u, #regularizer: penalize non-unimodal distribution of each cell's inferred location
        "lambda_cu": lambda_cu, #regularizer: penalize non-uniform distribution of the overall distribution of cells
        "d_source": d_source,
    }


    logging.info(
        "Begin training with {} genes and with density_prior...".format(
            len(training_genes)
        )
    )

    mapper = Mapper.Mapper(
        X=X, Y=Y, Z=d_hamming, Z_clone=Z_clone, v=v,v_gt=None, d=d, device=device, random_state=random_state, **hyperparameters
    )


    mapping_matrix, training_history = mapper.train(
        device=device,learning_rate=learning_rate, num_epochs=num_epochs, print_each=print_each
    )

    print ('Writing to '+ out_dir)
    mapping_matrix.tofile(out_dir + "mapping_matrix_new_density" + params_run + ".csv",sep = ",")

    print ("Writing training histories to " + out_dir)
    for loss in training_history:
        curve_vec = training_history[loss]
        numeric_values = []
        plt.figure()
        for string in curve_vec:
            if (len(string) > 5):
                sub1 = "tensor("
                sub2 = ", device"
                idx1 = string.index(sub1)
                idx2 = string.index(sub2)
                numeric_values.append(float(string[(idx1 + len(sub1)):idx2]))
        
        plt.plot(numeric_values, marker='o')
        plt.xlabel('Iteration')
        plt.ylabel(loss)
        plt.title('Training history of TemSOMap')
        plt.grid(True)
        plt.savefig(out_dir + loss + "_new_density" + params_run + '.pdf')


    v_gt = spatial_map[["x", "y"]]
    v_gt = v_gt.to_numpy()
    MSE_st,MSE_coords,N_mismatches,j_index_individual,j_index_whole,L1_coords= temso_utils.Benchmark_temso(X,Y,mapping_matrix,v_gt,v,method="mean")

    with open(out_dir + "metrics_new_density" + params_run + ".txt", "w") as text_file:
        text_file.write(str(MSE_st) + "," + str(MSE_coords) + "," + str(N_mismatches)+ "," + str(j_index_individual)+ "," + str(j_index_whole)+ "," + str(L1_coords))
    


    v_pred,M_pred = temso_utils.generate_spatial_map(X,mapping_matrix,v,method = "mean")

    plt.figure(figsize=(16,12))
    #plt.rcParams["figure.figsize"] = (160,120) 
    plt.subplot(2,2,1)
    xs = spatial_map.x.copy().values
    ys = spatial_map.y.copy().values
    plt.axis('off')
    color_map = "cell_type"
    plt.scatter(xs, ys, c= spatial_map[color_map].values,cmap=plt.cm.Spectral);
    plt.gca().invert_yaxis()
    plt.title("spatial map of cells colored by " + color_map);

    plt.subplot(2,2,2)
    xs = spatial_map.x.copy().values
    ys = spatial_map.y.copy().values
    plt.axis('off')
    color_map = "clonalid"
    plt.scatter(xs, ys, c= spatial_map[color_map].values,cmap=plt.cm.Spectral);
    plt.gca().invert_yaxis()
    plt.title("spatial map of cells colored by " + color_map);

    plt.subplot(2,2,3)
    xs = v_pred[:,0] + np.random.normal(0, 0.5, v_pred.shape[0])
    ys = v_pred[:,1] + np.random.normal(0, 0.5, v_pred.shape[0])
    plt.axis('off')
    color_map = "cell_type"
    plt.scatter(xs, ys, c= spatial_map[color_map].values,cmap=plt.cm.Spectral);
    plt.gca().invert_yaxis()
    plt.title("inferred spatial map of cells colored by " + color_map);

    plt.subplot(2,2,4)
    xs = v_pred[:,0] #+ np.random.normal(0, 0.5, v_pred.shape[0])
    ys = v_pred[:,1] #+ np.random.normal(0, 0.5, v_pred.shape[0])
    plt.axis('off')
    color_map = "clonalid"
    plt.scatter(xs, ys, c= spatial_map[color_map].values,cmap=plt.cm.Spectral);
    plt.gca().invert_yaxis()
    plt.title("inferred spatial map of cells colored by " + color_map);
    plt.tight_layout()
    plt.savefig(out_dir + "scatters_new_density" + params_run + '.pdf')


    sema.release()


if __name__ == "__main__":  # confirms that the code is under main function
    header = '/project/xpan78/sim_spatedsim_continuous_1024/'
    ncells_list = [1024]
    mr_list = [0]
    run_list = [7]
    rd_list = [3]
    rm_list = [4]
    lambda_l = 0.01
    procs = []
    #proc = Process(target=Cas_run)  # instantiating without any argument
    #procs.append(proc)
    #proc.start()
    args = []
    for ncells in ncells_list:
        for m_r in mr_list:
            for rd in rd_list:
                for run in run_list:
                    for rm in rm_list:
                        args.append((header,ncells,m_r,rd,rm,run,lambda_l))

    
    concurrency = 1
    sema = Semaphore(concurrency)
    for arg in args:
        # print(name)
        sema.acquire()
        proc = Process(target=Temso_run, args=arg)
        procs.append(proc)
        proc.start()
    
    
    for proc in procs:
        proc.join()
