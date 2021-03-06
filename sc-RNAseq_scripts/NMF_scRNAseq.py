
"""run NMF on sc-RNAseq data.

Usage:
    NMF_scRNAseq.py --min_rank=<int> --max_rank=<int> --max_depth=<int> --RHeatmap=<str> 
        [--method=<str>] [--seed=<str>] [--depth_level=<int>] 
        [--n_run=<int>] [--max_iter=<int>] [--cores=<int>] 
        [--dataSource=<str>] [--algor=<str>]

Options:
    -h --help              # show this screen
    --min_rank <int>       # minimum number of rank
    --max_rank <int>       # maximum number of rank
    --max_depth <int>      # maximum depth to run NMF
    --RHeatmap <str>       # full path to the heatmap R files
    --method <str>         # method implemented in nimfa. lsnmf, bd and nmf (standard NMF) are [default: lsnmf]
    --seed <str>           # seeding choice [default: nndsvd]
    --depth_level <int>    # depth_level, 1 for the start [default: 1]
    --n_run <int>          # number of run [default: 1]
    --max_iter <int>       # maximum number of iterations [default: 50000]
    --cores <int>          # cores to use (not implemented currently) [default: 3]
    --dataSource <str>     # name of the project [default: unknown]
    --algor <str>          # description of running method [default: unknown]

"""

import pandas as pd
import pdb
import subprocess
import glob
import os
import os.path
from docopt import docopt
import nimfa

from contextlib import contextmanager
import timeit
import time

@contextmanager
def time_elapsed(label):
    start_time = timeit.default_timer()
    try:
        yield
    finally:
        elapsed = float("%.2f" % (timeit.default_timer() - start_time))
        # end = time.clock()
        print "elapsed seconds by %s: %s" %(label, elapsed)
        print "\n\n"

################################################################################

def FUN_NMF_iterative(nimfa_method, 
        seed,
        depth_level, 
        min_rank,
        max_rank, 
        max_depth, 
        n_run, 
        max_iter, 
        # estimateRank = False, 
        cores,
        dataSource,
        algor,
        RHeatmap):

    # pdb.set_trace()
    if max_depth > 0:
        current_depth = max_depth - 1
        LI_inputFiles = glob.glob("*_input_*")
        # depth_level += 1
        for inputFile in LI_inputFiles:
            rank_ranges = range(min_rank, max_rank + 1)

            ## generate folders for different rank
            current_path = os.getcwd()
            if "metaCells" not in inputFile:
                LI_prefix = ["depth_" + str(depth_level) + "_r" + str(rank) for rank in rank_ranges]
            else:
                metaCells = inputFile.split(".")[1].rsplit("_", 2)[0]
                LI_prefix = [metaCells + "_depth_" + str(depth_level) + "_r" + str(rank)  for rank in rank_ranges]
            LI_folders = [os.path.join(current_path, folder_name) for folder_name in LI_prefix]
            print "(II) working with depth:", depth_level
            print "(II) input file:", inputFile
            ## results will be overwritten!
            ## read the data
            dat_path = os.path.join(current_path, inputFile)
            npDat= pd.read_csv(dat_path, index_col = 0, sep = '\t', header = 0)

            ## skip if less than 10 cells exit!
            if npDat.shape[0] < 10:
                print "(II) !Skip as there are less 10 cells in the cluster!"
            else:
                # FUN_estimateRank(npDat, nimfa_method, seed, rank_ranges)
                for current_folder in LI_folders: 
                    if not os.path.exists(current_folder):
                        # pdb.set_trace()
                        os.makedirs(current_folder)
                for current_folder in LI_folders:
                    ## folder name: depth_1_1, depth_1_2
                    os.chdir(current_folder)
                    ## actual run, will modify from existing codes
                    rank = int(os.path.basename(current_folder).split("_")[-1].split("r")[1])
                    FUN_NMF_iterative_run(npDat, 
                            nimfa_method,
                            seed, 
                            # estimateRank, 
                            rank,
                            cores,
                            n_run, 
                            max_iter,
                            dataSource, 
                            algor, 
                            RHeatmap)
                    # pdb.set_trace()
                    FUN_NMF_generate_files(npDat)
                    print "(II) start new run..."
                    FUN_NMF_iterative(nimfa_method,
                            seed, 
                            depth_level + 1, 
                            min_rank,
                            max_rank,
                            current_depth, 
                            n_run, 
                            max_iter,
                            # estimateRank, 
                            cores,
                            dataSource,
                            algor, 
                            RHeatmap)
                    os.chdir("..")
    else:
        print "(==) reach max depth\n"

################################################################################

def FUN_NMF_iterative_run(npDat, 
        nimfa_method, 
        seed, 
        # estimateRank, 
        rank, 
        cores, 
        n_run, 
        max_iter,
        datSource, 
        algor, 
        RHeatmap):
    ## process the path to the file;
    # pdb.set_trace()
    dat = npDat
    rowN = list(dat.index)
    colN = dat.columns.values.tolist()

    npDat = dat.as_matrix()
    n_run = 1 if seed == "nndsvd" else n_run
    # pdb.set_trace()
    input_prefix = os.path.basename(os.getcwd())
    fileName = input_prefix + "." + datSource + "." + algor + ".R" + str(rank)
    res = FUN_nmf_run(npDat, nimfa_method, seed, rank, n_run, max_iter)

    NMF_sparW, NMF_sparH = res.fit.sparseness()
    NMF_rss = res.fit.rss()
    NMF_evar = res.fit.evar()
    # print "(==) Algorithm: %s" %(algor)
    print "(==) rss: %.2f, evar: %.2f, spare_W: %.2f, spare_H: %.2f \n" %(NMF_rss, NMF_evar, NMF_sparW, NMF_sparH)
    ## to start with 1 instead of 0
    metaCells = ["metaCells_" + str(i + 1) for i in range(rank)]

    ## get W matrix
    try:
        resW = pd.DataFrame(res.basis().todense())
        # print " (II) scipy sparse matrixLI_cells"
    except:
        resW = pd.DataFrame(res.basis())
        # print "(II) numpy matrix"
    resW.index = rowN
    resW.columns = metaCells
    resW.to_csv(fileName + ".basis.csv", sep = "\t", quoting = False) 
    ## normalized by row;
    resWNor = resW.div(resW.sum(axis=1), axis=0)
    fileName_W = fileName + ".basis.nor.csv"
    resWNor.to_csv(fileName_W, sep = "\t", quoting = False) 
    subprocess.call(["Rscript", "--vanilla", RHeatmap, fileName_W])

    ## get H matrix
    try:
        resH = pd.DataFrame(res.coef().todense())
        # print "(II) scipy sparse matrix"
    except:
        resH = pd.DataFrame(res.coef())
        # print "(II) numpy matrix"
    resH.columns = colN 
    resH.index = metaCells
    resH = resH.T
    resH.to_csv(fileName + ".coefficient.csv", sep = "\t", quoting = False) 
    resHNor = resH.div(resH.sum(axis = 1), axis = 0)
    fileName_H = fileName + ".coefficient.nor.csv"
    resHNor.to_csv(fileName_H, sep = "\t", quoting = False) 

################################################################################
 
"""
subgroups cells with gene expression
"""
def FUN_NMF_generate_files(npDat):
    # pdb.set_trace()
    basisMat_file = glob.glob("*basis.nor.csv")[0]
    basis_mat = pd.read_csv(basisMat_file, index_col = 0, sep = '\t', header = 0)
    basis_cl = basis_mat.idxmax(axis = 1)
    df_anno = {"cells": basis_cl.index.tolist(), "cls": basis_cl}
    df_anno = pd.DataFrame(df_anno)
    ## convert the basis_cl to a dataFrame in pd
    cls = basis_cl.unique()
    cls.sort()
    for metaCell in cls:
        LI_cells = df_anno[(df_anno.cls == metaCell)]["cells"].tolist()
        # LI_cells = df_anno[(df_anno.cls == "metaCells_2")]["cells"].tolist()
        ## select cells only in each metacells
        dat_sel = npDat.loc[LI_cells]
        dat_sel = dat_sel[(dat_sel.T != 0).any()]
        ## get the basename of current path
        ## save the file
        input_prefix = os.path.basename(os.getcwd())
        # pdb.set_trace()
        dat_sel.to_csv(input_prefix + "." + metaCell + "_input_.csv", sep = "\t", quoting = False)

################################################################################

def FUN_nmf_run(npDat, 
        nimfa_method, 
        seed, 
        rank, 
        n_run, 
        max_iter):
        # **argsToMethod):
    nmfModel = nimfa_method(npDat,
        seed = seed, 
        rank = rank, 
        max_iter = max_iter, 
        n_run = n_run)
        # **argsToMethod)
    print "(II) rank: %d, method: %s, seed: %s, n_run: %d" %(nmfModel.rank, nmfModel.name, str(nmfModel.seed), nmfModel.n_run)
    print "(==) starting NMF calculation..."
    nmfFit = nmfModel()
    return nmfFit

################################################################################

if __name__ == '__main__':
    arguments = docopt(__doc__)
    # print arguments
    nimfa_methods = {"lsnmf":nimfa.Lsnmf, "bd":nimfa.Bd, "nmf":nimfa.Nmf}

    FUN_NMF_iterative(nimfa_methods[arguments["--method"]],
             seed = arguments["--seed"],
             depth_level = int(arguments["--depth_level"]),
             min_rank = int(arguments["--min_rank"]), 
             max_rank = int(arguments["--max_rank"]), 
             max_depth = int(arguments["--max_depth"]), 
             n_run = int(arguments["--n_run"]),
             cores = int(arguments["--cores"]),
             max_iter = int(arguments["--max_iter"]), 
             dataSource = arguments["--dataSource"], 
             algor = arguments["--algor"],
             RHeatmap = arguments["--RHeatmap"])


