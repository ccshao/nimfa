"""
File Name : PY_fun_snmf.py
Created By : Shao chunxuan
Creation Date : Wed Jan 14 15:22:56 2015
Last Modified : Wed 24 May 2017 04:10:00 PM CEST
Usage : python PY_fun_01.py
Description:
    - adapt to the latest nimfa version 1.3.2

"""

import nimfa
import pandas as pd
import pdb
# import matplotlib.pyplot as plt
import subprocess
import glob
import os
import os.path

RHeatmap = "/home/shao/Desktop/Projects/NMF/01_automate/test/nmfHeatmap.R"
# RHeatmap = "/home/shao/Desktop/Projects/NMF/basicScripts/nmfHeatmap.R"

"""
## pseudo-codes for iterative running
## the start input file and python script in the same folder
Func_:
    max_depth -= 1
    if max_depth > 0:
        list all input files
        for each input file:
            read the file
            get the rank ranges
            create subfolders in the same folder
            ## generate up to 20 ranks and calculate entroy and other things
            for each subfolder:
                cd to it
                NMF with desired rank number
                generate sub input files
                then run Func_
- npDat must contain "_input_" in the file name, e.g., Cortex_hippocampus_input_.csv
- do not estimate rank by dispersion
"""
def FUN_NMF_internative(nimfa_method, 
        seed = "nndsvd",
        depth_level = 1, 
        max_rank = 5, 
        max_depth = 2, 
        n_run = 50, 
        max_iter = 50000, 
        # estimateRank = False, 
        cores = 3,
        dataSource = "unknown",
        algor = "unknown"):
        # **argsToMethod):
    if max_depth > 0:
        current_depth = max_depth - 1
        LI_inputFiles = glob.glob("*_input_*")
        depth_level += 1
        for inputFile in LI_inputFiles:
            rank_ranges = range(2, max_rank + 1)

            ## generate folders for different rank
            current_path = os.getcwd()
            if "metaCells" not in inputFile:
                LI_prefix = ["depth_" + str(depth_level) + "_r" + str(rank)  for rank in rank_ranges]
            else:
                metaCells = inputFile.split(".")[1].rsplit("_", 2)[0]
                LI_prefix = [metaCells + "_depth_" + str(depth_level) + "_r" + str(rank)  for rank in rank_ranges]
            LI_folders = [os.path.join(current_path, folder_name) for folder_name in LI_prefix]
            print "(II) working with depth:", depth_level
            print "(II) input file is:", inputFile
            ## results will be overwritten!
            ## read the data
            dat_path = os.path.join(current_path, inputFile)
            npDat= pd.read_csv(dat_path, index_col = 0, sep = '\t', header = 0)

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
                FUN_NMF_internative_run(npDat, 
                        nimfa_method,
                        seed, 
                        # estimateRank, 
                        rank,
                        cores,
                        n_run, 
                        max_iter,
                        dataSource, 
                        algor)
                        # **argsToMethod)
                # pdb.set_trace()
                FUN_NMF_generate_files(npDat)
                FUN_NMF_internative(nimfa_method,
                        seed, 
                        depth_level, 
                        max_rank,
                        current_depth, 
                        n_run, 
                        max_iter,
                        # estimateRank, 
                        cores,
                        dataSource,
                        algor)
                        # **argsToMethod)
                os.chdir("..")
    else:
        print "(==) reach max depth"

################################################################################

"""
The wrapper function for every thing;
argsToMethod contains all parameters that are used by nimfa;
seed is an important parameter here, "nndsvd" will lead to run only once, and it is the default;
Not sure why multiple core computing does not work in interactive shells; but the results is correct;
"""
def FUN_NMF_internative_run(npDat, 
        nimfa_method, 
        seed, 
        # estimateRank, 
        rank, 
        cores, 
        n_run, 
        max_iter,
        datSource, 
        algor):
        # **argsToMethod):
    ## process the path to the file;
    # pdb.set_trace()
    dat = npDat
    rowN = list(dat.index)
    colN = dat.columns.values.tolist()

    npDat = dat.as_matrix()
    n_run = 1 if seed == "nndsvd" else n_run

    print "(==) Target matrix size:", dat.shape

    input_prefix = os.path.basename(os.getcwd())
    fileName = input_prefix + "." + datSource + "." + algor + ".R" + str(rank)
    # res = FUN_nmf_run(npDat, nimfa_method, seed, rank, n_run, max_iter, **argsToMethod)
    res = FUN_nmf_run(npDat, nimfa_method, seed, rank, n_run, max_iter)

    NMF_sparW, NMF_sparH = res.fit.sparseness()
    NMF_rss = res.fit.rss()
    NMF_evar = res.fit.evar()
    # print "(==) Algorithm: %s" %(algor)
    print "(==) NMF rss: %.4f, evar: %.4f, spare_W: %.4f, spare_H: %.4f" %(NMF_rss, NMF_evar, NMF_sparW, NMF_sparH)
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
## get the files for subgroups
## Two type of files to get:
##   - subgroups cells with gene expression
##   - marker genes in all cells (check the methods I used)
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

# """ 
# return the following value
#      sparseness
#      rss
#      evar
#      cophenetic
#      dispersion
# """
# def FUN_estimateRank(npDat, 
#         nimfa_method, 
#         seed, 
#         rank, 
#         n_run, 
#         max_iter,
#         **argsToMethod):
#     print ("(II) start to calcualte metrics for rank:", rank)
#     res = FUN_nmf_run(npDat, nimfa_method, seed, rank, n_run, **argsToMethod)
#     NMF_sparW, NMF_sparH = res.fit.sparseness()
#     NMF_rss = res.fit.rss()
#     NMF_evar = res.fit.evar()
#     NMF_cophenetic = res.fit.coph_cor2()
#     NMF_dispersion = res.fit.dispersion2()
#     NMF_quality = [rank, NMF_sparW, NMF_sparH, NMF_rss, NMF_evar, NMF_cophenetic, NMF_dispersion]
#     return NMF_quality

################################################################################
"""
call nimfa. Put additionial parameters like: min_residuals=1e-05
"""
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
# def FUN_nmf(datPath, nimfa_method, seed = "nndsvd", estimateRank = False, 
#         rankRange = range(2, 11), datSource = "notSpecified", cores = 3, 
#         rank = 1, n_run = 50, algor = None, **argsToMethod):
#     """
#     The wrapper function for every thing;
#     argsToMethod contains all parameters that are used by nimfa;
#     seed is an important parameter here, "nndsvd" will lead to run only once, and it is the default;
#     Not sure why multiple core computing does not work in interactive shells; but the results is correct;
#     """
#     ## process the path to the file;
#     # pdb.set_trace()
#     dat = pd.read_csv(datPath, index_col = 0, sep = '\t', header = 0)
#     rowN = list(dat.index)
#     colN = dat.columns.values.tolist()
# 
#     npDat = dat.as_matrix()
#     n_run = 1 if seed == "nndsvd" else n_run
# 
#     print "(II) Target matrix size:", dat.shape
#     if estimateRank == True:
#         ## parallel running;
#         print "(II) seed method in estimate rank:", seed
#         print "(II) number of runs in estimate rank:", n_run
#         res = joblib.Parallel(n_jobs=cores) (joblib.delayed(FUN_estimateRank) (
#             npDat, nimfa_method, seed = seed, rank = rank, n_run = n_run, **argsToMethod) for rank in rankRange)  
# 
#         # pdb.set_trace()
#         # fileName = datSource + "." + algor 
#         res = pd.DataFrame(res)
#         # res.columns = ["rank", "spare_W", "spare_H", "rss", "evar", "cophenetic", "dispersion"]
#         res.columns = ["rank", "spare_W", "spare_H", "rss", "evar", "cophenetic2", "dispersion2"]
#         # res.to_csv(fileName + ".metric.csv", sep = "\t", quoting = False) 
#         resPlot = res.plot(x = "rank", grid = False, subplots=True, layout=(3, 2))
#         plt.savefig(datSource + "." + algor + ".estimateRank.pdf")
#         res.to_csv(datSource + "." + algor + ".estimateRank.csv", sep = "\t", index = False, quoting = False)
#     else:
#         fileName = datSource + "." + algor + ".R" + str(rank)
#         # fileName_obj = fileName + ".pkl"
# 
#         ## save the data for later use;
#         res = FUN_nmf_run(npDat, nimfa_method, seed = seed, rank = rank, n_run = n_run, **argsToMethod)
#         # save_object(res, fileName_obj)
# 
#         NMF_sparW, NMF_sparH = res.fit.sparseness()
#         NMF_rss = res.fit.rss()
#         NMF_evar = res.fit.evar()
#         print ("(II) Algorithm: %s" %(algor))
#         print ("(II) NMF rss: %.4f, evar: %.4f, spare_W: %.4f, spare_H: %.4f" 
#                 %(NMF_rss, NMF_evar, NMF_sparW, NMF_sparH))
#         ## to start with 1 instead of 0
#         metaCells = ["metaCells_" + str(i + 1) for i in range(rank)]
# 
#         # pdb.set_trace()
#         # pdb.set_trace()
#         try:
#             resW = pd.DataFrame(res.basis().todense())
#             # print " (II) scipy sparse matrix"
#         except:
#             resW = pd.DataFrame(res.basis())
#             # print "(II) numpy matrix"
#         resW.index = rowN
#         resW.columns = metaCells
#         resW.to_csv(fileName + ".basis.csv", sep = "\t", quoting = False) 
#         ## normalized by row;
#         resWNor = resW.div(resW.sum(axis=1), axis=0)
#         fileName_W = fileName + ".basis.nor.csv"
#         resWNor.to_csv(fileName_W, sep = "\t", quoting = False) 
#         subprocess.call(["Rscript", RHeatmap, fileName_W])
# 
#         try:
#             resH = pd.DataFrame(res.coef().todense())
#             # print "(II) scipy sparse matrix"
#         except:
#             resH = pd.DataFrame(res.coef())
#             # print "(II) numpy matrix"
#         resH.columns = colN 
#         resH.index = metaCells
#         resH = resH.T
#         resH.to_csv(fileName + ".coefficient.csv", sep = "\t", quoting = False) 
#         resHNor = resH.div(resH.sum(axis = 1), axis = 0)
#         fileName_H = fileName + ".coefficient.nor.csv"
#         resHNor.to_csv(fileName_H, sep = "\t", quoting = False) 
# 
# 

## code for estimate ranks
#     if estimateRank == True:
#         ## parallel running;
#         print "(II) seed method in estimate rank:", seed
#         print "(II) number of runs in estimate rank:", n_run
#         rankRange = range(2, rank + 1)
#         res = joblib.Parallel(n_jobs=cores) (joblib.delayed(FUN_estimateRank) (
#             npDat, nimfa_method, seed, rank, n_run, max_iter, **argsToMethod) for rank in rankRange)  
# 
#         # pdb.set_trace()
#         res = pd.DataFrame(res)
#         res.columns = ["rank", "spare_W", "spare_H", "rss", "evar", "cophenetic2", "dispersion2"]
#         # res.to_csv(fileName + ".metric.csv", sep = "\t", quoting = False) 
#         resPlot = res.plot(x = "rank", grid = False, subplots=True, layout=(3, 2))
#         plt.savefig(datSource + "." + algor + ".estimateRank.pdf")
#         res.to_csv(datSource + "." + algor + ".estimateRank.csv", sep = "\t", index = False, quoting = False)


