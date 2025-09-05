from copy import deepcopy
import logging
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import matplotlib.cm as cm
import IPython.display as IPdisplay
import glob
import re
import gc
from numpy.random import randn
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm as normal
import scipy.stats
from scipy import sparse
from operator import itemgetter, attrgetter
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE #Barnes hut approximation
from sklearn.neighbors import NearestNeighbors
from numba import jit
import umap
import seaborn as sns
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
warnings.simplefilter("ignore", category=DeprecationWarning)
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter("ignore", category=NumbaPendingDeprecationWarning)
#pip install -U loompy
import loompy
from velocyto.neighbors import knn_distance_matrix, connectivity_to_weights, convolve_by_sparse_weights, BalancedKNN
from velocyto.estimation import fit_slope, fit_slope_offset, fit_slope_weighted, fit_slope_weighted_offset #if introduce own version of the model fit remove this line and add below
from velocyto.estimation import clusters_stats #if introduce own version of the model fit remove this line and add below
from velocyto.estimation import colDeltaCor, colDeltaCorSqrt, colDeltaCorLog10, colDeltaCorpartial, colDeltaCorSqrtpartial, colDeltaCorLog10partial #if introduce own version of the model fit remove this line and add below
#import sys
#sys.path.insert(0, '/home/cdunican/PhD/')
#import py_file_name_wo_py_at_end as estimation
#print (estimation.get_set1()) #add estimation. onto the names of: fit_slope, fit_slope_offset, fit_slope_weighted, fit_slope_weighted_offset, clusters_stats, colDeltaCor, colDeltaCorSqrt, colDeltaCorLog10, colDeltaCorpartial, colDeltaCorSqrtpartial, colDeltaCorLog10partial functions
from velocyto.diffusion import Diffusion
from velocyto.serialization import dump_hdf5, load_hdf5
from typing import *
import h5py
import pickle
import zlib
import os
from matplotlib import cm
from collections import OrderedDict
import sys
from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
from numpy import *
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from rpy2.rinterface import SexpVector
import seaborn as sns
from numpy_groupies import aggregate, aggregate_np
from rpy2.robjects import pandas2ri
import shutup; shutup.please()
gc.collect()

#warnings.filterwarnings("ignore") 
#import pandas.rpy.common as com
#Make sense of the last plot - email the authors because no indication of how it picks which point to circle and add an arrow to.

#low priority:
#look at the non-default plots in:
#https://github.com/velocyto-team/velocyto-notebooks/blob/master/python/DentateGyrus.ipynb
#plot in 3D tSNE, add arrows to 3D PCA, add arrows to 2D and 3D diffusion maps
#plot in 3d tSNE
#add arrows to 3d pca
#add arrows to 3d diffusion maps

class VelocytoLoom:
    def __init__(self, loom_filepath: str) -> None:
        self.loom_filepath = loom_filepath
        ds = loompy.connect(self.loom_filepath)
        ds.layers["spliced"] = spl     #has to be a matrix full of integers, rows: cells, columns: genes i.e.
        ds.layers["unspliced"] = unspl #np.array([[0,0,0,0], [0,0,0,0], [1,1,1,1]])
        self.S = ds.layer["spliced"][:, :]
        self.U = ds.layer["unspliced"][:, :]
        self.ca = dict(ds.col_attrs.items())
        self.ra = dict(ds.row_attrs.items())
        #ds.close()
        self.initial_cell_size = self.S.sum(0)
        self.initial_Ucell_size = self.U.sum(0) #sum(0) axis = 0 (sums all row values in each colum) [8. 8. 8. 8. 8. 8.] 8 genes, 6 samples, all 6 have all the spliced and unspliced counts
        try:
            if np.mean(self.ca["_Valid"]) < 1:
                logging.warn(f"fraction of _Valid cells is {np.mean(self.ca['_Valid'])} but all will be taken in consideration")
        except KeyError:
            pass
            # logging.debug("The file did not specify the _Valid column attribute")
        self.norm_factor = 1
        self.S_sz = self.norm_factor * self.S
        self.S_norm = self.S
        self.U_sz = self.norm_factor * self.U
        self.U_norm = self.U



    """A convenient object to store the data of a velocyto loom file.

    Data will be stored in memory

    Examples
    --------
    For usage examples consult the documentation

    Attributes
    ----------
    S: np.ndarray
        Expressed spliced molecules
    U: np.ndarray
        Unspliced molecule count
    ca: dict
        Column attributes of the loom file
    ra: dict
        Row attributes of the loom file
    loom_filepath: str
        The original path the loom files has been read from
    initial_cell_size: int
        The sum of spliced molecules
    initial_Ucell_size: int
        The sum of unspliced molecules
    """

    def plot_fractions(self, save2file: str=None) -> None: #to improve: yaxis: number of molecules, xaxis: cells/patients currently gives a weird not that useful graph so needs improving -cd2514/4/12/2018
        """Plots a barplot showing the abundance of spliced/unspliced molecules in the dataset

        Arguments
        ---------
        save2file: str (default: None)
            If not None specifies the file path to which plots get saved

        Returns
        -------
        Nothing, it plots a barplot
        """
        plt.figure(figsize=(3.2, 5))
        try:
            chips, chip_ix = np.unique(self.ca["Patients"], return_inverse=1)
        except KeyError:
            chips, chip_ix = np.unique([i.split(":")[0] for i in self.ca["Patients"]], return_inverse=1)
        n = len(chips)
        for i in np.unique(chip_ix):
            tot_mol_cell_submatrixes = [X[:, chip_ix == i].sum(0) for X in [self.S, self.U]]
            total = np.sum(tot_mol_cell_submatrixes, 0)
            _mean = [np.mean(j / total) for j in tot_mol_cell_submatrixes]
            _std = [np.std(j / total) for j in tot_mol_cell_submatrixes]
            plt.ylabel("Fraction")
            plt.bar(np.linspace(-0.2, 0.2, n)[i] + np.arange(3), _mean, 0.5 / (n * 1.05), label=chips[i])
            plt.errorbar(np.linspace(-0.2, 0.2, n)[i] + np.arange(3), _mean, _std, c="k", fmt="none", lw=1, capsize=2)
            # Hide the right and top spines
            #plt.gca().axison = False
            #plt.gca().stale = True
            #for k in plt.axes().spines.keys():
            #    plt.gca().spines['right'].set_visible(False)
            #    plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            # Only show ticks on the left and bottom spines
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('bottom')
            plt.agca().spines['left'].set_bounds(0, 0.8)
            plt.legend()

        plt.xticks(np.arange(2), ["spliced", "unspliced"])
        plt.tight_layout()
        if save2file:
            plt.savefig(save2file, bbox_inches="tight")
        plt.close()


    def set_clusters(self, cluster_labels: np.ndarray, cluster_colors_dict: Dict[str, List[float]]=None, colormap: Any=True) -> None: #this assigns the coloration for the points used later therefore this will later be changed for the phenotypes
        """Utility function to set cluster labels, names and colormapthis will lat

        Arguments
        ---------
        cluster_labels: np.ndarray
            A vector of strings containing the name of the cluster for each cells
        cluster_colors_dict: dict[str, List[float]]
            A mapping  cluster_name -> rgb_color_triplet for example "StemCell":[0.65,0.1,0.4]
        colormap:
            (optional)
            In alternative to cluster_colors_dict a colormap object (e.g. from matplotlib or similar callable) can be passed

        Returns
        -------
        Nothing, the attributes `cluster_labels, colorandum, cluster_ix, cluster_uid` are created.

        """
        self.cluster_labels = np.array(cluster_labels)
        if self.cluster_labels.dtype == "O":  # Fixes a bug when importing from pandas
            self.cluster_labels = self.cluster_labels.astype(np.bytes_)#AG has put np.bytes_ here as np.string_ is deprecated from numpy>2.0 and errors
            #print(self.cluster_labels)
            #for label in self.cluster_labels):
            #    str(label,'utf-8')
        if cluster_colors_dict:
            self.colorandum = np.array([cluster_colors_dict[i] for i in cluster_labels])
            self.cluster_colors_dict = cluster_colors_dict
            self.colormap = None
        else:
            if colormap is None:
                self.colorandum = colormap_fun[self.cluster_ix]
                cluster_uid = self.cluster_uid
                self.cluster_colors_dict = {cluster_uid[i]: colormap_fun(i) for i in range(len(cluster_uid))}
            else:
                self.colormap = colormap
                self.colorandum = self.colormap(self.cluster_ix)
                cluster_uid = self.cluster_uid
                self.cluster_colors_dict = {cluster_uid[i]: self.colormap(i) for i in range(len(cluster_uid))}

    @property
    def cluster_uid(self) -> np.ndarray:
        clusters_uid = np.unique(self.cluster_labels)
        return clusters_uid

    @property
    def cluster_ix(self) -> np.ndarray:
        _, cluster_ix = np.unique(self.cluster_labels, return_inverse=True)
        return cluster_ix

    def score_cluster_expression(self, min_avg_U: float=0.02, min_avg_S: float=0.08) -> np.ndarray:
        """Prepare filtering genes on the basis of cluster-wise expression threshold

        Arguments
        ---------
        min_avg_U: float
            Include genes that have unspliced average bigger than `min_avg_U` in at least one of the clusters
        min_avg_S: float
            Include genes that have spliced average bigger than `min_avg_U` in at least one of the clusters
        Note: the two conditions are combined by and "&" logical operator

        Returns
        -------
        Nothing but it creates the attribute
        clu_avg_selected: np.ndarray bool
            The gene cluster that is selected
        To perform the filtering use the method `filter_genes`
        """
        self.U_avgs, self.S_avgs = clusters_stats(self.U, self.S, self.cluster_uid, self.cluster_ix, size_limit=40)
        self.clu_avg_selected = (self.U_avgs.max(1) > min_avg_U) & (self.S_avgs.max(1) > min_avg_S)

    def custom_filter_attributes(self, attr_names: List[str], bool_filter: np.ndarray) -> None:
        """Filter attributes given a boolean array. attr_names can be dictionaries or numpy arrays

        Arguments
        ---------
        attr_names: List[str]
            a list of the attributes to be modified. The can be
            1d arrays, dictionary of 1d arrays, ndarrays, will be filtered by axis=0
            if .T is specified by axis=-1
        bool_filter:
            the boolean filter to be applied

        Returns
        -------
        Nothing it filters the specified attributes
        """
        transpose_flag = False
        for attr in attr_names:
            if attr[-2:] == ".T":
                obj = getattr(self, attr[:-2])
                transpose_flag = True
            else:
                obj = getattr(self, attr)
                transpose_flag = False
            if type(obj) is dict:
                setattr(self, attr, {k: v[bool_filter] for k, v in obj.items()})
            elif type(obj) is np.ndarray:
                if len(obj.shape) > 1:
                    if transpose_flag:
                        setattr(self, attr, obj[..., bool_filter])
                    else:
                        setattr(self, attr, obj[bool_filter, :])
                else:
                    setattr(self, attr, obj[bool_filter])
            else:
                raise NotImplementedError(f"The filtering of an object of type {type(obj)} is not defined")


    def perform_PCA(self, which: str="S_norm", n_components: int=None, div_by_std: bool=False) -> None:
        """Perform PCA (cells as samples) uses the self.pcs for the tSNE

        Arguments
        ---------
        which: str, default="S_norm"
            The name of the attribute to use for the calculation (e.g. S_norm or Sx_norm)
        n_components: int, default=None
            Number of components to keep. If None all the components will be kept.
        div_by_std: bool, default=False
            Wether to divide by standard deviation

        Returns
        -------
        Returns nothing but it creates the attributes:
        pca: np.ndarray
            a numpy array of shape (cells, npcs)

        """
        X = getattr(self, which)
        self.pca = PCA(n_components=n_components) #n_components experiment with 2 or 3, no matter the number of coordinates are the same as are the percentages
        if div_by_std:
            self.pcs = self.pca.fit_transform(X.T / X.std(0))
        else:
            self.pcs = self.pca.fit_transform(X.T)
        pcs_inte = []
        self.percentages = []
        self.percentagesRounded = []
        for per_variance in self.pca.explained_variance_ratio_:
            self.percentages.append(per_variance*100)
            self.percentagesRounded.append(round(per_variance*100, 2))


    def plot_pca(self, dim: List[int]=[0, 1, 2], elev: float=60, azim: float=-140, make_gify: bool=True, Filename: str="Original", typegroup: str="integers", changeShape: str="True", Colours: str="Paired") -> None:
        """Plot 3d PCA
        """
        PCNum = list(range(1, len(self.percentages)+1))
        Percentages = np.transpose(np.array([PCNum, self.percentages]))
        Percentages_df = pd.DataFrame(Percentages, columns=["PC", "%"])
        file_name = "PCs_percentages.csv"
        Percentages_df.to_csv(file_name, sep=",", header=True, index=False)
        if dim == [0,1]:
            if changeShape == "True":
                my_markers = self.ca["Markers"]
                mark = sorted(set(my_markers))
                corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
                count = 0
                my_types = {}
                for the_type in mark:
                    my_types[the_type] = corresp[count]
                    count = count + 1
                markers = []
                for mar in my_markers:
                    markers.append(my_types.get(mar))
                z = randn(10)
                fig, ax = plt.subplots(figsize=(8, 7)) 
                for i in range(0, len(markers)):
                    ax.scatter(self.pcs[i, dim[0]], self.pcs[i, dim[1]], marker=markers[i], color=np.array(self.colorandum[i]))
                hexy_codes = []
                my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
                decoded_my_keys = []
                for item in my_keys:
                    decoded_my_keys.append(str(item,'utf-8')) #int()
                cmap2 = matplotlib.colormaps[Colours]    # PiYG
                for i in range(cmap2.N):
                    rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                    hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
                legend_elements = []
                combined = list(zip(decoded_my_keys, hexy_codes))
                combined = sorted(combined, key=itemgetter(0))
                if (typegroup != "Strings"):
                    timepoints = [int(i[0]) for i in combined]
                    colys = [i[1] for i in combined]
                    tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                    tempDF = tempDF.sort_values('Group')
                    groupys2 = tempDF['Group'].tolist()
                    cols2 = tempDF['colours'].tolist()
                    combined = list(zip(groupys2, cols2))
                decoded_my_keys = [i[0] for i in combined]
                hexy_codes = [i[1] for i in combined]
                #get markers and names
                for sample_type in range(0, len(decoded_my_keys)):
                    to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                    legend_elements.append(to_add)
                for sample_type in range(0, len(mark)):
                    to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                    legend_elements.append(to_add)
            if changeShape == "False":
                z = randn(10)
                fig, ax = plt.subplots(figsize=(8, 7)) 
                ax.scatter(self.pcs[:,0], self.pcs[:,1], c=self.colorandum, marker="o")
                hexy_codes = []
                my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
                decoded_my_keys = []
                for item in my_keys:
                    decoded_my_keys.append(str(item,'utf-8')) #int()
                cmap2 = matplotlib.colormaps[Colours]    # PiYG
                for i in range(cmap2.N):
                    rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                    hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
                legend_elements = []
                combined = list(zip(decoded_my_keys, hexy_codes))
                combined = sorted(combined, key=itemgetter(0))
                if (typegroup != "Strings"):
                    timepoints = [int(i[0]) for i in combined]
                    colys = [i[1] for i in combined]
                    tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                    tempDF = tempDF.sort_values('Group')
                    groupys2 = tempDF['Group'].tolist()
                    cols2 = tempDF['colours'].tolist()
                    combined = list(zip(groupys2, cols2))
                decoded_my_keys = [i[0] for i in combined]
                hexy_codes = [i[1] for i in combined]
                for sample_type in range(0, len(decoded_my_keys)):
                    to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                    legend_elements.append(to_add)
            ax.legend(handles=legend_elements, title='Key') #Timepont \n(hours)
            ax.set_xlabel('PC1 '+'('+str(self.percentagesRounded[0])+'%)')
            ax.set_ylabel('PC2 '+'('+str(self.percentagesRounded[1])+'%)')
            static = "PCA_Plot_" + Filename + ".png"
            plt.savefig(static)
            plt.close()
            

        if dim == [0, 1, 2]:
            if changeShape == "True":
                my_markers = self.ca["Markers"]
                mark = sorted(set(my_markers))
                corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
                count = 0
                my_types = {}
                for the_type in mark:
                    my_types[the_type] = corresp[count]
                    count = count + 1
                markers = []
                for mar in my_markers:
                    markers.append(my_types.get(mar))
                z = randn(10)
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(111, projection='3d')
                for i in range(0, len(markers)):
                    ax.scatter(self.pcs[i, dim[0]], self.pcs[i, dim[1]], self.pcs[i, dim[2]], marker=markers[i], color=self.colorandum[i])
                hexy_codes = []
                my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
                decoded_my_keys = []
                for item in my_keys:
                    decoded_my_keys.append(str(item,'utf-8')) #int()
                cmap2 = matplotlib.colormaps[Colours]    # PiYG
                for i in range(cmap2.N):
                    rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                    hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
                legend_elements = []
                combined = list(zip(decoded_my_keys, hexy_codes))
                combined = sorted(combined, key=itemgetter(0))
                if (typegroup != "Strings"):
                    timepoints = [int(i[0]) for i in combined]
                    colys = [i[1] for i in combined]
                    tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                    tempDF = tempDF.sort_values('Group')
                    groupys2 = tempDF['Group'].tolist()
                    cols2 = tempDF['colours'].tolist()
                    combined = list(zip(groupys2, cols2))
                decoded_my_keys = [i[0] for i in combined]
                hexy_codes = [i[1] for i in combined]
                #get markers and names
                for sample_type in range(0, len(decoded_my_keys)):
                    to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                    legend_elements.append(to_add)
                for sample_type in range(0, len(mark)):
                    to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                    legend_elements.append(to_add)
            if changeShape == "False":
                z = randn(10)
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(self.pcs[:,0], self.pcs[:,1], self.pcs[:,2], c=self.colorandum, marker="o")
                hexy_codes = []
                my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
                decoded_my_keys = []
                for item in my_keys:
                    decoded_my_keys.append(str(item,'utf-8')) #int()
                cmap2 = matplotlib.colormaps[Colours]   # PiYG
                for i in range(cmap2.N):
                    rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                    hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
                legend_elements = []
                combined = list(zip(decoded_my_keys, hexy_codes))
                combined = sorted(combined, key=itemgetter(0))
                if (typegroup != "Strings"):
                    timepoints = [int(i[0]) for i in combined]
                    colys = [i[1] for i in combined]

                    tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                    tempDF = tempDF.sort_values('Group')
                    groupys2 = tempDF['Group'].tolist()
                    cols2 = tempDF['colours'].tolist()
                    combined = list(zip(groupys2, cols2))
                decoded_my_keys = [i[0] for i in combined]
                hexy_codes = [i[1] for i in combined]
                #get markers and names
                for sample_type in range(0, len(decoded_my_keys)):
                    to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                    legend_elements.append(to_add)
            ax.legend(handles=legend_elements, title='Key') #Timepont \n(hours)
            ax.set_xlabel('PC1 '+'('+str(self.percentagesRounded[0])+'%)')
            ax.set_ylabel('PC2 '+'('+str(self.percentagesRounded[1])+'%)')
            ax.set_zlabel('PC3 '+'('+str(self.percentagesRounded[2])+'%)')
            static = "PCA_" + Filename + ".png"
            plt.savefig(static)
            #plt.show()
            if (make_gify == True):
                def rotate(angle):
                    ax.view_init(azim=angle)
                rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
                name1 = "PCA_LeftRight" + Filename + ".gif"
                rot_animation.save(name1, dpi=80, writer='imagemagick')
                def init():
                    return fig,
                def animate(i):
                    ax.view_init(elev=(i-45)*4, azim=10)
                    return fig,
                ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
                name2 = "PCA_rotation_UpDown" + Filename + ".gif"
                ani_2.save(name2, dpi=80, writer='imagemagick')
            plt.close()
            """
            #animate it
            def make_views(ax,angles,elevation=None, width=8, height = 6, prefix='tmprot_',**kwargs):
                files = []
                ax.figure.set_size_inches(width,height)
                for i,angle in enumerate(angles):  
                    ax.view_init(elev = elevation, azim=angle)
                    fname = '%s%03d.png'%(prefix,i)
                    ax.figure.savefig(fname)
                    files.append(fname)
                return files 
            def make_movie(files,output, fps=10,bitrate=1800,**kwargs):
                output_name, output_ext = os.path.splitext(output)
                command = { '.mp4' : 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc lavc\
                             -lavcopts vcodec=msmpeg4v2:vbitrate=%d'
                             %(",".join(files),fps,output_name,bitrate)}                  
                command['.ogv'] = command['.mp4'] + '; ffmpeg -i %s.mp4 -r %d %s'%(output_name,fps,output)
                print(command[output_ext])
                output_ext = os.path.splitext(output)[1]
                os.system(command[output_ext])
            def make_gif(files,output,delay=100, repeat=True,**kwargs):     
                loop = -1 if repeat else 0
                os.system('convert -delay %d -loop %d %s %s'
                      %(delay,loop," ".join(files),output))
            def make_strip(files,output,**kwargs):     
                os.system('montage -tile 1x -geometry +0+0 %s %s'%(" ".join(files),output)) 
            def rotanimate(ax, angles, output, **kwargs):
                output_ext = os.path.splitext(output)[1]
                files = make_views(ax,angles, **kwargs) 
                D = { '.mp4' : make_movie,
                      '.ogv' : make_movie,
                      '.gif': make_gif ,
                      '.jpeg': make_strip,
                      '.png':make_strip}
                D[output_ext](files,output,**kwargs) 
                for f in files:
                    os.remove(f)   
            ax.get_legend().remove()
            angles = np.linspace(0,360,21)[:-1] # A list of 20 angles between 0 and 360
            # create an animated gif (20ms between frames)
            rotanimate(ax, angles,'pca.gif',delay=200) 
            # create a movie with 10 frames per seconds and 'quality' 2000
            rotanimate(ax, angles,'pca.mp4',fps=50,bitrate=4000)
            """

    def PCAVariationFS(self, n_pca_dims: int=3, spliced_file: str="Splicec.csv", resuts_directory: str="/filler/dir"):
          # import R's utility package
          utils = importr('utils')
          #utils.install_packages('foreign')
          #utils.install_packages('pbkrtest')
          #utils.install_packages('sp')
          #utils.install_packages('rio')
          #utils.install_packages('lattice')
          #utils.install_packages('maptools')
          #utils.install_packages('Matrix')
          #utils.install_packages('nlme')
          #utils.install_packages('car')
          ggplot2 = importr('ggplot2')
          car = importr('car')
          #utils.install_packages('cluster')
          car = importr('cluster')
          #utils.install_packages('FactoMineR')
          #utils.install_packages('factoextra')
          FactoMineR = importr('FactoMineR')
          factoextra = importr('factoextra')
          pandas2ri.activate()
          robjects.globalenv['spliced_file'] = spliced_file
          robjects.globalenv['resuts_directory'] = resuts_directory
          ro.r('spliced_data <- read.csv(spliced_file, sep=",", header=TRUE)') #mm10_pberghei_ercc_mouse_
          ro.r('rownames(spliced_data) <- spliced_data$ID')
          ro.r('spliced_data <- spliced_data[,-1]')
          ro.r('pcadata <- t(spliced_data)')
          ro.r('data.pca <- PCA(pcadata, scale.unit = FALSE, ncp = 3, graph = FALSE)')
          ro.r('setwd(resuts_directory)')
          ro.r('genes_contributions <- data.frame(data.pca$var$contrib)')
          ro.r('genes_contributions_pca_1_order <- genes_contributions[order(-genes_contributions$Dim.1),]')
          ro.r('genes_contributions_pca_2_order <- genes_contributions[order(-genes_contributions$Dim.2),]')
          #ro.r('print(genes_contributions_pca_1_order[1:10,])')
          #ro.r('print(genes_contributions_pca_2_order[1:10,])')
          ro.r('min_percent = 100/nrow(spliced_data)')
          ro.r('genes_contributionspc1 = genes_contributions[genes_contributions$Dim.1 > min_percent*6,]') #pvivax dat: pc1 loose as higher % contributors going wrong way so to ensure they are kept in when pc1 is more stringent need stringent pc2 and pc3 thresholds, as pc1 increase, pc2 and 3 decrease i.e. some kind of weighting is required
          ro.r('genes_contributionspc2 = genes_contributions[genes_contributions$Dim.2 > min_percent*30,]')
          ro.r('genes_contributionspc3 = genes_contributions[genes_contributions$Dim.3 > min_percent*30,]')
          ro.r('genes_contributionspc1a2 <- unique(c(rownames(genes_contributionspc1), rownames(genes_contributionspc2), rownames(genes_contributionspc3)))')
          ro.r('print(genes_contributionspc1a2)')
          ro.r('write.table(genes_contributionspc1, file = "PC1_Genes.csv", row.names=TRUE, sep=",")')
          ro.r('write.table(genes_contributionspc2, file = "PC2_Genes.csv", row.names=TRUE, sep=",")')
          ro.r('write.table(genes_contributionspc3, file = "PC3_Genes.csv", row.names=TRUE, sep=",")')
          ro.r('write.table(genes_contributionspc1a2, file = "selected_Genes.csv", row.names=TRUE, sep=",")')
          ro.r('write.table(min_percent, file = "PCAthreshold_used.csv", row.names=TRUE, sep=",")')
          important_genes = robjects.r.genes_contributionspc1a2
          genes = list(self.ra["Genes"])
          new_genes = []
          for gene in genes:
              new_genes.append(gene.replace('"',""))
          genes = new_genes
          important_genes = list(important_genes)
          save_genes = list(set(genes).intersection(important_genes))
          indexes = []
          for gene in save_genes:
              index = genes.index(gene)
              indexes.append(index)
          new_S = self.S
          new_S = new_S[indexes, ]
          self.S = new_S
          self.S_sz = new_S
          self.S_norm = new_S
          new_U = self.U
          new_U = new_U[indexes, ]
          self.U = new_U
          self.U_sz = new_U
          self.U_norm = new_U
          genesy = self.ra["Genes"]
          genesy = genesy[indexes]
          self.ra["Genes"] = genesy

    def TransitionProbabilityScatterplots(self, filename: str="SummaryTable_Mean_SelfNotRemoved_NN5_UMAP.csv", make_gify: bool=True, Colours: str="Paired"):
        Meandata = pd.read_csv(filename)
        GroupNames = sorted(set(list(Meandata.iloc[:,1])))
        NumberGroups = len(set(list(Meandata.iloc[:,1])))
        colours = ["#CC79A7", "#009E73", "#E69F00"]
        markers = ['o', '^', 'X']
        def rotate(angle):
            ax.view_init(azim=angle)
        def init():
            return fig,
        def animate(i):
            ax.view_init(elev=(i-45)*4, azim=10)
            return fig,
        if NumberGroups == 2:
            plotnames = []
            fig, ax = plt.subplots(figsize=(8, 7)) 
            count = 0
            for groupynam in GroupNames:
                Group1 = Meandata.loc[Meandata['Group'] == groupynam]
                Group1x = Group1[[str(GroupNames[-2])]]
                Group1y = Group1[[str(GroupNames[-1])]]
                myplot = ax.scatter(x=Group1x, y=Group1y, c=colours[count], marker=markers[count])
                plotnames.append(myplot)
                count = count + 1 
            xlab = 'Transition Probability to ' + str(GroupNames[-2])     
            ylab = 'Transition Probability to ' + str(GroupNames[-1])
            plt.xlabel(xlab) 
            plt.ylabel(ylab)
            ax.legend(plotnames, list(GroupNames), fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key', prop={'size':7})
            newfilename = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "",   re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_2DScatterplot.png"
            plt.savefig(newfilename)
        elif NumberGroups == 3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            Group1 = Meandata.loc[Meandata['Group'] == GroupNames[0]]
            Group1x = Group1.iloc[:,2]
            Group1y = Group1.iloc[:,3]
            Group1z = Group1.iloc[:,4]
            Group1plot = ax.scatter(Group1x, Group1y, Group1z, c=colours[0], marker=markers[0])
            Group2 = Meandata.loc[Meandata['Group'] == GroupNames[1]]
            Group2x = Group2.iloc[:,2]
            Group2y = Group2.iloc[:,3]
            Group2z = Group2.iloc[:,4]
            Group2plot = ax.scatter(Group2x, Group2y, Group2z, c=colours[1], marker=markers[1])
            Group3 = Meandata.loc[Meandata['Group'] == GroupNames[2]]
            Group3x = Group3.iloc[:,2]
            Group3y = Group3.iloc[:,3]
            Group3z = Group3.iloc[:,4]
            Group3plot = ax.scatter(Group3x, Group3y, Group3z, c=colours[2], marker=markers[2])
            xlab = 'Transition Probability to ' + Meandata.columns[2]
            ylab = 'Transition Probability to ' + Meandata.columns[3]      
            zlab = 'Transition Probability to ' + Meandata.columns[4]
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.set_zlabel(zlab)
            ax.legend((Group1plot, Group2plot, Group3plot), list(GroupNames), fancybox=True, framealpha=1, borderpad=1, loc=(0.85,0.5), title='Key')
            newfilename = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "",   re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_3DScatterplot.png"
            plt.savefig(newfilename)
            if (make_gify == True):
                gif_name_1 = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "", re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_3DScatterploLeftRight.gif"
                gif_name_2 = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "", re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_3DScatterploUpDown.gif"
                rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
                rot_animation.save(gif_name_1, dpi=80, writer='imagemagick')
                ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
                ani_2.save(gif_name_2, dpi=80, writer='imagemagick')
        elif NumberGroups > 3:
            numGraphstomake = NumberGroups/3
            numGraphstomake = int(numGraphstomake)
            cmap2 = matplotlib.colormaps[Colours]
            manycolours = []
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3]
                manycolours.append(matplotlib.colors.rgb2hex(rgb))
            if NumberGroups%3==0:
                for i in range(0, numGraphstomake):
                    Groupnums = list(range(i+i+i,i+i+i+2+1))
                    myGroups = GroupNames[Groupnums[0]: Groupnums[-1]+1]
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    plotnames = []
                    count = 0
                    for groupynam in GroupNames:
                        Group1 = Meandata.loc[Meandata['Group'] == groupynam]
                        Group1x = Group1[[str(myGroups[0])]]
                        Group1y = Group1[[str(myGroups[1])]]
                        Group1z = Group1[[str(myGroups[2])]]
                        myplot = ax.scatter(Group1x, Group1y, Group1z, c=manycolours[count], marker='o')
                        plotnames.append(myplot)
                        count = count + 1
                    xlab = 'Transition Probability to ' + str(myGroups[0])     
                    ylab = 'Transition Probability to ' + str(myGroups[1])     
                    zlab = 'Transition Probability to ' + str(myGroups[2])
                    ax.set_xlabel(xlab)
                    ax.set_ylabel(ylab)
                    ax.set_zlabel(zlab)
                    ax.legend(list(plotnames), list(GroupNames), title='Key')
                    newfilename = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "",   re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(i+1) + "_3DScatterplot.png"
                    plt.savefig(newfilename)
                    if (make_gify == True):
                        gif_name_1 = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "", re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(i+1) + "_3DScatterploLeftRight.gif"
                        gif_name_2 = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "", re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(i+1) + "_3DScatterploUpDown.gif"
                        rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
                        rot_animation.save(gif_name_1, dpi=80, writer='imagemagick')
                        ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
                        ani_2.save(gif_name_2, dpi=80, writer='imagemagick')
            else:
                for i in range(0, numGraphstomake):
                    Groupnums = list(range(i+i+i,i+i+i+2+1))
                    myGroups = GroupNames[Groupnums[0]: Groupnums[-1]+1]
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    plotnames = []
                    count = 0
                    for groupynam in GroupNames:
                        Group1 = Meandata.loc[Meandata['Group'] == groupynam]
                        Group1x = Group1[[str(myGroups[0])]]
                        Group1y = Group1[[str(myGroups[1])]]
                        Group1z = Group1[[str(myGroups[2])]]
                        myplot = ax.scatter(Group1x, Group1y, Group1z, c=manycolours[count], marker='o')
                        plotnames.append(myplot)
                        count = count + 1
                    xlab = 'Transition Probability to ' + str(myGroups[0])     
                    ylab = 'Transition Probability to ' + str(myGroups[1])     
                    zlab = 'Transition Probability to ' + str(myGroups[2])
                    ax.set_xlabel(xlab)
                    ax.set_ylabel(ylab)
                    ax.set_zlabel(zlab)
                    ax.legend(list(plotnames), list(GroupNames), title='Key')
                    newfilename = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "",   re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(i+1) + "_3DScatterplot.png"
                    plt.savefig(newfilename)
                    if (make_gify == True):
                        gif_name_1 = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "", re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(i+1) + "_3DScatterploLeftRight.gif"
                        gif_name_2 = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "", re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(i+1) + "_3DScatterploUpDown.gif"
                        rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
                        rot_animation.save(gif_name_1, dpi=80, writer='imagemagick')
                        ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
                        ani_2.save(gif_name_2, dpi=80, writer='imagemagick')
                #2d graph of last two groups:
                plotnames = []
                fig, ax = plt.subplots(figsize=(8, 7)) 
                count = 0
                for groupynam in GroupNames:
                    Group1 = Meandata.loc[Meandata['Group'] == groupynam]
                    Group1x = Group1[[str(GroupNames[-2])]]
                    Group1y = Group1[[str(GroupNames[-1])]]
                    myplot = ax.scatter(x=Group1x, y=Group1y, c=manycolours[count], marker='o')
                    plotnames.append(myplot)
                    count = count + 1 
                xlab = 'Transition Probability to ' + str(GroupNames[-2])     
                ylab = 'Transition Probability to ' + str(GroupNames[-1])
                plt.xlabel(xlab) 
                plt.ylabel(ylab)
                ax.legend(plotnames, list(GroupNames), fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key', prop={'size':7})
                newfilename = filename.rsplit('/', 1)[0] + '/' + re.sub("SummaryTable_", "",   re.sub(".csv", "", filename.rsplit('/', 1)[1])) + "_" + str(numGraphstomake+1) + "_2DScatterplot.png"
                plt.savefig(newfilename)
        plt.close("all")

    def CorrelationCoefficentFS(self, my_clusters):
       spliced = self.S
       unspliced = self.U
       genes = self.ra["Genes"]
       time_points = self.ca["ClusterName"]
       samples = self.ca["Patients"]
       #rows genes, columns: time-point
       spliced_df = pd.DataFrame.from_records(spliced)
       spliced_df = spliced_df.transpose()
       #spliced_df.rows = samples
       spliced_df.rename({1: 2, 2: 4}, axis='index')
       spliced_df.columns = genes
       spliced_df['Time_Point'] = list(map(int,time_points))
       spliced_df['Samples'] = samples
       spliced_df = pd.DataFrame(spliced_df).set_index('Samples')
       ccs = []
       for gene in genes:
           ccs.append(spliced_df[gene].corr(spliced_df["Time_Point"])) #default is pearson
       coefficient_dataframe = pd.DataFrame({'Genes': genes,'CorrelCoeff': ccs})
       increase = pd.DataFrame(coefficient_dataframe[coefficient_dataframe['CorrelCoeff'] > float(0.95)]) #222
       decrease = pd.DataFrame(coefficient_dataframe[coefficient_dataframe['CorrelCoeff'] < float(-0.95)]) #510
       save_genes = unique(increase['Genes'].tolist() + decrease['Genes'].tolist()) #concatenate with decrease genes
       print(len(save_genes))
       print(increase)
       print(decrease)
       increase.to_csv("Transcripts_Selected_positive_correlation.csv", index=False)
       decrease.to_csv("Transcripts_Selected_negative_correlation.csv", index=False)
       genes = list(self.ra["Genes"])
       indexes = []
       for gene in save_genes:
          index = genes.index(gene)
          indexes.append(index)
       new_S = self.S
       new_S = new_S[indexes, ]
       self.S = new_S
       self.S_sz = new_S
       self.S_norm = new_S
       new_U = self.U
       new_U = new_U[indexes, ]
       self.U = new_U
       self.U_sz = new_U
       self.U_norm = new_U
       genesy = self.ra["Genes"]
       genesy = genesy[indexes]
       self.ra["Genes"] = genesy

    def basicUMAP(self, num_nei, make_gify: bool=True, distance: float=0.1, metric_method: str="euclidean", typegroup: str="integers", do3D: str="True", UMAP_3D_loc: str="/location/", ChangeShape: str="True", Colours: str="Paired"):
       if ChangeShape == "True":
           my_markers = self.ca["Markers"]
           mark = sorted(set(my_markers))
           corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
           trans2comp = umap.UMAP(n_neighbors = num_nei, n_components=2).fit(self.S_norm.T)
           self.UMAPembed = trans2comp.embedding_
           if do3D == "True":
               os.chdir(UMAP_3D_loc)
               trans = umap.UMAP(n_neighbors = num_nei, min_dist=distance, n_components=3, metric=metric_method).fit(self.S_norm.T) #min_distance defualt is 0.1, varies 0 to 0.99
               self.UMAPembed3d = trans.embedding_
               count = 0
               my_types = {}
               for the_type in mark:
                   my_types[the_type] = corresp[count]
                   count = count + 1
               markers = []
               for mar in my_markers:
                   markers.append(my_types.get(mar))
               z = randn(10)
               fig = plt.figure(figsize=(8, 6))
               ax = fig.add_subplot(111, projection='3d')
               for i in range(0, len(markers)):
                   ax.scatter(trans.embedding_[i, 0], trans.embedding_[i, 1], trans.embedding_[i, 2], marker=markers[i], color=self.colorandum[i])
               hexy_codes = []
               my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
               decoded_my_keys = []
               for item in my_keys:
                   decoded_my_keys.append(str(item,'utf-8')) #int()
               cmap2 = matplotlib.colormaps[Colours]   # PiYG
               for i in range(cmap2.N):
                   rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                   hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
               legend_elements = []
               combined = list(zip(decoded_my_keys, hexy_codes))
               combined = sorted(combined, key=itemgetter(0))
               if (typegroup != "Strings"):
                   timepoints = [int(i[0]) for i in combined]
                   colys = [i[1] for i in combined]
                   tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                   tempDF = tempDF.sort_values('Group')
                   groupys2 = tempDF['Group'].tolist()
                   cols2 = tempDF['colours'].tolist()
                   combined = list(zip(groupys2, cols2))
               decoded_my_keys = [i[0] for i in combined]
               hexy_codes = [i[1] for i in combined]
               for sample_type in range(0, len(decoded_my_keys)):
                   to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                   legend_elements.append(to_add)
               for sample_type in range(0, len(mark)):
                   to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                   legend_elements.append(to_add)
               ax.legend(handles=legend_elements, title='Key') #Timepont \n(hours)
               ax.set_xlabel('UMAP1')
               ax.set_ylabel('UMAP2')
               ax.set_zlabel('UMAP3')
               figure_name3D = 'UMAP3D_NN' + str(num_nei) + '.png'
               plt.savefig(figure_name3D)
               #plt.title('UMAP projection', fontsize=10);
               #plt.show()
               if (make_gify == True):
                   gify_name_1 = 'UMAP_rotation_NN' + str(num_nei) + '.gif'
                   gify_name_2 = 'UMAP_rotation_up_down_NN' + str(num_nei) + '.gif'
                   def rotate(angle):
                       ax.view_init(azim=angle)
                   rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
                   rot_animation.save(gify_name_1, dpi=80, writer='imagemagick')
                   def init():
                       return fig,
                   def animate(i):
                       ax.view_init(elev=(i-45)*4, azim=10)
                       return fig,
                   ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
                   ani_2.save(gify_name_2, dpi=80, writer='imagemagick')
               plt.close()
       if ChangeShape == "False":
           trans2comp = umap.UMAP(n_neighbors = num_nei, n_components=2).fit(self.S_norm.T)
           self.UMAPembed = trans2comp.embedding_
           if do3D == "True":
               os.chdir(UMAP_3D_loc)
               trans = umap.UMAP(n_neighbors = num_nei, min_dist=distance, n_components=3, metric=metric_method).fit(self.S_norm.T) #min_distance defualt is 0.1, varies 0 to 0.99
               self.UMAPembed3d = trans.embedding_
               z = randn(10)
               fig = plt.figure(figsize=(8, 6))
               ax = fig.add_subplot(111, projection='3d')
               ax.scatter(trans.embedding_[:, 0], trans.embedding_[:, 1], trans.embedding_[:, 2], marker="o", color=self.colorandum)
               #for i in range(0, len(self.colorandum)):
                   #ax.scatter(trans.embedding_[i, 0], trans.embedding_[i, 1], trans.embedding_[i, 2], marker="o", color=self.colorandum[i])
               hexy_codes = []
               my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
               decoded_my_keys = []
               for item in my_keys:
                   decoded_my_keys.append(str(item,'utf-8')) #int()
               cmap2 = matplotlib.colormaps[Colours]    # PiYG
               for i in range(cmap2.N):
                   rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                   hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
               legend_elements = []
               combined = list(zip(decoded_my_keys, hexy_codes))
               combined = sorted(combined, key=itemgetter(0))
               if (typegroup != "Strings"):
                   timepoints = [int(i[0]) for i in combined]
                   colys = [i[1] for i in combined]
                   tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                   tempDF = tempDF.sort_values('Group')
                   groupys2 = tempDF['Group'].tolist()
                   cols2 = tempDF['colours'].tolist()
                   combined = list(zip(groupys2, cols2))
               decoded_my_keys = [i[0] for i in combined]
               hexy_codes = [i[1] for i in combined]
               for sample_type in range(0, len(decoded_my_keys)):
                   to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                   legend_elements.append(to_add)
               ax.legend(handles=legend_elements, title='Key') #Timepont \n(hours)
               ax.set_xlabel('UMAP1')
               ax.set_ylabel('UMAP2')
               ax.set_zlabel('UMAP3')
               figure_name3D = 'UMAP3D_NN' + str(num_nei) + '.png'
               plt.savefig(figure_name3D)
               #plt.title('UMAP projection', fontsize=10);
               #plt.show()
               if (make_gify == True):
                   gify_name_1 = 'UMAP_rotation_NN' + str(num_nei) + '.gif'
                   gify_name_2 = 'UMAP_rotation_up_down_NN' + str(num_nei) + '.gif'
                   def rotate(angle):
                       ax.view_init(azim=angle)
                   rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
                   rot_animation.save(gify_name_1, dpi=80, writer='imagemagick')
                   def init():
                       return fig,
                   def animate(i):
                       ax.view_init(elev=(i-45)*4, azim=10)
                       return fig,
                   ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
                   ani_2.save(gify_name_2, dpi=80, writer='imagemagick')
               plt.close()

    def estimate_transition_probUMAP(self, hidim: str="Sx_sz", embed: str="UMAPembed", transform: str="sqrt",ndims: int=None, n_sight: int=None, psc: float=None,knn_random: bool=False, sampled_fraction: float=0.3,sampling_probs: Tuple[float, float]=(0.5, 0.1), max_dist_embed: float=None, n_jobs: int=4, threads: int=None, calculate_randomized: bool=True, random_seed: int=15071990, **kwargs) -> None:
      numba_random_seed(random_seed)
      self.which_hidim = hidim #spliced states before timepoint same calculation used to calculate after just with different names for the same thing this is: self.Sx_sz
      if "n_neighbors" in kwargs:
          n_neighbors = kwargs.pop("n_neighbors")
          if len(kwargs) > 0:
              logging.warning(f"keyword arguments were passed but could not be interpreted {kwargs}")
      else:
          n_neighbors = None
      if n_sight is None and n_neighbors is None:
          n_neighbors = int(self.S.shape[1] / 5)
      if (n_sight is not None) and (n_neighbors is not None) and n_neighbors != n_sight:
          raise ValueError("n_sight and n_neighbors are different names for the same parameter, they cannot be set differently")
      if n_sight is not None and n_neighbors is None:
          n_neighbors = n_sight
      if psc is None:
          if transform == "log" or transform == "logratio":
              psc = 1.
          elif transform == "sqrt":
              psc = 1e-10  # for numerical stablity
          else:  # transform == "linear":
              psc = 0
      if knn_random:
          np.random.seed(random_seed)
          self.corr_calc = "knn_random"
          if "pcs" in hidim:  # sic
              hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
              hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
          else:
              if ndims is not None:
                  raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
              hi_dim = getattr(self, hidim)  # [:, :ndims] hi_dim is the
              hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
              if calculate_randomized:
                  self.delta_S_rndm = np.copy(self.delta_S)
                  permute_rows_nsign(self.delta_S_rndm)
                  hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
          embedding = getattr(self, embed) #
          self.embedding_UMAP = embedding
          logging.debug("Calculate KNN in the embedding space")
          nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
          nn.fit(embedding)  # NOTE should support knn in high dimensions
          self.embedding_UMAP_knn = nn.kneighbors_graph(mode="connectivity")
          # Pick random neighbours and prune the rest
          neigh_ixs = self.embedding_UMAP_knn.indices.reshape((-1, n_neighbors))
          p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
          p = p / p.sum()
          # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
          # resulting self.embedding_UMAPf with different number of nonzero entry per row.
          # Not updated yet not to break previous analyses
          # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
          # I change it here since I am doing some breaking changes
          sampling_ixs = np.stack((np.random.choice(neigh_ixs.shape[1],
                                                  size=(int(sampled_fraction * (n_neighbors)),),
                                                  replace=False,
                                                  p=p) for i in range(neigh_ixs.shape[0])), 0)
          self.embedding_UMAP = sampling_ixs
          neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
          nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
          self.embedding_UMAP_knn = sparse.csr_matrix((np.ones(nonzero),
                                                neigh_ixs.ravel(),
                                                np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                               shape=(neigh_ixs.shape[0],
                                                      neigh_ixs.shape[0]))
          logging.debug(f"Correlation Calculation '{self.corr_calc}'")
          if transform == "log":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
          elif transform == "logratio":
              log2hidim = np.log2(hi_dim + psc)
              delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
              self.corrcoef_UMAP = colDeltaCorpartial(log2hidim, delta_hi_dim, neigh_ixs, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + psc) - log2hidim
                  self.corrcoef_UMAP_random = colDeltaCorpartial(log2hidim, delta_hi_dim_rndm, neigh_ixs, threads=threads)
          elif transform == "linear":
              self.corrcoef_UMAP = colDeltaCorpartial(hi_dim, hi_dim_t - hi_dim, neigh_ixs, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  self.corrcoef_UMAP_random = colDeltaCorpartial(hi_dim, hi_dim_t_rndm - hi_dim, neigh_ixs, threads=threads)
          elif transform == "sqrt":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
          else:
              raise NotImplementedError(f"transform={transform} is not a valid parameter")
          np.fill_diagonal(self.embedding_UMAP, 0)
          if np.any(np.isnan(self.embedding_UMAP)):
              self.embedding_UMAP[np.isnan(self.embedding_UMAP)] = 1
              logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
          if calculate_randomized:
              np.fill_diagonal(self.corrcoef_UMAP_random, 0)
              if np.any(np.isnan(self.corrcoef_UMAP_random)):
                  self.corrcoef_UMAP_random[np.isnan(self.corrcoef_UMAP_random)] = 1
                  logging.warning("Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
          logging.debug(f"Done Correlation Calculation")
      else:
          self.corr_calc = "full"
          if "pcs" in hidim:  # sic
              hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
              hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
          else: #this route:
              if ndims is not None:
                  raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
              hi_dim = getattr(self, hidim)  # [:, :ndims]
              hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
              if calculate_randomized:
                  self.delta_S_rndm = np.copy(self.delta_S)
                  permute_rows_nsign(self.delta_S_rndm)
                  hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
          embedding = getattr(self, embed)
          self.embedding_UMAP = embedding
          logging.debug("Calculate KNN in the embedding space")
          nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
          nn.fit(embedding)
          self.embedding_UMAP_knn = nn.kneighbors_graph(mode="connectivity") #def of embedding_knn
          #print(transform) #sqrt
          logging.debug("Correlation Calculation 'full'")
          if transform == "log":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
          elif transform == "logratio":
              log2hidim = np.log2(hi_dim + psc)
              delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
              self.corrcoef_UMAP = colDeltaCor(log2hidim, delta_hi_dim, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + 1) - log2hidim
                  self.corrcoef_UMAP_random = colDeltaCor(log2hidim, delta_hi_dim_rndm, threads=threads)
          elif transform == "linear":
              self.corrcoef_UMAP = colDeltaCor(hi_dim, hi_dim_t - hi_dim, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  self.corrcoef_UMAP_random = colDeltaCor(hi_dim, hi_dim_t_rndm - hi_dim, threads=threads)
          elif transform == "sqrt":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
          else:
              raise NotImplementedError(f"transform={transform} is not a valid parameter")
          np.fill_diagonal(self.corrcoef_UMAP, 0)
          if calculate_randomized:
              np.fill_diagonal(self.corrcoef_UMAP_random, 0)
          colnames_cells = list(self.ca["Patients"]) #shape of all: (rows=genes, cols=samples)
          rownames_genes = list(self.ra["Genes"])
          tp = np.array(hi_dim)
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #emat= hi_dim #expression
          file_name = "emat_UMAP.csv"
          tp_df.to_csv(file_name, sep=",", header=True) 
          sqrtdmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim)
          tp = np.array(sqrtdmat)
          file_name = "sqrt_dmat_UMAP.csv"
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #dmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim) #dmat are the velocities
          tp_df.to_csv(file_name, sep=",", header=True)  
          alt_dmat = hi_dim_t - hi_dim
          tp = np.array(alt_dmat)
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  # alt dmat: hi_dim_t - hi_dim (for linear)
          file_name = "linear_dmat_UMAP.csv"
          tp_df.to_csv(file_name, sep=",", header=True)  
          tp = np.array(self.corrcoef_UMAP)
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)  # raw cc values withot any alt
          file_name = "CorrelationCoefficients_Raw_UMAP.csv"
          tp_df.to_csv(file_name, sep=",", header=True) 

    def calculate_embedding_shiftUMAP(self, sigma_corr: float=0.05, expression_scaling: bool=False, Myneighbours: float=None, scaling_penalty: float=1., metadataToUse: str="metadata.csv", pathToTransitionProbabilityCode: str="/", grouptype: str="Strings", secondgroup: str="metadataVersion.csv", manyfeatures: str="True", grouptype2: str="Strings", Randomsampling: str="False", extraNNvalue: int=5, DE_location: str="/location/", selfCC_folder: str="/location/", VecX_location: str="/location/", VecY_loc: str="/location/", TP_folder: str="/location/", UMAP_2D_GroupLevelTP: str="/location/", DoSecondPred: str="no", UserThreshold: float=0.5) -> None: #Use the transition probability to project the velocity direction on the embedding
        # Kernel evaluation
        logging.debug("Calculate transition probability")
        if self.corr_calc == "full" or self.corr_calc == "knn_random":
            # NOTE maybe sparse matrix here are slower than dense
            # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
            self.transition_prob_UMAP = np.exp(self.corrcoef_UMAP / sigma_corr) * self.embedding_UMAP_knn.toarray()  # naive, the exponential of the correlation coefficient/kernal scaling, AG aadded to.array()
            self.transition_prob_UMAP /= self.transition_prob_UMAP.sum(1)[:, None]
            if hasattr(self, "corrcoef_random"):
                logging.debug("Calculate transition probability for negative control")
                self.transition_prob_UMAP_random_UMAP = np.exp(self.self.corrcoef_random_UMAP / sigma_corr) * self.embedding_UMAP_knn.toarray()  # naive AG added to.array()
                self.transition_prob_UMAP_random_UMAP /= self.transition_prob_UMAP_random_UMAP.sum(1)[:, None]
            unitary_vectors = self.embedding_UMAP.T[:, None, :] - self.embedding_UMAP.T[:, :, None]  # shape (2,ncells,ncells) #for each cells coordinates - row1: cell1-cell1, cell1-cell2, cell1-cell3 cell1-cell4, cell1-cell1, cell2-cell1, cell3-cell1 etc.
            with np.errstate(divide='ignore', invalid='ignore'):
                unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2
                np.fill_diagonal(unitary_vectors[0, ...], 0)  # fix nans replace NAN with 0
                np.fill_diagonal(unitary_vectors[1, ...], 0) #embedding comes from the tSNE
            self.MY_ARROWS_umap = (self.transition_prob_UMAP * unitary_vectors).sum(2) #.A means change the data type from a matrix to an array #unit vector*transition probability, below: transition probab (no. samples x no.samples) is multiplied by each of the three (no.samplesxno.samples) matrices in the unitary vector individually. The result is then summed for each xm y and z across all number of samples rows for each number of sample columns
            self.MY_ARROWS_umap -= (self.embedding_UMAP_knn.toarray() * unitary_vectors).sum(2) / self.embedding_UMAP_knn.sum(1).toarray().T #AG aadded to.array()
            self.MY_ARROWS_umap = self.MY_ARROWS_umap.T #transposes the vector
            if expression_scaling:
                hi_dim = getattr(self, self.which_hidim)
                estim_delta = hi_dim.dot(self.transition_prob_UMAP.T) - hi_dim.dot((self.embedding_UMAP_knn.toarray() / self.embedding_UMAP_knn.sum(1).toarray()).T) #this is the same problem as above, AG aadded to.array()
                cos_proj = (self.delta_S * estim_delta).sum(0) / np.sqrt((estim_delta**2).sum(0))
                self.scaling_UMAP = np.clip(cos_proj / scaling_penalty, 0, 1)
                self.MY_ARROWS_umap = self.MY_ARROWS_umap * self.scaling_UMAP[:, None]
            if hasattr(self, "corrcoef_random"):
                self.MY_ARROWS_umap_random = (self.transition_prob_UMAP_random_UMAP * unitary_vectors).sum(2)
                self.MY_ARROWS_umap_random -= (self.embedding_UMAP_knn.toarray() * unitary_vectors).sum(2) / self.embedding_UMAP_knn.sum(1).toarray().T #AG aadded to.array()
                self.MY_ARROWS_umap_random = self.MY_ARROWS_umap_random.T
                if expression_scaling:
                    estim_delta_rndm = hi_dim.dot(self.transition_prob_UMAP_random_UMAP.T) - hi_dim.dot((self.embedding_UMAP_knn.toarray() / self.embedding_UMAP_knn.sum(1).toarray()).T) #AG aadded to.array()
                    cos_proj_rndm = (self.delta_S_rndm * estim_delta_rndm).sum(0) / np.sqrt((estim_delta_rndm**2).sum(0))
                    self.scaling_UMAP_rndm = np.clip(cos_proj_rndm / scaling_penalty, 0, 1)
                    self.MY_ARROWS_umap_random = self.MY_ARROWS_umap_random * self.scaling_UMAP_rndm[:, None]
        else:
            # NOTE should implement a version with cython
            raise NotImplementedError(f"Weird value self.corr_calc={self.corr_calc}")

        colnames_cells = list(self.ca["Patients"])
        genes_embedding = np.array(self.MY_ARROWS_umap)
        spliced_df = pd.DataFrame(genes_embedding, columns=["x", "y"], index=colnames_cells)
        file_name = "Delta_embedding_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP.csv"
        os.chdir(DE_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(self.corrcoef_UMAP)
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "selfCC_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP.csv"
        os.chdir(selfCC_folder)
        spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[0])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecX_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP.csv"
        os.chdir(VecX_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[1])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecY_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP.csv"
        os.chdir(VecY_loc)
        spliced_df.to_csv(file_name, sep=",", header=True)
        #print(self.embedding_knn) just a bunch of 1.0s
        tp = np.array(self.transition_prob_UMAP) #number of samples x number of samples
        tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)
        file_name = "Transition_Probability_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP.csv"
        os.chdir(TP_folder)
        tp_df.to_csv(file_name, sep=",", header=True)
        #####RUN TP code here################################
        #note: the following requires these R packages to be installed:  'gtools', 'reshape2', ggplot2', grid' and 'tidyverse'
        #to code to install them has been commented out but is below:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        utils = importr('utils')
        #utils.install_packages('gtools')
        #utils.install_packages('reshape2')
        #utils.install_packages('ggplot2')        
        #utils.install_packages('grid')
        #utils.install_packages('tidyverse')
        gtools = importr("gtools")
        reshape2 = importr("reshape2")
        ggplot2 = importr("ggplot2")
        grid = importr("grid")
        tidyverse = importr("tidyverse")
        path=pathToTransitionProbabilityCode
        version = "SampleGroup1"
        def RunTransitionProbabilityCode(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold):
            r=ro.r
            r.source(path+"TransitionProbability.R")
            p=r.TransitionProbability(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold)
            return p
        #a=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, metadataToUse, UMAP_2D_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold) 
        b=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, metadataToUse, UMAP_2D_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold)
        #if second metadata feature type == True
        if manyfeatures == "True" and DoSecondPred == "yes":
            version = "SampleGroup2" #secondgroup
            #c=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, secondgroup, UMAP_2D_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold) 
            d=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, secondgroup, UMAP_2D_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold)

    def UMAPwithArrows2D(self, UMAPembedding: str="UMAPembed", UMAP_Arrows: str="MY_ARROWS_umap", numberNeighb: float=2., typegroup: str="integers", extraNNvalue: int=5, ChangeShape: str="True", Colours: str="Paired"):
       #trans = umap.UMAP(n_neighbors = num_nei, n_components=3).fit(self.S_norm.T)
       UMAP_samplesembedding = getattr(self, UMAPembedding)
       self.UMAPsampleembedding = UMAP_samplesembedding
       umap_arrows_embedding_stuff = getattr(self, UMAP_Arrows)
       self.arrowUMAPembed = umap_arrows_embedding_stuff
       if ChangeShape == "True":
           my_markers = self.ca["Markers"]
           mark = sorted(set(my_markers))
           corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
           count = 0
           my_types = {}
           for the_type in mark:
               my_types[the_type] = corresp[count]
               count = count + 1
           markers = []
           for mar in my_markers:
               markers.append(my_types.get(mar))
           z = randn(10)
           fig = plt.figure(figsize=(15, 10))
           quiver_scale = 60
           for i in range(0, len(markers)):
               plt.scatter(self.UMAPsampleembedding[i, 0], self.UMAPsampleembedding[i, 1], marker=markers[i], color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
               quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
               plt.quiver(self.UMAPsampleembedding[i, 0], self.UMAPsampleembedding[i, 1], self.arrowUMAPembed[i, 0], self.arrowUMAPembed[i, 1], scale=quiver_scale, **quiver_kwargs)
           hexy_codes = []
           my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
           decoded_my_keys = []
           for item in my_keys:
               decoded_my_keys.append(str(item,'utf-8')) #int()
           cmap2 = matplotlib.colormaps[Colours]    # PiYG
           for i in range(cmap2.N):
               rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
               hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
           legend_elements = []
           combined = list(zip(decoded_my_keys, hexy_codes))
           combined = sorted(combined, key=itemgetter(0))
           if (typegroup != "Strings"):
               timepoints = [int(i[0]) for i in combined]
               colys = [i[1] for i in combined]
               tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
               tempDF = tempDF.sort_values('Group')
               groupys2 = tempDF['Group'].tolist()
               cols2 = tempDF['colours'].tolist()
               combined = list(zip(groupys2, cols2))
           decoded_my_keys = [i[0] for i in combined]
           hexy_codes = [i[1] for i in combined]
           for sample_type in range(0, len(decoded_my_keys)):
               to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
               legend_elements.append(to_add)
           for sample_type in range(0, len(mark)):
               to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
               legend_elements.append(to_add)
       if ChangeShape == "False":
           z = randn(10)
           fig = plt.figure(figsize=(15, 10))
           quiver_scale = 60
           plt.scatter(self.UMAPsampleembedding[:, 0], self.UMAPsampleembedding[:, 1], marker="o", color=self.colorandum, alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
           quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
           plt.quiver(self.UMAPsampleembedding[:, 0], self.UMAPsampleembedding[:, 1], self.arrowUMAPembed[:, 0], self.arrowUMAPembed[:, 1], scale=quiver_scale, **quiver_kwargs)
           #for i in range(0, len(self.colorandum)):
           #    plt.scatter(self.UMAPsampleembedding[i, 0], self.UMAPsampleembedding[i, 1], marker="o", color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
           #    quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
           #    plt.quiver(self.UMAPsampleembedding[i, 0], self.UMAPsampleembedding[i, 1], self.arrowUMAPembed[i, 0], self.arrowUMAPembed[i, 1], scale=quiver_scale, **quiver_kwargs)
           hexy_codes = []
           my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
           decoded_my_keys = []
           for item in my_keys:
               decoded_my_keys.append(str(item,'utf-8')) #int()
           cmap2 = matplotlib.colormaps[Colours]    # PiYG
           for i in range(cmap2.N):
               rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
               hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
           legend_elements = []
           combined = list(zip(decoded_my_keys, hexy_codes))
           combined = sorted(combined, key=itemgetter(0))
           if (typegroup != "Strings"):
               timepoints = [int(i[0]) for i in combined]
               colys = [i[1] for i in combined]
               tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
               tempDF = tempDF.sort_values('Group')
               groupys2 = tempDF['Group'].tolist()
               cols2 = tempDF['colours'].tolist()
               combined = list(zip(groupys2, cols2))
           decoded_my_keys = [i[0] for i in combined]
           hexy_codes = [i[1] for i in combined]
           for sample_type in range(0, len(decoded_my_keys)):
               to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
               legend_elements.append(to_add)
       plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
       plt.xlabel('UMAP1')
       plt.ylabel('UMAP2')
       figure_name = 'UMAP2D_Arrow_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
       plt.savefig(figure_name)
       #plt.title('UMAP projection', fontsize=10);
       #plt.show()
       plt.close()

    def calculate_embedding_shiftUMAP3D(self, sigma_corr: float=0.05, expression_scaling: bool=False, Myneighbours: float=None, scaling_penalty: float=1., metadataToUse: str="metadata.csv", pathToTransitionProbabilityCode: str="/", grouptype: str="Strings", secondgroup: str="metadataVersion.csv", manyfeatures: str="True", grouptype2: str="Strings", Randomsampling: str="False", extraNNvalue: int=5, DE_location: str="/location/", selfCC_folder: str="/location/", VecX_location: str="/location/", VecY_loc: str="/location/", VecZ_loc: str="/location/", TP_folder: str="/location/", UMAP_3D_GroupLevelTP: str="/location/", DoSecondPred: str="no", UserThreshold: float=0.5) -> None: #Use the transition probability to project the velocity direction on the embedding
        # Kernel evaluation
        logging.debug("Calculate transition probability")
        if self.corr_calc == "full" or self.corr_calc == "knn_random": #thit way:
            # NOTE maybe sparse matrix here are slower than dense
            # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
            self.transition_prob_UMAP_3d = np.exp(self.corrcoef_UMAP_3d / sigma_corr) * self.embedding_UMAP_knn_3d.A  # naive, the exponential of the correlation coefficient/kernal scaling
            self.transition_prob_UMAP_3d /= self.transition_prob_UMAP_3d.sum(1)[:, None]
            if hasattr(self, "corrcoef_random"):
                logging.debug("Calculate transition probability for negative control")
                self.transition_prob_UMAP_random_UMAP_3d = np.exp(self.self.corrcoef_random_UMAP_3d / sigma_corr) * self.embedding_UMAP_knn_3d.A  # naive
                self.transition_prob_UMAP_random_UMAP_3d /= self.transition_prob_UMAP_random_UMAP_3d.sum(1)[:, None]
                #print(len(embedding_UMAP_knn_3d))
                #print(len(embedding_UMAP_knn_3d[0]))
            unitary_vectors = self.embedding_UMAP_3d.T[:, None, :] - self.embedding_UMAP_3d.T[:, :, None]  # shape (3,ncells,ncells) #for each cells coordinates - row1: cell1-cell1, cell2-cell1, cell3-cell1 cell4-cell1
            with np.errstate(divide='ignore', invalid='ignore'):
                unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2, divides by the vector normed version, 15x15
                np.fill_diagonal(unitary_vectors[0, ...], 0)  # fix nans replace NAN with 0
                np.fill_diagonal(unitary_vectors[1, ...], 0) #embedding comes from the UMAP
                np.fill_diagonal(unitary_vectors[2, ...], 0) #embedding comes from the UMAP
            self.MY_ARROWS_umap_3d = (self.transition_prob_UMAP_3d * unitary_vectors).sum(2) #.A means change the data type from a matrix to an array #unit vector*transition probability
            self.MY_ARROWS_umap_3d -= (self.embedding_UMAP_knn_3d.A * unitary_vectors).sum(2) / self.embedding_UMAP_knn_3d.sum(1).A.T
            self.MY_ARROWS_umap_3d = self.MY_ARROWS_umap_3d.T #transposes the vector
            if expression_scaling:
                hi_dim = getattr(self, self.which_hidim)
                estim_delta = hi_dim.dot(self.transition_prob_UMAP_3d.T) - hi_dim.dot((self.embedding_UMAP_knn_3d.A / self.embedding_UMAP_knn_3d.sum(1).A).T)
                cos_proj = (self.delta_S * estim_delta).sum(0) / np.sqrt((estim_delta**2).sum(0))
                self.scaling_UMAP_3d = np.clip(cos_proj / scaling_penalty, 0, 1)
                self.MY_ARROWS_umap_3d = self.MY_ARROWS_umap_3d * self.scaling_UMAP_3d[:, None]
            if hasattr(self, "corrcoef_random"):
                self.MY_ARROWS_umap_random_3d = (self.transition_prob_UMAP_random_UMAP_3d * unitary_vectors).sum(2)
                self.MY_ARROWS_umap_random_3d -= (self.embedding_UMAP_knn_3d.A * unitary_vectors).sum(2) / self.embedding_UMAP_knn_3d.sum(1).A.T
                self.MY_ARROWS_umap_random_3d = self.MY_ARROWS_umap_random_3d.T
                if expression_scaling:
                    estim_delta_rndm = hi_dim.dot(self.transition_prob_UMAP_random_UMAP_3d.T) - hi_dim.dot((self.embedding_UMAP_knn_3d.A / self.embedding_UMAP_knn_3d.sum(1).A).T)
                    cos_proj_rndm = (self.delta_S_rndm * estim_delta_rndm).sum(0) / np.sqrt((estim_delta_rndm**2).sum(0))
                    self.scaling_UMAP_rndm_3d = np.clip(cos_proj_rndm / scaling_penalty, 0, 1)
                    self.MY_ARROWS_umap_random_3d = self.MY_ARROWS_umap_random_3d * self.scaling_UMAP_rndm_3d[:, None]
        else:
            # NOTE should implement a version with cython
            raise NotImplementedError(f"Weird value self.corr_calc={self.corr_calc}")
        colnames_cells = list(self.ca["Patients"])
        genes_embedding = np.array(self.MY_ARROWS_umap_3d)
        spliced_df = pd.DataFrame(genes_embedding, columns=["x", "y", "z"], index=colnames_cells)
        file_name = "Delta_embedding_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP3D.csv"
        os.chdir(DE_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(self.corrcoef_UMAP_3d)
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "selfCC_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP3D.csv"
        os.chdir(selfCC_folder)
        spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[0])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecX_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP3D.csv"
        os.chdir(VecX_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[1])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecY_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP3D.csv"
        os.chdir(VecY_loc)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[2])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecZ_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP3D.csv"
        os.chdir(VecZ_loc)
        #spliced_df.to_csv(file_name, sep=",", header=True)        
        tp = np.array(self.transition_prob_UMAP_3d) #number of samples x number of samples
        tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)
        file_name = "Transition_Probability_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_UMAP3D.csv"
        os.chdir(TP_folder)
        tp_df.to_csv(file_name, sep=",", header=True)
        #####RUN TP code here################################
        #note: the following requires these R packages to be installed:  'gtools', 'reshape2', ggplot2', grid' and 'tidyverse'
        #to code to install them has been commented out but is below:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        utils = importr('utils')
        #utils.install_packages('gtools')
        #utils.install_packages('reshape2')
        #utils.install_packages('ggplot2')        
        #utils.install_packages('grid')
        #utils.install_packages('tidyverse')
        gtools = importr("gtools")
        reshape2 = importr("reshape2")
        ggplot2 = importr("ggplot2")
        grid = importr("grid")
        tidyverse = importr("tidyverse")
        path=pathToTransitionProbabilityCode
        version = "SampleGroup1"
        def RunTransitionProbabilityCode(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold):
            r=ro.r
            r.source(path+"TransitionProbability.R")
            p=r.TransitionProbability(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold)
            return p
        #a=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, metadataToUse, UMAP_3D_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold) 
        b=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, metadataToUse, UMAP_3D_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold)
        #if second metadata feature type == True
        if manyfeatures == "True" and DoSecondPred == "yes":
            version = "SampleGroup2" #secondgroup
            #c=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, secondgroup, UMAP_3D_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold) 
            d=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, secondgroup, UMAP_3D_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold)

    def estimate_transition_probUMAP3D(self, hidim: str="Sx_sz", embed: str="UMAPembed3d", transform: str="sqrt",ndims: int=None, n_sight: int=None, psc: float=None,knn_random: bool=False, sampled_fraction: float=0.3, sampling_probs: Tuple[float, float]=(0.5, 0.1), max_dist_embed: float=None, n_jobs: int=4, threads: int=None, calculate_randomized: bool=True, random_seed: int=15071990, **kwargs) -> None:
      numba_random_seed(random_seed)
      self.which_hidim = hidim #spliced states before timepoint same calculation used to calculate after just with different names for the same thing this is: self.Sx_sz
      if "n_neighbors" in kwargs:
          n_neighbors = kwargs.pop("n_neighbors")
          if len(kwargs) > 0:
              logging.warning(f"keyword arguments were passed but could not be interpreted {kwargs}")
      else:
          n_neighbors = None
      if n_sight is None and n_neighbors is None:
          n_neighbors = int(self.S.shape[1] / 5)
      if (n_sight is not None) and (n_neighbors is not None) and n_neighbors != n_sight:
          raise ValueError("n_sight and n_neighbors are different names for the same parameter, they cannot be set differently")
      if n_sight is not None and n_neighbors is None:
          n_neighbors = n_sight
      if psc is None:
          if transform == "log" or transform == "logratio":
              psc = 1.
          elif transform == "sqrt":
              psc = 1e-10  # for numerical stablity
          else:  # transform == "linear":
              psc = 0
      if knn_random:
          np.random.seed(random_seed)
          self.corr_calc = "knn_random"
          if "pcs" in hidim:  # sic
              hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
              hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
          else:
              if ndims is not None:
                  raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
              hi_dim = getattr(self, hidim)  # [:, :ndims] hi_dim is the
              hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
              if calculate_randomized:
                  self.delta_S_rndm = np.copy(self.delta_S)
                  permute_rows_nsign(self.delta_S_rndm)
                  hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
          embedding = getattr(self, embed)
          self.embedding_UMAP_3d = embedding
          logging.debug("Calculate KNN in the embedding space")
          nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
          nn.fit(embedding)  # NOTE should support knn in high dimensions
          self.embedding_UMAP_knn_3d = nn.kneighbors_graph(mode="connectivity")
          # Pick random neighbours and prune the rest
          neigh_ixs = self.embedding_UMAP_knn_3d.indices.reshape((-1, n_neighbors))
          p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
          p = p / p.sum()
          # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
          # resulting self.embedding_UMAPf with different number of nonzero entry per row.
          # Not updated yet not to break previous analyses
          # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
          # I change it here since I am doing some breaking changes
          sampling_ixs = np.stack((np.random.choice(neigh_ixs.shape[1],
                                                  size=(int(sampled_fraction * (n_neighbors)),),
                                                  replace=False,
                                                  p=p) for i in range(neigh_ixs.shape[0])), 0)
          self.embedding_umap_3d = sampling_ixs
          neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
          nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
          self.embedding_UMAP_knn_3d = sparse.csr_matrix((np.ones(nonzero),
                                                neigh_ixs.ravel(),
                                                np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                               shape=(neigh_ixs.shape[0],
                                                      neigh_ixs.shape[0]))
          logging.debug(f"Correlation Calculation '{self.corr_calc}'")
          if transform == "log":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP_3d = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random_3d = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
          elif transform == "logratio":
              log2hidim = np.log2(hi_dim + psc)
              delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
              self.corrcoef_UMAP_3d = colDeltaCorpartial(log2hidim, delta_hi_dim, neigh_ixs, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + psc) - log2hidim
                  self.corrcoef_UMAP_random_3d = colDeltaCorpartial(log2hidim, delta_hi_dim_rndm, neigh_ixs, threads=threads)
          elif transform == "linear":
              self.corrcoef_UMAP_3d = colDeltaCorpartial(hi_dim, hi_dim_t - hi_dim, neigh_ixs, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  self.corrcoef_UMAP_random_3d = colDeltaCorpartial(hi_dim, hi_dim_t_rndm - hi_dim, neigh_ixs, threads=threads)
          elif transform == "sqrt":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP_3d = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random_3d = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
          else:
              raise NotImplementedError(f"transform={transform} is not a valid parameter")
          np.fill_diagonal(self.embedding_UMAP_3d, 0)
          if np.any(np.isnan(self.embedding_UMAP_3d)):
              self.embedding_UMAP_3d[np.isnan(self.embedding_UMAP_3d)] = 1
              logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
          if calculate_randomized:
              np.fill_diagonal(self.corrcoef_UMAP_random_3d, 0)
              if np.any(np.isnan(self.corrcoef_UMAP_random_3d)):
                  self.corrcoef_UMAP_random_3d[np.isnan(self.corrcoef_UMAP_random_3d)] = 1
                  logging.warning("Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
          logging.debug(f"Done Correlation Calculation")
      else:
          self.corr_calc = "full"
          if "pcs" in hidim:  # sic
              hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
              hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
          else: #this route:
              if ndims is not None:
                  raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
              hi_dim = getattr(self, hidim)  # [:, :ndims]
              hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
              if calculate_randomized:
                  self.delta_S_rndm = np.copy(self.delta_S)
                  permute_rows_nsign(self.delta_S_rndm)
                  hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
          embedding = getattr(self, embed)
          self.embedding_UMAP_3d = embedding
          logging.debug("Calculate KNN in the embedding space")
          nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs) #auto algorithm means it will tyr and determine the best method itself.000
          nn.fit(embedding)  # NOTE should support knn in high dimensions
          self.embedding_UMAP_knn_3d = nn.kneighbors_graph(mode="connectivity") #def of embedding_knn
          #print(transform) #sqrt
          logging.debug("Correlation Calculation 'full'")
          if transform == "log":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP_3d = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random_3d = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
          elif transform == "logratio":
              log2hidim = np.log2(hi_dim + psc)
              delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
              self.corrcoef_UMAP_3d = colDeltaCor(log2hidim, delta_hi_dim, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + 1) - log2hidim
                  self.corrcoef_UMAP_random_3d = colDeltaCor(log2hidim, delta_hi_dim_rndm, threads=threads)
          elif transform == "linear":
              self.corrcoef_UMAP_3d = colDeltaCor(hi_dim, hi_dim_t - hi_dim, threads=threads)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  self.corrcoef_UMAP_random_3d = colDeltaCor(hi_dim, hi_dim_t_rndm - hi_dim, threads=threads)
          elif transform == "sqrt":
              delta_hi_dim = hi_dim_t - hi_dim
              self.corrcoef_UMAP_3d = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
              if calculate_randomized:
                  logging.debug(f"Correlation Calculation for negative control")
                  delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                  self.corrcoef_UMAP_random_3d = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc) #30 * number of genes matrix of submatrices
          else:
              raise NotImplementedError(f"transform={transform} is not a valid parameter")
          np.fill_diagonal(self.corrcoef_UMAP_3d, 0)
          if calculate_randomized:
              np.fill_diagonal(self.corrcoef_UMAP_random_3d, 0)
          colnames_cells = list(self.ca["Patients"]) #shape of all: (rows=genes, cols=samples)
          rownames_genes = list(self.ra["Genes"])
          tp = np.array(hi_dim)
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #emat= hi_dim #expression
          file_name = "emat_UMAP_3D.csv"
          tp_df.to_csv(file_name, sep=",", header=True) 
          sqrtdmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim)
          tp = np.array(sqrtdmat)
          file_name = "sqrt_dmat_UMAP_3D.csv"
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #dmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim) #dmat are the velocities
          tp_df.to_csv(file_name, sep=",", header=True)  
          alt_dmat = hi_dim_t - hi_dim
          tp = np.array(alt_dmat)
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  # alt dmat: hi_dim_t - hi_dim (for linear)
          file_name = "linear_dmat_UMAP_3D.csv"
          tp_df.to_csv(file_name, sep=",", header=True)  
          tp = np.array(self.corrcoef_UMAP_3d)
          tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)  # raw cc values withot any alt
          file_name = "CorrelationCoefficients_Raw_UMAP_3D.csv"
          tp_df.to_csv(file_name, sep=",", header=True) 


    def UMAPwithArrows3D(self, UMAPembedding: str="UMAPembed3d", UMAP_Arrows: str="MY_ARROWS_umap_3d", make_gify: bool=True, numberNeighb: float=2., typegroup: str="integers", extraNNvalue: int=5, ChangeShape: str="True", Colours: str="Paired"):
        #trans = umap.UMAP(n_neighbors = num_nei, n_components=3).fit(self.S_norm.T)
        UMAP_samplesembedding = getattr(self, UMAPembedding)
        self.UMAPsampleembedding3d = UMAP_samplesembedding
        umap_arrows_embedding_stuff = getattr(self, UMAP_Arrows)
        self.arrowUMAPembed3d = umap_arrows_embedding_stuff
        if ChangeShape == "True":
            my_markers = self.ca["Markers"]
            mark = sorted(set(my_markers))
            corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
            count = 0
            my_types = {}
            for the_type in mark:
                my_types[the_type] = corresp[count]
                count = count + 1
            markers = []
            for mar in my_markers:
                markers.append(my_types.get(mar))
            z = randn(10)
            fig = plt.figure(figsize=(8, 6))
            quiver_scale = 60
            #ax = fig.add_subplot(111, projection='3d')
            ax = fig.add_subplot(projection='3d')
            for i in range(0, len(markers)):
                ax.scatter(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], marker=markers[i], color=self.colorandum[i], s=35)
                ax.quiver(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], self.arrowUMAPembed3d[i, 0], self.arrowUMAPembed3d[i, 1], self.arrowUMAPembed3d[i, 2], arrow_length_ratio=0.35, normalize=False, color="black") #arrow_length_ratio=1.1, color="black") 
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            for sample_type in range(0, len(mark)):
                to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                legend_elements.append(to_add)
        if ChangeShape == "False":
            z = randn(10)
            fig = plt.figure(figsize=(8, 6))
            quiver_scale = 60
            #ax = fig.add_subplot(111, projection='3d')
            ax = fig.add_subplot(projection='3d')
            ax.scatter(self.UMAPsampleembedding3d[:, 0], self.UMAPsampleembedding3d[:, 1], self.UMAPsampleembedding3d[:, 2], marker="o", color=self.colorandum, s=35)
            ax.quiver(self.UMAPsampleembedding3d[:, 0], self.UMAPsampleembedding3d[:, 1], self.UMAPsampleembedding3d[:, 2], self.arrowUMAPembed3d[:, 0], self.arrowUMAPembed3d[:, 1], self.arrowUMAPembed3d[:, 2], arrow_length_ratio=0.35, normalize=False, color="black") #arrow_length_ratio=1.1, color="black") 
            #for i in range(0, len(self.colorandum)):
            #    ax.scatter(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], marker="o", color=self.colorandum[i], s=35)
            #    ax.quiver(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], self.arrowUMAPembed3d[i, 0], self.arrowUMAPembed3d[i, 1], self.arrowUMAPembed3d[i, 2], arrow_length_ratio=0.35, normalize=False, color="black") #arrow_length_ratio=1.1, color="black") 
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
        #ax.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
        ax.legend(handles=legend_elements, title='Key')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_zlabel('UMAP3')
        figure_name3D = 'UMAP3D_Arrows_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
        plt.savefig(figure_name3D)
        #plt.title('UMAP 3D projection', fontsize=10);
        #plt.show()
        gif_name_1 = 'UMAP3D_Arrows_rotation_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.gif'
        gif_name_2 = 'UMAP3D_Arrows_rotation_up_down_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.gif'
        if (make_gify == True):
            def rotate(angle):
                ax.view_init(azim=angle)
            rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
            rot_animation.save(gif_name_1, dpi=80, writer='imagemagick')
            def init():
                return fig,
            def animate(i):
                ax.view_init(elev=(i-45)*4, azim=10)
                return fig,
            ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
            ani_2.save(gif_name_2, dpi=80, writer='imagemagick')
        if ChangeShape == "True":
            #UMAP1 vs UMAP2
            fig = plt.figure(figsize=(15, 10))
            quiver_scale = 60
            for i in range(0, len(markers)):
                plt.scatter(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], marker=markers[i], color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                plt.quiver(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], self.arrowUMAPembed3d[i, 0], self.arrowUMAPembed3d[i, 1], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            for sample_type in range(0, len(mark)):
                to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP2')
            figure_name = 'UMAP3D_Arrow_UMAP1vsUMAP2_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #UMAP1 vs UMAP3
            fig = plt.figure(figsize=(15, 10))
            quiver_scale = 60
            for i in range(0, len(markers)):
                plt.scatter(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 2], marker=markers[i], color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                plt.quiver(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 2], self.arrowUMAPembed3d[i, 0], self.arrowUMAPembed3d[i, 2], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            for sample_type in range(0, len(mark)):
                to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP3')
            figure_name = 'UMAP3D_Arrow_UMAP1vsUMAP3_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #UMAP2 vs UMAP3
            fig = plt.figure(figsize=(15, 10))
            quiver_scale = 60
            for i in range(0, len(markers)):
                plt.scatter(self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], marker=markers[i], color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                plt.quiver(self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], self.arrowUMAPembed3d[i, 1], self.arrowUMAPembed3d[i, 2], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            decoded_my_keys = [i[0] for i in combined]
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            for sample_type in range(0, len(mark)):
                to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            plt.xlabel('UMAP2')
            plt.ylabel('UMAP3')
            figure_name = 'UMAP3D_Arrow_UMAP2vsUMAP3_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            plt.close()
        if ChangeShape == "False":
            #UMAP1 vs UMAP2
            fig = plt.figure(figsize=(15, 10))
            quiver_scale = 60
            plt.scatter(self.UMAPsampleembedding3d[:, 0], self.UMAPsampleembedding3d[:, 1], marker="o", color=self.colorandum, alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
            quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
            plt.quiver(self.UMAPsampleembedding3d[:, 0], self.UMAPsampleembedding3d[:, 1], self.arrowUMAPembed3d[:, 0], self.arrowUMAPembed3d[:, 1], scale=quiver_scale, **quiver_kwargs)
            #for i in range(0, len(self.colorandum)):
                #plt.scatter(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], marker="o", color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                #quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                #plt.quiver(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 1], self.arrowUMAPembed3d[i, 0], self.arrowUMAPembed3d[i, 1], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP2')
            figure_name = 'UMAP3D_Arrow_UMAP1vsUMAP2_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #UMAP1 vs UMAP3
            fig = plt.figure(figsize=(15, 10))
            quiver_scale = 60
            plt.scatter(self.UMAPsampleembedding3d[:, 0], self.UMAPsampleembedding3d[:, 2], marker="o", color=self.colorandum, alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
            quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
            plt.quiver(self.UMAPsampleembedding3d[:, 0], self.UMAPsampleembedding3d[:, 2], self.arrowUMAPembed3d[:, 0], self.arrowUMAPembed3d[:, 2], scale=quiver_scale, **quiver_kwargs)
            #for i in range(0, len(self.colorandum)):
                #plt.scatter(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 2], marker="o", color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                #quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                #plt.quiver(self.UMAPsampleembedding3d[i, 0], self.UMAPsampleembedding3d[i, 2], self.arrowUMAPembed3d[i, 0], self.arrowUMAPembed3d[i, 2], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP3')
            figure_name = 'UMAP3D_Arrow_UMAP1vsUMAP3_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #UMAP2 vs UMAP3
            fig = plt.figure(figsize=(15, 10))
            quiver_scale = 60
            plt.scatter(self.UMAPsampleembedding3d[:, 1], self.UMAPsampleembedding3d[:, 2], marker="o", color=self.colorandum, alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
            quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
            plt.quiver(self.UMAPsampleembedding3d[:, 1], self.UMAPsampleembedding3d[:, 2], self.arrowUMAPembed3d[:, 1], self.arrowUMAPembed3d[:, 2], scale=quiver_scale, **quiver_kwargs)
            #for i in range(0, len(self.colorandum)):
                #plt.scatter(self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], marker="o", color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                #quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                #plt.quiver(self.UMAPsampleembedding3d[i, 1], self.UMAPsampleembedding3d[i, 2], self.arrowUMAPembed3d[i, 1], self.arrowUMAPembed3d[i, 2], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            decoded_my_keys = [i[0] for i in combined]
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            plt.xlabel('UMAP2')
            plt.ylabel('UMAP3')
            figure_name = 'UMAP3D_Arrow_UMAP2vsUMAP3_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            plt.close()

    def estimate_transition_probPCA(self, hidim: str="Sx_sz", embed: str="pcs", transform: str="sqrt",
                                 ndims: int=None, n_sight: int=None, psc: float=None,
                                 knn_random: bool=False, sampled_fraction: float=0.3,
                                 sampling_probs: Tuple[float, float]=(0.5, 0.1), max_dist_embed: float=None,
                                 n_jobs: int=4, threads: int=None, calculate_randomized: bool=True,
                                 random_seed: int=15071990, AdjustPCs: str="False", **kwargs) -> None:
        """Use correlation to estimate transition probabilities for every cells to its embedding neighborhood

        Arguments
        ---------
        hidim: str, default="Sx_sz"
            The name of the attribute containing the high dimensional space. It will be retrieved as getattr(self, hidim)
            The updated vector at time t is assumed to be getattr(self, hidim + "_t")
            Appending .T to the string will transpose the matrix (useful in case we want to use S or Sx)
        embed: str, default="ts"
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
        transform: str, default="sqrt"
            The transformation that is applies on the high dimensional space.
            If None the raw data will be used
        ndims: int, default=None
            The number of dimensions of the high dimensional space to work with. If None all will be considered
            It makes sense only when using principal components
        n_sight: int, default=None (also n_neighbors)
            The number of neighbors to take into account when performing the projection
        psc: float, default=None
            pseudocount added in variance normalizing transform
            If None, 1 would be used for log, 0 otherwise
        knn_random: bool, default=True
            whether to random sample the neighborhoods to speedup calculation
        sampling_probs: Tuple, default=(0.5, 1)
        max_dist_embed: float, default=None
            CURRENTLY NOT USED
            The maximum distance allowed
            If None it will be set to 0.25 * average_distance_two_points_taken_at_random
        n_jobs: int, default=4
            number of jobs to calculate knn
            this only applies to the knn search, for the more time consuming correlation computation see threads
        threads: int, default=None
            The threads will be used for the actual correlation computation by default half of the total.
        calculate_randomized: bool, default=True
            Calculate the transition probabilities with randomized residuals.
            This can be plotted downstream as a negative control and can be used to adjust the visualization scale of the velocity field.
        random_seed: int, default=15071990
            Random seed to make knn_random mode reproducible

        Returns
        -------
        """
        numba_random_seed(random_seed)
        self.which_hidim = hidim #spliced states before timepoint same calculation used to calculate after just with different names for the same thing this is: self.Sx_sz
        if "n_neighbors" in kwargs:
            n_neighbors = kwargs.pop("n_neighbors")
            if len(kwargs) > 0:
                logging.warning(f"keyword arguments were passed but could not be interpreted {kwargs}")
        else:
            n_neighbors = None

        if n_sight is None and n_neighbors is None:
            n_neighbors = int(self.S.shape[1] / 5)

        if (n_sight is not None) and (n_neighbors is not None) and n_neighbors != n_sight:
            raise ValueError("n_sight and n_neighbors are different names for the same parameter, they cannot be set differently")

        if n_sight is not None and n_neighbors is None:
            n_neighbors = n_sight

        if psc is None:
            if transform == "log" or transform == "logratio":
                psc = 1.
            elif transform == "sqrt":
                psc = 1e-10  # for numerical stablity
            else:  # transform == "linear":
                psc = 0

        if knn_random:
            np.random.seed(random_seed)
            self.corr_calc = "knn_random"
            if "pcs" in hidim:  # sic
                hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
                hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
            else:
                if ndims is not None:
                    raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
                hi_dim = getattr(self, hidim)  # [:, :ndims] hi_dim is the
                hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
                if calculate_randomized:
                    self.delta_S_rndm = np.copy(self.delta_S)
                    permute_rows_nsign(self.delta_S_rndm)
                    hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm

            embedding = getattr(self, embed)
            self.embedding_PCA = embedding
            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs) #+1?
            if AdjustPCs == "True":
                pcs_edited = np.empty((np.shape(self.pcs)[0], 0), dtype=float)
                pcs_edited = np.append(pcs_edited, np.array([self.pcs[:, 0]]).transpose(), axis=1)
                myRange = range(np.shape(self.pcs)[1]-1)
                for n in myRange:
                    ScalingFactor = self.percentages[0]/self.percentages[n+1]
                    pcs_edited = np.append(pcs_edited, np.array([ScalingFactor * self.pcs[:, n+1]]).transpose(), axis=1)
                nn.fit(pcs_edited)  # NOTE should support knn
            else:
                nn.fit(embedding)
            self.embedding_PCA_knn = nn.kneighbors_graph(mode="connectivity")
            # Pick random neighbours and prune the rest
            neigh_ixs = self.embedding_PCA_knn.indices.reshape((-1, n_neighbors))
            p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
            p = p / p.sum()

            # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
            # resulting self.corrcoef_PCAf with different number of nonzero entry per row.
            # Not updated yet not to break previous analyses
            # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
            # I change it here since I am doing some breaking changes
            sampling_ixs = np.stack((np.random.choice(neigh_ixs.shape[1],
                                                      size=(int(sampled_fraction * (n_neighbors)),),
                                                      replace=False,
                                                      p=p) for i in range(neigh_ixs.shape[0])), 0)
            self.sampling_ixs = sampling_ixs
            neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
            nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
            self.embedding_PCA_knn = sparse.csr_matrix((np.ones(nonzero),
                                                    neigh_ixs.ravel(),
                                                    np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                                   shape=(neigh_ixs.shape[0],
                                                          neigh_ixs.shape[0]))
            logging.debug(f"Correlation Calculation '{self.corr_calc}'")
            if transform == "log":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef_PCA = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_PCA_random = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
            elif transform == "logratio":
                log2hidim = np.log2(hi_dim + psc)
                delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
                self.corrcoef_PCA = colDeltaCorpartial(log2hidim, delta_hi_dim, neigh_ixs, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + psc) - log2hidim
                    self.corrcoef_PCA_random = colDeltaCorpartial(log2hidim, delta_hi_dim_rndm, neigh_ixs, threads=threads)
            elif transform == "linear":
                self.corrcoef_PCA = colDeltaCorpartial(hi_dim, hi_dim_t - hi_dim, neigh_ixs, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    self.corrcoef_PCA_random = colDeltaCorpartial(hi_dim, hi_dim_t_rndm - hi_dim, neigh_ixs, threads=threads)
            elif transform == "sqrt":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef_PCA = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_PCA_random = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
            else:
                raise NotImplementedError(f"transform={transform} is not a valid parameter")
            np.fill_diagonal(self.corrcoef_PCA, 0)
            if np.any(np.isnan(self.corrcoef_PCA)):
                self.corrcoef_PCA[np.isnan(self.corrcoef_PCA)] = 1
                logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_PCA_random, 0)
                if np.any(np.isnan(self.corrcoef_PCA_random)):
                    self.corrcoef_PCA_random[np.isnan(self.corrcoef_PCA_random)] = 1
                    logging.warning("Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            logging.debug(f"Done Correlation Calculation")
        else:
            self.corr_calc = "full"
            if "pcs" in hidim:  # sic
                hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
                hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
            else:
                if ndims is not None:
                    raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
                hi_dim = getattr(self, hidim)  # [:, :ndims]
                hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
                if calculate_randomized:
                    self.delta_S_rndm = np.copy(self.delta_S)
                    permute_rows_nsign(self.delta_S_rndm)
                    hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
            embedding = getattr(self, embed)
            self.embedding_PCA = embedding
            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
            if AdjustPCs == "True":
                pcs_edited = np.empty((np.shape(self.pcs)[0], 0), dtype=float)
                pcs_edited = np.append(pcs_edited, np.array([self.pcs[:, 0]]).transpose(), axis=1)
                myRange = range(np.shape(self.pcs)[1]-1)
                for n in myRange:
                    ScalingFactor = self.percentages[0]/self.percentages[n+1]
                    pcs_edited = np.append(pcs_edited, np.array([ScalingFactor * self.pcs[:, n+1]]).transpose(), axis=1)
                nn.fit(pcs_edited)  # NOTE should support knn
            else:
                nn.fit(embedding)
            self.embedding_PCA_knn = nn.kneighbors_graph(mode="connectivity") #def of embedding_knn
            logging.debug("Correlation Calculation 'full'")
            if transform == "log":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef_PCA = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_PCA_random = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
            elif transform == "logratio":
                log2hidim = np.log2(hi_dim + psc)
                delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
                self.corrcoef_PCA = colDeltaCor(log2hidim, delta_hi_dim, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + 1) - log2hidim
                    self.corrcoef_PCA_random = colDeltaCor(log2hidim, delta_hi_dim_rndm, threads=threads)
            elif transform == "linear":
                self.corrcoef_PCA = colDeltaCor(hi_dim, hi_dim_t - hi_dim, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    self.corrcoef_PCA_random = colDeltaCor(hi_dim, hi_dim_t_rndm - hi_dim, threads=threads)
            elif transform == "sqrt":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef_PCA = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)

                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_PCA_random = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
            else:
                raise NotImplementedError(f"transform={transform} is not a valid parameter")
            np.fill_diagonal(self.corrcoef_PCA, 0)
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_PCA_random, 0)
            colnames_cells = list(self.ca["Patients"]) #shape of all: (rows=genes, cols=samples)
            rownames_genes = list(self.ra["Genes"])
            tp = np.array(hi_dim)
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #emat= hi_dim #expression
            file_name = "emat_PCA.csv"
            tp_df.to_csv(file_name, sep=",", header=True) 
            sqrtdmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim)
            tp = np.array(sqrtdmat)
            file_name = "sqrt_dmat_PCA.csv"
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #dmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim) #dmat are the ve,locities
            tp_df.to_csv(file_name, sep=",", header=True)  
            alt_dmat = hi_dim_t - hi_dim
            tp = np.array(alt_dmat)
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  # alt dmat: hi_dim_t - hi_dim (for linear)
            file_name = "linear_dmat_PCA.csv"
            tp_df.to_csv(file_name, sep=",", header=True)
            tp = np.array(self.corrcoef_PCA)
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)  # raw cc values withot any alt
            file_name = "CorrelationCoefficients_Raw_PCA.csv"
            tp_df.to_csv(file_name, sep=",", header=True)  
 	#AG added the definitions of TPMin, TPMax TPStep below and incorported a loop so that they are directly imported from the shell script used to execute veloCD. These can then be flexible to whatever the user enters.
    import sys

    TPMin = int(sys.argv[-12])
    TPMax = int(sys.argv[-11])
    TPStep = int(sys.argv[-10])
    def calculate_embedding_shiftPCA(self, sigma_corr: float=0.05, expression_scaling: bool=False, scaling_penalty: float=1., metadataToUse: str="metadata.csv", 
                                 pathToTransitionProbabilityCode: str="/", grouptype: str="Strings", secondgroup: str="metadataVersion.csv", manyfeatures: str="True", 
                                 grouptype2: str="Strings", Randomsampling: str="False", extraNNvalue: int=5,TPMin: int =2, TPMax: int =4, TPStep: int =1, DE_location: str="/location/", selfCC_folder: str="/location/", 
                                 VecX_location: str="/location/", VecY_loc: str="/location/", VecZloc: str="/location/", TP_folder: str="/location/", 
                                 PCA_GroupLevelTP: str="/location/", 
                                 ThreeDVersion: str="yes", DoSecondPred: str="no", UserThreshold: float=0.5) -> None:
        logging.debug("Calculate transition probability")
        # --------- AG MADE THIS A LOOP FOR RANGE OF TPNNs ---------
        for nn in range(TPMin, TPMax, TPStep):  # TPMin, TPMax, TPStep should come from function arguments or parsed CLI
            extraNNvalue = nn
            self.embedding_PCA_Red = self.embedding_PCA[:,:2]
            if self.corr_calc == "full" or self.corr_calc == "knn_random":
            # NOTE maybe sparse matrix here are slower than dense
            # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
                self.transition_prob_PCA_Red = np.exp(self.corrcoef_PCA / sigma_corr) * self.embedding_PCA_knn.toarray()  # naive, the exponential of the correlation coefficient/kernal scaling, AG changed A to to.array()
                self.transition_prob_PCA_Red /= self.transition_prob_PCA_Red.sum(1)[:, None]
                if hasattr(self, "corrcoef_random"):
                    logging.debug("Calculate transition probability for negative control")
                    self.transition_prob_PCA_Red_random = np.exp(self.corrcoef_PCA_random / sigma_corr) * self.embedding_PCA_knn.toarray() # naive, AG changed A to to.array()
                    self.transition_prob_PCA_Red_random /= self.transition_prob_PCA_Red_random.sum(1)[:, None]
                unitary_vectors_PCA_Red = self.embedding_PCA_Red.T[:, None, :] - self.embedding_PCA_Red.T[:, :, None]
                with np.errstate(divide='ignore', invalid='ignore'):
                    unitary_vectors_PCA_Red /= np.linalg.norm(unitary_vectors_PCA_Red, ord=2, axis=0)  # divide by L2
                    np.fill_diagonal(unitary_vectors_PCA_Red[0, ...], 0)  # fix nans replace NAN with 0
                    np.fill_diagonal(unitary_vectors_PCA_Red[1, ...], 0) #embedding comes from the PCA
                self.delta_embedding_PCA_Red = (self.transition_prob_PCA_Red * unitary_vectors_PCA_Red).sum(2) #.A means change the data type from a matrix to an array #unit vector*transition probability
                self.delta_embedding_PCA_Red -= (self.embedding_PCA_knn.toarray() * unitary_vectors_PCA_Red).sum(2) / self.embedding_PCA_knn.toarray().sum(1).T #AG changed A to to.array()
                self.delta_embedding_PCA_Red = self.delta_embedding_PCA_Red.T #transposes the vector
                if expression_scaling:
                    hi_dim = getattr(self, self.which_hidim)
                    estim_delta = hi_dim.dot(self.transition_prob_PCA_Red.T) - hi_dim.dot((self.embedding_PCA_knn.toarray() / self.embedding_PCA_knn.toarray().sum(1))).T #this is the same problem as above, AG changed A to to.array()
                    cos_proj = (self.delta_S * estim_delta).sum(0) / np.sqrt((estim_delta**2).sum(0))
                    self.scaling_PCA= np.clip(cos_proj / scaling_penalty, 0, 1)
                    self.delta_embedding_PCA_Red = self.delta_embedding_PCA_Red * self.scaling[:, None]
                if hasattr(self, "corrcoef_random"):
                    self.delta_embedding_PCA_Red_random = (self.transition_prob_PCA_Red_random * unitary_vectors_PCA_Red).sum(2)
                    self.delta_embedding_PCA_Red_random -= (self.embedding_PCA_knn.toarray() * unitary_vectors_PCA_Red).sum(2) / self.embedding_PCA_knn.toarray().sum(1).T #AG changed A to to.array()
                    self.delta_embedding_PCA_Red_random = self.delta_embedding_PCA_Red_random.T
                    if expression_scaling:
                        estim_delta_rndm_Red = hi_dim.dot(self.transition_prob_PCA_Red_random.T) - hi_dim.dot((self.embedding_PCA_knn.toarray() / self.embedding_PCA_knn.toarray().sum(1))).T #AG changed A to to.array()
                        cos_proj_rndm = (self.delta_S_rndm * estim_delta_rndm_Red).sum(0) / np.sqrt((estim_delta_rndm_Red**2).sum(0))
                        self.scaling_rndm = np.clip(cos_proj_rndm / scaling_penalty, 0, 1)
                        self.delta_embedding_PCA_Red_random = self.delta_embedding_PCA_Red_random * self.scaling_rndm[:, None]
            else:
                # NOTE should implement a version with cython
                raise NotImplementedError(f"Weird value self.corr_calc={self.corr_calc}")
        colnames_cells = list(self.ca["Patients"])
        genes_embedding = np.array(self.delta_embedding_PCA_Red)
        spliced_df = pd.DataFrame(genes_embedding, columns=["x", "y"], index=colnames_cells)
        file_name = "Delta_embedding_TPNN" + str(extraNNvalue) + "_PCA.csv"
        os.chdir(DE_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(self.corrcoef_PCA)
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "selfCC_TPNN" + str(extraNNvalue) + "_PCA.csv"
        os.chdir(selfCC_folder)
        spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors_PCA_Red[0])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecX_TPNN" + str(extraNNvalue) + "_PCA.csv"
        os.chdir(VecX_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors_PCA_Red[1])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecY_TPNN" + str(extraNNvalue) + ".csv"
        #os.chdir(VecY_loc)
        #spliced_df.to_csv(file_name, sep=",", header=True)

        tp = np.array(self.transition_prob_PCA_Red) #number of samples x number of samples
        tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)
        file_name = "Transition_Probability_PCA_2D_TPNN" + str(extraNNvalue) + ".csv"
        os.chdir(TP_folder)
        tp_df.to_csv(file_name, sep=",", header=True)

        #####RUN TP code here################################
        #note: the following requires these R packages to be installed:  'gtools', 'reshape2', ggplot2', grid' and 'tidyverse'
        #to code to install them has been commented out but is below:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        utils = importr('utils')
        #utils.install_packages('gtools')
        #utils.install_packages('reshape2')
        #utils.install_packages('ggplot2')        
        #utils.install_packages('grid')
        #utils.install_packages('tidyverse')
        gtools = importr("gtools")
        reshape2 = importr("reshape2")
        ggplot2 = importr("ggplot2")
        grid = importr("grid")
        tidyverse = importr("tidyverse")
        path=pathToTransitionProbabilityCode
        version = "SampleGroup1"
        def RunTransitionProbabilityCode(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold):
            r=ro.r
            r.source(path+"TransitionProbability.R")
            p=r.TransitionProbability(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold)
            return p
        #a=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, metadataToUse, PCA_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold) 
        b=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, metadataToUse, PCA_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold)
        if manyfeatures == "True" and DoSecondPred == "yes":
            version = "SampleGroup2" #secondgroup
            #c=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, secondgroup, PCA_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold) 
            d=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, secondgroup, PCA_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold)
        if ThreeDVersion == "yes":
            self.embedding_PCA_Three = self.embedding_PCA[:,:3]
            if self.corr_calc == "full" or self.corr_calc == "knn_random":
                # NOTE maybe sparse matrix here are slower than dense
                # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
                self.transition_prob_PCA_Three = np.exp(self.corrcoef_PCA / sigma_corr) * self.embedding_PCA_knn.toarray()  # naive, the exponential of the correlation coefficient/kernal scaling, #AG changed A to to.array()
                self.transition_prob_PCA_Three /= self.transition_prob_PCA_Three.sum(1)[:, None]
                if hasattr(self, "corrcoef_random"):
                    logging.debug("Calculate transition probability for negative control")
                    self.transition_prob_PCA_Three_random = np.exp(self.corrcoef_PCA_random / sigma_corr) * self.embedding_PCA_knn.toarray() # naive, #AG changed A to to.array()
                    self.transition_prob_PCA_Three_random /= self.transition_prob_PCA_Three_random.sum(1)[:, None]
                unitary_vectors_PCA = self.embedding_PCA_Three.T[:, None, :] - self.embedding_PCA_Three.T[:, :, None]  # shape (2,ncells,ncells) #for each cells coordinates - row1: cell1-cell1, cell2-cell1, cell3-cell1 cell4-cell1
                with np.errstate(divide='ignore', invalid='ignore'):
                    unitary_vectors_PCA /= np.linalg.norm(unitary_vectors_PCA, ord=2, axis=0)  # divide by L2
                    np.fill_diagonal(unitary_vectors_PCA[0, ...], 0)  # fix nans replace NAN with 0
                    np.fill_diagonal(unitary_vectors_PCA[1, ...], 0) #embedding comes from the PCA
                    np.fill_diagonal(unitary_vectors_PCA[2, ...], 0) 
                self.delta_embedding_PCA_Three = (self.transition_prob_PCA_Three * unitary_vectors_PCA).sum(2) #.A means change the data type from a matrix to an array #unit vector*transition probabilit
                self.delta_embedding_PCA_Three -= (self.embedding_PCA_knn.toarray() * unitary_vectors_PCA).sum(2) / self.embedding_PCA_knn.toarray().sum(1) ##AG changed A to to.array(), remove .T
                self.delta_embedding_PCA_Three = self.delta_embedding_PCA_Three.T #transposes the vector
                if expression_scaling:
                    hi_dim = getattr(self, self.which_hidim)
                    estim_delta = hi_dim.dot(self.transition_prob_PCA_Three.T) - hi_dim.dot((self.embedding_PCA_knn.toarray() / self.embedding_PCA_knn.toarray().sum(1)).T) #this is the same problem as above, #AG changed A to toarray()
                    cos_proj = (self.delta_S * estim_delta).sum(0) / np.sqrt((estim_delta**2).sum(0))
                    self.scaling_PCA= np.clip(cos_proj / scaling_penalty, 0, 1)
                    self.delta_embedding_PCA_Three = self.delta_embedding_PCA_Three * self.scaling[:, None]
                if hasattr(self, "corrcoef_random"):
                    self.delta_embedding_PCA_Three_random = (self.transition_prob_PCA_Three_random * unitary_vectors_PCA).sum(2)
                    self.delta_embedding_PCA_Three_random -= (self.embedding_PCA_knn.toarray() * unitary_vectors_PCA).sum(2) / self.embedding_PCA_knn.toarray().sum(1).T ##AG changed A to to.array()
                    self.delta_embedding_PCA_Three_random = self.delta_embedding_PCA_Three_random.T
                    if expression_scaling:
                        estim_delta_rndm = hi_dim.dot(self.transition_prob_PCA_Three_random.T) - hi_dim.dot((self.embedding_PCA_knn.toarray() / self.embedding_PCA_knn.toarray().sum(1)).T) ##AG changed A to to.array()
                        cos_proj_rndm = (self.delta_S_rndm * estim_delta_rndm).sum(0) / np.sqrt((estim_delta_rndm**2).sum(0))
                        self.scaling_rndm = np.clip(cos_proj_rndm / scaling_penalty, 0, 1)
                        self.delta_embedding_PCA_Three_random = self.delta_embedding_PCA_Three_random * self.scaling_rndm[:, None]
            else:
                # NOTE should implement a version with cython
                raise NotImplementedError(f"Weird value self.corr_calc={self.corr_calc}")
            colnames_cells = list(self.ca["Patients"])
            genes_embedding = np.array(self.delta_embedding_PCA_Three)
            spliced_df = pd.DataFrame(genes_embedding, columns=["x", "y", "Z"], index=colnames_cells)
            file_name = "PCA_3D_Delta_embedding_TPNN" + str(extraNNvalue) + "_PCA.csv"
            os.chdir(DE_location)
            #spliced_df.to_csv(file_name, sep=",", header=True)
            #genes_embedding = np.array(self.corrcoef_PCA)
            #spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
            #file_name = "PCA_3D_selfCC_TPNN" + str(extraNNvalue) + "_PCA.csv"
            #os.chdir(selfCC_folder)
            #spliced_df.to_csv(file_name, sep=",", header=True)
            genes_embedding = np.array(unitary_vectors_PCA[0])
            spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
            file_name = "PCA_3D_UnVecX_TPNN" + str(extraNNvalue) + "_PCA.csv"
            os.chdir(VecX_location)
            #spliced_df.to_csv(file_name, sep=",", header=True)
            genes_embedding = np.array(unitary_vectors_PCA[1])
            spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
            file_name = "PCA_3D_UnVecY_TPNN" + str(extraNNvalue) + "_PCA.csv"
            os.chdir(VecY_loc)
            #spliced_df.to_csv(file_name, sep=",", header=True)
            genes_embeddingZ = np.array(unitary_vectors_PCA[2])
            spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
            file_name = "PCA_3D_UnVecZ_TPNN" + str(extraNNvalue) + "_PCA.csv"
            os.chdir(VecZloc)
            #spliced_df.to_csv(file_name, sep=",", header=True)
            tp = np.array(self.transition_prob_PCA_Three) #number of samples x number of samples
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)
            file_name = "Transition_Probability_PCA_3D_TPNN" + str(extraNNvalue) + ".csv"
            os.chdir(TP_folder)
            tp_df.to_csv(file_name, sep=",", header=True)
            #####RUN TP code here################################
            #note: the following requires these R packages to be installed:  'gtools', 'reshape2', ggplot2', grid' and 'tidyverse'
            #to code to install them has been commented out but is below:
            import rpy2.robjects as ro
            from rpy2.robjects.packages import importr
            utils = importr('utils')
            #utils.install_packages('gtools')
            #utils.install_packages('reshape2')
            #utils.install_packages('ggplot2')        
            #utils.install_packages('grid')
            #utils.install_packages('tidyverse')
            gtools = importr("gtools")
            reshape2 = importr("reshape2")
            ggplot2 = importr("ggplot2")
            grid = importr("grid")
            tidyverse = importr("tidyverse")
            path=pathToTransitionProbabilityCode
            version = "SampleGroup1"
            def RunTransitionProbabilityCode(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold):
                r=ro.r
                r.source(path+"TransitionProbability.R")
                p=r.TransitionProbability(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold)
                return p
            #a=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, metadataToUse, PCA_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold) 
            b=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, metadataToUse, PCA_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold)
            if manyfeatures == "True" and DoSecondPred == "yes":
                version = "SampleGroup2" #secondgroup
                #c=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, secondgroup, PCA_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold) 
                d=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, secondgroup, PCA_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold)    


    def one_arrowPCA(self, typegroup: str="integers", extraNNvalue: int=5, ChangeShape: str="True", Colours: str="Paired", ThreeDVersion: str="yes", make_gify: bool=True):
        if ChangeShape == "True":
            plt.figure(None,(15,12))
            quiver_scale = 60
            my_markers = self.ca["Markers"]
            mark = sorted(set(my_markers))
            corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
            count = 0
            my_types = {}
            for the_type in mark:
                my_types[the_type] = corresp[count]
                count = count + 1
            markers = []
            for mar in my_markers:
                markers.append(my_types.get(mar))
            z = randn(10)
            for i in range(0, len(markers)):
                plt.scatter(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1], marker=markers[i], color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                plt.quiver(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1], self.delta_embedding_PCA_Red[i, 0], self.delta_embedding_PCA_Red[i, 1], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            for sample_type in range(0, len(mark)):
                to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            #ax.set_xlim(0, 10)
            #plt.axis("off")
            plt.xlabel('PC1 '+'('+str(self.percentagesRounded[0])+'%)', fontsize=12)
            plt.ylabel('PC2 '+'('+str(self.percentagesRounded[1])+'%)', fontsize=12)
            figure_name = 'PCA_2D_Arrows_TPNN_' + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #plt.show()
            plt.close()
        if ChangeShape == "False":
            plt.figure(None,(15,12))
            quiver_scale = 60
            z = randn(10)
            #ax = fig.add_subplot(111, projection='3d')
            #plt.scatter(self.embedding_PCA[:, 0], self.embedding_PCA[:, 1],c="0.8", alpha=1, s=120, edgecolor="")
            plt.scatter(self.embedding_PCA[:, 0], self.embedding_PCA[:, 1], marker="o", color=self.colorandum, alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
            quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
            plt.quiver(self.embedding_PCA[:, 0], self.embedding_PCA[:, 1], self.delta_embedding_PCA_Red[:, 0], self.delta_embedding_PCA_Red[:, 1], scale=quiver_scale, **quiver_kwargs)
            #for i in range(0, len(self.colorandum)):
                #plt.scatter(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1],c="0.8", alpha=1, s=120, edgecolor="")
                #plt.scatter(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1], marker="o", color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                #quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                #plt.quiver(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1], self.delta_embedding_PCA_Red[i, 0], self.delta_embedding_PCA_Red[i, 1], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            #ax.set_xlim(0, 10)
            #plt.axis("off")
            plt.xlabel('PC1 '+'('+str(self.percentagesRounded[0])+'%)', fontsize=12)
            plt.ylabel('PC2 '+'('+str(self.percentagesRounded[1])+'%)', fontsize=12)
            figure_name = 'PCA_2D_Arrows_TPNN_' + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #plt.show()
            plt.close()
        if ThreeDVersion == "yes":
            if ChangeShape == "True":
                my_markers = self.ca["Markers"]
                mark = sorted(set(my_markers))
                corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
                count = 0
                my_types = {}
                for the_type in mark:
                    my_types[the_type] = corresp[count]
                    count = count + 1
                markers = []
                for mar in my_markers:
                    markers.append(my_types.get(mar))
                z = randn(10)
                fig = plt.figure(figsize=(8, 6))
                quiver_scale = 60
                ax = fig.add_subplot(111, projection='3d')
                for i in range(0, len(markers)):
                    ax.scatter(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1], self.embedding_PCA[i, 2], marker=markers[i], color=self.colorandum[i], s=35)
                    ax.quiver(self.embedding_PCA[i, 0], self.embedding_PCA[i, 1], self.embedding_PCA[i, 2], self.delta_embedding_PCA_Three[i, 0], self.delta_embedding_PCA_Three[i, 1], self.delta_embedding_PCA_Three[i, 2], arrow_length_ratio=0.35, normalize=False, color="black")
                hexy_codes = []
                my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
                decoded_my_keys = []
                for item in my_keys:
                    decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
                cmap2 = matplotlib.colormaps[Colours]    # PiYG
                for i in range(cmap2.N):
                    rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                    hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
                legend_elements = []
                combined = list(zip(decoded_my_keys, hexy_codes))
                combined = sorted(combined, key=itemgetter(0))
                if (typegroup != "Strings"):
                    timepoints = [int(i[0]) for i in combined]
                    colys = [i[1] for i in combined]
                    tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                    tempDF = tempDF.sort_values('Group')
                    groupys2 = tempDF['Group'].tolist()
                    cols2 = tempDF['colours'].tolist()
                    combined = list(zip(groupys2, cols2))
                decoded_my_keys = [i[0] for i in combined]
                hexy_codes = [i[1] for i in combined]
                for sample_type in range(0, len(decoded_my_keys)):
                    to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                    legend_elements.append(to_add)
                for sample_type in range(0, len(mark)):
                    to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                    legend_elements.append(to_add)
                #ax.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
                ax.legend(handles=legend_elements, title='Key') #Timepont \n(hours)
                #ax.set_xlim(0, 10)
                #plt.axis("off")
                ax.set_xlabel('PC1 '+'('+str(self.percentagesRounded[0])+'%)', fontsize=12)
                ax.set_ylabel('PC2 '+'('+str(self.percentagesRounded[1])+'%)', fontsize=12)
                ax.set_zlabel('PC3 '+'('+str(self.percentagesRounded[2])+'%)', fontsize=12)
                figure_name = 'PCA_Arrows_3D_TPNN_' + str(extraNNvalue) + '.png'
                plt.savefig(figure_name)
                #plt.show()
                plt.close()
            if ChangeShape == "False":
                z = randn(10)
                fig = plt.figure(figsize=(8, 6))
                quiver_scale = 60
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(self.embedding_PCA[:, 0], self.embedding_PCA[:, 1], self.embedding_PCA[:, 2], marker="o", color=self.colorandum, s=35)
                ax.quiver(self.embedding_PCA[:, 0], self.embedding_PCA[:, 1], self.embedding_PCA[:, 2], self.delta_embedding_PCA_Three[:, 0], self.delta_embedding_PCA_Three[:, 1], self.delta_embedding_PCA_Three[:, 2], arrow_length_ratio=0.35, normalize=False, color="black")
                hexy_codes = []
                my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
                decoded_my_keys = []
                for item in my_keys:
                    decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
                cmap2 = matplotlib.colormaps[Colours]    # PiYG
                for i in range(cmap2.N):
                    rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                    hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
                legend_elements = []
                combined = list(zip(decoded_my_keys, hexy_codes))
                combined = sorted(combined, key=itemgetter(0))
                if (typegroup != "Strings"):
                    timepoints = [int(i[0]) for i in combined]
                    colys = [i[1] for i in combined]
                    tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                    tempDF = tempDF.sort_values('Group')
                    groupys2 = tempDF['Group'].tolist()
                    cols2 = tempDF['colours'].tolist()
                    combined = list(zip(groupys2, cols2))
                decoded_my_keys = [i[0] for i in combined]
                hexy_codes = [i[1] for i in combined]
                for sample_type in range(0, len(decoded_my_keys)):
                    to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                    legend_elements.append(to_add)
                #ax.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
                ax.legend(handles=legend_elements, title='Key') #Timepont \n(hours)
                #ax.set_xlim(0, 10)
                #plt.axis("off")
                ax.set_xlabel('PC1 '+'('+str(self.percentagesRounded[0])+'%)', fontsize=12)
                ax.set_ylabel('PC2 '+'('+str(self.percentagesRounded[1])+'%)', fontsize=12)
                ax.set_zlabel('PC3 '+'('+str(self.percentagesRounded[2])+'%)', fontsize=12)
                figure_name = 'PCA_Arrows_3D_TPNN' + str(extraNNvalue) + '.png'
                plt.savefig(figure_name)
                #plt.show()
                plt.close()
        if (make_gify == True):
            def rotate(angle):
                ax.view_init(azim=angle)
            rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
            name1 = "PCA_arrows_3D_LeftRight_TPNN" + str(extraNNvalue) + ".gif"
            rot_animation.save(name1, dpi=80, writer='imagemagick')
            def init():
                return fig,
            def animate(i):
                ax.view_init(elev=(i-45)*4, azim=10)
                return fig,
            ani_2 = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=100, blit=True)
            name2 = "PCA_arrows_3D_rotation_UpDown_TPNN" + str(extraNNvalue) + ".gif"
            ani_2.save(name2, dpi=80, writer='imagemagick')
        plt.close()



    def _perform_PCA_imputed(self, n_components: int=None) -> None:
        """Simply performs PCA of `Sx_norm` and save the result as  `pcax`"""
        self.pcax = PCA(n_components=n_components)
        self.pcsx = self.pcax.fit_transform(self.Sx_norm.T)

    def _plot_pca_imputed(self, dim: List[int]=[0, 1, 2], elev: float=60, azim: float=-140) -> None:
        """Plot 3d PCA of the smoothed data
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.pcsx[:, dim[0]],
                   self.pcsx[:, dim[1]],
                   self.pcsx[:, dim[2]],
                   c=self.colorandum)
        ax.view_init(elev=elev, azim=azim)



    def fit_gammas(self, steady_state_bool: np.ndarray=None, use_imputed_data: bool=True, use_size_norm: bool=False,
                   fit_offset: bool=True, fixperc_q: bool=False, weighted: bool=True, weights: np.ndarray = "maxmin_diag",
                   limit_gamma: bool=False, maxmin_perc: List[float]=[2, 98], maxmin_weighted_pow: float=15) -> None:
        """Fit gamma using spliced and unspliced data

        Arguments
        ---------
        steady_state_bool: np.ndarray, default=None
            if a boolean array is specified, gamma is fitted using only the corresponding samples
        use_imputed_data: bool, default=True
            use knn smoothed data
        use_size_norm: bool, default=False
            use size normalized data for the fit
        fit_offset: bool, default=True
            Fit with offset
        fixperc_q: bool, default=False
            (when fit_offset==False) Wether to fix the offset to a lower percentile of the unspliced
        weighted: bool, default=True
            use weights for the least squares fit
        weights: string or np.ndarray, default="maxmin_diag"
            the method to determine the weights of the least squares fit.
            "maxmin_diag", "maxmin", "sum", "prod", "maxmin_weighted" are supported
            if a 2d np.ndarray is provided the entry (i,j) is the weight of the cell j when fitting gamma to gene i
        limit_gamma: np.ndarray, default=True
            whether to limit gamma when unspliced is much higher than spliced
        maxmin_perc: List[float], default=[2,98]
            the percentile to use if weights = "maxmin" or "maxmin_diag"

        Returns
        -------
        Nothing it just creates the attributes:
        gammas: np.ndarray
            the vector of the gammas fit to each gene
        q: np.ndarray
            the vector of offsets of the fit
        R2: np.ndarray (optional)
            The vector of squared coefficient of determination

        """
        self.Sx = np.copy(self.S)
        self.Ux = np.copy(self.U)
        self.Sx_sz = np.copy(self.S)
        self.Ux_sz = np.copy(self.U)
        if steady_state_bool:
            self.steady_state = steady_state_bool
        else:
            self.steady_state = np.ones(self.S.shape[1], dtype=bool)

        if use_imputed_data:
            if use_size_norm:
                tmpS = self.Sx_sz
                tmpU = self.Ux_sz
            else:
                tmpS = self.Sx
                tmpU = self.Ux
        else:
            if use_size_norm:
                tmpS = self.S_sz
                tmpU = self.U_sz
            else:
                tmpS = self.S
                tmpU = self.U

        if weighted:
            if type(weights) is np.ndarray:
                W = weights
            elif weights == "sum":
                W = (tmpS / np.percentile(tmpS, 99, 1)[:, None]) + (tmpU / np.percentile(tmpU, 99, 1)[:, None])
            elif weights == "prod":
                W = (tmpS / np.percentile(tmpS, 99, 1)[:, None]) * (tmpU / np.percentile(tmpU, 99, 1)[:, None])
            elif weights == "maxmin_weighted":
                # Slightly smoother than just takin top and bottom percentile
                down, up = np.percentile(tmpS, maxmin_perc, 1)  # Do this asymmetrically, data is sparse!
                Srange = np.clip(tmpS, down[:, None], up[:, None])
                Srange -= Srange.min(1)[:, None]
                Srange /= Srange.max(1)[:, None]
                W = 0.5 * (Srange**maxmin_weighted_pow + (1 - Srange)**maxmin_weighted_pow)
            elif weights == "maxmin":
                down, up = np.percentile(tmpS, maxmin_perc, 1)  # Do this asymmetrically, data is sparse!
                W = ((tmpS <= down[:, None]) | (tmpS >= up[:, None])).astype(float)
            elif weights == "maxmin_diag": #default goes down this route
                denom_Sx = np.percentile(self.Sx, 99.9, 1)
                if np.sum(denom_Sx == 0):
                    denom_Sx[denom_Sx == 0] = np.maximum(np.max(self.Sx[denom_Sx == 0, :], 1), 0.001)
                denom_Ux = np.percentile(self.Ux, 99.9, 1)
                if np.sum(denom_Ux == 0):
                    denom_Ux[denom_Ux == 0] = np.maximum(np.max(self.Ux[denom_Ux == 0, :], 1), 0.001)
                Sx_maxnorm = self.Sx / denom_Sx[:, None]
                Ux_maxnorm = self.Ux / denom_Ux[:, None]
                X = Sx_maxnorm + Ux_maxnorm
                down, up = np.percentile(X, maxmin_perc, axis=1)
                W = ((X <= down[:, None]) | (X >= up[:, None])).astype(float)
            elif weights == "maxmin_double":
                denom_Sx = np.percentile(self.Sx, 99.9, 1)
                denom_Sx[denom_Sx == 0] = np.maximum(np.max(self.Sx[denom_Sx == 0, :], 1), 0.001)
                denom_Ux = np.percentile(self.Ux, 99.9, 1)
                denom_Ux[denom_Ux == 0] = np.maximum(np.max(self.Ux[denom_Ux == 0, :], 1), 0.001)
                Sx_maxnorm = self.Sx / denom_Sx[:, None]
                Ux_maxnorm = self.Ux / denom_Ux[:, None]
                X = Sx_maxnorm + Ux_maxnorm
                down, up = np.percentile(X, maxmin_perc, axis=1)
                W = ((X <= down[:, None]) | (X >= up[:, None])).astype(float)
                down, up = np.percentile(self.Sx, maxmin_perc, 1)
                W += ((self.Sx <= down[:, None]) | (self.Sx >= up[:, None])).astype(float)

        if fit_offset:
            if weighted:
                self.gammas, self.q, self.R2 = fit_slope_weighted_offset(tmpU[:, self.steady_state], #loops through the genes and fits it to the slope, gammas = slopes, Y: np.ndarray, shape=(genes, cells)the dependent variable (unspliced), X: np.ndarray, shape=(genes, cells): the independent variable (spliced)
                                                                         tmpS[:, self.steady_state],
                                                                         W,
                                                                         return_R2=True, #add third parameter: np.array of the svp's
                                                                         limit_gamma=limit_gamma) #Performs the _fit1_slope_weighted_offset function from the estimation source file this fits a weighted linear regression model with intercept
            else:
                if limit_gamma:
                    logging.warning("limit_gamma not implemented with this settings")
                self.gammas, self.q = fit_slope_offset(tmpU[:, self.steady_state],
                                                       tmpS[:, self.steady_state]) #add third parameter: np.array of the svp's
        elif fixperc_q:
            if weighted:
                self.gammas, self.q = fit_slope_weighted_offset(tmpU[:, self.steady_state],
                                                                tmpS[:, self.steady_state],
                                                                W, fixperc_q=True, limit_gamma=limit_gamma) #add third parameter: np.array of the svp's
            else:
                if limit_gamma:
                    logging.warning("limit_gamma not implemented with this settings")
                self.gammas, self.q = fit_slope_offset(tmpU[:, self.steady_state],
                                                       tmpS[:, self.steady_state],
                                                       fixperc_q=True) #add third parameter: np.array of the svp's
        else:
            if weighted:
                self.gammas, self.R2 = fit_slope_weighted(tmpU[:, self.steady_state],
                                                          tmpS[:, self.steady_state],
                                                          W,
                                                          return_R2=True,
                                                          limit_gamma=limit_gamma) #add third parameter: np.array of the svp's
                self.q = np.zeros_like(self.gammas)
            else:
                if limit_gamma:
                    logging.warning("limit_gamma not implemented with this settings")
                self.gammas = fit_slope(tmpU[:, self.steady_state],
                                        tmpS[:, self.steady_state]) #add third parameter: np.array of the svp's
                self.q = np.zeros_like(self.gammas)

        # Fix gammas
        self.gammas[~np.isfinite(self.gammas)] = 0


    #def filter_genes_good_fit(self, minR: float=0.1, min_gamma: float=0.01) -> None:
    #    """For backwards compatibility a wrapper around filter_genes_by_phase_portrait
    #    """
    #    return self.filter_genes_by_phase_portrait(minR2=minR, min_gamma=min_gamma, minCorr=None)


    def filter_genes_by_phase_portrait(self, minR2: float=0.1, min_gamma: float=0.01, minCorr: float=0.1) -> None:
        """Use the coefficient of determination to filter away genes that have an irregular/complex phase portrait

        Arguments
        ---------
        minR2: float, default=0.1
            Filter away low coefficient of determination fits. If None this filtering will be skipped
        min_gamma: float, default=0.01
            Filter away low gammas. If None this filtering will be skipped
        minCorr: flaot, default=0.2
            Filter away low spliced-usnpliced correlation. If None this filtering will be skipped

        Returns
        -------
        Nothing but modifies it filters out the genes that do not satisfy the conditions
        This affects: "U", "U_sz", "U_norm", "Ux", "Ux_sz", "Ux_norm", "S", "S_sz", "S_norm", "Sx", "Sx_sz", "Sx_norm", "gammas", "q", "R2"
        """

        def paired_correlation_rows(A: np.array, B: np.array) -> np.array:
            A_m = A - A.mean(1)[:, None]
            B_m = B - B.mean(1)[:, None]
            return (A_m * B_m).sum(1) / (np.linalg.norm(A_m, 2, 1) * np.linalg.norm(B_m, 2, 1))

        # @numba.njit()
        # def paired_correlation_rows(A, B):
        #     res = np.zeros(A.shape[0])
        #     for i in range(A.shape[0]):
        #         a = A[i,:] - np.sum(A[i,:]) / A.shape[1]
        #         b = B[i,:] - np.sum(B[i,:]) / B.shape[1]
        #         res[i] = np.sum(a * b) / (np.sqrt(np.sum(a*a)) * np.sqrt(np.sum(b*b)))
        #     return res

        tmp_filter = np.ones(self.gammas.shape, dtype=bool)
        if minR2 is not None:
            # NOTE Should be: tmp_filter = np.sqrt(self.R2) > minR but since the fit is weighted and constrained R2 can be negative
            R2_corrected = np.sqrt(np.abs(self.R2)) * np.sign(self.R2)
            tmp_filter = tmp_filter & (R2_corrected > minR2)
        if min_gamma is not None:
            tmp_filter = tmp_filter & (self.gammas > min_gamma)
        if minCorr is not None:
            Corr = paired_correlation_rows(self.Sx_sz, self.Ux_sz)
            tmp_filter = tmp_filter & (Corr > minCorr)
        # Perform the filtering
        self.ra = {k: v[tmp_filter] for k, v in self.ra.items()}
        matrixes2filter = ["U", "U_sz", "U_norm", "Ux", "Ux_sz", "Ux_norm",
                           "S", "S_sz", "S_norm", "Sx", "Sx_sz", "Sx_norm"]
        vectors2filter = ["gammas", "q", "R2"]
        for name_attr in matrixes2filter:
            if hasattr(self, name_attr):
                setattr(self, name_attr, getattr(self, name_attr)[tmp_filter, :])
        for name_attr in vectors2filter:
            if hasattr(self, name_attr):
                setattr(self, name_attr, getattr(self, name_attr)[tmp_filter])


    def predict_U(self, which_gamma: str="gammas", which_S: str= "Sx_sz", which_offset: str="q") -> None:
        """Predict U (gamma * S) given the gamma model fit

        Arguments
        ---------
        which_gamma: str, default="gammas"
            name of the attribute to use as gamma
        which_S: str, default="Sx_sz"
            name of the attribute to use as S
        which_offset: str, default="q"
            name of the attribute containing the offset

        Returns
        ------
        Noting but it creates the attribute
        Upred: np.ndarray
           unspliced estimated as `gamma * S`
        """
        self.which_S_for_pred = which_S
        if which_offset is None:
            if hasattr(self, "q_W") or hasattr(self, "q"):
                logging.warn("Predicting U without intercept but intercept was previously fit! Set which_offset='q' or 'q_W' ")
            self.Upred = getattr(self, which_gamma)[:, None] * getattr(self, which_S)
            # self.Upred = self.gammas[:, None] * self.Sx_sz
        else:
            self.Upred = getattr(self, which_gamma)[:, None] * getattr(self, which_S) + getattr(self, which_offset)[:, None]


    def calculate_velocity(self, kind: str="residual", eps: float=None, Myneighbours: float=None, VelocityFolder: str="/location/", UnsplicedFolder: str="/location/") -> None:
        """Calculate velocity
        Arguments
        ---------
        kind: str, default="residual"
            "residual" calculates the velocity as U_measured - U_predicted

        eps: float, default=None
            if specified, velocities with absolute value smaller than eps * range_of_U will be set to 0
            if None this step will be skipped
        Results
        -------
        Nothing but it creates the attribute:
        velocity: np.ndarray
            U_measured - U_predicted
        """
        if kind == "residual":
            if self.which_S_for_pred == "Sx_sz":
                self.velocity = self.Ux_sz - self.Upred
            elif self.which_S_for_pred == "Sx":
                self.velocity = self.Ux - self.Upred
            else:
                NotImplementedError(f"Not implemented with which_S = {self.which_S_for_pred}")
        else:
            raise NotImplementedError(f"Velocity calculation kind={kind} is not implemented")

        if eps:
            minimal_signed_res = self.Upred.max(1) * eps
            self.velocity[np.abs(self.velocity) < minimal_signed_res[:, None]] = 0
        rownames_genes = list(self.ra["Genes"])
        colnames_cells = list(self.ca["Patients"])
        genes_velocity = np.array(self.velocity)
        spliced_df = pd.DataFrame(genes_velocity, columns=colnames_cells, index=rownames_genes)
        file_name = "VelocityValues_EmbedNN" + str(Myneighbours) + ".csv"
        os.chdir(VelocityFolder)
        spliced_df.to_csv(file_name, sep=",", header=True)
        #self.Upred
        genes_PrUnsp = np.array(self.Upred)
        preunspdf = pd.DataFrame(genes_PrUnsp, columns=colnames_cells, index=rownames_genes)
        file_name = "PredictedUnspiced_EmbedNN" + str(Myneighbours) + ".csv"
        os.chdir(UnsplicedFolder)
        preunspdf.to_csv(file_name, sep=",", header=True)


    def calculate_shift(self, assumption: str="constant_velocity", delta_t: float=1) -> None:
        """Find the change (deltaS) in gene expression for every cell

        Arguments
        ---------
        assumption: str, default="constant_velocity"
            constant_velocity (described in the paper as Model I)
            constant_unspliced (described in the paper as Model II)
        delta_t: float, default=1
            the time step for extrapolation

        Returns
        -------
        Nothing it only creates the following attributes
        delta_S: np.ndarray
            The variation in gene expression
        """
        if assumption == "constant_velocity":
            self.delta_S = delta_t * self.velocity #where the timestep is incorporated into the calculation #here!!!!
        elif assumption == "constant_unspliced":
            # Ux_sz = self.Ux_sz - offset; Ux_sz[Ux_sz<0] = 0
            # maybe I should say ratio see below
            Ux_szo = self.Ux_sz - self.q[:, None]
            Ux_szo[Ux_szo < 0] = 0
            egt = np.exp(-self.gammas * delta_t)[:, None]
            self.delta_S = self.Sx_sz * egt + (1 - egt) * Ux_szo / self.gammas[:, None] - self.Sx_sz
        else:
            raise NotImplementedError(f"Assumption {assumption} is not implemented")



    def extrapolate_cell_at_t(self, delta_t: float=1, clip: bool=True, Myneighbours: float=1, location: str="/location/") -> None:
        """Extrapolate the gene expression profile for each cell after delta_t

        Arguments
        ---------
        delta_t: float, default=1
            the time step considered for the extrapolation
        clip: bool, default=True
            If True negative values are clipped to zero

        Returns
        -------
        Nothing but it creates the attributes:
        Sx_sz_t: np.ndarray
            the extrapolated expression profile
        used_delta_t: float
            stores delta_t for future usage
        """
        if self.which_S_for_pred == "Sx_sz":
            self.Sx_sz_t = self.Sx_sz + delta_t * self.delta_S #extrapolated expression at the timestep = expression profile + timestep x variation in expression
            if clip:
                self.Sx_sz_t = np.clip(self.Sx_sz_t, 0, None)
                self.used_delta_t = delta_t
        elif self.which_S_for_pred == "Sx":
            self.Sx_t = self.Sx + delta_t * self.delta_S
            if clip:
                self.Sx_t = np.clip(self.Sx_t, 0, None)
                self.used_delta_t = delta_t
        else:
            NotImplementedError("not implemented for other situations other than Sx or Sx_sz")
        #print("genes future expression level: ", self.Sx_sz_t) #matrix of genes in each cell exprapolated future state
        os.chdir(location)
        rownames_genes = list(self.ra["Genes"])
        colnames_cells = list(self.ca["Patients"])
        genes_spliced = np.array(self.Sx_sz_t)
        spliced_df = pd.DataFrame(genes_spliced, columns=colnames_cells, index=rownames_genes)
        file_name = "Future_gene_spliced_state_EmbedNN" + str(Myneighbours) + ".csv"
        spliced_df.to_csv(file_name, sep=",", header=True)


    def perform_TSNE(self, n_dims: int=2, perplexity: float=30, initial_pos: np.ndarray=None, #not used
                     theta: float=0.5, n_pca_dim: int=None, max_iter: int=1000) -> None:
        """Perform TSNE on the PCA using barnes hut approximation
        """
        # Perform TSNE
        logging.debug("Running bhtsne")
        if initial_pos is None:
            initial_pos = "random"
        bh_tsne = TSNE(n_components=n_dims, perplexity=perplexity, angle=theta, init=initial_pos, n_iter=max_iter)
        self.ts = bh_tsne.fit_transform(self.pcs[:, :n_pca_dim])

    def estimate_transition_prob(self, hidim: str="Sx_sz", embed: str="ts", transform: str="sqrt", #uses the principle components from the pca to calculate embedding coordinates for each point (embedding)
                                 ndims: int=None, n_sight: int=None, psc: float=None,
                                 knn_random: bool=False, sampled_fraction: float=0.3,
                                 sampling_probs: Tuple[float, float]=(0.5, 0.1), max_dist_embed: float=None,
                                 n_jobs: int=4, threads: int=None, calculate_randomized: bool=True,
                                 random_seed: int=15071990, **kwargs) -> None:
        """Use correlation to estimate transition probabilities for every cells to its embedding neighborhood

        Arguments
        ---------
        hidim: str, default="Sx_sz"
            The name of the attribute containing the high dimensional space. It will be retrieved as getattr(self, hidim)
            The updated vector at time t is assumed to be getattr(self, hidim + "_t")
            Appending .T to the string will transpose the matrix (useful in case we want to use S or Sx)
        embed: str, default="ts"
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
        transform: str, default="sqrt"
            The transformation that is applies on the high dimensional space.
            If None the raw data will be used
        ndims: int, default=None
            The number of dimensions of the high dimensional space to work with. If None all will be considered
            It makes sense only when using principal components
        n_sight: int, default=None (also n_neighbors)
            The number of neighbors to take into account when performing the projection
        psc: float, default=None
            pseudocount added in variance normalizing transform
            If None, 1 would be used for log, 0 otherwise
        knn_random: bool, default=True
            whether to random sample the neighborhoods to speedup calculation
        sampling_probs: Tuple, default=(0.5, 1)
        max_dist_embed: float, default=None
            CURRENTLY NOT USED
            The maximum distance allowed
            If None it will be set to 0.25 * average_distance_two_points_taken_at_random
        n_jobs: int, default=4
            number of jobs to calculate knn
            this only applies to the knn search, for the more time consuming correlation computation see threads
        threads: int, default=None
            The threads will be used for the actual correlation computation by default half of the total.
        calculate_randomized: bool, default=True
            Calculate the transition probabilities with randomized residuals.
            This can be plotted downstream as a negative control and can be used to adjust the visualization scale of the velocity field.
        random_seed: int, default=15071990
            Random seed to make knn_random mode reproducible

        Returns
        -------
        """
        numba_random_seed(random_seed)
        self.which_hidim = hidim #spliced states before timepoint same calculation used to calculate after just with different names for the same thing this is: self.Sx_sz
        if "n_neighbors" in kwargs:
            n_neighbors = kwargs.pop("n_neighbors")
            if len(kwargs) > 0:
                logging.warning(f"keyword arguments were passed but could not be interpreted {kwargs}")
        else:
            n_neighbors = None

        if n_sight is None and n_neighbors is None:
            n_neighbors = int(self.S.shape[1] / 5)

        if (n_sight is not None) and (n_neighbors is not None) and n_neighbors != n_sight:
            raise ValueError("n_sight and n_neighbors are different names for the same parameter, they cannot be set differently")

        if n_sight is not None and n_neighbors is None:
            n_neighbors = n_sight

        if psc is None:
            if transform == "log" or transform == "logratio":
                psc = 1.
            elif transform == "sqrt":
                psc = 1e-10  # for numerical stablity
            else:  # transform == "linear":
                psc = 0

        if knn_random:
            np.random.seed(random_seed)
            self.corr_calc = "knn_random"
            if "pcs" in hidim:  # sic
                hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
                hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
            else:
                if ndims is not None:
                    raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
                hi_dim = getattr(self, hidim)  # [:, :ndims] hi_dim is the
                hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
                if calculate_randomized:
                    self.delta_S_rndm = np.copy(self.delta_S)
                    permute_rows_nsign(self.delta_S_rndm)
                    hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm

            embedding = getattr(self, embed)
            self.embedding = embedding
            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs) #+1?
            nn.fit(embedding)  # NOTE should support knn in high dimensions
            self.embedding_knn = nn.kneighbors_graph(mode="connectivity")
            # Pick random neighbours and prune the rest
            neigh_ixs = self.embedding_knn.indices.reshape((-1, n_neighbors))
            p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
            p = p / p.sum()

            # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
            # resulting self.corrcoeff with different number of nonzero entry per row.
            # Not updated yet not to break previous analyses
            # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
            # I change it here since I am doing some breaking changes
            sampling_ixs = np.stack((np.random.choice(neigh_ixs.shape[1],
                                                      size=(int(sampled_fraction * (n_neighbors)),),
                                                      replace=False,
                                                      p=p) for i in range(neigh_ixs.shape[0])), 0)
            self.sampling_ixs = sampling_ixs
            neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
            nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
            self.embedding_knn = sparse.csr_matrix((np.ones(nonzero),
                                                    neigh_ixs.ravel(),
                                                    np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                                   shape=(neigh_ixs.shape[0],
                                                          neigh_ixs.shape[0]))
            logging.debug(f"Correlation Calculation '{self.corr_calc}'")
            if transform == "log":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
            elif transform == "logratio":
                log2hidim = np.log2(hi_dim + psc)
                delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
                self.corrcoef = colDeltaCorpartial(log2hidim, delta_hi_dim, neigh_ixs, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + psc) - log2hidim
                    self.corrcoef_random = colDeltaCorpartial(log2hidim, delta_hi_dim_rndm, neigh_ixs, threads=threads)
            elif transform == "linear":
                self.corrcoef = colDeltaCorpartial(hi_dim, hi_dim_t - hi_dim, neigh_ixs, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    self.corrcoef_random = colDeltaCorpartial(hi_dim, hi_dim_t_rndm - hi_dim, neigh_ixs, threads=threads)
            elif transform == "sqrt":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
            else:
                raise NotImplementedError(f"transform={transform} is not a valid parameter")
            np.fill_diagonal(self.corrcoef, 0)
            if np.any(np.isnan(self.corrcoef)):
                self.corrcoef[np.isnan(self.corrcoef)] = 1
                logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_random, 0)
                if np.any(np.isnan(self.corrcoef_random)):
                    self.corrcoef_random[np.isnan(self.corrcoef_random)] = 1
                    logging.warning("Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            logging.debug(f"Done Correlation Calculation")
        else:
            self.corr_calc = "full"
            if "pcs" in hidim:  # sic
                hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
                hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
            else:
                if ndims is not None:
                    raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
                hi_dim = getattr(self, hidim)  # [:, :ndims]
                hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
                if calculate_randomized:
                    self.delta_S_rndm = np.copy(self.delta_S)
                    permute_rows_nsign(self.delta_S_rndm)
                    hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
            embedding = getattr(self, embed)
            self.embedding = embedding
            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
            nn.fit(embedding)
            self.embedding_knn = nn.kneighbors_graph(mode="connectivity").toarray() #def of embedding_knn, AG has added in toarray()
            logging.debug("Correlation Calculation 'full'")
            if transform == "log":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
            elif transform == "logratio":
                log2hidim = np.log2(hi_dim + psc)
                delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
                self.corrcoef = colDeltaCor(log2hidim, delta_hi_dim, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + 1) - log2hidim
                    self.corrcoef_random = colDeltaCor(log2hidim, delta_hi_dim_rndm, threads=threads)
            elif transform == "linear":
                self.corrcoef = colDeltaCor(hi_dim, hi_dim_t - hi_dim, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    self.corrcoef_random = colDeltaCor(hi_dim, hi_dim_t_rndm - hi_dim, threads=threads)
            elif transform == "sqrt":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)

                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
            else:
                raise NotImplementedError(f"transform={transform} is not a valid parameter")
            np.fill_diagonal(self.corrcoef, 0)
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_random, 0)

            
            colnames_cells = list(self.ca["Patients"]) #shape of all: (rows=genes, cols=samples)
            rownames_genes = list(self.ra["Genes"])
            tp = np.array(hi_dim)
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #emat= hi_dim #expression
            file_name = "emat_tSNE.csv"
            tp_df.to_csv(file_name, sep=",", header=True) 
 
            sqrtdmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim)
            tp = np.array(sqrtdmat)
            file_name = "sqrt_dmat_tSNE.csv"
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  #dmat = np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim) #dmat are the ve,locities
            tp_df.to_csv(file_name, sep=",", header=True)  

            alt_dmat = hi_dim_t - hi_dim
            tp = np.array(alt_dmat)
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=rownames_genes)  # alt dmat: hi_dim_t - hi_dim (for linear)
            file_name = "linear_dmat_tSNE.csv"
            tp_df.to_csv(file_name, sep=",", header=True)  
   
            tp = np.array(self.corrcoef)
            tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)  # raw cc values withot any alt
            file_name = "CorrelationCoefficients_Raw_tSNE.csv"
            tp_df.to_csv(file_name, sep=",", header=True)  


    def calculate_embedding_shift(self, sigma_corr: float=0.05, expression_scaling: bool=False, Myneighbours: float=None, scaling_penalty: float=1., metadataToUse: str="metadata.csv", pathToTransitionProbabilityCode: str="/", grouptype: str="Strings", secondgroup: str="metadataVersion.csv", manyfeatures: str="True", grouptype2: str="Strings", Randomsampling: str="False", extraNNvalue: int=5, DE_location: str="/location/", selfCC_folder: str="/location/", VecX_location: str="/location/", VecY_loc: str="/location/", TP_folder: str="/location/", tSNE_GroupLevelTP: str="/location/", DoSecondPred: str="no", UserThreshold: float=0.5) -> None: #Use the transition probability to project the velocity direction on the embedding
        """Use the transition probability to project the velocity direction on the embedding

        Arguments
        ---------
        sigma_corr: float, default=0.05
            the kernel scaling
        expression_scaling: bool, default=True
            rescale arrow intensity penalizing arrows that explain very small amount of expression differences
        scaling_penalty: float, default=1
            Higher values correspond to a stronger penalty


        Returns
        -------
        Nothing but it creates the following attributes:
        transition_prob: np.ndarray
            the transition probability calculated using the exponential kernel on the correlation coefficient
        delta_embedding: np.ndarray
            The resulting vector
        """
        # Kernel evaluation
        logging.debug("Calculate transition probability")
        if self.corr_calc == "full" or self.corr_calc == "knn_random":
            # NOTE maybe sparse matrix here are slower than dense
            # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
            self.transition_prob = np.exp(self.corrcoef / sigma_corr) * self.embedding_knn  # naive, the exponential of the correlation coefficient/kernal scaling, AG removed .A 
            self.transition_prob /= self.transition_prob.sum(1)[:, None]
            if hasattr(self, "corrcoef_random"):
                logging.debug("Calculate transition probability for negative control")
                self.transition_prob_random = np.exp(self.corrcoef_random / sigma_corr) * self.embedding_knn  # naive, AG removed .A 
                self.transition_prob_random /= self.transition_prob_random.sum(1)[:, None]
            unitary_vectors = self.embedding.T[:, None, :] - self.embedding.T[:, :, None]  # shape (2,ncells,ncells) #for each cells coordinates - row1: cell1-cell1, cell2-cell1, cell3-cell1 cell4-cell1
            with np.errstate(divide='ignore', invalid='ignore'):
                unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2
                np.fill_diagonal(unitary_vectors[0, ...], 0)  # fix nans replace NAN with 0
                np.fill_diagonal(unitary_vectors[1, ...], 0) #embedding comes from the tSNE
            self.delta_embedding = (self.transition_prob * unitary_vectors).sum(2) #.A means change the data type from a matrix to an array #unit vector*transition probability
            self.delta_embedding -= (self.embedding_knn * unitary_vectors).sum(2) / self.embedding_knn.sum(1).T #AG removed .A 
            self.delta_embedding = self.delta_embedding.T #transposes the vector
            if expression_scaling:
                hi_dim = getattr(self, self.which_hidim)
                estim_delta = hi_dim.dot(self.transition_prob.T) - hi_dim.dot((self.embedding_knn / self.embedding_knn.sum(1)).T) #this is the same problem as above, AG removed .A 
                cos_proj = (self.delta_S * estim_delta).sum(0) / np.sqrt((estim_delta**2).sum(0))
                self.scaling = np.clip(cos_proj / scaling_penalty, 0, 1)
                self.delta_embedding = self.delta_embedding * self.scaling[:, None]

            if hasattr(self, "corrcoef_random"):
                self.delta_embedding_random = (self.transition_prob_random * unitary_vectors).sum(2)
                self.delta_embedding_random -= (self.embedding_knn * unitary_vectors).sum(2) / self.embedding_knn.sum(1).T #AG removed .A 
                self.delta_embedding_random = self.delta_embedding_random.T
                if expression_scaling:
                    estim_delta_rndm = hi_dim.dot(self.transition_prob_random.T) - hi_dim.dot((self.embedding_knn / self.embedding_knn.sum(1)()).T) #AG removed .A 
                    cos_proj_rndm = (self.delta_S_rndm * estim_delta_rndm).sum(0) / np.sqrt((estim_delta_rndm**2).sum(0))
                    self.scaling_rndm = np.clip(cos_proj_rndm / scaling_penalty, 0, 1)
                    self.delta_embedding_random = self.delta_embedding_random * self.scaling_rndm[:, None]
        else:
            # NOTE should implement a version with cython
            raise NotImplementedError(f"Weird value self.corr_calc={self.corr_calc}")
        colnames_cells = list(self.ca["Patients"])
        genes_embedding = np.array(self.delta_embedding)
        spliced_df = pd.DataFrame(genes_embedding, columns=["x", "y"], index=colnames_cells)
        file_name = "Delta_embedding_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_tSNE.csv"
        os.chdir(DE_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(self.corrcoef)
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "selfCC_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_tSNE.csv"
        os.chdir(selfCC_folder)
        spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[0])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecX_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_tSNE.csv"
        os.chdir(VecX_location)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        genes_embedding = np.array(unitary_vectors[1])
        spliced_df = pd.DataFrame(genes_embedding, columns=colnames_cells, index=colnames_cells)
        file_name = "UnVecY_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_tSNE.csv"
        os.chdir(VecY_loc)
        #spliced_df.to_csv(file_name, sep=",", header=True)
        tp = np.array(self.transition_prob) #number of samples x number of samples
        tp_df = pd.DataFrame(tp, columns=colnames_cells, index=colnames_cells)
        file_name = "Transition_Probability_NN" + str(Myneighbours) + "_" + str(extraNNvalue) + "_tSNE.csv"
        os.chdir(TP_folder)
        tp_df.to_csv(file_name, sep=",", header=True)
        #####RUN TP code here################################
        #note: the following requires these R packages to be installed:  'gtools', 'reshape2', ggplot2', grid' and 'tidyverse'
        #to code to install them has been commented out but is below:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        utils = importr('utils')
        #utils.install_packages('gtools')
        #utils.install_packages('reshape2')
        #utils.install_packages('ggplot2')        
        #utils.install_packages('grid')
        #utils.install_packages('tidyverse')
        gtools = importr("gtools")
        reshape2 = importr("reshape2")
        ggplot2 = importr("ggplot2")
        grid = importr("grid")
        tidyverse = importr("tidyverse")
        path=pathToTransitionProbabilityCode
        version = "SampleGroup1"
        def RunTransitionProbabilityCode(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold):
            r=ro.r
            r.source(path+"TransitionProbability.R")
            p=r.TransitionProbability(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold)
            return p
        #a=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, metadataToUse, tSNE_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold) 
        b=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, metadataToUse, tSNE_GroupLevelTP, grouptype, version, Randomsampling, UserThreshold)
        if manyfeatures == "True" and DoSecondPred == "yes":
            version = "SampleGroup2" #secondgroup
            #c=RunTransitionProbabilityCode(file_name, "TRUE", TP_folder, secondgroup, tSNE_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold) 
            d=RunTransitionProbabilityCode(file_name, "FALSE", TP_folder, secondgroup, tSNE_GroupLevelTP, grouptype2, version, Randomsampling, UserThreshold)
       

 
    def calculate_grid_arrows(self, embed: str="embedding", smooth: float=0.5, steps: Tuple=(40, 40),
                              n_neighbors: int=100, n_jobs: int=4) -> None:
        """Calculate the velocity using a points on a regular grid and a gaussian kernel
        Note: the function should work also for n-dimensional grid
        Arguments
        ---------
        embed: str, default=embedding
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
            The difference vector is getattr(self, 'delta' + '_' + embed)
        smooth: float, smooth=0.5
            Higher value correspond to taking in consideration further points
            the standard deviation of the gaussian kernel is smooth * stepsize
        steps: tuple, default
            the number of steps in the grid for each axis
        n_neighbors:
            number of neighbors to use in the calculation, bigger number should not change too much the results..
            ...as soon as smooth is small
            Higher value correspond to slower execution time
        n_jobs:
            number of processes for parallel computing
        Returns
        -------
        Nothing but it sets the attributes:
        flow_embedding: np.ndarray
            the coordinates of the embedding
        flow_grid: np.ndarray
            the gridpoints
        flow: np.ndarray
            vector field coordinates
        flow_magnitude: np.ndarray
            magnitude of each vector on the grid
        total_p_mass: np.ndarray
            density at each point of the grid
        """
        embedding = getattr(self, embed)
        if hasattr(self, f"delta_{embed}"):
            delta_embedding = getattr(self, f"delta_{embed}")
            if hasattr(self, "corrcoef_random"):
                delta_embedding_random = getattr(self, f"delta_{embed}_random")
        else:
            raise KeyError("This embedding does not have a delta_*")
        # Prepare the grid
        grs = []
        for dim_i in range(embedding.shape[1]):
            m, M = np.min(embedding[:, dim_i]), np.max(embedding[:, dim_i])
            m = m - 0.025 * np.abs(M - m)
            M = M + 0.025 * np.abs(M - m)
            gr = np.linspace(m, M, steps[dim_i])
            grs.append(gr)

        meshes_tuple = np.meshgrid(*grs)
        gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

        nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
        nn.fit(embedding)
        dists, neighs = nn.kneighbors(gridpoints_coordinates)

        std = np.mean([(g[1] - g[0]) for g in grs])
        #isotropic gaussian kernel
        gaussian_w = normal.pdf(loc=0, scale=smooth * std, x=dists)
        self.total_p_mass = gaussian_w.sum(1)
        UZ = (delta_embedding[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, self.total_p_mass)[:, None]  # weighed average, delta_embedding is calculated by calculate embedding shift
        magnitude = np.linalg.norm(UZ, axis=1)
        #Assign attributes
        self.flow_embedding = embedding
        self.flow_grid = gridpoints_coordinates
        self.flow = UZ
        self.flow_norm = UZ / np.percentile(magnitude, 99.5)
        self.flow_norm_magnitude = np.linalg.norm(self.flow_norm, axis=1) #1600 length
        if hasattr(self, "corrcoef_random"):
            UZ_rndm = (delta_embedding_random[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, self.total_p_mass)[:, None]  # weighed average
            magnitude_rndm = np.linalg.norm(UZ, axis=1)
            #Assign attributes
            self.flow_rndm = UZ_rndm
            self.flow_norm_rndm = UZ_rndm / np.percentile(magnitude_rndm, 99.5)
            self.flow_norm_magnitude_rndm = np.linalg.norm(self.flow_norm_rndm, axis=1)

    def prepare_markov(self, sigma_D: np.ndarray, sigma_W: np.ndarray, direction: str="forward", cells_ixs: np.ndarray=None) -> None:
        """Prepare a transition probability for Markov process

        Arguments
        ---------
        sigma_D: float
            the standard deviation used on the locality-limiting component
        sigma_W: float
            the standard deviation used on the noise component
        direction: str, default="backwards"
            whether to diffuse forward of backwards
        cells_ixs: np.ndarray, default=None
            Cells to use, if None all the cells will be considered.
        Returns
        -------
        Nothing but it creates the following attributes:
        tr: np.ndarray
            the transition probability matrix

        """
        if cells_ixs is None:
            cells_ixs = np.arange(self.transition_prob.shape[0])

        # NOTE: This implementation is not speed optimized to improve the speed of the implementation:
        # - the C/Fortran contiguity of the transition matrix should be taken into account
        # - a knn implementation would reduce computation
        # - should avoid transformation to and from dense-sparse formats
        if direction == "forward":
            self.tr = np.array(self.transition_prob[cells_ixs, :][:, cells_ixs])
        elif direction == "backwards":
            self.tr = np.array((self.transition_prob[cells_ixs, :][:, cells_ixs]).T, order="C")
        else:
            raise NotImplementedError(f"{direction} is not an implemented direction")
        dist_matrix = squareform(pdist(self.embedding[cells_ixs, :]))
        K_D = gaussian_kernel(dist_matrix, sigma=sigma_D)
        self.tr = self.tr * K_D
        # Fill diagonal with max or the row and sum=1 normalize
        np.fill_diagonal(self.tr, self.tr.max(1))
        self.tr = self.tr / self.tr.sum(1)[:, None]

        K_W = gaussian_kernel(dist_matrix, sigma=sigma_W)
        K_W = K_W / K_W.sum(1)[:, None]
        self.tr = 0.8 * self.tr + 0.2 * K_W
        self.tr = self.tr / self.tr.sum(1)[:, None]
        self.tr = scipy.sparse.csr_matrix(self.tr)


    def run_markov(self, starting_p: np.ndarray=None, n_steps: int=2500, mode: str="time_evolution") -> None: #cd2514 not currently being used but could be useful
    #    Run a Markov process

    #    Arguments
    #    ---------
    #    starting_p: np.ndarray, default=None
    #        specifies the starting density
    #        if None is passed an array of 1/self.tr.shape[0] will be created
    #    n_steps: np.ndarray, default=2500
    #        Numbers of steps to be performed
    #    mode: str, default="time_evolution"
    #        this argument is passed to the Diffusion.diffuse call

    #    Returns
    #    -------
    #    Nothing but it creates the attribute:
    #    diffused: np.ndarray
    #        The probability to be found at any of the states
        if starting_p is None:
            starting_p = np.ones(self.tr.shape[0]) / self.tr.shape[0]
        diffusor = Diffusion()
        self.diffused = diffusor.diffuse(starting_p, self.tr, n_steps=n_steps, mode=mode)[0]

    def _plot_phase_portrait(self, gene: str, gs_i: Any=None) -> None:
        """Plot spliced-unspliced scatterplot resembling phase portrait
        """
        if gene is None:
            plt.subplot(111)
        else:
            plt.subplot(gs_i)
        ix = np.where(self.ra["Genes"] == gene)[0][0]
        scatter_viz(self.Sx_sz[ix, :], self.Ux_sz[ix, :], c=self.colorandum, s=5, alpha=0.4)
        plt.title(gene)
        xnew = np.linspace(0, self.Sx_sz[ix, :].max())
        plt.plot(xnew, self.gammas[ix] * xnew + self.q[ix], c="k")

    def plot_phase_portraits(self, genes: List[str]) -> None:
        """Plot spliced-unspliced scatterplots resembling phase portraits

        Arguments
        ---------
        genes: List[str]
            A list of gene symbols.
        """
        genes = self.ra["Genes"]
        n = len(genes)
        sqrtn = int(np.ceil(np.sqrt(n)))
        gs = plt.GridSpec(sqrtn, int(np.ceil(n / sqrtn)))
        for i, gn in enumerate(genes):
            self._plot_phase_portrait(gn, gs[i])
        plt.savefig("phase_portraits.png")
        plt.close()

    def one_arrow(self, numberNeighb, typegroup: str="integers", extraNNvalue: int=5, ChangeShape: str="True", Colours: str="Paired"):
        if ChangeShape == "True":
            plt.figure(None,(16,10))
            quiver_scale = 60
            my_markers = self.ca["Markers"]
            mark = sorted(set(my_markers))
            corresp = ["o", "s", "^", "P", "D", "X", "*", "p", "+", ">", "<"]
            count = 0
            my_types = {}
            for the_type in mark:
                my_types[the_type] = corresp[count]
                count = count + 1
            markers = []
            for mar in my_markers:
                markers.append(my_types.get(mar))
            z = randn(10)
            #ax = fig.add_subplot(111, projection='3d')
            for i in range(0, len(markers)):
                #plt.scatter(self.embedding[i, 0], self.embedding[i, 1],c="0.8", alpha=1, s=120, edgecolor="")
                plt.scatter(self.embedding[i, 0], self.embedding[i, 1], marker=markers[i], color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                plt.quiver(self.embedding[i, 0], self.embedding[i, 1], self.delta_embedding[i, 0], self.delta_embedding[i, 1], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]   # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            for sample_type in range(0, len(mark)):
                to_add = Line2D([0], [0], marker=[*my_types.values()][sample_type], color='w', label=str([*my_types.keys()][sample_type]).strip('"'), markerfacecolor='black', markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            #ax.set_xlim(0, 10)
            #plt.axis("off")
            plt.xlabel(r"$t$-SNE1", fontsize=12)
            plt.ylabel(r"$t$-SNE2", fontsize=12)
            figure_name = 'tSNE_Arrows_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #plt.show()
            plt.close()
        if ChangeShape == "False":
            plt.figure(None,(16,10))
            quiver_scale = 60
            z = randn(10)
            #ax = fig.add_subplot(111, projection='3d')
            #plt.scatter(self.embedding[:, 0], self.embedding[:, 1],c="0.8", alpha=1, s=120, edgecolor="")
            plt.scatter(self.embedding[:, 0], self.embedding[:, 1], marker="o", color=self.colorandum, alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
            quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
            plt.quiver(self.embedding[:, 0], self.embedding[:, 1], self.delta_embedding[:, 0], self.delta_embedding[:, 1], scale=quiver_scale, **quiver_kwargs)
            #for i in range(0, len(self.colorandum)):
                #plt.scatter(self.embedding[i, 0], self.embedding[i, 1],c="0.8", alpha=1, s=120, edgecolor="")
                #plt.scatter(self.embedding[i, 0], self.embedding[i, 1], marker="o", color=self.colorandum[i], alpha=1, s=120, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)
                #quiver_kwargs=dict(headaxislength=7, headlength=13, headwidth=10,linewidths=0.5, width=0.005,edgecolors="k", color="black", alpha=1, angles='xy')
                #plt.quiver(self.embedding[i, 0], self.embedding[i, 1], self.delta_embedding[i, 0], self.delta_embedding[i, 1], scale=quiver_scale, **quiver_kwargs)
            hexy_codes = []
            my_keys = [k  for  k in  self.cluster_colors_dict.keys()]
            decoded_my_keys = []
            for item in my_keys:
                decoded_my_keys.append(str(item,'utf-8'))#,'utf-8')) #int()
            cmap2 = matplotlib.colormaps[Colours]    # PiYG
            for i in range(cmap2.N):
                rgb = cmap2(i)[:3] # will return rgba, we take only first 3 so we get rgb
                hexy_codes.append(matplotlib.colors.rgb2hex(rgb))
            legend_elements = []
            combined = list(zip(decoded_my_keys, hexy_codes))
            combined = sorted(combined, key=itemgetter(0))
            if (typegroup != "Strings"):
                timepoints = [int(i[0]) for i in combined]
                colys = [i[1] for i in combined]
                tempDF = pd.DataFrame({'Group':timepoints, 'colours':colys})
                tempDF = tempDF.sort_values('Group')
                groupys2 = tempDF['Group'].tolist()
                cols2 = tempDF['colours'].tolist()
                combined = list(zip(groupys2, cols2))
            decoded_my_keys = [i[0] for i in combined]
            hexy_codes = [i[1] for i in combined]
            for sample_type in range(0, len(decoded_my_keys)):
                to_add = Line2D([0], [0], marker='s', color='w', label=str(decoded_my_keys[sample_type]).strip('"'), markerfacecolor=str(hexy_codes[sample_type]), markersize=7)
                legend_elements.append(to_add)
            plt.legend(handles=legend_elements, fancybox=True, framealpha=1, borderpad=1, loc=(1.01,0.5), title='Key')
            #ax.set_xlim(0, 10)
            #plt.axis("off")
            plt.xlabel(r"$t$-SNE1", fontsize=12)
            plt.ylabel(r"$t$-SNE2", fontsize=12)
            figure_name = 'tSNE_Arrows_NN' + str(numberNeighb) + "_" + str(extraNNvalue) + '.png'
            plt.savefig(figure_name)
            #plt.show()
            plt.close()


    def plot_grid_arrows(self, quiver_scale: Union[str, float]="auto", scale_type: str= "relative", min_mass: float=0.0, min_magnitude: float=0.0,
                         scatter_kwargs_dict: Dict= None, plot_dots: bool=False, plot_random: bool=False, **quiver_kwargs: Any) -> None:
        """Plots vector field averaging velocity vectors on a grid
        Arguments
        ---------
        quiver_scale: float, default="auto"
            Rescaling factor applied to the arrow field to enhance visibility
            If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
            If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
            NOTE: In the situation where "auto" is resulting in very small or big velocities, pass a float to this parameter
            The float will be interpreted as a scaling, importantly both the data and the control will be scaled
            in this way you can rescale the velocity arbitrarily without the risk of observing just an overfit of the noise
        scale_type: str, default="relative"
            How to interpret `quiver_scale`:
            If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
            (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
            If "absolute" the value will be passed to the matplotlib quiver function (not recommended if you are not sure what this implies)
        min_mass: float, default=1
            the minimum density around a grid point for it to be considered and plotted
        min_magnitude: float, default=None
            the minimum magnitude of the velocity for it to be considered and plotted
        scatter_kwargs_dict: dict, default=None
            a dictionary of keyword arguments to pass to scatter
            by default the following are passed: s=20, zorder=-1, alpha=0.2, lw=0, c=self.colorandum. But they can be overridden.
        plot_dots: bool, default=True
            whether to plot dots in correspondence of all low velocity grid points
        plot_random: bool, default=True
            whether to plot the randomized control next to the plot
        **quiver_kwargs: dict
            keyword arguments to pass to quiver
            By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.
        """
        # plt.figure(figsize=(10, 10))
        _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
        _quiver_kwargs.update(quiver_kwargs)

        scatter_dict = {"s": 20, "zorder": -1, "alpha": 0.2, "lw": 0, "c": self.colorandum}
        if scatter_kwargs_dict is not None:
            scatter_dict.update(scatter_kwargs_dict)

        # Determine quiver scale
        if scale_type == "relative":
            if hasattr(self, "flow_rndm"):
                plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0), 2)  # Diagonal of the plot
                arrows_scale = np.percentile(np.linalg.norm(self.flow_rndm[self.total_p_mass >= min_mass, :], 2, 1), 90)  # Typical lenght of an arrow
                if quiver_scale == "auto":
                    quiver_scale = arrows_scale / (plot_scale * 0.0025)
                else:
                    quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.0025)
            else:
                raise ValueError(""""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
                Please run estimate_transition_prob or set `scale_type` to `absolute`""")
        else:
            logging.warning("The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")
        mass_filter = self.total_p_mass < min_mass
        if min_magnitude is None:
            XY, UV = np.copy(self.flow_grid), np.copy(self.flow) #flow grid has current x, y coordinates, flow or flow_norm has the thing to be added on
            if not plot_dots:
                UV = UV[~mass_filter, :]
                XY = XY[~mass_filter, :]
            else:
                UV[mass_filter, :] = 0
        else:
            XY, UV = np.copy(self.flow_grid), np.copy(self.flow_norm)
            if not plot_dots:
                UV = UV[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
                XY = XY[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
            else:
                UV[mass_filter | (self.flow_norm_magnitude < min_magnitude), :] = 0

        if plot_random:
            if min_magnitude is None:
                XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_rndm)
                if not plot_dots:
                    UV_rndm = UV_rndm[~mass_filter, :]
                    XY = XY[~mass_filter, :]
                else:
                    UV_rndm[mass_filter, :] = 0
            else:
                XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_norm_rndm)
                if not plot_dots:
                    UV_rndm = UV_rndm[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
                    XY = XY[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
                else:
                    UV_rndm[mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude), :] = 0

            plt.subplot(122)
            plt.title("Randomized")
            plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
            plt.quiver(XY[:, 0], XY[:, 1], UV_rndm[:, 0], UV_rndm[:, 1],
                       scale=quiver_scale, zorder=20000, **_quiver_kwargs) #XY[:, 0] + UV_rndm[:, 0] = end of arrow x coordinate, XY[:, 1] + UV_rndm[:, 1] = end of arrow y coordinate
            plt.axis("off")
            plt.subplot(121)
            plt.title("Data")
        #will need third coordinate to be added to embedding
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter()
        #ax.set_xlabel('x-tSNE')
        #ax.set_ylabel('y-tSNE')
        #ax.set_zlabel('z-tSNE')
        plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
        plt.quiver(XY[:, 0], XY[:, 1], UV[:, 0], UV[:, 1],
                   scale=quiver_scale, zorder=20000, **_quiver_kwargs)  #XY[:, 0] + UV[:, 0] = end of arrow x coordinate, XY[:, 1] + UV[:, 1] = end of arrow y coordinate scale: scale : None, float, optional #Number of data units per arrow length unit, e.g., m/s per plot width; a smaller scale parameter makes the arrow longer. Default is None.
        plt.axis("off")
        plt.savefig("2D_tSNE.png")
        #plt.show()


    def plot_arrows_embedding(self, choice: Union[str, int]=10, quiver_scale: Union[str, float]='auto', scale_type: str="relative",
                              plot_scatter: bool=False, scatter_kwargs: Dict={}, color_arrow: str="cluster",
                              new_fig: bool=False, plot_random: bool=True, **quiver_kwargs: Any) -> None:
        """Plots velocity on the embedding cell-wise

        Arguments
        ---------
        choice: int, default = "auto"
            the number of cells to randomly pick to plot the arrows (To avoid overcrowding)
        quiver_scale: float, default="auto"
            Rescaling factor applied to the arrow field to enhance visibility
            If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
            If a float is provided the interpretation of the value depends on the parameter `plt.shoscale_type`, see below.
            NOTE: Despite a similar option than plot_grid_arrows, here there is no strong motivation to calculate the scale relative to the randomized control
            This is because the randomized doesn't have to have smaller velocity cell-wise, there might be for example
            scattered cells that will have strong velocity but they will, correctly just average out when calculating the average velocity field.
        scale_type: str, default="relative"
            How to interpret `quiver_scale`:
            If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
            (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
            If "absolute" the value will be passed to the matplotlib quiver function
        plot_scatter: bool, default = False
            whether to plot the points
        scatter_kwargs: Dict
            A dictionary containing all the keywords arguments to pass to matplotlib scatter
            by default the following are passed: c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3. But they can be overridden.
        color_arrow: str, default = "cluster"
            the color of the arrows, if "cluster" the arrows are colored the same as the cluster
        new_fig: bool, default=False
            whether to create a new figure
        plot_random: bool, default=True
            whether to plot the randomized control next to the plot
        **quiver_kwargs: dict
            keyword arguments to pass to quiver
            By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.

        Returns
        -------
        Nothing, just plots the tsne with arrows
        """
        if choice == "auto":
            choice = int(self.S.shape[1] / 3)
            logging.warning(f"Only {choice} arrows will be shown to avoid overcrowding, you can choose the exact number setting the `choice` argument")
        _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
        _scatter_kwargs = dict(c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3)
        _scatter_kwargs.update(scatter_kwargs)
        if new_fig:
            if plot_random and hasattr(self, "delta_embedding_random"):
                plt.figure(figsize=(22, 12))
            else:
                plt.figure(figsize=(14, 14))

        #ix_choice = np.random.choice(self.embedding.shape[0], size=choice, replace=False)
        ix_choice = np.arange(0,choice)*1+0
        # Determine quiver scale
        if scale_type == "relative":
            if hasattr(self, "delta_embedding_random"):
                plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0), 2)  # Diagonal of the plot
                arrows_scale = np.percentile(np.linalg.norm(self.delta_embedding_random, 2, 1), 80)  # Typical length of an arrow
                if quiver_scale == "auto":
                    quiver_scale = arrows_scale / (plot_scale * 0.005)
                else:
                    quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.005)
            else:
                raise ValueError("""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
                Please run estimate_transition_prob or set `scale_type` to `absolute`""")
        else:
            logging.warning("The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")

        if color_arrow == "cluster":
            colorandum = self.colorandum[ix_choice, :]
        else:
            colorandum = color_arrow
        _quiver_kwargs.update({"color": colorandum})
        _quiver_kwargs.update(quiver_kwargs)



        if plot_random and hasattr(self, "delta_embedding_random"):
            plt.subplot(122)
            plt.title("Randomized")
            if plot_scatter:
                plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)
            plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
                       self.delta_embedding_random[ix_choice, 0], self.delta_embedding_random[ix_choice, 1],
                       scale=quiver_scale, **_quiver_kwargs)
            plt.axis("off")
            plt.subplot(121)
            plt.title("Data")
        if plot_scatter:
            plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)

        plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1], #self.embedding[ix_choice, 0] + self.delta_embedding[ix_choice, 0] = new coordinate (end of arrow)
                   self.delta_embedding[ix_choice, 0], self.delta_embedding[ix_choice, 1], #self.embedding[ix_choice, 1] + self.delta_embedding[ix_choice, 1] = new coordinate
                   scale=quiver_scale, **_quiver_kwargs)
        #plt.savefig("with_arrows.png")
        plt.ylabel('y-tSNE')
        plt.xlabel('x-tSNE')
        plt.axis("off")
        plt.savefig("tSNE.pbg==ng")
        #plt.show()
        #plt.close()


    def plot_cell_transitions(self, cell_ix: int=0, alpha: float=0.1, alpha_neigh: float=0.2,
                              cmap_name: str="RdBu_r", plot_arrow: bool=True,
                              mark_cell: bool=True, head_width: int=1) -> None:
        """Plot the probability of a cell to transition to any other cell

        This function is untested
        """
        cmap = matplotlib.colormaps[cmap_name]
        #colorandum = np.ones((self.embedding.shape[0], 4))
        #colorandum *= 0.3
        #colorandum[:, -1] = alpha
        ix_choice = np.random.choice(self.embedding.shape[0], size=int(self.embedding.shape[0]/1.), replace=False)
        colorandum = self.colorandum[ix_choice, :]
        plt.scatter(self.embedding[:, 0], self.embedding[:, 1], c=colorandum, s=50, edgecolor="")#self.embedding[:, 0], self.embedding[:, 1] are the x and y coordinates of each point on the tsne graph
        plt.axis("off")
        #plt.xlabel("x-tSNE")
        #plt.ylabel("y-tSNE")
        if mark_cell:
            plt.scatter(self.embedding[cell_ix, 0], self.embedding[cell_ix, 1],
                        facecolor="none", s=100, edgecolor="k")
        if plot_arrow:
            plt.arrow(self.embedding[cell_ix, 0], self.embedding[cell_ix, 1],
                      self.delta_embedding[cell_ix, 0], self.delta_embedding[cell_ix, 1],
                      head_width=head_width, length_includes_head=True) #x,y coordinates of arrow
        plt.savefig("tSNE_cell_transitions.png")
        #plt.show()
        plt.close()

    def plot_velocity_as_color(self, gene_name: str=None, cmap: Any= plt.cm.RdBu_r,
                               gs: Any=None, which_tsne: str="ts", **kwargs: Dict) -> None:
        """Plot velocity as color on the Tsne embedding

        Arguments
        ---------
        gene_name: str
            The name of the gene, should be present in self.S
        cmap: maplotlib.cm.Colormap, default=maplotlib.cm.RdBu_r
            Colormap to use, divergent ones are better, RdBu_r is default
            Notice that 0 will be always set as the center of the colormap. (e.g. white in RdBu_r)
        gs: Gridspec subplot
            Gridspec subplot to plot on.
        which_tsne: str, default="ts"
            the name of the attributed where the desired embedding is stored
        **kwargs: dict
            other keywords arguments will be passed to the plt.scatter call

        Returns
        -------
        Nothing
        """
        ix = np.where(self.ra["Genes"] == gene_name)[0][0]
        kwarg_plot = {"alpha": 0.5, "s": 8, "edgecolor": "0.8", "lw": 0.15}
        kwarg_plot.update(kwargs)
        if gs is None:
            fig = plt.figure(figsize=(10, 10))
            plt.subplot(111)
        else:
            plt.subplot(gs)

        tsne = getattr(self, which_tsne)
        if self.which_S_for_pred == "Sx_sz":
            tmp_colorandum = self.Sx_sz_t[ix, :] - self.Sx_sz[ix, :]
        else:
            tmp_colorandum = self.Sx_t[ix, :] - self.Sx[ix, :]
        if (np.abs(tmp_colorandum) > 0.00005).sum() < 10:  # If S vs U scatterplot it is flat
            print("S vs U scatterplot it is flat")
            return
        limit = np.max(np.abs(np.percentile(tmp_colorandum, [1, 99])))  # upper and lowe limit / saturation
        tmp_colorandum = tmp_colorandum + limit  # that is: tmp_colorandum - (-limit)
        tmp_colorandum = tmp_colorandum / (2 * limit)  # that is: tmp_colorandum / (limit - (-limit))
        tmp_colorandum = np.clip(tmp_colorandum, 0, 1)

        scatter_viz(tsne[:, 0], tsne[:, 1],
                    c=cmap(tmp_colorandum), **kwarg_plot)
        plt.axis("off")
        plt.title(f"{gene_name}")


def array_to_rmatrix(X):
    import rpy2.robjects as robj
    nr, nc = X.shape
    xvec = robj.FloatVector(X.transpose().reshape((X.size)))
    xr = robj.r.matrix(xvec, nrow=nr, ncol=nc)
    return xr

def principal_curve(X, pca=True):
    """
    input : numpy.array
    returns:
    Result::Object
        Methods:
        projections - the matrix of the projectiond
        ixsort - the order ot the points (as in argsort)
        arclength - the lenght of the arc from the beginning to the point
    """
    # convert array to R matrix
    xr = array_to_rmatrix(X)
    if pca:
        #perform pca
        import rpy2.robjects as robj
        t = robj.r.prcomp(xr)
        #determine dimensionality reduction
        usedcomp = max( sum( np.array(t[t.names.index('sdev')]) > 1.1) , 4)
        usedcomp = min([usedcomp, sum( np.array(t[t.names.index('sdev')]) > 0.25), X.shape[0]])
        Xpc = np.array(t[t.names.index('x')])[:,:usedcomp]
        # convert array to R matrix
        xr = array_to_rmatrix(Xpc)
    #import the correct namespace
    #utils = importr('utils')
    #utils.install_packages('princurve')
    princurve = importr("princurve", on_conflict="warn")
    #call the function
    fit1 = princurve.principal_curve(xr)
    #extract the outputs
    class Results:
        pass
    results = Results()
    results.projections = np.array( fit1[0] )
    results.ixsort = np.array( fit1[1] ) - 1 # R is 1 indexed
    results.arclength = np.array( fit1[2] )
    results.dist = np.array( fit1[3] )
    if pca:
        results.PCs = np.array(xr) #only the used components

    return results

def scatter_viz(x: np.ndarray, y: np.ndarray, *args: Any, **kwargs: Any) -> Any:
    """A wrapper of scatter plot that guarantees that every point is visible in a very crowded scatterplot

    Args
    ----
    x: np.ndarray
        x axis coordinates
    y: np.ndarray
        y axis coordinates
    args and kwargs:
        positional and keyword arguments as in matplotplib.pyplot.scatter

    Returns
    -------
    Plots the graph and returns the axes object
    """
    ix_x_sort = np.argsort(x, kind="mergesort")
    ix_yx_sort = np.argsort(y[ix_x_sort], kind="mergesort")
    args_new = []
    kwargs_new = {}
    for arg in args:
        if type(arg) is np.ndarray:
            args_new.append(arg[ix_x_sort][ix_yx_sort])
        else:
            args_new.append(arg)
    for karg, varg in kwargs.items():
        if type(varg) is np.ndarray:
            kwargs_new[karg] = varg[ix_x_sort][ix_yx_sort]
        else:
            kwargs_new[karg] = varg
    ax = plt.scatter(x[ix_x_sort][ix_yx_sort], y[ix_x_sort][ix_yx_sort], *args_new, **kwargs_new)
    return ax

def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool=True) -> np.ndarray:
    "This is super duper magic sauce to make the order of one list to be like another"
    if check_content:
        assert len(np.intersect1d(a, b)) == len(a), f"The two arrays are not matching"
    return np.argsort(a)[np.argsort(np.argsort(b))]


colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))


def colormap_fun(x: np.ndarray) -> np.ndarray:
    return colors20[np.mod(x, 20)]



@jit("float64[:](float64[:], int32[:], int32[:], float64[:])", nopython=True)
def _scale_to_match_median(data: np.ndarray, indices: np.ndarray,
                           indptr: np.ndarray, genes_total: np.ndarray) -> np.ndarray:
    # Helper function that operates directly on the .data array of a sparse matrix object
    new_data = np.zeros(data.shape)
    # Loop through the columns
    for i in range(genes_total.shape[0]):
        # Retrieve the values
        non_zero_genes_total = genes_total[indices[indptr[i]:indptr[i + 1]]]
        # Find the normalization factor
        w = np.minimum(1, np.median(non_zero_genes_total) / non_zero_genes_total)
        new_data[indptr[i]:indptr[i + 1]] = w * data[indptr[i]:indptr[i + 1]]
    return new_data


@jit(nopython=True)
def numba_random_seed(value: int) -> None:
    """Same as np.random.seed but for numba"""
    np.random.seed(value)



@jit(nopython=True)
def permute_rows_nsign(A: np.ndarray) -> None:
    """Permute in place the entries and randomly switch the sign for each row of a matrix independently.
    """
    plmi = np.array([+1, -1])
    for i in range(A.shape[0]):
        np.random.shuffle(A[i, :])
        A[i, :] = A[i, :] * np.random.choice(plmi, size=A.shape[1])

def scale_to_match_median(sparse_matrix: sparse.csr_matrix, genes_total: np.ndarray) -> sparse.csr_matrix:
    """Normalize contribution of different neighbor genes to match the median totals

    Arguments
    ---------
    sparse_matrix: sparse.csr_matrix
        weights matrix

    genes_total: sparse.csr_matrix shape=(sparse_matrix.shape[0])
        array of the total molecules detected for each gene

    Returns
    -------
    knn_weights: sparse.csr_matrix
        sparse_matrix after the normalization

    # NOTE, since the use I made of this later I could have changed sparse_matrix in place
    """
    newdata = _scale_to_match_median(sparse_matrix.data, sparse_matrix.indices, sparse_matrix.indptr, genes_total)
    return sparse.csc_matrix((newdata,
                              sparse_matrix.indices,
                              sparse_matrix.indptr),
                             shape=sparse_matrix.shape,
                             copy=True)

def gaussian_kernel(X: np.ndarray, mu: float=0, sigma: float=1) -> np.ndarray:
    """Compute gaussian kernel"""
    return np.exp(-(X - mu)**2 / (2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)

print("VeloCD is starting!")
print("Reading in your data......")

#read in the csv file of RNA-Seq data

spliced_reads = np.genfromtxt(sys.argv[3], delimiter=',', dtype=None, encoding=None)

col_gene_name_attribute = spliced_reads[:,[0]]
genes = [] #row attribute, extract the gene names from the file, delete and add back as a row attribute (ra) later
for gene_num in range(1, len(col_gene_name_attribute)):
	gene = str(col_gene_name_attribute[gene_num])
	gene = gene.replace('"', '')
	genes.append(gene)

genes = [s.replace("[", "") for s in genes]
genes = [s.replace("]", "") for s in genes]
genes = [s.replace("'", "") for s in genes]

spliced_reads = np.delete(spliced_reads, 0, 1)  #sample names #1 = col, 0 = row

row_patient_name_attribute = spliced_reads[[0], :]
patients = []
for patient in row_patient_name_attribute: #column attribute, extract the patient/sample IDs from the file, delete and add back as a column attribute (ca) later when the file is made
	for word in patient:
		decoded = str(word)
		decoded_final = decoded.replace('"', "")
		patients.append(decoded_final)

spliced_reads = np.delete(spliced_reads, 0, 0)  #gene names

spl = []
for row in spliced_reads:
	rows = []
	for value in row: #remove the encoding added and ensure all the matrix values are integers
		rows.append(float(str(value)))
	spl.append(rows)

my_shape = np.array(spl).shape
ncols = my_shape[0] #number of genes
nrows = my_shape[1] #number of cells

spl = np.array(spl).reshape(ncols,nrows)

print("Your data contains", ncols, "genes and", nrows, "samples")
#get the clinical complications as the clusters

metadata = np.genfromtxt(sys.argv[5], delimiter=',', dtype=None, encoding=None)
metadatanames = list(metadata[:,0])
metadatanames.pop(0)
metadata = list(metadata[:,1])

metadata.pop(0)
clusters = [] #,'utf-8')
for meta in metadata:
    clusters.append(str(meta))

metadata2 = np.genfromtxt(sys.argv[6], delimiter=',', dtype=None, encoding=None)

metadatanames2 = list(metadata2[:,0])
metadatanames2.pop(0)
metadata2 = list(metadata2[:,1])

metadata2.pop(0)
MarkerData = [] #,'utf-8')
for meta in metadata2:
    MarkerData.append(str(meta))

#check both metadata first column matches first column of spliced 
metadatanames_s = []
for sample in metadatanames: 
	decoded = str(sample)
	decoded_final = decoded.replace('"', "")
	metadatanames_s.append(decoded_final)

metadatanames2_s = []
for patienty in metadatanames2: 
	decoded = str(patienty)
	decoded_final = decoded.replace('"', "")
	metadatanames2_s.append(decoded_final)

if (metadatanames_s == patients):
    print("Congratulations - The sample names of the expression data match that of the first metadata!")
else:
    print("The sample names of the expression data DO NOT match that of the first metadata - please make sure they do and restart the analysis!")
    print("The sample names of the expression data:")
    print(patients)
    print("....do not match the first metadata file sample names:")
    print(metadatanames_s)
    sys.exit()

if sys.argv[29] == "True":
    if (metadatanames2_s == patients):
        print("Congratulations - The sample names of the expression data match that of the second metadata file!")
    else:
        print("The sample names of the expression data DO NOT match that of the second metadata - please make sure they do and restart the analysis!")
        print("The sample names of the expression data:")
        print(patients)
        print("....do not match the second metadata file sample names:")
        print(metadatanames2_s)
        sys.exit()

#clusters = ["cluster1"]*nrows #assign all to the same cluster
#for patient in range(1, nrows+1): #assign each to their own unique cluster:
#    clusters.append("cluster"+str(patient))
#print(loompy_file[:, :]) #returns the matrix

unspliced_reads = np.genfromtxt(sys.argv[4], delimiter=',', dtype=None, encoding=None)
unspliced_reads = np.delete(unspliced_reads, 0, 1)  #1 = col, 0 = row
unspliced_reads = np.delete(unspliced_reads, 0, 0)
unspl = []
for row in unspliced_reads:
	rows = []
	for value in row: #remove the encoding added and ensure all the matrix values are integers
		rows.append(float(str(value)))
	unspl.append(rows)
unspl = np.array(unspl).reshape(ncols,nrows)

#feature_selection_practise(spl, 0.1, "variance", unspl, my_matrix, genes)

filename = "test.loom"
row_attrs = { "Genes": np.array(genes) } #genes are the rows, cells are the columns
col_attrs = { "Patients": np.array(patients), "ClusterName": np.array(clusters), "Markers": np.array(MarkerData) }

my_matrix = np.add(spl, unspl) #sum spliced and unspliced

loompy.create(filename, my_matrix, row_attrs, col_attrs)
#loompy_file = loompy.connect("test.loom")

#Make subdirectories for high-throughput analysis (many sample runs for a single dataset like FEAST)
PCA_Folder = str(sys.argv[2] + sys.argv[8] + "Fate_Maps/PCA/" + str(sys.argv[28]) + "/")
if str(sys.argv[22]) == "True":
    UMAP_Folder = str(sys.argv[2] + sys.argv[8] + "Fate_Maps/UMAP_2D/" + str(sys.argv[28]) + "/")
    if str(sys.argv[23]) == "True":
        UMAP_3D_Folder = str(sys.argv[2] + sys.argv[8] + "Fate_Maps/UMAP_3D/" + str(sys.argv[28]) + "/")
    else:
        UMAP_3D_Folder = str(sys.argv[2] + sys.argv[8] + "Fate_Maps/UMAP_2D/" + str(sys.argv[28]) + "/")

if str(sys.argv[21]) == "True":
    tSNE_Folder = str(sys.argv[2] + sys.argv[8] + "Fate_Maps/tSNE/" + str(sys.argv[28]) + "/")

VelocityValues_Folder = str(sys.argv[2] + sys.argv[8] + "Files/Velocity_Values/" + str(sys.argv[28])+ "/")
Predicted_Unspiced_Folder = str(sys.argv[2] + sys.argv[8] + "Files/Predicted_Unspiced/" + str(sys.argv[28]) + "/")
FutureSpliced_Folder = str(sys.argv[2] + sys.argv[8] + "Files/Future_Spliced/" + str(sys.argv[28]) + "/")

if str(sys.argv[32]) == "yes":
    PCA_Delta_embedding_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Delta_Embedding/" + str(sys.argv[28]) + "/")
    PCA_selfCC_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Self_Correlation_Coefficients/" + str(sys.argv[28]) + "/")
    PCA_UnVecX_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Unitary_Vector_X/" + str(sys.argv[28]) + "/")
    PCA_UnVecY_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Unitary_Vector_Y/" + str(sys.argv[28]) + "/")
    PCA_Transition_Probability_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Transition_Probability/Sample_Level/" + str(sys.argv[28]) + "/")
    PCA_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/")
    if str(sys.argv[33]) == "yes":
        PCA_UnVecZ_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Unitary_Vector_Z/" + str(sys.argv[28])+ "/")
    else:
        PCA_UnVecZ_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Unitary_Vector_Y/" + str(sys.argv[28])+ "/")

if str(sys.argv[23]) == "True":
    UMAP_3D_Delta_embedding_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Delta_Embedding/" + str(sys.argv[28]) + "/")
    UMAP_3D_selfCC_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Self_Correlation_Coefficients/" + str(sys.argv[28]) + "/")
    UMAP_3D_UnVecX_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Unitary_Vector_X/" + str(sys.argv[28]) + "/")
    UMAP_3D_UnVecY_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Unitary_Vector_Y/" + str(sys.argv[28]) + "/")
    UMAP_3D_UnVecZ_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Unitary_Vector_Z/" + str(sys.argv[28])+ "/")
    UMAP_3D_Transition_Probability_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Transition_Probability/Sample_Level/" + str(sys.argv[28]) + "/")
    UMAP_3D_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/")

if str(sys.argv[22]) == "True":
    UMAP_2D_Delta_embedding_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Delta_Embedding/" + str(sys.argv[28]) + "/")
    UMAP_2D_selfCC_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Self_Correlation_Coefficients/" + str(sys.argv[28]) + "/")
    UMAP_2D_UnVecX_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Unitary_Vector_X/" + str(sys.argv[28]) + "/")
    UMAP_2D_UnVecY_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Unitary_Vector_Y/" + str(sys.argv[28]) + "/")
    UMAP_2D_Transition_Probability_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Transition_Probability/Sample_Level/" + str(sys.argv[28]) + "/")
    UMAP_2D_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/")

if str(sys.argv[21]) == "True":
    tSNE_Delta_embedding_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Delta_Embedding/" + str(sys.argv[28]) + "/")
    tSNE_selfCC_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Self_Correlation_Coefficients/" + str(sys.argv[28]) + "/")
    tSNE_UnVecX_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Unitary_Vector_X/" + str(sys.argv[28]) + "/")
    tSNE_UnVecY_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Unitary_Vector_Y/" + str(sys.argv[28]) + "/")
    tSNE_Transition_Probability_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Transition_Probability/Sample_Level/" + str(sys.argv[28]) + "/")
    tSNE_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/")


if str(sys.argv[19]) == "True":
    if str(sys.argv[22]) == "True":
        CI_UMAP_2D_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Confidence_Interval/" 
    if str(sys.argv[23]) == "True":
        CI_UMAP_3D_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Confidence_Interval/"
    if str(sys.argv[21]) == "True":
        CI_tSNE_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Confidence_Interval/" 
    if str(sys.argv[32]) == "yes":
        CI_PCA_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Confidence_Interval/" 

if str(sys.argv[18]) == "True":
    if str(sys.argv[22]) == "True":
        BS_UMAP_2D_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_2D/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Bootstrapping/" 
    if str(sys.argv[23]) == "True":
        BS_UMAP_3D_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/UMAP_3D/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Bootstrapping/"
    if str(sys.argv[21]) == "True":
        BS_tSNE_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/tSNE/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Bootstrapping/" 
    if str(sys.argv[32]) == "yes":
        BS_PCA_GroupLevelTP_Folder = str(sys.argv[2] + sys.argv[8] + "Files/PCA/Transition_Probability/Group_Level/" + str(sys.argv[28]) + "/") + "Bootstrapping/" 

#choose the colour scheme:
import matplotlib as mpl
import matplotlib.cm as cm
#norm = mpl.colors.Normalize(vmin=-20, vmax=10) #convert to greyscake

cmap = matplotlib.colormaps[sys.argv[31]] #cm.Paired #https://matplotlib.org/users/colormaps.html

#from pylab import *

#set1: #e41a1c #377eb8 #4daf4a #984ea3 #ff7f00 #ffff33 #a65628 #f781bf #999999

#Running the pipeline
vlm = VelocytoLoom("test.loom") #still requires a loom file as input, The input bam file needs to be sorted by position: samtools sort mybam.bam -o sorted_bam.bam
location = sys.argv[2] + sys.argv[8]
os.chdir(location)

##vlm.plot_fractions("plot.png")
vlm.set_clusters(vlm.ca["ClusterName"], colormap=cmap)

if sys.argv[10] == "False":
    scaling_to_DO = False
else:
    scaling_to_DO = True

if sys.argv[15] == "False":
    makegifs = False
else:
    makegifs = True

if int(sys.argv[30]) == 3:
    make = [0,1,2]
if int(sys.argv[30]) == 2:
    make = [0,1]

#performs PCA: smoothing the data using kNN neighbours pooling
PCA_FileName1 = str(sys.argv[28])
PCA_FileName2 = str(sys.argv[28]) + "_FeatureSelection"
grouptype = sys.argv[14]
os.chdir(PCA_Folder)
print("Performing PCA...")
vlm.perform_PCA(n_components=int(sys.argv[24])) #must be performed on normalised data, Barnes hut approximation
vlm.plot_pca(dim=make, make_gify=makegifs, Filename=PCA_FileName1, typegroup=grouptype, changeShape=sys.argv[29], Colours=sys.argv[31]) #plots 3D PCA, need to change
if sys.argv[7] == "PCAVar":
    location2 = sys.argv[2]
    os.chdir(location2)
    vlm.PCAVariationFS(n_pca_dims=int(sys.argv[24]), spliced_file=sys.argv[3], resuts_directory=PCA_Folder)
    os.chdir(PCA_Folder)
    vlm.perform_PCA(n_components=int(sys.argv[24])) #must be performed on normalised data, Barnes hut approximation
    vlm.plot_pca(dim=make, make_gify=makegifs, Filename=PCA_FileName2, typegroup=grouptype, changeShape=sys.argv[29], Colours=sys.argv[31])
    print("Performing PCA-based feature selection")

if sys.argv[9] == "CoorelationCoeff":
    os.chdir(location2)
    vlm.CorrelationCoefficentFS(clusters)
    os.chdir(PCA_Folder)
    vlm.perform_PCA(n_components=int(sys.argv[24])) #must be performed on normalised data, Barnes hut approximation
    vlm.plot_pca(dim=make, make_gify=makegifs, Filename="CCFS", typegroup=grouptype, changeShape=sys.argv[29], Colours=sys.argv[31])
    print("Performing correlation coefficient feature selection")

if sys.argv[1] != "NotApplicable":
    number_neighbous = int(sys.argv[1]) #has to be 2 or greater for UMAP
else:
    number_neighbous = sys.argv[1] #has to be 2 or greater for UMAP

#Estimating Velocity:
print("Calculating RNA Velocities...")
vlm.fit_gammas(weighted=True) #apply a gamma fit
my_genes = []
#RNA velocity calculation
vlm.predict_U()
vlm.calculate_velocity(Myneighbours=number_neighbous, VelocityFolder=VelocityValues_Folder, UnsplicedFolder=Predicted_Unspiced_Folder)
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1., Myneighbours=number_neighbous, location=FutureSpliced_Folder) #delta_t=1???
#or under constant unspliced assumption:
#vlm.calculate_shift(assumption="constant_unspliced", delta_t=1.0)
#vlm.extrapolate_cell_at_t(delta_t=1., Myneighbours=number_neighbous, location=FutureSpliced_Folder)

metaloc = sys.argv[5]
metadataGroup2 = sys.argv[6]
R_Transition_Probability_Code = sys.argv[13]

#Embedding under UMAP
#UMAP
if str(sys.argv[22]) == "True":
    print("Performing UMAP...")
    os.chdir(UMAP_Folder)
    vlm.basicUMAP(number_neighbous, make_gify=makegifs, distance=float(sys.argv[11]), metric_method=sys.argv[12], typegroup=grouptype, do3D=str(sys.argv[23]), UMAP_3D_loc=UMAP_3D_Folder, ChangeShape=sys.argv[29], Colours=sys.argv[31])
    #Fit to 2 dimensions:
    for i in range(int(sys.argv[25]), int(sys.argv[26]), int(sys.argv[27])):
        os.chdir(UMAP_Folder)
        vlm.estimate_transition_probUMAP(hidim="Sx_sz", embed="UMAPembed", transform="sqrt", psc=None, n_neighbors=i, knn_random=False, sampled_fraction=0.5)
        vlm.calculate_embedding_shiftUMAP(sigma_corr = 0.05, expression_scaling=scaling_to_DO, Myneighbours=number_neighbous, metadataToUse=metaloc, pathToTransitionProbabilityCode = R_Transition_Probability_Code, grouptype=grouptype, secondgroup=metadataGroup2, manyfeatures=str(sys.argv[16]), grouptype2=str(sys.argv[17]), Randomsampling=str(sys.argv[18]), extraNNvalue=i, DE_location=UMAP_2D_Delta_embedding_Folder, selfCC_folder=UMAP_2D_selfCC_Folder, VecX_location=UMAP_2D_UnVecX_Folder, VecY_loc=UMAP_2D_UnVecY_Folder, TP_folder=UMAP_2D_Transition_Probability_Folder, UMAP_2D_GroupLevelTP=UMAP_2D_GroupLevelTP_Folder, DoSecondPred=str(sys.argv[34]), UserThreshold=float(sys.argv[35]))
        os.chdir(UMAP_Folder)
        vlm.UMAPwithArrows2D(UMAPembedding="UMAPembed", UMAP_Arrows="MY_ARROWS_umap", numberNeighb=number_neighbous, typegroup=grouptype, extraNNvalue=i, ChangeShape=sys.argv[29], Colours=sys.argv[31])
        plt.close("all")
        #Fit to 3 dimensions
        if str(sys.argv[23]) == "True":
            os.chdir(UMAP_3D_Folder)
            vlm.estimate_transition_probUMAP3D(hidim="Sx_sz", embed="UMAPembed3d", transform="sqrt", psc=None, n_neighbors=i, knn_random=False, sampled_fraction=0.5) 
            vlm.calculate_embedding_shiftUMAP3D(sigma_corr = 0.05, expression_scaling=scaling_to_DO, Myneighbours=number_neighbous, metadataToUse=metaloc, pathToTransitionProbabilityCode = R_Transition_Probability_Code, grouptype=grouptype, secondgroup=metadataGroup2, manyfeatures=str(sys.argv[16]), grouptype2=str(sys.argv[17]), Randomsampling=str(sys.argv[18]), extraNNvalue=i, DE_location=UMAP_3D_Delta_embedding_Folder, selfCC_folder=UMAP_3D_selfCC_Folder, VecX_location=UMAP_3D_UnVecX_Folder, VecY_loc=UMAP_3D_UnVecY_Folder, VecZ_loc=UMAP_3D_UnVecZ_Folder, TP_folder=UMAP_3D_Transition_Probability_Folder, UMAP_3D_GroupLevelTP=UMAP_3D_GroupLevelTP_Folder, DoSecondPred=str(sys.argv[34]), UserThreshold=float(sys.argv[35]))
            os.chdir(UMAP_3D_Folder)
            vlm.UMAPwithArrows3D(UMAPembedding="UMAPembed3d", UMAP_Arrows="MY_ARROWS_umap_3d", make_gify=makegifs, numberNeighb=number_neighbous, typegroup=grouptype, extraNNvalue=i, ChangeShape=sys.argv[29], Colours=sys.argv[31])
            plt.close("all")
#tSNE
if str(sys.argv[21]) == "True":
    print("Performing tSNE...")
    os.chdir(tSNE_Folder)
    bh_tsne = TSNE(perplexity=number_neighbous)
    vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :int(sys.argv[24])]) #this is the speed-limiting step
    for i in range(int(sys.argv[25]), int(sys.argv[26]), int(sys.argv[27])):
        vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=None, n_neighbors=i, knn_random=False, sampled_fraction=0.5)
        vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=scaling_to_DO, Myneighbours=number_neighbous, metadataToUse=metaloc, pathToTransitionProbabilityCode = R_Transition_Probability_Code, grouptype=grouptype, secondgroup=metadataGroup2, manyfeatures=str(sys.argv[16]), grouptype2=str(sys.argv[17]), Randomsampling=str(sys.argv[18]), extraNNvalue=i, DE_location=tSNE_Delta_embedding_Folder, selfCC_folder=tSNE_selfCC_Folder, VecX_location=tSNE_UnVecX_Folder, VecY_loc=tSNE_UnVecY_Folder, TP_folder=tSNE_Transition_Probability_Folder, tSNE_GroupLevelTP=tSNE_GroupLevelTP_Folder, DoSecondPred=str(sys.argv[34]), UserThreshold=float(sys.argv[35]))
        os.chdir(tSNE_Folder)
        vlm.one_arrow(numberNeighb=number_neighbous, typegroup=grouptype, extraNNvalue=i, ChangeShape=sys.argv[29], Colours=sys.argv[31])
        plt.close("all")


import os.path

if str(sys.argv[32]) == "yes": #don't want to run this for every embedding NN
    MyPath = str(sys.argv[2] + sys.argv[8] + "Fate_Maps/PCA/"  + str(sys.argv[28]) + "/" + "emat_PCA.csv")
    file_exists = os.path.exists(MyPath)
    if file_exists == False:
        print("Generating PCA-based fate maps...")
        os.chdir(PCA_Folder)
        for i in range(int(sys.argv[25]), int(sys.argv[26]), int(sys.argv[27])):
            vlm.estimate_transition_probPCA(hidim="Sx_sz", embed="pcs", transform="sqrt", psc=None, n_neighbors=i, knn_random=False, sampled_fraction=0.5, AdjustPCs=str(sys.argv[36]))
            vlm.calculate_embedding_shiftPCA(sigma_corr = 0.05, expression_scaling=scaling_to_DO, metadataToUse=metaloc, pathToTransitionProbabilityCode = R_Transition_Probability_Code, grouptype=grouptype, secondgroup=metadataGroup2, manyfeatures=str(sys.argv[16]), grouptype2=str(sys.argv[17]), Randomsampling=str(sys.argv[18]), extraNNvalue=i, DE_location=PCA_Delta_embedding_Folder, selfCC_folder=PCA_selfCC_Folder, VecX_location=PCA_UnVecX_Folder, VecY_loc=PCA_UnVecY_Folder, VecZloc=PCA_UnVecZ_Folder, TP_folder=PCA_Transition_Probability_Folder, PCA_GroupLevelTP=PCA_GroupLevelTP_Folder, ThreeDVersion=str(sys.argv[33]), DoSecondPred=str(sys.argv[34]), UserThreshold=float(sys.argv[35]))
            os.chdir(PCA_Folder)
            vlm.one_arrowPCA(typegroup=grouptype, extraNNvalue=i, ChangeShape=sys.argv[29], Colours=sys.argv[31], ThreeDVersion=sys.argv[33], make_gify=makegifs)
            plt.close("all")

print("VeloCD is finished!")
print("----------------------------------------------------------")


os.chdir(sys.argv[2])
os.remove("test.loom")

"""
mylocytSNE = tSNE_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + ""

for Summaryfile in glob.glob(mylocytSNE):
    vlm.TransitionProbabilityScatterplots(filename=Summaryfile, make_gify=makegifs, Colours=sys.argv[31]) #makes nice presentations of the transition probabilities

mylocyUMAP2D = UMAP_2D_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + ".csv"
for Summaryfile in glob.glob(mylocyUMAP2D):
    vlm.TransitionProbabilityScatterplots(filename=Summaryfile, make_gify=makegifs, Colours=sys.argv[31]) #makes nice presentations of the transition probabilities

mylocyUMAP3D = UMAP_3D_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + ".csv"
for Summaryfile in glob.glob(mylocyUMAP3D):
    vlm.TransitionProbabilityScatterplots(filename=Summaryfile, make_gify=makegifs, Colours=sys.argv[31]) #makes nice presentations of the transition probabilities

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
utils = importr('utils')
path = R_Transition_Probability_Code

def RunTransitionTypeEstimates(filename, location): #estimates the % future, self and past estimates for linear data
    r=ro.r
    r.source(path+"TransitionTypeEstimates.R")
    p=r.TransitionTypeEstimates(filename, location)
    return p

mylocytsNE1 = tSNE_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + "SampleGroup1" + "*" + ".csv"
for Summaryfile in glob.glob(mylocytsNE1):
    if grouptype == "integers": #and integers like time-course
        f=RunTransitionTypeEstimates(Summaryfile, tSNE_GroupLevelTP_Folder)

mylocyUMAP2D1 = UMAP_2D_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + "SampleGroup1" + "*" + ".csv"
for Summaryfile in glob.glob(mylocyUMAP2D1):
    if grouptype == "integers": #and integers like time-course
        f=RunTransitionTypeEstimates(Summaryfile, UMAP_2D_GroupLevelTP_Folder)

mylocyUMAP3D1 = UMAP_3D_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + "SampleGroup1" + "*" + ".csv"
for Summaryfile in glob.glob(mylocyUMAP3D1):
    if grouptype == "integers": #and integers like time-course
        f=RunTransitionTypeEstimates(Summaryfile, UMAP_3D_GroupLevelTP_Folder)

mylocyPCA = PCA_GroupLevelTP_Folder + "PCA_2D*SummaryTable*" + "SampleGroup1" + "*" + ".csv"
for Summaryfile in glob.glob(mylocyPCA):
    if grouptype == "integers": #and integers like time-course
        f=RunTransitionTypeEstimates(Summaryfile, PCA_GroupLevelTP_Folder)

mylocyPCA3D = PCA_GroupLevelTP_Folder + "*PCA_3D*SummaryTable*" + "PCASampleGroup1" + "*" + ".csv"
for Summaryfile in glob.glob(mylocyPCA3D):
    if grouptype == "integers": #and integers like time-course
        f=RunTransitionTypeEstimates(Summaryfile, PCA_Transition_Probability_Folder)

if str(sys.argv[16]) == "True":
    mylocytSNE2 = tSNE_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + "SampleGroup2" + "*" + ".csv"
    for Summaryfile in glob.glob(mylocytSNE2):
        if str(sys.argv[17]) == "integers": #and integers like time-course
            f=RunTransitionTypeEstimates(Summaryfile, tSNE_GroupLevelTP_Folder)
            plt.close("all")
    mylocyUMAP2D2 = UMAP_2D_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + "SampleGroup2" + "*" + ".csv"
    for Summaryfile in glob.glob(mylocyUMAP2D2):
        if str(sys.argv[17]) == "integers": #and integers like time-course
            f=RunTransitionTypeEstimates(Summaryfile, UMAP_2D_GroupLevelTP_Folder)
            plt.close("all")
    mylocyUMAP3D2 = UMAP_3D_GroupLevelTP_Folder + "SummaryTable*" + "NN" + str(number_neighbous) + "*" + "SampleGroup2" + "*" + ".csv"
    for Summaryfile in glob.glob(mylocyUMAP3D2):
        if str(sys.argv[17]) == "integers": #and integers like time-course
            f=RunTransitionTypeEstimates(Summaryfile, UMAP_3D_GroupLevelTP_Folder)
            plt.close("all")

#Run Confidence Interval calculation for each sample to each groups transition probability values across TP NN values!
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
utils = importr('utils')
boot = importr('ggplot2')
boot = importr('Rmisc')
print("Performing Confidence Interval Calculations across the transition probability NN values!")

def ConfidenceIntervalsAcrossTPNN(wd, MainPattern, NN, embedType, outdir):
    r=ro.r
    r.source(path+"ConfidenceIntervalsAcrossTPNN.R")
    l=r.ConfidenceIntervalsAcrossTPNN(wd, MainPattern, NN, embedType, outdir)
    return l

NNval = "_NN"+ str(number_neighbous) + "_"
my_Filelist = ["SummaryTable_Sum_SelfRemoved", "SummaryTable_Mean_SelfNotRemoved", "SummaryTable_Mean_SelfRemoved"]
for Summaryfile in my_Filelist:
    j=ConfidenceIntervalsAcrossTPNN(UMAP_2D_GroupLevelTP_Folder, Summaryfile, NNval, "UMAPSample", CI_UMAP_2D_GroupLevelTP_Folder)
    plt.close("all")
    k=ConfidenceIntervalsAcrossTPNN(UMAP_3D_GroupLevelTP_Folder, Summaryfile, NNval, "UMAP3DSample", CI_UMAP_3D_GroupLevelTP_Folder)
    plt.close("all")
    l=ConfidenceIntervalsAcrossTPNN(tSNE_GroupLevelTP_Folder, Summaryfile, NNval, "tSNESample", CI_tSNE_GroupLevelTP_Folder)
    plt.close("all")

#Confidence Interval via bootstrapping Calculation 
if str(sys.argv[19]) == "True":
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    utils = importr('utils')
    #utils.install_packages('boot')
    #utils.install_packages('tidyverse') #library(foreach); library(doParallel); library(parallel)
    #utils.install_packages('reshape2')
    #utils.install_packages('Rmisc')
    #utils.install_packages('gridExtra')
    boot = importr('boot')
    boot = importr('tidyverse')
    boot = importr('reshape2')
    boot = importr('Rmisc')
    boot = importr('gridExtra')
    def RunInternalBootStrappingTransitonProbabilities(filename, selfremoved, location, metaloc, outputdir, grouptype, numberofIterations):
        r=ro.r
        r.source(path+"InternalBootStrappingTransitionProbabilitiesOriginal.R")
        l=r.InternalBootStrappingTransitonProbabilities(filename, selfremoved, location, metaloc, outputdir, grouptype, numberofIterations)
        return l
    mylocy3 = "Transition_Probability_NN" + str(number_neighbous) + "_*.csv"
    print("Performing Confidence Interval Calculations")
    for Summaryfile in glob.glob(mylocy3):
        j=RunInternalBootStrappingTransitonProbabilities(Summaryfile, "FALSE", BS_UMAP_2D_GroupLevelTP_Folder, metaloc, UMAP_2D_GroupLevelTP_Folder, grouptype, int(sys.argv[20]))
        plt.close("all")
        k=RunInternalBootStrappingTransitonProbabilities(Summaryfile, "TRUE", BS_UMAP_2D_GroupLevelTP_Folder, metaloc, UMAP_2D_GroupLevelTP_Folder, grouptype, int(sys.argv[20]))
        plt.close("all")
        l=RunInternalBootStrappingTransitonProbabilities(Summaryfile, "FALSE", BS_UMAP_3D_GroupLevelTP_Folder, metaloc, UMAP_3D_GroupLevelTP_Folder, grouptype, int(sys.argv[20]))
        plt.close("all")
        m=RunInternalBootStrappingTransitonProbabilities(Summaryfile, "TRUE", BS_UMAP_3D_GroupLevelTP_Folder, metaloc, UMAP_3D_GroupLevelTP_Folder, grouptype, int(sys.argv[20]))
        plt.close("all")
        n=RunInternalBootStrappingTransitonProbabilities(Summaryfile, "FALSE", BS_tSNE_GroupLevelTP_Folder, metaloc, tSNE_GroupLevelTP_Folder, grouptype, int(sys.argv[20]))
        plt.close("all")
        o=RunInternalBootStrappingTransitonProbabilities(Summaryfile, "TRUE", BS_tSNE_GroupLevelTP_Folder, metaloc, tSNE_GroupLevelTP_Folder, grouptype, int(sys.argv[20]))
        plt.close("all")

"""

