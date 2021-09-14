# Fitting functional data splines to raw toe-in gait patterns
# nrokh 2021

#########################################################################
#   input: KJC_x.mat, KJC_y.mat, COP_x.mat, COP_y_mat 
#   output: LOOCV patterns for validation, combined patterns from all 12 subjects
#   dependencies: must have scikit-fda installed
#########################################################################

# INSTALLING SCIKIT-FDA:
# pip install scikit-fda --user

# SCIKIT-FDA DOCUMENTATION: 
# https://fda.readthedocs.io/en/latest/auto_tutorial/plot_getting_data.html

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import skfda
from itertools import chain

# LOAD DATA
KJCx = sio.loadmat('KJC_x.mat')
KJCy = sio.loadmat('KJC_y.mat')
COPx = sio.loadmat('COP_x.mat')
COPy = sio.loadmat('COP_y.mat')

################################### FDA LOOCV: ############################################
for i in range(12):
    # COPy
    COPy_arr = COPy['all_COP_t_y']
    COPy_arr_r = np.copy(COPy_arr)
    COPy_arr_r = np.delete(COPy_arr_r, (i), axis=0)
    COPy_arr_r = np.reshape(COPy_arr_r, (110))
    COPy_arr_r = np.array(list(chain(*COPy_arr_r)))
    
    # compute basis for all FPA
    COPy_fpa = np.empty((10,100))
    for j in range(10):
        COPy_fpa[j,:] = np.mean(COPy_arr_r[COPy_arr_r[:,0]==-j,1:],0)
          
    fd = skfda.FDataGrid(
        data_matrix = COPy_fpa)

    basis = skfda.representation.basis.BSpline(n_basis = 12)
    X_basis_spline = fd.to_basis(basis)
    # save 
    splines_COPy = X_basis_spline.to_grid().data_matrix[:,1::5,0]
    
    
    # COPx
    COPx_arr = COPx['all_COP_t_x']
    COPx_arr_r = np.copy(COPx_arr)
    COPx_arr_r = np.delete(COPx_arr_r, (i), axis=0)
    COPx_arr_r = np.reshape(COPx_arr_r, (110))
    COPx_arr_r = np.array(list(chain(*COPx_arr_r)))
    
    # compute basis for all FPA
    COPx_fpa = np.empty((10,100))
    for j in range(10):
        COPx_fpa[j,:] = np.mean(COPx_arr_r[COPx_arr_r[:,0]==-j,1:],0)
          
    fd = skfda.FDataGrid(
        data_matrix = COPx_fpa)

    basis = skfda.representation.basis.BSpline(n_basis = 12)
    X_basis_spline = fd.to_basis(basis)
    # save 
    splines_COPx = X_basis_spline.to_grid().data_matrix[:,1::5,0]
    
    
    # KJCx
    KJCx_arr = KJCx['all_KJC_t_x']
    KJCx_arr_r = np.copy(KJCx_arr)
    KJCx_arr_r = np.delete(KJCx_arr_r, (i), axis=0)
    KJCx_arr_r = np.reshape(KJCx_arr_r, (110))
    KJCx_arr_r = np.array(list(chain(*KJCx_arr_r)))
    
    # compute basis for all FPA
    KJCx_fpa = np.empty((10,100))
    for j in range(10):
        KJCx_fpa[j,:] = np.mean(KJCx_arr_r[KJCx_arr_r[:,0]==-j,1:],0)
          
    fd = skfda.FDataGrid(
        data_matrix = KJCx_fpa)

    basis = skfda.representation.basis.BSpline(n_basis = 12)
    X_basis_spline = fd.to_basis(basis)
    # save 
    splines_KJCx = X_basis_spline.to_grid().data_matrix[:,1::5,0]
    
    
    # KJCY
    KJCy_arr = KJCy['all_KJC_t_y']
    KJCy_arr_r = np.copy(KJCy_arr)
    KJCy_arr_r = np.delete(KJCy_arr_r, (i), axis=0)
    KJCy_arr_r = np.reshape(KJCy_arr_r, (110))
    KJCy_arr_r = np.array(list(chain(*KJCy_arr_r)))
    
    # compute basis for all FPA
    KJCy_fpa = np.empty((10,100))
    for j in range(10):
        KJCy_fpa[j,:] = np.mean(KJCy_arr_r[KJCy_arr_r[:,0]==-j,1:],0)
          
    fd = skfda.FDataGrid(
        data_matrix = KJCy_fpa)

    basis = skfda.representation.basis.BSpline(n_basis = 12)
    X_basis_spline = fd.to_basis(basis)
    # save 
    splines_KJCy = X_basis_spline.to_grid().data_matrix[:,1::5,0]
    

    # save that subject's LOOCV
    sio.savemat("COPy_FDA" + str(i+1) + ".mat", {'COPy_offset': splines_COPy})
    sio.savemat("COPx_FDA" + str(i+1) + ".mat", {'COPx_offset': splines_COPx})
    sio.savemat("KJCy_FDA" + str(i+1) + ".mat", {'KJCy_offset': splines_KJCy})
    sio.savemat("KJCx_FDA" + str(i+1) + ".mat", {'KJCx_offset': splines_KJCx})


################################### FDA COMBINED ############################################
# KJC_y
fd = skfda.FDataGrid(
        data_matrix = KJCy_fpa)

basis = skfda.representation.basis.BSpline(n_basis = 5)
X_basis_spline = fd.to_basis(basis)
# save 
splines_KJCy = X_basis_spline.to_grid().data_matrix[:,1::5,0]

# KJC_x
fd = skfda.FDataGrid(
        data_matrix = KJCx_fpa)

basis = skfda.representation.basis.BSpline(n_basis = 5)
X_basis_spline = fd.to_basis(basis)
# save 
splines_KJCx = X_basis_spline.to_grid().data_matrix[:,1::5,0]

# COP_y
fd = skfda.FDataGrid(
        data_matrix = COPy_fpa)

basis = skfda.representation.basis.BSpline(n_basis = 5)
X_basis_spline = fd.to_basis(basis)
# save 
splines_COPy = X_basis_spline.to_grid().data_matrix[:,1::5,0]

# COP_y
fd = skfda.FDataGrid(
        data_matrix = COPx_fpa)

basis = skfda.representation.basis.BSpline(n_basis = 5)
X_basis_spline = fd.to_basis(basis)
# save 
splines_COPx = X_basis_spline.to_grid().data_matrix[:,1::5,0]

sio.savemat("offsets/y_FDA_COP.mat", {'COPy_offset': splines_COPy})
sio.savemat('offsets/x_FDA_COP.mat', {'COPx_offset': splines_COPx})
sio.savemat('offsets/y_FDA_KJC.mat', {'KJCy_offset': splines_KJCy})
sio.savemat('offsets/x_FDA_KJC.mat', {'KJCx_offset': splines_KJCx})