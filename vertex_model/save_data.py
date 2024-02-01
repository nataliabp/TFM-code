
from scipy.spatial import Delaunay
import vertex_model.characterization as crt
import numpy as np
import os
# to save data from simulation
def save_simulation(I, data, mutant_cells, id_ecad_cells, date, size, noise): # data is going to be the simulation (only one) 
    """
    I : number of the simulation
    data: simulation data
    name_folder: name of the folder where data is kept. It has to include the percentage of mutant
                cells of the current simulation 
    """
    # main folder 
    path=r'C:\Users\natal\OneDrive\Documentos\MASTER\TFM\Simulations' 
    name_folder = 'sim_'+'%0.3f'%mutant_cells +f'_{I}'+ f'_{size}'+f'_{noise}_' + date

    new_path = os.path.join(path, name_folder)
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    # last meshes folder 
    path_final_meshes = os.path.join(new_path, 'final_meshes')
    if not os.path.exists(path_final_meshes):
        os.makedirs(path_final_meshes)
        
    save_last_mesh(data, path_final_meshes, id_ecad_cells)    

## Characterization of the final mesh of each simulation 
def save_last_mesh(simulation, path, id_ecad_cells):
# simulation will contain all time steps meshes. We only need the last one 
    area = []
    length_distortion = []
    angle_distortion = []
    edges_length = []
    si =[]
    te = []
    last_mesh = simulation[-1].mesh
    centres_x, centres_y, centres = mesh_centres(last_mesh)
    tri = Delaunay(centres)
    triangles = centres[tri.simplices]
    non_obtuse, non_obtuse_indices, obtuse, obtuse_indices = crt.triangle_sieve(triangles, np.pi / 10, np.pi - (np.pi / 10))
    trisieve = tri.simplices[non_obtuse_indices]
    length_distortion.append(np.mean(crt.length_distortion(non_obtuse)))
    angle_distortion.append(np.mean(crt.angle_distortion(non_obtuse)))
    edges_length.append(crt.mean_perimeter(non_obtuse))
    ar, std = crt.mean_area_triangles(non_obtuse)
    area.append(ar)
    shape_index = crt.shape_index(last_mesh)
    si.append(shape_index)
    time_extrusion = t_extrusion(simulation, id_ecad_cells)
    te.append(time_extrusion)
    file_ld = os.path.join(path,  'length_distortion.txt')
    np.savetxt(file_ld, length_distortion, delimiter=',')

    file_ad = os.path.join(path, 'angle_distortion.txt')
    np.savetxt(file_ad, angle_distortion, delimiter=',')

    file_el = os.path.join(path, 'edges_length.txt')
    np.savetxt(file_el, edges_length, delimiter=',')

    file_a = os.path.join(path, 'area.txt')
    np.savetxt(file_a, area, delimiter=',')

    file_si = os.path.join(path, 'shape_index.txt')
    np.savetxt(file_si, si, delimiter=',')

    file_te = os.path.join(path, 'time_extrusion.txt')
    np.savetxt(file_te, te, delimiter=',')

def cell_vertices(mesh, id):    
    cell_i = []
    for i in range(len(mesh.face_id_by_edge)):
        if mesh.face_id_by_edge[i] == id:
            cell_i.append(i)
    return cell_i
def centroid(vertexes_x, vertexes_y):   
     _len = len(vertexes_x)
     _x = sum(vertexes_x) / _len
     _y = sum(vertexes_y) / _len
     return(_x, _y)
def mesh_centres(mesh):
    '''
    Compute the centres of the cells of a given mesh 
    '''
    centres_x =[]
    centres_y = [] 
    centres =np.zeros((mesh.n_face, 2))
    for i in range(mesh.n_face):       
       cell_i = cell_vertices(mesh, i) #vertices of the ith cell 
       if len(mesh.vertices[0][cell_i]) != 0:
        cell_center_i = centroid(mesh.vertices[0][cell_i], mesh.vertices[1][cell_i])
        centres_x.append(cell_center_i[0])
        centres_y.append(cell_center_i[1])
        #centres.append(np.array([cell_center_i[0],cell_center_i[1] ]))
        centres[i, 0] = cell_center_i[0]
        centres[i, 1] = cell_center_i[1]
    return centres_x, centres_y, centres


## Extrusion Time 
def t_extrusion(simulation, id_ecad_cells): 
    t = 0
    for time_step in simulation:
        extrused_cells = 0
        t = t+1
        area = time_step.mesh.area
        for a in area:
            if a == 0:
                extrused_cells= extrused_cells+1
        if extrused_cells == len(id_ecad_cells): 
            break
    return t

