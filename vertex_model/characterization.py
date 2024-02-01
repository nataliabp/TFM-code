"""This class aims to characterize cell networks at a time step. To characterize a network we propose 
    the study of four parameters:
        i) Edges length d_ij 
        ii) Triangle (Polygon) area A_i
        iii) Length distortion gamma_k
        iv) Angle distortion phi_k

    The two last metrics are computed to assess network regularity. gamma_k measeures the length variance 
    of the kth triangle when compared to the ideal equilateral triangle (regular polygon when correspond)
    for which gamma_k is supposed to be 0. phi_k measures the variance of the angles in comparison with 
    an equilateral shape. 
    They are calculated as follows: 
    - gamma_k  = sum(d_ki-mu_k)^2/(3mu_k^2)
    - phi_ k = (3/pi ^ 2) sum(thetha_ki -pi/3)^2 

    Both sums are over i from i = 1 to i = 3 bc it is assumed that the polygons are triangles in 
    this case. If the polygon is not a triangle the summatory must be adapted. The index k goes from 
    1 to the total number of triangles (polygons) in the mesh. 
    
    The value mu_k is the mean length of the triangle edges, i.e, 1/3*sum(d_ki) and theta_ki is the 
    angle bewteen the edges k and i of the triangle. It can be computed using the well known formula
    of the cosine. 

    Reference: Scientific Reports, 
    Quantification of topological features in cell meshes to explore E-cadherin dysfunction
"""

## IMPORTS
import numpy as np
from scipy.spatial import Delaunay

def mean_area(self):
    # return the mean of all cells areas 
    area_vec = self.mesh.area()
    area = sum(area_vec) /  area_vec.size
    return area
def one_cell_area(self, cell):
    return self.mesh.area[cell]

def triangle_area(t): 
    '''
        Inputs: 
        triangle: array of 3 arrays, the vertices of the triangle,
        each one containing the x and the y coordinate of each vertex
    Output: area of the triangle 
    '''
    a = t[0][0]*(t[1][1]-t[2][1])
    b = t[1][0]*(t[2][1]-t[0][1])
    c = t[2][0]*(t[0][1]-t[1][1])
    area_t =0.5*(a+b+c)
    return area_t
    

## other functions to characterize a mesh 
def mean_area_triangles( triangles):
        """Function to calculate the mean of the areas of the Delaunay triangulation of a set of points"""
        area_vec = []
        for t in triangles: 
            a = t[0][0]*(t[1][1]-t[2][1])
            b = t[1][0]*(t[2][1]-t[0][1])
            c = t[2][0]*(t[0][1]-t[1][1])
            area_t =0.5*(a+b+c)
            area_vec.append(area_t)
        area = np.mean(area_vec)
        std = np.std(area_vec)
        return area, std

def edges_length_triangles(triangles): 
    """
    Computes the length of the edges of each triangle
    edges_length is an array of arrays of 3 elements, each one containing the three edges of a triangle
    eij are the vectors (for posterior calculations of angles, etc)
    """
    edges_length = []
    eij= []
    for t in triangles: #recall that t is an array of 3 elements, each of them with the x and y coordinate of each vertex 
        e01 = np.linalg.norm(t[0]-t[1])
        e02 = np.linalg.norm(t[0]-t[2])
        e12 = np.linalg.norm(t[1]-t[2])
        edges_length.append([e01, e02, e12])
        eij.append([t[0]-t[1],t[0]-t[2],t[1]-t[2]] )
    return edges_length, eij

def mean_perimeter(triangles):
    edges_length, eij = edges_length_triangles(triangles)
    perimeter = [sum(i) for i in edges_length]
    p_mean = np.mean(perimeter)
    return p_mean

def total_length(triangles):
    edges_length, eij = edges_length_triangles(triangles)
    perimeter = [sum(i) for i in edges_length]
    total_p =sum(perimeter)
    return total_p


def mean_length_single_t(triangles): 
    """
    Computes the mean of the length of the edges of each triangle in the Delaunay triangulation of the given points
    """
    edges_length, eij = edges_length_triangles(triangles)
    return np.mean(edges_length, 1)

def angles_t(triangles): 
    edges_length, eij = edges_length_triangles(triangles)
    angles = []
    for t in eij:
        theta01 = angle_between(t[0], t[1])
        theta02 = angle_between(t[0], -t[2])
        theta21 = angle_between(t[2], t[1])
        angles.append([theta01, theta02, theta21])
    return angles
    
def angle_between(v1, v2): 
    num = np.dot(v1, v2)
    den = np.linalg.norm(v1)*np.linalg.norm(v2)
    angle = np.arccos(num/den)
    return angle

def length_distortion(triangles): 
    t, eij = edges_length_triangles(triangles)
    mu = mean_length_single_t(triangles)
    l_distortion = []
    for i in np.arange(len(t)): 
        m = mu[i]
        n =  np.sum((t[i]-m)**2)
        d = 3*(m**2)
        gamma_t = n/d
        l_distortion.append(gamma_t)
    return l_distortion

def angle_distortion(triangles): 
    angles = angles_t(triangles)
    angle_distortion = []
    r = np.pi / 3
    for  i in np.arange(len(angles)): 
        n = 3* (((angles[i][0]-r)**2) +((angles[i][1]-r)**2 )+((angles[i][2]-r)**2))
        d = np.pi ** 2
        phi_k = n/d
        #dif = (t[0]-r)**2 +(t[1]-r)**2 +(t[2]-r)**2 
        angle_distortion.append(phi_k)
    return angle_distortion

def triangle_sieve(triangles, min, max):
    '''
    Remove triangles of the triangulation such that they are highly obtuse triangles,
    containing very large and small angles that are generated by the Delaunay 
    triangulation. They represent outliers so they have to be removed and not considered
    for analysis. 
    For the moment the triangles that will be removed are those such that one of its angles
    alpha is either >9/10pi or > pi/10
    '''
    obtuse =[]
    non_obtuse = []
    obtuse_indices =[]
    non_obtuse_indices = []
    angles =  angles_t(triangles) #each element of this array contains the angles of a triangle 
    for  i in np.arange(len(angles)): 
            if ((angles[i][0] <min) | (angles[i][0] > max)
            |(angles[i][1] < min) | (angles[i][1] >max)|
            (angles[i][2] < min) | (angles[i][2] > max)):
                obtuse_indices.append(i)
                obtuse.append(triangles[i])
            else:
                non_obtuse_indices.append(i)
                non_obtuse.append(triangles[i])
    return non_obtuse, non_obtuse_indices, obtuse, obtuse_indices
def triangle_sieve_zscore(triangles,threshold):
    '''
    Remove triangles of the triangulation such that they are highly obtuse triangles,
    containing very large and small angles that are generated by the Delaunay 
    triangulation. They represent outliers so they have to be removed and not considered
    for analysis. 
    For the moment the triangles that will be removed are those such that one of its angles
    alpha is either >9/10pi or > pi/10
    '''
    obtuse =[]
    non_obtuse = []
    obtuse_indices =[]
    non_obtuse_indices = []
    angles =  angles_t(triangles) #each element of this array contains the angles of a triangle 
    for  i in np.arange(len(angles)): 
            if ( (angles[i][0] > threshold)
            | (angles[i][1] >threshold)|
             (angles[i][2] > threshold)):
                obtuse_indices.append(i)
                obtuse.append(triangles[i])
            else:
                non_obtuse_indices.append(i)
                non_obtuse.append(triangles[i])
    return non_obtuse, non_obtuse_indices, obtuse, obtuse_indices
def cell_vertices(mesh, id):    
    '''
    Compute the (x,y) coordinates of the vertices of the given cell
        mesh: mesh object where the cell belongs to 
        id: id of the cell whose vertices we want to compute 
    '''
    cell_i = []
    for i in range(len(mesh.face_id_by_edge)):
        if mesh.face_id_by_edge[i] == id:
            cell_i.append(i)
    return cell_i
def centroid(vertexes_x, vertexes_y):   
     '''
     Computes the centroid of a polygon of the given coordinates
     '''
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
    for i in range(len(mesh.n_face)):       
       cell_i = cell_vertices(mesh, i) #vertices of the ith cell 
       cell_center_i = centroid(mesh.vertices[0][cell_i], mesh.vertices[1][cell_i])
       centres_x.append(cell_center_i[0])
       centres_y.append(cell_center_i[1])
    return centres_x, centres_y

def mean_area_hex(mesh):
    areas_non_zero = []
    for a in mesh.area:
        if a != 0:
            areas_non_zero.append[a]
    return np.mean(areas_non_zero)

def mean_perimeter_hex(mesh): 
    per_non_zero = []
    for p in mesh.perimeter:
        if p != 0:
            per_non_zero.append[p]
    return np.mean(per_non_zero)



def shape_index(mesh): 
    """Returns the mean of all cells shape index"""
    shape_index =[]
    for j in range(len(mesh.area)):
        l = mesh.perimeter[j]
        a = mesh.area[j]
        if (a != 0) & (l !=0):
            shape_index.append(l/np.sqrt(a))
    return np.mean(shape_index)




    
    






    


