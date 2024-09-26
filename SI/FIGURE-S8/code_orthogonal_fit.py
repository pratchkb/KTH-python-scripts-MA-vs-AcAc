import numpy as np
import copy



# function definitions
def normalize(vec):

    evec = copy.copy(vec)
    lvec = np.linalg.norm(evec)
    evec = evec/lvec

    return evec
 
def orthogonal_plane_fit(X):
    """Takes an array of points (3D or higher)
    and computes the least-squares orthogonal line fit

    Parametric form: (X - X0)^T * n = 0

    where X0 is a point in the plane and n the unit normal vector orthogonal to the plane

    Args:
        X - (np.ndarray): (npoints, dim)-array with coordinates for all points

    Returns:
        (np.ndarray, np.ndarray): point on the fit plane, normal vector of the fit plane
    """

    dim = np.shape(X)[1]
    if dim < 3:
        raise ValueError("The dimension is less than 3.")

    npoints = np.shape(X)[0]
    X0 = np.zeros(dim)

    for idim in range(dim):
        X0[idim] = np.sum(X[:,idim])/npoints

    X_X0 = X - np.outer(np.ones(npoints), X0)

    M = np.zeros((dim, dim))
    H = np.zeros((dim, dim))

    #print(M.shape)

    for idim in range(dim):
        for jdim in range(idim, dim):
            #print(idim, jdim)
            M[idim, jdim] = np.sum(X_X0[:,idim]*X_X0[:,jdim])
            if idim != jdim:
                M[jdim, idim] = M[idim, jdim]

    u, S, vh = np.linalg.svd(M)

    idx_min = np.argmin(S)
    n = vh[idx_min, :]

    return X0, n

# input x,y,z coordinates of N points for fitting to plane
def coord(xyz,list_3d_points):
    """
    Takes the xyz coordinates and the list of 3d point from which to fit plane
    """
    xs = np.zeros(shape=(len(list_3d_points), 3))
      
    for i in range(len(list_3d_points)):
        if len(list_3d_points) >= 3:
            b_i = list_3d_points[i]
            xs[i] = xyz[b_i]

        else:
            #xs[i] = [0,0,0]
            raise ValueError("Need more than 3 points in 3D space")
    return xs

#def normal_direction_fix():

            

def compute_torsion_orthogonal_plane(
    xyz,
    triple,
    list_3d_points
    ):
    """ Dihedral angle routine with respect to a fitted plane which minimizes orthogonal distances 
    The angle is defined by two planes from:
        [a1, a2, b1]  and
        [b1, b2, b3, b4, b5, b6, ...] to be fitted
                                                                   
    Arguments:                                                     
        xyz - list/ndarray with coordinates or full xyz list       
        triple - list(int) - dices : [a1, a2, b1]                        
        list_3d_points - list(int) -  [b1, b2, b3, b4, b5, b6, ...] at least 3 points                         

                            
        \                b3 
         \              /   
          \            /    
        a1 \__________/ b1  
           /   v2->   \  (plane 2)
     v1   / (plane 1)  \   Out
         /    Out       \  
       a2                b2

    Return:                                                        
        float - torsion angle in degrees                           
    """

    a1 = triple[0]
    a2 = triple[1]
    b1 = triple[2]

    b = list_3d_points[0]
    b2 = list_3d_points[1]
    b3 = list_3d_points[2]

    # vectors to define 1st plane by a2, a1, b1
    v1 = xyz[a2,:] - xyz[a1,:]
    v2 = xyz[b1,:] - xyz[a1,:]

    # vectors to define fix direction by b1(b) here, b2, b3
    vec1 = xyz[b2,:] - xyz[b,:]
    vec2 = xyz[b3,:] - xyz[b,:]

    # input x,y,z coordinates of N points for fitting to plane
    xs = np.zeros(shape=(len(list_3d_points), 3))
    #print(xs)
    for i in range(len(list_3d_points)):
        if len(list_3d_points) >= 3:
            b_i = list_3d_points[i]
            xs[i] = xyz[b_i]
            #print(xyz[b_i])
            
        else:
            print("Error: Provide minimum 3 points for plane fitting")

    
    fix_direc = normalize(np.cross(vec1, vec2))
    #print(fix_direc)

    normal_to_ortho_fit_plane = np.squeeze(np.asarray(orthogonal_plane_fit(xs)[1]))
    #print(normal_to_ortho_fit_plane)

    angle = np.arccos(np.dot(fix_direc,normal_to_ortho_fit_plane)) * 180. / np.pi
    #print(angle)

    if angle > 90:
        normal_to_ortho_fit_plane = normal_to_ortho_fit_plane
    else:
        normal_to_ortho_fit_plane = (-1)*normal_to_ortho_fit_plane
    #print(normal_to_ortho_fit_plane)

    n1 = normalize(np.cross(v1, v2))
    #print(n1)
    n2 = normalize(normal_to_ortho_fit_plane)
    #print(n2)
    n1dotn2 = np.dot(n1, n2)
    #print(n1dotn2)
    theta = np.arccos(n1dotn2) * 180. / np.pi

    # Adjust sign to get a +/-180 range
    sign = np.sign(np.dot(np.cross(n1, n2), v2))
    theta *= sign

    return theta, n1, normal_to_ortho_fit_plane




if __name__ == "__main__":
  
    import get_optim
    #from aimsrelated.tcutil.code.geom_param import geom_param
    #from aimsrelated.tcutil.code.get_optim.get_optim import read_xyz_as_np
    data = get_optim.read_xyz_as_np('example.xyz')

    print(compute_torsion_orthogonal_plane(data, [5,6,0], [1,0,3,4]))
    print(compute_torsion_orthogonal_plane(data, [5,7,0], [1,0,3,4]))
    print(compute_torsion_orthogonal_plane(data, [5,8,0], [1,0,3,4]))

    print(compute_torsion_orthogonal_plane(data, [9,3,10], [1,0,3,4]))
    print(compute_torsion_orthogonal_plane(data, [9,3,11], [1,0,3,4]))
    print(compute_torsion_orthogonal_plane(data, [9,3,12], [1,0,3,4]))


    #print(compute_torsion_plane(data, [6,5,0], [1,0,3,4]))
    #print(compute_torsion_plane(data, [7,5,0], [1,0,3,4]))
    #print(compute_torsion_plane(data, [8,5,0], [1,0,3,4]))

    #print(compute_torsion_plane(data, [10,9,3], [1,0,3,4]))
    #print(compute_torsion_plane(data, [11,9,3], [1,0,3,4]))
    #print(compute_torsion_plane(data, [12,9,3], [1,0,3,4]))


#    print(compute_torsion_orthogonal_plane(data, [6,5,0], [1,0,3,4]))
#    print(compute_torsion_orthogonal_plane(data, [7,5,0], [1,0,3,4]))
#    print(compute_torsion_orthogonal_plane(data, [8,5,0], [1,0,3,4]))

#    print(compute_torsion_orthogonal_plane(data, [10,9,3], [1,0,3,4]))
#    print(compute_torsion_orthogonal_plane(data, [11,9,3], [1,0,3,4]))
#    print(compute_torsion_orthogonal_plane(data, [12,9,3], [1,0,3,4]))

#    ndotn = np.dot(compute_torsion_plane(data, [6,5,0], [1,0,3,4])[1], compute_torsion_orthogonal_plane(data, [6,5,0], [1,0,3,4])[1])
#    ntheta = np.arccos(ndotn) * 180. / np.pi

#    print(ntheta)
