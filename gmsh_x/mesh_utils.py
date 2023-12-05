import gmsh
import sys
import numpy as np

def copyAndRotate(geometry, startAngle, N):
    tags = [geometry[0]]
    jumpAngle = 2*np.pi/N
    angle = startAngle
    for i in range(1, N):
        copiedGeometry = gmsh.model.occ.copy(geometry)
        tags.append(copiedGeometry[0])
        angle += jumpAngle
        gmsh.model.occ.rotate(copiedGeometry, 0, 0, 0, 0, 0, 1, angle)
    gmsh.model.occ.synchronize()
    return tags

def mirror_mesh_x_axis(reflection_axis = [1,-1,1]):
    # get the mesh data
    m = {}
    for e in gmsh.model.getEntities():
        # print(e)
        bnd = gmsh.model.getBoundary([e])
        nod = gmsh.model.mesh.getNodes(e[0], e[1])
        ele = gmsh.model.mesh.getElements(e[0], e[1])
        m[e] = (bnd, nod, ele)

    # transform the mesh and create new discrete entities to store it
    def transform(m, offset_entity, offset_node, offset_element, tx, ty, tz):
        for e in sorted(m):
            gmsh.model.addDiscreteEntity(e[0], e[1] + offset_entity,
                                        [abs(b[1]) + offset_entity for b in m[e][0]])

            coord = []
            for i in range(0, len(m[e][1][1]), 3):
                x = m[e][1][1][i] * tx
                y = m[e][1][1][i + 1] * ty
                z = m[e][1][1][i + 2] * tz
                coord.append(x)
                coord.append(y)
                coord.append(z)
            gmsh.model.mesh.addNodes(e[0], e[1] + offset_entity,
                                        m[e][1][0] + offset_node, coord)
            gmsh.model.mesh.addElements(e[0], e[1] + offset_entity, m[e][2][0],
                                        [t + offset_element for t in m[e][2][1]],
                                        [n + offset_node for n in m[e][2][2]])
            if (tx * ty * tz) < 0: # reverse the orientation
                gmsh.model.mesh.reverse([(e[0], e[1] + offset_entity)])

    mx,my,mz = reflection_axis[0],reflection_axis[1], reflection_axis[2]
    max_node_tag = gmsh.model.mesh.getMaxNodeTag() 
    max_element_tag = gmsh.model.mesh.getMaxElementTag() 
    transform(m, 1000, 100000, 100, mx, my, mz)
    
    gmsh.model.mesh.removeDuplicateNodes()

def yaw(angle_deg):
    angle = np.deg2rad(angle_deg)

    yaw = np.array([[np.cos(angle), -np.sin(angle), 0],
              [np.sin(angle),  np.cos(angle), 0],
              [0,0,1]])
    return yaw

def getMeshData(gmsh_model):
    # get the mesh data
    m = {}
    for e in gmsh_model.getEntities():
        # print(e)
        bnd = gmsh_model.getBoundary([e])
        nod = gmsh_model.mesh.getNodes(e[0], e[1])
        ele = gmsh_model.mesh.getElements(e[0], e[1])
        m[e] = (bnd, nod, ele)
    
    return m

# rotate the mesh and create new discrete entities to store it
def rotate_z(m, offset_entity, offset_node, offset_element, yaw_angle):
    for e in sorted(m):
        gmsh.model.addDiscreteEntity(e[0], e[1] + offset_entity,
                                    [abs(b[1]) + offset_entity for b in m[e][0]])

        coord = []
        for i in range(0, len(m[e][1][1]), 3):

            point_old = m[e][1][1][i:i+3]
            point_new = point_old @ yaw(yaw_angle)
            x = point_new[0]
            y = point_new[1]
            z = point_new[2]
            coord.append(x)
            coord.append(y)
            coord.append(z)
        gmsh.model.mesh.addNodes(e[0], e[1] + offset_entity,
                                    m[e][1][0] + offset_node, coord)
        gmsh.model.mesh.addElements(e[0], e[1] + offset_entity, m[e][2][0],
                                    [t + offset_element for t in m[e][2][1]],
                                    [n + offset_node for n in m[e][2][2]])
        gmsh.model.mesh.reverse([(e[0], e[1] + offset_entity)])

def cart2cyl(x, y, z):
    # cylindrical to Cartesian
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    zeta = z
    return rho, phi, zeta 

def cyl2cart(rho, phi, zeta):
    # cylindrical to Cartesian
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    z = zeta
    return x, y, z

def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 2)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 2)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 0)

    gmsh.option.setNumber("Mesh.Lines", 0)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 2)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 2) # CHANGE THIS FLAG TO 0 TO SEE LABELS

    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    gmsh.option.setNumber("Mesh.VolumeFaces", 0)