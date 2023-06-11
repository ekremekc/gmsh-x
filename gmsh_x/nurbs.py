import gmsh
import numpy as np

def generate_radius_dict(radius, N_u, N_v):
        radius_dict = {}
        for key_u in range(0,N_u-1):
            radius_dict[key_u] = {}                
            for key_v in range(0,N_v):
                radius_dict[key_u][key_v] = radius
        
        # Write the same radii values for the last dict
        # This is extremely important to ensure the continuity on the BSPline surface 
        radius_dict[N_u-1] = radius_dict[0]

        return radius_dict

def calculate_knot_v(knots_v, multiplicities_v, deg=2):
    
    multiplicities_v[0] = deg+1
    multiplicities_v[-1] = deg+1
    knots_array_v = []
    for i in range(len(knots_v)):
        knots_array_v += [knots_v[i]] * int(multiplicities_v[i])

    return knots_array_v

def nurbs_cylinder(z_start, z_end, radius, N_cpt_v, lc, center=[0,0,0], deg=2):

    N_u = 9
    
    knots_u= [0, 1/4, 1/2, 3/4, 1]
    multiplicities_u = [3, 2, 2, 2, 3]
    knots_array_u =[0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]

    knots_v = np.linspace(0,1,N_cpt_v-1)
    multiplicities_v = np.ones(len(knots_v))
    multiplicities_v[0] = deg+1
    multiplicities_v[-1] = deg+1
    knots_array_v = []
    for i in range(len(knots_v)):
        knots_array_v += [knots_v[i]] * int(multiplicities_v[i])

    weights_surface = [1, 2**0.5/2, 1, 2**0.5/2, 1, 2**0.5/2, 1, 2**0.5/2, 1]
    weights = np.tile(weights_surface,N_cpt_v)
    z_wall = np.linspace(z_start,z_end,N_cpt_v)

    vertex_list = []

    for r_index, z in enumerate(z_wall):

        vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[0][r_index], center[1], z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[1][r_index],center[1]+radius[1][r_index],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0],center[1]+radius[2][r_index],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0]-radius[3][r_index],center[1]+radius[3][r_index],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0]-radius[4][r_index],center[1],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0]-radius[5][r_index],center[1]-radius[5][r_index],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0],center[1]-radius[6][r_index],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[7][r_index],center[1]-radius[7][r_index],z,lc))
        vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[8][r_index], center[1], z,lc))

    # print(vertex_list)
    # for i in range(len(z_wall)):
    #     print(i*(N_u-1)+1)
    outer_surface = gmsh.model.occ.addBSplineSurface(vertex_list, N_u, degreeU=2, degreeV=2,
                                    weights=weights, knotsU=knots_u, multiplicitiesU=multiplicities_u,
                                    knotsV=knots_v, multiplicitiesV=multiplicities_v)

    gmsh.model.occ.synchronize()
    curve_entities = gmsh.model.getEntities(1)
    curve_tags = []
    for tag, curve in enumerate(curve_entities):
        mass = gmsh.model.occ.getMass(curve[0], curve[1])
        if mass - 2*np.pi*radius[0][0] < 1e-6 or mass - 2*np.pi*radius[0][N_cpt_v-1] < 1e-6 :
            curve_tags.append(curve[1])

    for tag in curve_tags:
        com = gmsh.model.occ.getCenterOfMass(1, tag)
        if np.isclose(com[2], z_start, rtol=1e-6):
            inlet_wire_tag = tag
        elif np.isclose(com[2], z_end, rtol=1e-6):
            outlet_wire_tag = tag

    # inlet surface
    W_inlet = gmsh.model.occ.addWire([inlet_wire_tag])
    inlet = gmsh.model.occ.addSurfaceFilling(W_inlet)
    gmsh.model.occ.synchronize()

    # outlet surface
    W_outlet = gmsh.model.occ.addWire([outlet_wire_tag])
    outlet = gmsh.model.occ.addSurfaceFilling(W_outlet)
    gmsh.model.occ.synchronize()

    surface_tags = [inlet, outlet, outer_surface]

    surface_loop = gmsh.model.occ.addSurfaceLoop(surface_tags)
    cylinder = gmsh.model.occ.addVolume([surface_loop])

    gmsh.model.occ.synchronize()

    return cylinder

# def nurbs_cylinder(z_start, z_end, radius, N_cpt_v, lc, center=[0,0,0], deg=2):

#     N_u = 9
    
#     knots_u= [0, 1/4, 1/2, 3/4, 1]
#     multiplicities_u = [3, 2, 2, 2, 3]
#     knots_array_u =[0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]

#     knots_v = np.linspace(0,1,N_cpt_v-1)
#     multiplicities_v = np.ones(len(knots_v))
#     knots_array_v = calculate_knot_v(knots_v, multiplicities_v, deg=2)

#     weights_surface = [1, 2**0.5/2, 1, 2**0.5/2, 1, 2**0.5/2, 1, 2**0.5/2, 1]
#     weights = np.tile(weights_surface,N_cpt_v)
#     z_wall = np.linspace(z_start,z_end,N_cpt_v)
    
#     vertex_list = []

#     for r_index, z in enumerate(z_wall):

#         vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[r_index], center[1], z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[r_index],center[1]+radius[r_index],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0],center[1]+radius[r_index],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0]-radius[r_index],center[1]+radius[r_index],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0]-radius[r_index],center[1],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0]-radius[r_index],center[1]-radius[r_index],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0],center[1]-radius[r_index],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[r_index],center[1]-radius[r_index],z,lc))
#         vertex_list.append(gmsh.model.occ.addPoint(center[0]+radius[r_index], center[1], z,lc))

#     outer_surface = gmsh.model.occ.addBSplineSurface(vertex_list, N_u, degreeU=2, degreeV=2,
#                                     weights=weights, knotsU=knots_u, multiplicitiesU=multiplicities_u,
#                                     knotsV=knots_v, multiplicitiesV=multiplicities_v)

#     gmsh.model.occ.synchronize()
#     curve_entities = gmsh.model.getEntities(1)
#     curve_tags = []
#     for tag, curve in enumerate(curve_entities):
#         mass = gmsh.model.occ.getMass(curve[0], curve[1])
#         if mass - 2*np.pi*radius[0] < 1e-6 or mass - 2*np.pi*list(radius)[-1] < 1e-6 :
#             curve_tags.append(curve[1])

#     for tag in curve_tags:
#         com = gmsh.model.occ.getCenterOfMass(1, tag)
#         if np.isclose(com[2], z_start, rtol=1e-6):
#             inlet_wire_tag = tag
#         elif np.isclose(com[2], z_end, rtol=1e-6):
#             outlet_wire_tag = tag

#     # inlet surface
#     W_inlet = gmsh.model.occ.addWire([inlet_wire_tag])
#     inlet = gmsh.model.occ.addSurfaceFilling(W_inlet)
#     gmsh.model.occ.synchronize()

#     # outlet surface
#     W_outlet = gmsh.model.occ.addWire([outlet_wire_tag])
#     outlet = gmsh.model.occ.addSurfaceFilling(W_outlet)
#     gmsh.model.occ.synchronize()

#     surface_tags = [inlet, outlet, outer_surface]

#     surface_loop = gmsh.model.occ.addSurfaceLoop(surface_tags)
#     outer_cylinder = gmsh.model.occ.addVolume([surface_loop])

#     gmsh.model.occ.synchronize()

#     return outer_cylinder

def nurbs_annulus(z_start, z_end, radius_outer, radius_inner, N_cpt_v, lc, deg=2):
    
    outer_plenum = nurbs_cylinder(z_start, z_end, radius_outer, N_cpt_v, lc, deg=deg)

    inner_plenum = nurbs_cylinder(z_start, z_end, radius_inner, N_cpt_v, lc, deg=deg)

    # Generate annulus
    annulus = gmsh.model.occ.cut([(3,outer_plenum)],[(3,inner_plenum)])

    gmsh.model.occ.synchronize()
    
    return annulus




if __name__ == "__main__":
    
    import sys

    # Mesh generation
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    N_u = 9
    N_cpt_v = 5

    radius = 0.2
    # Radius Dictionary
    radius_dict = {}
    for key_u in range(0,N_u):
        radius_dict[key_u] = {}                
        for key_v in range(0,N_cpt_v):
            radius_dict[key_u][key_v] = radius
    
    # perturbation index
    u_index = 2
    v_index = 2
    perturbation_magnitude = 0.05

    radius_dict[u_index][v_index] += perturbation_magnitude 

    z_plenum_start = 0.
    l_plenum = 0.07
    z_plenum_end = z_plenum_start + l_plenum

    lc_plenum = l_plenum/10

    N_cpt_plenum = 5 # number of control points in the +z direction 

    cylinder = nurbs_cylinder(z_plenum_start, z_plenum_end, radius_dict, N_cpt_v, lc_plenum)
    print("worked")
    gmsh.model.mesh.generate(3)

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()