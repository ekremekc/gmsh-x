from gmsh_x.nurbs import nurbs_cylinder, nurbs_annulus, generate_radius_dict
from gmsh_x.mesh_utils import copyAndRotate
from helmholtz_x.dolfinx_utils import dict_loader
import numpy as np
import sys
import gmsh

def micca_parametric(mesh_name):

    N_sector = 16
    N_u = 9

    # Mesh generation
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Plenum geometry
    r_plenum_outer_base = 0.21
    r_plenum_inner_base = 0.14

    z_plenum_start = 0.
    l_plenum = 0.07
    z_plenum_end = z_plenum_start + l_plenum

    lc_plenum = l_plenum/10 

    N_cpt_plenum = 5 # number of control points in the +z direction 

    r_plenum_outer = generate_radius_dict(r_plenum_outer_base, N_u, N_cpt_plenum)
    r_plenum_inner = generate_radius_dict(r_plenum_inner_base, N_u, N_cpt_plenum)

    plenum = nurbs_annulus(z_plenum_start, z_plenum_end, r_plenum_outer, r_plenum_inner, N_cpt_plenum, lc_plenum, deg=2)

    import math
    # def apply_fragment(volume):
    def generate_slicers():
        # Create the first cutting plane:
        L = 0.5
        H = 0.5
        xmin, ymin, zmin = 0,0,0

        slicer1 = []
        slicer1.append((2, gmsh.model.occ.addRectangle(xmin, ymin, zmin, L, H)))
        gmsh.model.occ.rotate([slicer1[0]], xmin, ymin, zmin, 1, 0, 0, math.pi/2)

        slicer2 = []
        slicer2.append((2, gmsh.model.occ.addRectangle(xmin, ymin, zmin, L, H)))
        gmsh.model.occ.rotate([slicer2[0]], xmin, ymin, zmin, 1, 0, 0, math.pi/2)
        gmsh.model.occ.rotate([slicer2[0]], xmin, ymin, zmin, 0, 0, 1, math.pi/16)

        gmsh.model.occ.synchronize()
        
        return slicer1, slicer2

    def clean_partition(three_dim_tag, index=-1):
        gmsh.model.occ.fragment([(3,three_dim_tag)], slicer1)
        gmsh.model.occ.fragment([(3,three_dim_tag)], slicer2)
        gmsh.model.occ.synchronize()

        gmsh.model.occ.remove(gmsh.model.occ.getEntities(2), True)
        
        gmsh.model.occ.synchronize()
        to_be_removed = gmsh.model.getEntities(3)[index]
        gmsh.model.occ.remove([to_be_removed], True)
        gmsh.model.occ.synchronize()
    
    slicer1, slicer2 = generate_slicers()

    clean_partition(plenum[0][0][1])
    
    # Burner geometry
    z_burner_start = z_plenum_end
    l_burner = 0.014 
    z_burner_end = z_burner_start + l_burner
    r_burner_base = 0.033/2
    r_center = (r_plenum_outer_base+r_plenum_inner_base)/2 
    burner_center=[r_center,0,z_burner_start]

    N_cpt_burner = 3

    r_burner = generate_radius_dict(r_burner_base, N_u, N_cpt_burner)

    lc_burner = l_burner/3

    single_burner = nurbs_cylinder(z_burner_start, z_burner_end, r_burner, N_cpt_burner, lc_burner, center=burner_center)

    slicer1, slicer2 = generate_slicers()
    clean_partition(single_burner)

    # Perforated Plate Geometry
    z_perforated_plate_start = z_burner_end 
    l_perforated_plate = 0.006
    z_perforated_plate_end = z_perforated_plate_start + l_perforated_plate
    r_perforated_plate_base = 0.0189/2

    N_cpt_perforated_plate = 3

    r_perforated_plate = generate_radius_dict(r_perforated_plate_base, N_u, N_cpt_perforated_plate)

    lc_perforated_plate = l_perforated_plate/3

    single_perforated_plate = nurbs_cylinder(z_perforated_plate_start, z_perforated_plate_end, 
                            r_perforated_plate, N_cpt_perforated_plate, 
                            lc_perforated_plate, center=burner_center)

    slicer1, slicer2 = generate_slicers()
    clean_partition(single_perforated_plate)

    # Fusion the upstream
    upstream = gmsh.model.occ.fuse(plenum[0], [(3,single_burner),(3,single_perforated_plate)])
    gmsh.model.occ.synchronize()

    # # # Combustion Chamber Geometry

    z_cc_start = z_perforated_plate_end
    l_cc = 0.2
    l_ec = 0.041
    z_cc_end = z_cc_start + l_cc + l_ec
    r_cc_outer_base = 0.20
    r_cc_inner_base = 0.15

    N_wall_cc = 5

    r_cc_outer = generate_radius_dict(r_cc_outer_base, N_u, N_wall_cc)

    r_cc_inner = generate_radius_dict(r_cc_inner_base, N_u, N_wall_cc)

    lc_cc = l_cc/20

    combustion_chamber = nurbs_annulus(z_cc_start, z_cc_end, r_cc_outer, r_cc_inner, N_wall_cc, lc_cc)

    slicer1, slicer2 = generate_slicers()
    clean_partition(combustion_chamber[0][0][1])

    # Flame geometry
    z_flame_start = z_cc_start
    l_flame = 0.006 
    z_flame_end = z_flame_start + l_flame
    r_flame_base = 0.036/2

    N_cpt_flame = 3

    r_flame = generate_radius_dict(r_flame_base, N_u, N_cpt_flame)

    lc_flame = l_flame/2

    single_flame = nurbs_cylinder(z_flame_start, z_flame_end, r_flame, N_cpt_flame, lc_flame, center=burner_center)
    slicer1, slicer2 = generate_slicers()
    clean_partition(single_flame, index=-2)
    gmsh.model.occ.synchronize()

    # Seperate CC and Flame 
    downstream = gmsh.model.occ.cut(combustion_chamber[0], [(3,4)],removeTool=False, removeObject=True)
    gmsh.model.occ.synchronize()
    print(gmsh.model.occ.getEntities(3))

    # Fuse everything
    half_sector = gmsh.model.occ.fuse([(3,1),(3,2)], [(3,4)], removeTool=False, removeObject=True)
    gmsh.model.occ.synchronize()

    gmsh.model.occ.removeAllDuplicates()
    # gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.occ.synchronize()

    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    # gmsh.option.setNumber('Mesh.MeshSizeMin', 0.002) 
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.02) # 0.004

    gmsh.model.mesh.generate(3)
    
    # Add surfaces to the physical group

    # # Elementary tags
    plenum_inlet_tag_el = 1
    plenum_outlet_tag_el = 6
    plenum_inner_lateral_tag_el = 3
    plenum_outer_lateral_tag_el = 5
    cc_inlet_tag_el = 10
    cc_inner_lateral_tag_el = 11
    cc_outer_lateral_tag_el = 13
    cc_outlet_tag_el = 12

    # # Physical tags 
    plenum_inlet_tag = 2
    plenum_outlet_tag = 3
    plenum_outer_lateral_tag = 1
    plenum_inner_lateral_tag = 4
    cc_inlet_tag = 5 
    cc_inner_lateral_tag = 6 
    cc_outer_lateral_tag = 7
    cc_outlet_tag = 11

    # gmsh.model.addPhysicalGroup(2,[plenum_inlet_tag_el], plenum_inlet_tag)
    
    # gmsh.model.addPhysicalGroup(2,[cc_outlet_tag_el], cc_outlet_tag)


    from gmsh_x.mesh_utils import getMeshData, rotate_z, mirror_mesh_x_axis
    
    # Mirror
    mirror_mesh_x_axis([1,-1,1])

    # Copy and rotate
    m = getMeshData(gmsh.model)

    # angles = np.linspace(360/16, 15*360/16, 15)
    angles = [360/16, 360/8, 360/4, 360/2]
    # angles = [360/16]
    for angle in angles:
        print(angle)
        m = getMeshData(gmsh.model)

        max_node_tag = gmsh.model.mesh.getMaxNodeTag() + 1000
        max_element_tag = gmsh.model.mesh.getMaxElementTag() + 1000
        print(max_node_tag)
        rotate_z(m, max_node_tag, max_node_tag, max_element_tag, angle)

    gmsh.model.occ.synchronize()

    # gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.occ.removeAllDuplicates()
    
    gmsh.model.occ.synchronize()


    # TAGGING
    plenum_inlet_tags = []
    cc_outlet_tags = []

    # Surface tags
    boundaries_entities = [i[1] for i in gmsh.model.getEntities(2)]
    boundary_tags = np.array(boundaries_entities).reshape(-1,20)

    gmsh.model.addPhysicalGroup(2,boundary_tags[:,15], plenum_inlet_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,13], plenum_outlet_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,14], plenum_outer_lateral_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,16], plenum_inner_lateral_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,0], cc_inlet_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,3], cc_inner_lateral_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,5], cc_outer_lateral_tag)
    gmsh.model.addPhysicalGroup(2,boundary_tags[:,7], cc_outlet_tag)


    # Volume tags
    vol_entities = [i[1] for i in gmsh.model.getEntities(3)] # get just tags
    vol_tags = np.array(vol_entities).reshape(-1,3)          # reshape tags
    non_flame_tags = np.delete(vol_tags,0,axis=1).flatten()  # clean flame tags
    flame_tags = vol_tags[:,0].reshape(-1,2)                 # get only flame tags

    gmsh.model.addPhysicalGroup(3,non_flame_tags.tolist(), 16)
    for tag, flame in enumerate(flame_tags):
        gmsh.model.addPhysicalGroup(3,flame, tag)


    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.write("{}.msh".format(mesh_name))

    return gmsh.model

if __name__=="__main__":

    import datetime
    start_time = datetime.datetime.now()
    
    mesh_name = "MeshDir/MICCA_parametric"
    model = micca_parametric(mesh_name)

    from helmholtz_x.dolfinx_utils import write_xdmf_mesh, XDMFReader

    write_xdmf_mesh(mesh_name,dimension=3)
    Micca = XDMFReader(mesh_name)
    Micca.getInfo()

    print("Total Execution Time for Mesh Generation: ", datetime.datetime.now()-start_time)
