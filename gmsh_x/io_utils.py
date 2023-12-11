import meshio

def create_mesh(mesh, cell_type, prune_z):
    """Subroutine for mesh creation by using meshio library

    Args:
        mesh (meshio._mesh.Mesh): mesh to be converted into Dolfinx mesh
        cell_type ('str'): type of cell (it becomes tetrahedral most of the time)
        prune_z ('bool'): whether consider the 3th dimension's coordinate or not, (it should be False for 2D cases)

    Returns:
        meshio._mesh.Mesh: converted xdmf mesh
    """

    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:,:2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    
    if cell_type=="tetra":
        print("Number of 3D cells:  {:,}".format(len(cells)))
    elif cell_type=="triangle":
        print("Number of 2D cells (facets):  {:,}".format(len(cells)))

    return out_mesh

def write_xdmf_mesh(name, dimension, write_edge=False):
    """Writes gmsh (.msh) mesh as an .xdmf mesh

    Args:
        name ('str'): filename
        dimension ('int'): Dimension of the mesh (2 or 3)
    """

    # Read in mesh
    msh_name = name + ".msh"
    msh = meshio.read(msh_name)

    if dimension == 1:
        prune_z = True
        volume_mesh = create_mesh(msh, "line",prune_z)
        tag_mesh = create_mesh(msh, "vertex",prune_z)

    if dimension == 2:
        prune_z = True
        volume_mesh = create_mesh(msh, "triangle",prune_z)
        tag_mesh = create_mesh(msh, "line",prune_z)

    elif dimension == 3:
        prune_z = False
        volume_mesh = create_mesh(msh, "tetra",prune_z)
        tag_mesh = create_mesh(msh, "triangle",prune_z)
        if write_edge:
            edge_mesh = create_mesh(msh, "line",prune_z)
            xdmf_edge_name = name + "_edgetags.xdmf"
            meshio.write(xdmf_edge_name, edge_mesh)
        
    # Create and save one file for the mesh, and one file for the facets and one file for the edges
    xdmf_name = name + ".xdmf"
    xdmf_tags_name = name + "_tags.xdmf"
    
    meshio.write(xdmf_name, volume_mesh)
    meshio.write(xdmf_tags_name, tag_mesh)

    print(str(dimension)+"D XDMF mesh is generated.")

def animate_pvds(pvdfile, filename, num_timesteps):

    file1 = open(pvdfile, "w") 

    header = '<?xml version="1.0"?>'
    subheader1 = '<VTKFile type="Collection" version="0.1"'
    subheader2 = '         byte_order="LittleEndian"'
    subheader3 = '         compressor="vtkZLibDataCompressor">'

    file1.write(header)
    file1.write("\n")
    file1.write(subheader1)
    file1.write("\n")
    file1.write(subheader2)
    file1.write("\n")
    file1.write(subheader3)
    file1.write("\n")

    start_collection = '  <Collection>'
    file1.write(start_collection)
    file1.write("\n")
    for i in range(num_timesteps):
        timestep = str(i)
        name = filename+timestep+'_p0_000000.vtu'

        dataset = '    <DataSet timestep="'+timestep+'" group="" part="0"'
        file1.write(dataset)
        file1.write("\n")
        datafile = '             file="'+name+'"/>'
        file1.write(datafile)
        file1.write("\n")

    # name2 = '../VTK/CylinderDeformed1_p0_000000.vtu'

    # dataset = '    <DataSet timestep="'+str(1)+'" group="" part="0"'
    # file1.write(dataset)
    # file1.write("\n")
    # datafile = '             file="'+name2+'"/>'
    # file1.write(datafile)
    # file1.write("\n")

    upfooter = '  </Collection>'
    file1.write(upfooter)
    file1.write("\n")
    footer = '</VTKFile>'
    file1.write(footer)

    file1.close()