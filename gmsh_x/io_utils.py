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