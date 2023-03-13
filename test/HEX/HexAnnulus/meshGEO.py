# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 6
#
#  Transfinite meshes
#
# ------------------------------------------------------------------------------

import gmsh
import math
import sys

mesh_dir = "MeshDir/"
geom_name = "AnnulusGEOHEX"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("Structured")

lc = 1e-2
r_outer = 0.25
r_inner = 0.2 
L = 1

p0 = gmsh.model.geo.addPoint(0, 0, 0, lc,99)
# Outer Circle

p1 = gmsh.model.geo.addPoint(r_outer, 0, 0, lc,1)
p2 = gmsh.model.geo.addPoint(0, r_outer, 0, lc,2)

p3 = gmsh.model.geo.addPoint(r_inner, 0, 0, lc,3)
p4 = gmsh.model.geo.addPoint(0, r_inner, 0, lc,4)

l1 = gmsh.model.geo.addLine(3, 1,1)
l2 = gmsh.model.geo.addCircleArc(p1, p0, p2)
l3 = gmsh.model.geo.addLine(2, 4,3)
l4 = gmsh.model.geo.addCircleArc(p4, p0, p3)
# gmsh.model.geo.addLine(4, 3,4)

gmsh.model.geo.addCurveLoop([l1, l2, l3, l4], 1)
gmsh.model.geo.addPlaneSurface([1],1)


N_radius = 5
gmsh.model.geo.mesh.setTransfiniteCurve(1, N_radius)
gmsh.model.geo.mesh.setTransfiniteCurve(3, N_radius)
N_theta = 20
gmsh.model.geo.mesh.setTransfiniteCurve(2, N_theta)
gmsh.model.geo.mesh.setTransfiniteCurve(4, N_theta)
N_zet = 50

gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])


# To create quadrangles instead of triangles, one can use the `setRecombine'
# constraint:

# ov2 = gmsh.model.geo.extrude([(2,1)], 0, 0, 0.12)

gmsh.model.geo.mesh.setRecombine(2, 1)

ov = gmsh.model.geo.extrude([(2, 1)], 0, 0, L, [N_zet], recombine=True)

ov = gmsh.model.geo.revolve([(2, 21)], 0, 0, 0, 0, 0, 1, +math.pi/2 ,[N_theta], recombine=True)


gmsh.model.geo.mesh.setRecombine(3, 2)

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(3)

def mirror_mesh_x_axis():
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

    transform(m, 1000, 100000, 100, 1, -1, 1)
    gmsh.model.mesh.removeDuplicateNodes()

mirror_mesh_x_axis()

gmsh.model.geo.synchronize()

# Physical tags

vols = gmsh.model.getEntities(dim=3)
vol_tags = [el[1] for el in vols]
print("Vol tags: ", vols)
gmsh.model.addPhysicalGroup(3, vol_tags, tag=1) 

surfaces = gmsh.model.getEntities(dim=2)

for surface in surfaces:

    gmsh.model.addPhysicalGroup(surface[0], [surface[1]])

gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write("{}.msh".format(mesh_dir+geom_name))

# Launch the GUI to see the results:
# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()

gmsh.finalize()

from helmholtz_x.dolfinx_utils import write_xdmf_hexmesh, XDMFReader
write_xdmf_hexmesh(mesh_dir+geom_name,dimension=3)

geom = XDMFReader(mesh_dir+geom_name)
geom.getInfo()

