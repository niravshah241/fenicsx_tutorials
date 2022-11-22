datfile = input("Input .dat file name (without .dat extension): \n")

""" next, we create the new geo file"""
mesh = open(datfile + '.py', 'r+')
mesh.write("import gmsh \n \n")
mesh.write("gmsh.initialize('',False) \n \n")
aero = open(datfile + '.dat', 'r')
mesh.write("lc = 0.001 \n \n")
""" this section creates the points for the mesh and set up the spline"""
linecnt = 0
spline = []

for line in aero:
        linecnt +=1
        new = line.split(' ')
        print(new)
        if len(new) == 16 and linecnt > 4:
            mesh.write("gmsh.model.occ.addPoint(" + str(new[-11]) + "," + str(new[-6]) + "," + str(new[-1]).replace('\n','') + "," + " lc" + "," + str(linecnt-4) + ")\n")
            spline.append(linecnt-4)
        elif len(new) == 15 and linecnt > 4:
            mesh.write("gmsh.model.occ.addPoint(" + str(new[-10]) + "," + str(new[-6]) + "," + str(new[-1]).replace('\n','') + "," + " lc" + "," + str(linecnt-4) + ")\n")
            spline.append(linecnt-4)

"""this section determines whether an extra line is required to close off the aerofoil"""
'''extraline = False
with open(datfile + '.dat','r')as aero:
        lines = aero.readlines()
        if lines[3] == lines[int(linecnt/2)+2]:
                extraline = True'''
"""time to close the original .dat file"""
aero.close()

""" this section created the spline"""
'''spline[:int((linecnt/2)-2)] = reversed(spline[:int(linecnt/2)-2])
splines ='Spline(1) ={'+(str(spline))+'}'
splines = splines.replace('[','')
splines = splines.replace(']','')
mesh.write(splines)'''
   
mesh.close()



'''file_name_dat = "aerofoil_data.geo" #input("Enter .dat file name: \n")
f = open(str(file_name_dat),"r")

f.readline()
f.readline()
f.readline()

g = open(str(file_name_dat)+"_gmsh.py","r+")
step_size = 0.0003
point_number = 1
for line in f:
    if True:#line.rstrip('\n'):  # Check if it's empty after removing EOL character
        #g.write("Point(" + str(point_number) + ") = gmsh.model.occ.addPoint(" + str(line[4:9]) + "," + str(line[15:20]) + ",0.," + str(step_size) + ")" + "\n")
        g.write("gmsh.model.occ.addPoint(" + str(line[4:9]) + "," + str(line[15:20]) + ",0.," + str(step_size) + "," + str(point_number) + ")" + "\n")
        point_number += 1
    else:
        print("Empty line")

f.close()'''


'''TODO
1. this will only add points in the file and nothing more. Cretaing spline and surfaces etc. is extra.
'''

#4:12 15:23

'''
spline_line = gmsh.model.occ.addSpline([1000,1010,1020,1030,1040,1000])
spline_curve = gmsh.model.occ.addCurveLoop([spline_curve],1)
spline_surface = gmsh.model.occ.addPlaneSurface([spline_curve],1)

gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.Algorithm", 8) # 8=Frontal-Delaunay for Quads (See section 7.4,  https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2) # 2=simple full-quad (See section 7.4,  https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options)
gmsh.option.setNumber("Mesh.RecombineAll", 1) # Apply recombination algorithm to all surfaces, ignoring per-surface spec
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) # Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3: barycentric)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1) # Minimum characteristic element size
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.3) # Maximum characteristic element size
gmsh.model.mesh.generate(gdim) # Mesh generation
gmsh.model.mesh.setOrder(1) # Mesh order
gmsh.model.mesh.optimize("Netgen") # Mesh optimisation or improving quality of mesh

gdim = 2

surfaces = gmsh.model.getEntities(dim=gdim)
print(surfaces)
'''
