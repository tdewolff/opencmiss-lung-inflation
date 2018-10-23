#!/usr/bin/env python

import os
from medpy.io import load
from medpy.io import header
import numpy as np
from opencmiss.iron import iron
from morphic.utils import convert_hermite_lagrange
import mesh_tools

path = '/hpc/tdew803/Lung/Data/Pig_PE_Study_HRC/'
subject = 'AP00149'
pressure = '10cmH2O'
registration = '10_25'

c10 = 5000.0  # in Pa
c01 = 2000.0  # in Pa
k = 6000.0  # in Pa
density = 1000.0  # in kg m^-3
gravity = [0.0, 0.0, 0.0]  # in m s^-2

interpolation = iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE

##################################################################

exelemFile = path + subject + '/' + pressure + '/Lung/FEMesh/Left_RefittedNoVersion.exelem'
exnodeFile = path + subject + '/' + pressure + '/Lung/FEMesh/Left_RefittedNoVersion.exnode'
coordinatesField = 'coordinates'
displacementFiles = path + subject + '/reg_' + registration + '_d%s.nii'
#displacementFiles = './displacement/AP00149_d%s.nii'

#exelemFile = './cube_hermite.exelem'
#exnodeFile = './cube_hermite.exnode'
#displacementFiles = './displacement/cube_d%s.nii'

def Pix2Coord(pixdim, i, j, k):
    x = i * pixdim[0]
    y = 256 - (j * pixdim[1])
    z = -k * pixdim[2]
    return x, y, z

def Coord2Pix(pixdim, x, y, z):
    i = x / pixdim[0]
    j = (256-y) / pixdim[1]
    k = -z / pixdim[2]
    return i, j, k

def TrilinearInterpolation(displacement, pixdim, X, Y, Z):
    (x,y,z) = Coord2Pix(pixdim,X,Y,Z)
    x0 = np.floor(x)
    y0 = np.floor(y)
    z0 = np.floor(z)
    x1 = np.ceil(x)
    y1 = np.ceil(y)
    z1 = np.ceil(z)
    if x0 < 0 or y0 < 0 or z0 < 0 or x1 >= displacement.shape[0] or y1 >= displacement.shape[1] or z1 >= displacement.shape[2]:
        print(x0,y0,z0, x1,y1,z1)
        raise ValueError('TrilinearInterpolation: probe is outside displacement field')
    xd = 0.0
    yd = 0.0
    zd = 0.0
    if x1 > x0:
        xd = (x - x0) / (x1 - x0)
    if y1 > y0:
        yd = (y - y0) / (y1 - y0)
    if z1 > z0:
        zd = (z - z0) / (z1 - z0)
    c00 = displacement[x0, y0, z0] * (1 - xd) + displacement[x1, y0, z0] * xd
    c01 = displacement[x0, y0, z1] * (1 - xd) + displacement[x1, y0, z1] * xd
    c10 = displacement[x0, y1, z0] * (1 - xd) + displacement[x1, y1, z0] * xd
    c11 = displacement[x0, y1, z1] * (1 - xd) + displacement[x1, y1, z1] * xd
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    return c0 * (1 - zd) + c1 * zd

def Run(cubic_lagrange_morphic_mesh, interpolation, displacement, pixdim, c10, c01, k, density, gravity):
    # Get the number of computational nodes and this computational node number
    numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()

    # Set all diagnostic levels on for testing
    # iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL, [1, 2, 3, 4, 5], "Diagnostics",
    #                      ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

    coordinateSystemUserNumber = 1
    regionUserNumber = 1
    basisUserNumber = 1
    meshUserNumber = 1
    decompositionUserNumber = 1
    geometricFieldUserNumber = 1
    fibreFieldUserNumber = 2
    materialFieldUserNumber = 3
    dependentFieldUserNumber = 4
    sourceFieldUserNumber = 5
    equationsSetFieldUserNumber = 6
    equationsSetUserNumber = 1
    problemUserNumber = 1

    # Create a 3D rectangular cartesian coordinate system
    coordinateSystem = iron.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.DimensionSet(3)
    coordinateSystem.CreateFinish()

    # Create a region and assign the coordinate system to the region
    region = iron.Region()
    region.CreateStart(regionUserNumber, iron.WorldRegion)
    region.LabelSet("Region")
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()

    # Define basis
    basis = iron.Basis()
    basis.CreateStart(basisUserNumber)
    basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    basis.numberOfXi = 3
    basis.interpolationXi = [interpolation] * 3
    if interpolation == iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE:
        basis.quadratureNumberOfGaussXi = [2] * 3
    elif interpolation == iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE:
        basis.quadratureNumberOfGaussXi = [3] * 3
    elif interpolation == iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE:
        basis.quadratureNumberOfGaussXi = [4] * 3
    elif interpolation == iron.BasisInterpolationSpecifications.CUBIC_HERMITE:
        basis.quadratureNumberOfGaussXi = [4] * 3
    basis.CreateFinish()

    mesh, coordinates, node_nums, element_nums = mesh_tools.morphic_to_OpenCMISS(cubic_lagrange_morphic_mesh, region,
                                                                                 basis, meshUserNumber,
                                                                                 dimension=3,
                                                                                 interpolation='cubic')

    print('Number of nodes:', len(node_nums))
    print('Number of elements:', len(element_nums))

    # Create a decomposition for the mesh
    decomposition = iron.Decomposition()
    decomposition.CreateStart(decompositionUserNumber, mesh)
    decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
    decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    decomposition.CreateFinish()

    # Create a field for the geometry
    geometricField = iron.Field()
    geometricField.CreateStart(geometricFieldUserNumber, region)
    geometricField.MeshDecompositionSet(decomposition)
    geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
    geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    geometricField.CreateFinish()

    # Update the geometric field parameters
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    for node_idx, node in enumerate(node_nums):
        for component_idx, component in enumerate([1, 2, 3]):
            for derivative_idx, derivative in enumerate(range(1, coordinates.shape[2] + 1)):
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                        1, derivative, node, component,
                                                        coordinates[node_idx, component_idx, derivative_idx])

    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

    # Create a fibre field and attach it to the geometric field
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber, region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
    fibreField.CreateFinish()

    # stressField = iron.Field()
    # stressField.CreateStart(stressFieldUserNumber, region)
    # stressField.MeshDecompositionSet(decomposition)
    # stressField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.ELEMENT_BASED)
    # stressField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.ELEMENT_BASED)
    # stressField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 3, iron.FieldInterpolationTypes.ELEMENT_BASED)
    # stressField.VariableLabelSet(iron.FieldVariableTypes.U, "Cauchy Stress")
    # stressField.CreateFinish()

    ##################################################################
    # Setup Mooney-Rivlin equations
    ##################################################################

    equationsSetField = iron.Field()
    equationsSet = iron.EquationsSet()
    equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                 iron.EquationsSetTypes.FINITE_ELASTICITY,
                                 iron.EquationsSetSubtypes.COMPRESSIBLE_FINITE_ELASTICITY]
    equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, equationsSetSpecification,
                             equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()

    # Setup material field
    materialField = iron.Field()
    equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
    materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
    materialField.VariableLabelSet(iron.FieldVariableTypes.V, "Density")
    equationsSet.MaterialsCreateFinish()

    print("Material field number of components (U,V):", materialField.NumberOfComponentsGet(iron.FieldVariableTypes.U),
          materialField.NumberOfComponentsGet(iron.FieldVariableTypes.V))

    # Set Mooney-Rivlin constants c10 and c01 respectively.
    materialField.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, c10)
    materialField.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, c01)
    materialField.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, k)
    materialField.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.V, iron.FieldParameterSetTypes.VALUES, 1, density)

    # Setup dependent field
    dependentField = iron.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
    equationsSet.DependentCreateFinish()

    print("Dependent field number of components (U,DELUDELN):",
          dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.U),
          dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.DELUDELN))

    if gravity != [0.0, 0.0, 0.0]:
        # Setup gravity source field
        sourceField = iron.Field()
        equationsSet.SourceCreateStart(sourceFieldUserNumber, sourceField)
        sourceField.fieldScalingType = iron.FieldScalingTypes.UNIT
        equationsSet.SourceCreateFinish()

        # Set the gravity vector component values
        sourceField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, gravity[0])
        sourceField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, gravity[1])
        sourceField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, gravity[2])

    # Create equations
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    equations.outputType = iron.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    # Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
        dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
        dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3,
        dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3)

    ##################################################################
    # Define the problem
    ##################################################################

    problem = iron.Problem()
    problemSpecification = [iron.ProblemClasses.ELASTICITY,
                            iron.ProblemTypes.FINITE_ELASTICITY,
                            iron.ProblemSubtypes.NONE]
    problem.CreateStart(problemUserNumber, problemSpecification)
    problem.CreateFinish()

    # Create control loops
    problem.ControlLoopCreateStart()
    problem.ControlLoopCreateFinish()

    # Create problem solver
    nonLinearSolver = iron.Solver()
    linearSolver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
    nonLinearSolver.outputType = iron.SolverOutputTypes.MONITOR
    nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
    nonLinearSolver.NewtonLinearSolverGet(linearSolver)
    nonLinearSolver.NewtonAbsoluteToleranceSet(1e-11)
    nonLinearSolver.NewtonSolutionToleranceSet(1e-11)
    nonLinearSolver.NewtonRelativeToleranceSet(1e-11)
    linearSolver.linearType = iron.LinearSolverTypes.DIRECT
    # linearSolver.libraryType = iron.SolverLibraries.LAPACK
    problem.SolversCreateFinish()

    # Create solver equations and add equations set to solver equations
    solver = iron.Solver()
    solverEquations = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
    solver.SolverEquationsGet(solverEquations)
    solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    # Prescribe boundary conditions (absolute nodal parameters)
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    #dim = displacement['x'].shape
    #avgRadius = 0  # in voxels
    nodes = iron.MeshNodes()
    mesh.NodesGet(1, nodes)
    for nid in node_nums:
        if nodes.NodeOnBoundaryGet(nid) == iron.MeshBoundaryTypes.ON:
            X = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                     nid, 1)
            Y = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                     nid, 2)
            Z = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                     nid, 3)

            # if avgRadius > 0:
            #     dx = 0
            #     dy = 0
            #     dz = 0
            #     n = 0
            #     for k in range(max(0, int(np.floor(K - avgRadius))), min(dim[2], int(np.ceil(K + avgRadius) + 1))):
            #         for j in range(max(0, int(np.floor(J - avgRadius))), min(dim[1], int(np.ceil(J + avgRadius) + 1))):
            #             for i in range(max(0, int(np.floor(I - avgRadius))), min(dim[0], int(np.ceil(I + avgRadius) + 1))):
            #                 if (i - I) ** 2 + (j - J) ** 2 + (k - K) ** 2 <= avgRadius ** 2:
            #                     dx += displacement['x'][i, j, k]
            #                     dy += displacement['y'][i, j, k]
            #                     dz += displacement['z'][i, j, k]
            #                     n += 1
            #
            #     if n == 0:
            #         print('No displacement field defined at: ', [X, Y, Z])
            #         continue
            #     dx /= n
            #     dy /= n
            #     dz /= n
            # else:

            (dx, dy, dz) = TrilinearInterpolation(displacement, pixdim, X, Y, Z)

            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
                                       iron.BoundaryConditionsTypes.FIXED, -dx)
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
                                       iron.BoundaryConditionsTypes.FIXED, -dy)
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
                                       iron.BoundaryConditionsTypes.FIXED, -dz)

            print('BC set on node %d at %.2f %.2f %.2f += %f %f %f' % (nid, X,Y,Z, dx,dy,dz))

    solverEquations.BoundaryConditionsCreateFinish()

    ##################################################################
    # Solve!
    ##################################################################
    problem.Solve()

    maskFile = open('./results/mask.exdata', 'w')
    maskFile.write(' Group name: Mask\n')
    maskFile.write(' #Fields=2\n')
    maskFile.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    maskFile.write('   x.  Value index= 1, #Derivatives=0\n')
    maskFile.write('   y.  Value index= 2, #Derivatives=0\n')
    maskFile.write('   z.  Value index= 3, #Derivatives=0\n')
    maskFile.write(' 2) MaskValue, field, rectangular cartesian, #Components=1\n')
    maskFile.write('   1.  Value index= 4, #Derivatives=0\n')
    mask, _ = load(path + subject + '/' + pressure + '/Lung/' + subject + '.Lung.mask.img')
    nid = 0
    for k in range(0, mask.shape[2], 10):
        for j in range(0, mask.shape[1], 10):
            for i in range(0, mask.shape[0], 10):
                if mask[i,j,k] == 5:
                    (x,y,z) = Pix2Coord(pixdim,i,j,k)

                    nid += 1
                    maskFile.write(' Node: %d\n' % nid)
                    maskFile.write('  %E %E %E %E\n' % (x,y,z,mask[i,j,k]))
    maskFile.close()

    minX = float('inf')
    minY = float('inf')
    minZ = float('inf')
    maxX = -float('inf')
    maxY = -float('inf')
    maxZ = -float('inf')
    for nid in node_nums:
        x = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 1)
        y = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 2)
        z = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 3)
        minX = min(minX, x)
        minY = min(minY, y)
        minZ = min(minZ, z)
        maxX = max(maxX, x)
        maxY = max(maxY, y)
        maxZ = max(maxZ, z)
    print('x = (%f,%f)' % (minX, maxX))
    print('y = (%f,%f)' % (minY, maxY))
    print('z = (%f,%f)' % (minZ, maxZ))

    #
    # Expected values for strain tensor in the cube example are:
    # 0.22 0 0
    # 0    0 0
    # 0    0 0
    #
    # Expected values for stress tensor in the cube example are:
    # 24000     0     0
    #     0 18867     0
    #     0     0 18867
    #

    valuesFile = open('./results/values.exdata', 'w')
    valuesFile.write(' Group name: Values\n')
    valuesFile.write(' #Fields=3\n')
    valuesFile.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    valuesFile.write('   x.  Value index= 1, #Derivatives=0\n')
    valuesFile.write('   y.  Value index= 2, #Derivatives=0\n')
    valuesFile.write('   z.  Value index= 3, #Derivatives=0\n')
    valuesFile.write(' 2) Strain, field, rectangular cartesian, #Components=1\n')
    valuesFile.write('   1.  Value index= 4, #Derivatives=0\n')
    valuesFile.write(' 3) Stress, field, rectangular cartesian, #Components=1\n')
    valuesFile.write('   1.  Value index= 5, #Derivatives=0\n')

    elements = iron.MeshElements()
    mesh.ElementsGet(1, elements)

    valuesSizes = (3, 3)
    valuesInterp = 2
    nid = 0
    evaluatedCoords = set()
    print(element_nums)
    for eid in element_nums:
        for xiZ in np.linspace(0.0, 1.0, valuesInterp):
            for xiY in np.linspace(0.0, 1.0, valuesInterp):
                for xiX in np.linspace(0.0, 1.0, valuesInterp):
                    coords = dependentField.ParameterSetInterpolateSingleXiDP(iron.FieldVariableTypes.U,
                                                                              iron.FieldParameterSetTypes.VALUES,
                                                                              1, eid + 1, [xiX, xiY, xiZ], 3)
                    if tuple(coords) not in evaluatedCoords:
                        evaluatedCoords.add(tuple(coords))
                        strain = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.GREEN_LAGRANGE_STRAIN,
                                                                  eid, [xiX, xiY, xiZ], valuesSizes)
                        stress = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS, eid,
                                                                  [xiX, xiY, xiZ], valuesSizes)

                        #print('Coord:', [xiX, xiY, xiZ], coords)
                        if np.isfinite(strain).all() and np.isfinite(stress).all():
                            eigStrain = np.linalg.eigvals(strain)
                            eigStress = np.linalg.eigvals(stress)

                            avgStrain = np.mean(eigStrain)
                            avgStress = np.mean(eigStress)

                            nid += 1
                            valuesFile.write(' Node: %d\n' % nid)
                            valuesFile.write('  %E %E %E %E %E\n' % (coords[0], coords[1], coords[2], avgStrain, avgStress))

                            #print('Strain:', strain, eigStrain, avgStrain)
                            #print('Stress:', stress, eigStress, avgStress)
                        else:
                            print('Strain or stress had non-finite values at element ID %d and coords %E %E %E' % (eid, coords[0], coords[1], coords[2]))
                        #print()

    valuesFile.close()

    fields = iron.Fields()
    fields.CreateRegion(region)
    return fields


if not os.path.exists("./results"):
    os.makedirs("./results")

print('Load dx...')
dx, image_header = load(displacementFiles % ('x'))
displacement = np.zeros((dx.shape[0], dx.shape[1], dx.shape[2], 3))
displacement[:, :, :, 0] = dx
print('Load dy...')
displacement[:, :, :, 1], _ = load(displacementFiles % ('y'))
print('Load dz...')
displacement[:, :, :, 2], _ = load(displacementFiles % ('z'))
pixdim = header.get_pixel_spacing(image_header)
print('Image dimensions:', displacement.shape)
print('Voxel dimensions:', pixdim)

cubic_hermite_morphic_mesh = mesh_tools.exfile_to_morphic(exnodeFile, exelemFile, coordinatesField, dimension=3,
                                                          interpolation='hermite')
cubic_lagrange_morphic_mesh = convert_hermite_lagrange(cubic_hermite_morphic_mesh, tol=1e-9)

(x,y,z) = Pix2Coord(pixdim,100,100,100)
(i,j,k) = Coord2Pix(pixdim,x,y,z)
if i != 100 or j != 100 or k != 100:
    raise ValueError('Coord2Pix(Pix2Coord(100,100,100)) != (100,100,100)')

fields = Run(cubic_lagrange_morphic_mesh, interpolation, displacement, pixdim, c10, c01, k, density, gravity)
fields.NodesExport("./results/out", "FORTRAN")
fields.ElementsExport("./results/out", "FORTRAN")
fields.Finalise()
