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
pressure = '5cmH2O'
registration = '5_10'

c10 = 5000.0  # in Pa
c01 = 2000.0  # in Pa
k = 6000.0  # in Pa
density = 1000.0  # in kg m^-3
gravity = [0.0, 0.0, 0.0]  # in m s^-2

stretch = 1.2

lengthX = 0.1
lengthY = 0.1
lengthZ = 0.1

interpolation = iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE

##################################################################

exelemFile = path + subject + '/' + pressure + '/Lung/FEMesh/Left_RefittedNoVersion.exelem'
exnodeFile = path + subject + '/' + pressure + '/Lung/FEMesh/Left_RefittedNoVersion.exnode'

exelemFile = './cube.exelem'
exnodeFile = './cube.exnode'

displacementFiles = path + subject + '/reg_' + registration + '_d%s.nii'
displacement = {}

print('Load dx...')
displacement['x'], image_header = load(displacementFiles % ('x'))
print('Load dy...')
displacement['y'], _ = load(displacementFiles % ('y'))
print('Load dz...')
displacement['z'], _ = load(displacementFiles % ('z'))

displacementPixdim = header.get_pixel_spacing(image_header)

print('Image dimensions:', displacement['x'].shape)
print('Voxel dimensions:', displacementPixdim)

def Run(exelemFile, exnodeFile, interpolation, displacement, displacementPixdim, c10, c01, k, density, gravity):

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

    #mesh, coordinates, node_nums = mesh_tools.exfile_to_OpenCMISS(exnodeFile, exelemFile, 'Geometry', region, basis, meshUserNumber, dimension=3,
    #                                                          interpolation='hermite')

    cubic_hermite_morphic_mesh = mesh_tools.exfile_to_morphic(exnodeFile, exelemFile, 'Geometry', dimension=3,
                                                              interpolation='hermite')
    print('Number of hermite nodes: ', len(cubic_hermite_morphic_mesh.get_nodes()))
    cubic_lagrange_morphic_mesh = convert_hermite_lagrange(cubic_hermite_morphic_mesh, tol=1e-9)
    print('Number of lagrange nodes: ', len(cubic_lagrange_morphic_mesh.get_nodes()))
    mesh, coordinates, node_nums, element_nums = mesh_tools.morphic_to_OpenCMISS(cubic_lagrange_morphic_mesh, region,
                                                                                 basis, meshUserNumber,
                                                                                 dimension=3,
                                                                                 interpolation='cubicLagrange')
    #mesh, coordinates, node_nums, element_nums = mesh_tools.morphic_to_OpenCMISS(cubic_hermite_morphic_mesh, region,
    #                                                                             basis, meshUserNumber,
    #                                                                             dimension=3,
    #                                                                             interpolation='hermite')

    #mesh = iron.Mesh()
    #generatedMesh = iron.GeneratedMesh()
    #generatedMesh.CreateStart(generatedMeshUserNumber, region)
    #generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
    #generatedMesh.basis = [basis]
    #generatedMesh.extent = [lengthX, lengthY, lengthZ]
    #generatedMesh.numberOfElements = [numberXElements, numberYElements, numberZElements]
    #generatedMesh.CreateFinish(meshUserNumber, mesh)

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
    #geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
    geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    geometricField.CreateFinish()
    #generatedMesh.GeometricParametersCalculate(geometricField)

    # Update the geometric field parameters
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    for node_idx, node in enumerate(node_nums):
        for component_idx, component in enumerate([1, 2, 3]):
            for derivative_idx, derivative in enumerate(range(1, coordinates.shape[2] + 1)):
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                        1, derivative, node, component,
                                                        coordinates[node_idx, component_idx, derivative_idx])

    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


    fields = iron.Fields()
    fields.CreateRegion(region)
    return fields

    # Create a fibre field and attach it to the geometric field
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber, region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
    fibreField.CreateFinish()

    #stressField = iron.Field()
    #stressField.CreateStart(stressFieldUserNumber, region)
    #stressField.MeshDecompositionSet(decomposition)
    #stressField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.ELEMENT_BASED)
    #stressField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.ELEMENT_BASED)
    #stressField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 3, iron.FieldInterpolationTypes.ELEMENT_BASED)
    #stressField.VariableLabelSet(iron.FieldVariableTypes.U, "Cauchy Stress")
    #stressField.CreateFinish()

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

    print("Material field number of components (U,V):", materialField.NumberOfComponentsGet(iron.FieldVariableTypes.U), materialField.NumberOfComponentsGet(iron.FieldVariableTypes.V))

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

    print("Dependent field number of components (U,DELUDELN):", dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.U), dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.DELUDELN))

    if gravity != [0.0, 0.0, 0.0]:
        # Setup gravity source field
        sourceField = iron.Field()
        equationsSet.SourceCreateStart(sourceFieldUserNumber, sourceField)
        sourceField.fieldScalingType = iron.FieldScalingTypes.UNIT
        equationsSet.SourceCreateFinish()

        #Set the gravity vector component values
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
    #nonLinearSolver.NewtonAbsoluteToleranceSet(1e-11)
    #nonLinearSolver.NewtonSolutionToleranceSet(1e-11)
    #nonLinearSolver.NewtonRelativeToleranceSet(1e-11)
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
    equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    # Prescribe boundary conditions (absolute nodal parameters)
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    # N = 8  # number of nodes
    # avgRadius = 4  # in voxels
    # for nid in range(1,N+1):
    #     X = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 1)
    #     Y = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 2)
    #     Z = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 3)
    #
    #     dx = 0
    #     dy = 0
    #     dz = 0
    #     n = 0
    #     for k in range(max(0, int(np.floor(Z-avgRadius))), min(dim[2], int(np.ceil(Z+avgRadius)+1))):
    #         for j in range(max(0, int(np.floor(Y-avgRadius))), min(dim[1], int(np.ceil(Y+avgRadius)+1))):
    #             for i in range(max(0, int(np.floor(X-avgRadius))), min(dim[0], int(np.ceil(X+avgRadius)+1))):
    #                 if (i-X)**2 + (j-Y)**2 + (k-Z)**2 <= avgRadius**2:
    #                     dx += image_dx[i,j,k]
    #                     dy += image_dy[i,j,k]
    #                     dz += image_dz[i,j,k]
    #                     n += 1
    #     dx /= n
    #     dy /= n
    #     dz /= n
    #     x = X + dx
    #     y = Y + dy
    #     z = Z + dz
    #
    #     #dx = 0.0
    #     #if X == 1.0:
    #         #dx = 1.0
    #         #boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
    #                                #iron.BoundaryConditionsTypes.FIXED, dx)
    #     #if nid == 1:
    #     boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
    #                                iron.BoundaryConditionsTypes.FIXED, dx)
    #     boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
    #                                iron.BoundaryConditionsTypes.FIXED, dy)
    #     boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
    #                                iron.BoundaryConditionsTypes.FIXED, dz)
    #
    #     print(nid, X,Y,Z, '=>', x,y,z)

    nodes = iron.Nodes()
    region.NodesGet(nodes)
    print('Number of nodes:', nodes.numberOfNodes)

    for nid in range(1,nodes.numberOfNodes+1):
        X = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 1)
        Y = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 2)
        Z = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nid, 3)

        if X == 0.0:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
                                       iron.BoundaryConditionsTypes.FIXED, 0.0)
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
                                       iron.BoundaryConditionsTypes.FIXED, 0.0)
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
                                       iron.BoundaryConditionsTypes.FIXED, 0.0)
        elif X == lengthX:
            pass
            #boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
            #                           iron.BoundaryConditionsTypes.FIXED, lengthX*(stretch-1.0))
            #boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
            #                           iron.BoundaryConditionsTypes.FIXED, 0.0)
            #boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
            #                           iron.BoundaryConditionsTypes.FIXED, 0.0)

        elif Y == 0.0 or Z == 0.0 or Y == lengthY or Z == lengthZ:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
                                       iron.BoundaryConditionsTypes.FIXED, 0.0)
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
                                       iron.BoundaryConditionsTypes.FIXED, 0.0)

        else:
            print('Node unconstrained:', nid, X,Y,Z)

    solverEquations.BoundaryConditionsCreateFinish()

    ##################################################################
    # Solve!
    ##################################################################
    problem.Solve()

    #gaussPointNumber = 1
    #userElementNumber = 1
    #cauchy = equationsSet.TensorInterpolateGaussPoint(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS, gaussPointNumber,
    #                                         userElementNumber, valuesSizes)
    #print(cauchy)

    ssFile = open('./results/ss.exdata', 'w')
    ssFile.write(' Group name: Strain-Stress\n')
    ssFile.write(' #Fields=3\n')
    ssFile.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    ssFile.write('   x.  Value index= 1, #Derivatives=0\n')
    ssFile.write('   y.  Value index= 2, #Derivatives=0\n')
    ssFile.write('   z.  Value index= 3, #Derivatives=0\n')
    ssFile.write(' 2) Strain, field, rectangular cartesian, #Components=1\n')
    ssFile.write('   1.  Value index= 4, #Derivatives=0\n')
    ssFile.write(' 3) Stress, field, rectangular cartesian, #Components=1\n')
    ssFile.write('   1.  Value index= 5, #Derivatives=0\n')

    #
    # Expected values for strain tensor are:
    # 0.22 0 0
    # 0    0 0
    # 0    0 0
    #
    # Expected values for stress tensor are:
    # 24000     0     0
    #     0 18867     0
    #     0     0 18867
    #

    valuesSizes = (3,3)
    ssInterp = 2
    for eid in [0]:
        for xiZ in np.linspace(0.0, 1.0, ssInterp):
            for xiY in np.linspace(0.0, 1.0, ssInterp):
                for xiX in np.linspace(0.0, 1.0, ssInterp):
                    coords = dependentField.ParameterSetInterpolateSingleXiDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                                              1, eid+1, [xiX, xiY, xiZ], 3)
                    strain = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.GREEN_LAGRANGE_STRAIN, eid,
                                                              [xiX, xiY, xiZ], valuesSizes)
                    stress = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS, eid,
                                                              [xiX, xiY, xiZ], valuesSizes)

                    print(coords)
                    print(strain)
                    print(stress)
                    print()

        #print(eid, cauchy)

    #stressField.ParameterSetAddElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
    #                                     userElementNumber, 1, 9.0)

    fields = iron.Fields()
    fields.CreateRegion(region)
    return fields


if not os.path.exists("./results"):
    os.makedirs("./results")

fields = Run(exelemFile, exnodeFile, interpolation, displacement, displacementPixdim, c10, c01, k, density, gravity)
fields.NodesExport("./results/out", "FORTRAN")
fields.ElementsExport("./results/out", "FORTRAN")
fields.Finalise()
