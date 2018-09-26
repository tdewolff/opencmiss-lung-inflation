#!/usr/bin/env python

import os
from opencmiss.iron import iron
import mesh_io
from medpy.io import load
from medpy.io import header
import numpy as np

path = '/hpc/tdew803/Lung/Data/Pig_PE_Study_HRC/'
subject = 'AP00149'
pressure = '5cmH2O'
registration = '5_10'

c01 = 1.0
c10 = 0.0
c11 = 1.0e9
density = 9.0e-4  # in g mm^-3
gravity = [0.0, 0.0, 9.81]  # in m s^-2

height = 1.0
width = 1.0
length = 1.0

NumberOfGaussXi = 2
numberGlobalXElements = 2
numberGlobalYElements = 2
numberGlobalZElements = 2
numberOfXi = 3

##################################################################

#exelemFile = path + subject + '/' + pressure + '/Lung/FEMesh/Left_Refitted.exelem'
#exnodeFile = path + subject + '/' + pressure + '/Lung/FEMesh/Left_Refitted.exnode'
deformationFiles = path + subject + '/reg_' + registration + '_d%s.nii'

print('Load dx...')
image_dx, image_header = load(deformationFiles % ('x'))
print('Load dy...')
image_dy, _ = load(deformationFiles % ('y'))
print('Load dz...')
image_dz, _ = load(deformationFiles % ('z'))

dim = image_dx.shape
pixdim = header.get_pixel_spacing(image_header)

print('Image dimensions:', dim)
print('Voxel dimensions:', pixdim)

####

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Set all diagnostic levels on for testing
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL, [1, 2, 3, 4, 5], "Diagnostics",
#                      ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

##################################################################

coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
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
basis.numberOfXi = numberOfXi
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numberOfXi
if (NumberOfGaussXi > 0):
    basis.quadratureNumberOfGaussXi = [NumberOfGaussXi] * numberOfXi
basis.CreateFinish()

#mesh, coordinates, _ = mesh_io.exfile_to_OpenCMISS('/hpc/tdew803/Lung/Data/Pig_PE_Study_HRC/AP00149/5cmH2O/Lung/FEMesh/Left_Refitted.exnode',
#                            '/hpc/tdew803/Lung/Data/Pig_PE_Study_HRC/AP00149/5cmH2O/Lung/FEMesh/Left_Refitted.exelem',
#                            coordinateSystem, region, basis, meshUserNumber)

mesh = iron.Mesh()
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width, length, height]
generatedMesh.numberOfElements = [numberGlobalXElements, numberGlobalYElements, numberGlobalZElements]
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
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
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
fibreField.CreateFinish()


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
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, c11)
materialField.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.V, iron.FieldParameterSetTypes.VALUES, 1, density)

# Setup dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
equationsSet.DependentCreateFinish()

print("Dependent field number of components (U,DELUDELN):", dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.U), dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.DELUDELN))

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
nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
nonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
nonLinearSolver.NewtonSolutionToleranceSet(1e-14)
nonLinearSolver.NewtonRelativeToleranceSet(1e-14)
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

for nid in [1,7,19,25, 13,22,16,4,10]:
    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)

for nid in [3,9,21,27, 15,18,6,12,24]:
    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
                               iron.BoundaryConditionsTypes.FIXED, 0.2)
    #boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
    #                           iron.BoundaryConditionsTypes.FIXED, 0.0)
    #boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
    #                           iron.BoundaryConditionsTypes.FIXED, 0.0)

#for nid in [26,17,8,5,2,11,20,23]:
#    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
#                               iron.BoundaryConditionsTypes.FIXED, 0.0)
#    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
#                               iron.BoundaryConditionsTypes.FIXED, 0.0)

solverEquations.BoundaryConditionsCreateFinish()


##################################################################
# Solve!
##################################################################
problem.Solve()

if not os.path.exists("./results"):
    os.makedirs("./results")

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/out", "FORTRAN")
fields.ElementsExport("./results/out", "FORTRAN")
fields.Finalise()
