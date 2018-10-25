#!/usr/bin/env python

from medpy.io import load
from medpy.io import header
import numpy as np
from opencmiss.iron import iron
from morphic.utils import convert_hermite_lagrange
import mesh_tools
from enum import Enum

class FittingVariableTypes(Enum):
    CAUCHY_STRESS = 1
    GREEN_LAGRANGE_STRAIN = 2
    DEFORMATION_GRADIENT = 3
    PRINCIPAL_STRESS = 4
    PRINCIPAL_STRAIN = 5
    PRINCIPAL_DEFORMATION = 6
    AVERAGE_STRESS = 7
    AVERAGE_STRAIN = 8
    JACOBIAN = 9

class LungInflation:
    def __init__(self):
        self.c10 = 5000.0  # in Pa
        self.c01 = 2000.0  # in Pa
        self.k = 6000.0  # in Pa
        self.density = 1000.0  # in kg m^-3
        self.gravity = [0.0, 0.0, 0.0]  # in m s^-2

        self.interpolation = iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE

        self.transformLungCoordinates = False
        self.mesh = None
        self.coordinates = None
        self.node_nums = None
        self.element_nums = None
        self.displacement = None
        self.pixdim = None

        self.decomposition = None
        self.geometricField = None
        self.fibreField = None
        self.dependentField = None
        self.equationsSet = None

        self.fittingUserNumberCounter = 20
        self.fittingFinalisers = []

        coordinateSystemUserNumber = 1
        regionUserNumber = 1
        basisUserNumber = 1

        # Create a 3D rectangular cartesian coordinate system
        self.coordinateSystem = iron.CoordinateSystem()
        self.coordinateSystem.CreateStart(coordinateSystemUserNumber)
        self.coordinateSystem.DimensionSet(3)
        self.coordinateSystem.CreateFinish()

        # Create a region and assign the coordinate system to the region
        self.region = iron.Region()
        self.region.CreateStart(regionUserNumber, iron.WorldRegion)
        self.region.LabelSet("Region")
        self.region.coordinateSystem = self.coordinateSystem
        self.region.CreateFinish()

        # Define basis
        self.basis = iron.Basis()
        self.basis.CreateStart(basisUserNumber)
        self.basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        self.basis.numberOfXi = 3
        self.basis.interpolationXi = [self.interpolation] * 3
        if self.interpolation == iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE:
            self.basis.quadratureNumberOfGaussXi = [2] * 3
        elif self.interpolation == iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE:
            self.basis.quadratureNumberOfGaussXi = [3] * 3
        elif self.interpolation == iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE:
            self.basis.quadratureNumberOfGaussXi = [4] * 3
        elif self.interpolation == iron.BasisInterpolationSpecifications.CUBIC_HERMITE:
            self.basis.quadratureNumberOfGaussXi = [4] * 3
        self.basis.CreateFinish()

    def SetMaterialParameters(self, c10, c01, k, density):
        self.c10 = c10
        self.c01 = c01
        self.k = k
        self.density = density

    def SetGravity(self, zForce):
        self.gravity = [0.0, 0.0, -zForce]  # in m s^-2

    def LoadHermiteMesh(self, exelemFilename, exnodeFilename, coordinatesFieldName, transformLungCoordinates=False):
        meshUserNumber = 1

        self.transformLungCoordinates = transformLungCoordinates
        cubic_hermite_morphic_mesh = mesh_tools.exfile_to_morphic(exnodeFilename, exelemFilename, coordinatesFieldName,
                                                                  dimension=3, interpolation='hermite')
        cubic_lagrange_morphic_mesh = convert_hermite_lagrange(cubic_hermite_morphic_mesh, tol=1e-9)
        self.mesh, self.coordinates, self.node_nums, self.element_nums = mesh_tools.morphic_to_OpenCMISS(
            cubic_lagrange_morphic_mesh,
            self.region, self.basis,
            meshUserNumber,
            dimension=3,
            interpolation='cubic')

    def LoadDisplacements(self, dxFilename, dyFilename, dzFilename):
        print('Load dx...')
        dx, image_header = load(dxFilename)
        self.displacement = np.zeros((dx.shape[0], dx.shape[1], dx.shape[2], 3))
        self.displacement[:, :, :, 0] = dx
        print('Load dy...')
        self.displacement[:, :, :, 1], _ = load(dyFilename)
        print('Load dz...')
        self.displacement[:, :, :, 2], _ = load(dzFilename)
        self.pixdim = header.get_pixel_spacing(image_header)
        print('Image dimensions:', self.displacement.shape)
        print('Voxel dimensions:', self.pixdim)

    def TrilinearInterpolation(self, X, Y, Z):
        # lung mesh coordinates (in mm) have:
        #  x increasing from right to left
        #  y increasing from ventral to dorsal
        #  z increasing from caudal to cradal
        # pixel coordinates (in px) have:
        #  x increasing from right to left
        #  y increasing from dorsal to ventral
        #  z increasing from cradal to caudal
        # where y has been shifted (by half its size in pixels?)
        if self.transformLungCoordinates:
            X = X
            Y = 256 - Y
            Z = -Z
        x = X / self.pixdim[0]
        y = Y / self.pixdim[1]
        z = Z / self.pixdim[2]

        x0 = np.floor(x)
        y0 = np.floor(y)
        z0 = np.floor(z)
        x1 = np.ceil(x)
        y1 = np.ceil(y)
        z1 = np.ceil(z)
        if x0 < 0 or y0 < 0 or z0 < 0 or x1 >= self.displacement.shape[0] or y1 >= self.displacement.shape[1] or z1 >= \
                self.displacement.shape[2]:
            print(x0, y0, z0, x1, y1, z1)
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
        c00 = self.displacement[x0, y0, z0] * (1 - xd) + self.displacement[x1, y0, z0] * xd
        c01 = self.displacement[x0, y0, z1] * (1 - xd) + self.displacement[x1, y0, z1] * xd
        c10 = self.displacement[x0, y1, z0] * (1 - xd) + self.displacement[x1, y1, z0] * xd
        c11 = self.displacement[x0, y1, z1] * (1 - xd) + self.displacement[x1, y1, z1] * xd
        c0 = c00 * (1 - yd) + c10 * yd
        c1 = c01 * (1 - yd) + c11 * yd
        (dx, dy, dz) = c0 * (1 - zd) + c1 * zd

        # it appears that displacement data wth RAI orientation gives:
        #  x increasing from left to right
        #  y increasing from dorsal to ventral
        #  z increasing from cradal to caudal
        # but our mesh coordinate system has all of these reversed
        if self.transformLungCoordinates:
            dx = -dx
            dy = -dy
            dz = -dz
        return dx, dy, dz

    def Setup(self):
        decompositionUserNumber = 1
        geometricFieldUserNumber = 1
        fibreFieldUserNumber = 2
        materialFieldUserNumber = 3
        dependentFieldUserNumber = 4
        sourceFieldUserNumber = 5
        equationsSetFieldUserNumber = 6
        equationsSetUserNumber = 1

        # Create a decomposition for the mesh
        numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
        self.decomposition = iron.Decomposition()
        self.decomposition.CreateStart(decompositionUserNumber, self.mesh)
        self.decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
        self.decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
        self.decomposition.CreateFinish()

        # Create a field for the geometry
        self.geometricField = iron.Field()
        self.geometricField.CreateStart(geometricFieldUserNumber, self.region)
        self.geometricField.MeshDecompositionSet(self.decomposition)
        self.geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        self.geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
        self.geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
        self.geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
        self.geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
        self.geometricField.CreateFinish()

        # Update the geometric field parameters
        self.geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
        for node_idx, node in enumerate(self.node_nums):
            for component_idx, component in enumerate([1, 2, 3]):
                for derivative_idx, derivative in enumerate(range(1, self.coordinates.shape[2] + 1)):
                    self.geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES,
                                                                 1, derivative, node, component,
                                                                 self.coordinates[node_idx, component_idx, derivative_idx])

        self.geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES)

        # Create a fibre field and attach it to the geometric field
        self.fibreField = iron.Field()
        self.fibreField.CreateStart(fibreFieldUserNumber, self.region)
        self.fibreField.TypeSet(iron.FieldTypes.FIBRE)
        self.fibreField.MeshDecompositionSet(self.decomposition)
        self.fibreField.GeometricFieldSet(self.geometricField)
        self.fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
        self.fibreField.CreateFinish()

        ##################################################################
        # Setup Mooney-Rivlin equations
        ##################################################################

        equationsSetField = iron.Field()
        self.equationsSet = iron.EquationsSet()
        equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                     iron.EquationsSetTypes.FINITE_ELASTICITY,
                                     iron.EquationsSetSubtypes.COMPRESSIBLE_FINITE_ELASTICITY]
        self.equationsSet.CreateStart(equationsSetUserNumber, self.region, self.fibreField, equationsSetSpecification,
                                      equationsSetFieldUserNumber, equationsSetField)
        self.equationsSet.CreateFinish()

        # Setup material field
        materialField = iron.Field()
        self.equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
        materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
        materialField.VariableLabelSet(iron.FieldVariableTypes.V, "Density")
        self.equationsSet.MaterialsCreateFinish()

        print("Material field number of components (U,V):",
              materialField.NumberOfComponentsGet(iron.FieldVariableTypes.U),
              materialField.NumberOfComponentsGet(iron.FieldVariableTypes.V))

        # Set Mooney-Rivlin constants c10 and c01 respectively.
        materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, self.c10)
        materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, self.c01)
        materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, self.k)
        materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.V, iron.FieldParameterSetTypes.VALUES, 1, self.density)

        # Setup dependent field
        self.dependentField = iron.Field()
        self.equationsSet.DependentCreateStart(dependentFieldUserNumber, self.dependentField)
        self.dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
        self.dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN, "Dependent_delU_delN")
        self.equationsSet.DependentCreateFinish()

        print("Dependent field number of components (U,DELUDELN):",
              self.dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.U),
              self.dependentField.NumberOfComponentsGet(iron.FieldVariableTypes.DELUDELN))

        if self.gravity != [0.0, 0.0, 0.0]:
            # Setup gravity source field
            sourceField = iron.Field()
            self.equationsSet.SourceCreateStart(sourceFieldUserNumber, sourceField)
            sourceField.fieldScalingType = iron.FieldScalingTypes.UNIT
            self.equationsSet.SourceCreateFinish()

            # Set the gravity vector component values
            sourceField.ComponentValuesInitialiseDP(
                iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, self.gravity[0])
            sourceField.ComponentValuesInitialiseDP(
                iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, self.gravity[1])
            sourceField.ComponentValuesInitialiseDP(
                iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, self.gravity[2])

        # Create equations
        equations = iron.Equations()
        self.equationsSet.EquationsCreateStart(equations)
        equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
        equations.outputType = iron.EquationsOutputTypes.NONE
        self.equationsSet.EquationsCreateFinish()

        # Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
        iron.Field.ParametersToFieldParametersComponentCopy(
            self.geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
            self.dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1)
        iron.Field.ParametersToFieldParametersComponentCopy(
            self.geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
            self.dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2)
        iron.Field.ParametersToFieldParametersComponentCopy(
            self.geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3,
            self.dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3)

    def Solve(self, outputFilename, timeSteps):
        problemUserNumber = 1

        self.FitField("CauchyStress", FittingVariableTypes.CAUCHY_STRESS)
        self.FitField("GreenLagrangeStrain", FittingVariableTypes.GREEN_LAGRANGE_STRAIN)
        self.FitField("AverageCauchyStress", FittingVariableTypes.AVERAGE_STRESS)
        self.FitField("AverageGreenLagrangeStrain", FittingVariableTypes.AVERAGE_STRAIN)
        self.FitField("Jacobian", FittingVariableTypes.JACOBIAN)

        fields = iron.Fields()
        fields.CreateRegion(self.region)
        fields.NodesExport(outputFilename + "_0", "FORTRAN")
        fields.ElementsExport(outputFilename, "FORTRAN")
        fields.Finalise()

        for timeStep in range(1, timeSteps + 1):
            loadRatio = timeStep / timeSteps

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
            linearSolver = iron.Solver()
            nonLinearSolver = iron.Solver()
            problem.SolversCreateStart()
            problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
            nonLinearSolver.outputType = iron.SolverOutputTypes.MONITOR
            nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
            nonLinearSolver.NewtonLinearSolverGet(linearSolver)
            nonLinearSolver.NewtonAbsoluteToleranceSet(1e-9)
            nonLinearSolver.NewtonSolutionToleranceSet(1e-9)
            nonLinearSolver.NewtonRelativeToleranceSet(1e-9)
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
            solverEquations.EquationsSetAdd(self.equationsSet)
            problem.SolverEquationsCreateFinish()

            # Prescribe boundary conditions (absolute nodal parameters)
            boundaryConditions = iron.BoundaryConditions()
            solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

            nodes = iron.MeshNodes()
            self.mesh.NodesGet(1, nodes)
            for nid in self.node_nums:
                if nodes.NodeOnBoundaryGet(nid) == iron.MeshBoundaryTypes.ON:
                    X = self.geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                  iron.FieldParameterSetTypes.VALUES, 1,
                                                                  1,
                                                                  nid, 1)
                    Y = self.geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                  iron.FieldParameterSetTypes.VALUES, 1,
                                                                  1,
                                                                  nid, 2)
                    Z = self.geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                  iron.FieldParameterSetTypes.VALUES, 1,
                                                                  1,
                                                                  nid, 3)

                    (dx, dy, dz) = self.TrilinearInterpolation(X, Y, Z)

                    boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 1,
                                               iron.BoundaryConditionsTypes.FIXED, dx * loadRatio)
                    boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 2,
                                               iron.BoundaryConditionsTypes.FIXED, dy * loadRatio)
                    boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, nid, 3,
                                               iron.BoundaryConditionsTypes.FIXED, dz * loadRatio)

                    print('BC set on node %d at %.2f %.2f %.2f += %f %f %f' % (nid, X, Y, Z, dx, dy, dz))

            solverEquations.BoundaryConditionsCreateFinish()

            ##################################################################
            # Solve!
            ##################################################################
            problem.Solve()

            self.FitField("CauchyStress", FittingVariableTypes.CAUCHY_STRESS)
            self.FitField("GreenLagrangeStrain", FittingVariableTypes.GREEN_LAGRANGE_STRAIN)
            self.FitField("AverageCauchyStress", FittingVariableTypes.AVERAGE_STRESS)
            self.FitField("AverageGreenLagrangeStrain", FittingVariableTypes.AVERAGE_STRAIN)
            self.FitField("Jacobian", FittingVariableTypes.JACOBIAN)

            fields = iron.Fields()
            fields.CreateRegion(self.region)
            fields.NodesExport(outputFilename + "_%d" % timeStep, "FORTRAN")
            fields.Finalise()

            solverEquations.Finalise()
            problem.Finalise()

        for finaliser in self.fittingFinalisers:
            finaliser.Finalise()
        self.basis.Finalise()
        self.region.Finalise()
        self.coordinateSystem.Finalise()

    def FitField(self, name, variable, tau=0.01, kappa=0.0005):
        fittingFieldUserNumber = self.fittingUserNumberCounter + 0
        fittingEquationsSetUserNumber = self.fittingUserNumberCounter + 1
        fittingEquationsSetFieldUserNumber = self.fittingUserNumberCounter + 2
        fittingDependentFieldUserNumber = self.fittingUserNumberCounter + 3
        fittingMaterialsFieldUserNumber = self.fittingUserNumberCounter + 4
        fittingProblemUserNumber = self.fittingUserNumberCounter + 5
        self.fittingUserNumberCounter += 6

        numComponents = 1
        if variable in [FittingVariableTypes.CAUCHY_STRESS, FittingVariableTypes.GREEN_LAGRANGE_STRAIN,
                        FittingVariableTypes.DEFORMATION_GRADIENT]:
            numComponents = 6
        elif variable in [FittingVariableTypes.PRINCIPAL_STRESS, FittingVariableTypes.PRINCIPAL_STRAIN,
                          FittingVariableTypes.PRINCIPAL_DEFORMATION]:
            numComponents = 3

        fittingField = iron.Field()
        fittingField.CreateStart(fittingFieldUserNumber, self.region)
        fittingField.TypeSet(iron.FieldTypes.GENERAL)
        fittingField.MeshDecompositionSet(self.decomposition)
        fittingField.GeometricFieldSet(self.geometricField)
        fittingField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
        fittingField.NumberOfVariablesSet(2)
        fittingField.VariableTypesSet([iron.FieldVariableTypes.U, iron.FieldVariableTypes.V])
        fittingField.VariableLabelSet(iron.FieldVariableTypes.U, name + "_GaussValue")
        fittingField.VariableLabelSet(iron.FieldVariableTypes.V, name + "_GaussWeight")
        fittingField.NumberOfComponentsSet(iron.FieldVariableTypes.U, numComponents)
        fittingField.NumberOfComponentsSet(iron.FieldVariableTypes.V, numComponents)
        for component in range(1, numComponents + 1):
            fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, component, 1)
            fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, component, 1)
            fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                                   iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
            fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, component,
                                                   iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
        fittingField.CreateFinish()

        # Update the geometric field parameters
        numElements = self.mesh.NumberOfElementsGet()
        numGaussPoints = int(np.prod(self.basis.QuadratureNumberOfGaussXiGet(3)))
        fittingField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
        for eid in range(0, numElements):
            for gid in range(1, numGaussPoints + 1):
                xi = self.basis.QuadratureSingleGaussXiGet(iron.BasisQuadratureSchemes.DEFAULT, gid, 3)
                if variable in [FittingVariableTypes.CAUCHY_STRESS, FittingVariableTypes.PRINCIPAL_STRESS,
                                FittingVariableTypes.AVERAGE_STRESS]:
                    tensor = self.equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,
                                                                   eid, xi, (3, 3))
                elif variable in [FittingVariableTypes.GREEN_LAGRANGE_STRAIN, FittingVariableTypes.PRINCIPAL_STRAIN,
                                  FittingVariableTypes.AVERAGE_STRAIN]:
                    tensor = self.equationsSet.TensorInterpolateXi(
                        iron.EquationsSetDerivedTensorTypes.GREEN_LAGRANGE_STRAIN, eid, xi, (3, 3))
                elif variable in [FittingVariableTypes.DEFORMATION_GRADIENT, FittingVariableTypes.PRINCIPAL_DEFORMATION,
                                  FittingVariableTypes.JACOBIAN]:
                    tensor = self.equationsSet.TensorInterpolateXi(
                        iron.EquationsSetDerivedTensorTypes.DEFORMATION_GRADIENT, eid, xi, (3, 3))
                else:
                    raise ValueError("Fitting field variable doesn't exist")

                if numComponents == 6:
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 1,
                                                                tensor[0, 0])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 2,
                                                                tensor[0, 1])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 3,
                                                                tensor[0, 2])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 4,
                                                                tensor[1, 1])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 5,
                                                                tensor[1, 2])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 6,
                                                                tensor[2, 2])
                elif numComponents == 3:
                    eigs = np.linalg.eigvals(tensor)
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 1,
                                                                eigs[0])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 2,
                                                                eigs[1])
                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 3,
                                                                eigs[2])
                elif numComponents == 1:
                    eigs = np.linalg.eigvals(tensor)
                    if FittingVariableTypes.JACOBIAN:
                        value = np.prod(eigs)
                    else:
                        value = np.mean(eigs)

                    if isinstance(value, complex):
                        if abs(value.imag) > 1e-9:
                            print('WARNING: tensor interpolation in element %d gives complex number' % eid)
                        value = value.real

                    fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, gid, eid, 1, value)

        # Initialise Gauss point weight field to 1.0
        for component in range(1, numComponents + 1):
            fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V, iron.FieldParameterSetTypes.VALUES,
                                                     component, 1.0)

        # Create the fitting equations set
        fittingEquationsSetField = iron.Field()
        fittingEquationsSet = iron.EquationsSet()
        fittingEquationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                                            iron.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                                            iron.EquationsSetSubtypes.GAUSS_POINT_FITTING,
                                            iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
        fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber, self.region, self.geometricField,
                                        fittingEquationsSetSpecification, fittingEquationsSetFieldUserNumber,
                                        fittingEquationsSetField)
        fittingEquationsSet.CreateFinish()

        # Create the fitting dependent field
        fittingDependentField = iron.Field()
        fittingEquationsSet.DependentCreateStart(fittingDependentFieldUserNumber, fittingDependentField)
        fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.U, name)
        fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN, name + "_delU_delN")
        # Set the number of components to 2
        fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U, numComponents)
        fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, numComponents)
        # Set the field variables to be triquadratic Lagrange
        for component in range(1, numComponents + 1):
            fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, component, 1)
            fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, component, 1)
        fittingEquationsSet.DependentCreateFinish()

        # Create the fitting independent field
        fittingEquationsSet.IndependentCreateStart(fittingFieldUserNumber, fittingField)
        fittingEquationsSet.IndependentCreateFinish()

        # Create material field (Sobolev parameters)
        fittingMaterialField = iron.Field()
        fittingEquationsSet.MaterialsCreateStart(fittingMaterialsFieldUserNumber, fittingMaterialField)
        fittingMaterialField.VariableLabelSet(iron.FieldVariableTypes.U, name + "_SmoothingParameters")
        fittingEquationsSet.MaterialsCreateFinish()
        fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                         1, tau)
        fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                         2, kappa)

        # Create the fitting equations
        fittingEquations = iron.Equations()
        fittingEquationsSet.EquationsCreateStart(fittingEquations)
        fittingEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
        fittingEquations.outputType = iron.EquationsOutputTypes.NONE
        fittingEquationsSet.EquationsCreateFinish()

        # Create fitting problem
        fittingProblem = iron.Problem()
        fittingProblemSpecification = [iron.ProblemClasses.FITTING,
                                       iron.ProblemTypes.DATA_FITTING,
                                       iron.ProblemSubtypes.STATIC_FITTING]
        fittingProblem.CreateStart(fittingProblemUserNumber, fittingProblemSpecification)
        fittingProblem.CreateFinish()

        # Create control loops
        fittingProblem.ControlLoopCreateStart()
        fittingProblem.ControlLoopCreateFinish()

        # Create problem solver
        fittingSolver = iron.Solver()
        fittingProblem.SolversCreateStart()
        fittingProblem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, fittingSolver)
        fittingSolver.outputType = iron.SolverOutputTypes.MONITOR
        fittingProblem.SolversCreateFinish()

        # Create fitting solver equations and add fitting equations set to solver equations
        fittingSolverEquations = iron.SolverEquations()
        fittingProblem.SolverEquationsCreateStart()
        # Get the solver equations
        fittingSolver.SolverEquationsGet(fittingSolverEquations)
        fittingSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
        fittingSolverEquations.EquationsSetAdd(fittingEquationsSet)
        fittingProblem.SolverEquationsCreateFinish()

        # Prescribe boundary conditions for the fitting problem
        fittingBoundaryConditions = iron.BoundaryConditions()
        fittingSolverEquations.BoundaryConditionsCreateStart(fittingBoundaryConditions)
        fittingSolverEquations.BoundaryConditionsCreateFinish()

        # Solve the fitting problem
        fittingProblem.Solve()

        self.fittingFinalisers.append(fittingProblem)
        self.fittingFinalisers.append(fittingMaterialField)
        self.fittingFinalisers.append(fittingDependentField)
        self.fittingFinalisers.append(fittingEquationsSet)
        self.fittingFinalisers.append(fittingField)
