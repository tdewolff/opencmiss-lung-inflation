from opencmiss.iron import iron
import numpy as np
from enum import Enum

globalFittingUserNumber = 100


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


def FitField(name, variable, basis, region, mesh, decomposition, geometricField, equationsSet, tau=0.01, kappa=0.0005):
    global globalFittingUserNumber
    fittingFieldUserNumber = globalFittingUserNumber + 0
    fittingEquationsSetUserNumber = globalFittingUserNumber + 1
    fittingEquationsSetFieldUserNumber = globalFittingUserNumber + 2
    fittingDependentFieldUserNumber = globalFittingUserNumber + 3
    fittingMaterialsFieldUserNumber = globalFittingUserNumber + 4
    fittingProblemUserNumber = globalFittingUserNumber + 5
    globalFittingUserNumber += 10

    numComponents = 1
    if variable in [FittingVariableTypes.CAUCHY_STRESS, FittingVariableTypes.GREEN_LAGRANGE_STRAIN,
                    FittingVariableTypes.DEFORMATION_GRADIENT]:
        numComponents = 6
    elif variable in [FittingVariableTypes.PRINCIPAL_STRESS, FittingVariableTypes.PRINCIPAL_STRAIN,
                    FittingVariableTypes.PRINCIPAL_DEFORMATION]:
        numComponents = 3

    fittingField = iron.Field()
    fittingField.CreateStart(fittingFieldUserNumber, region)
    fittingField.TypeSet(iron.FieldTypes.GENERAL)
    fittingField.MeshDecompositionSet(decomposition)
    fittingField.GeometricFieldSet(geometricField)
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
    numElements = mesh.NumberOfElementsGet()
    numGaussPoints = np.prod(basis.QuadratureNumberOfGaussXiGet(3))
    fittingField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    for eid in range(0, numElements):
        for gid in range(1, numGaussPoints + 1):
            xi = basis.QuadratureSingleGaussXiGet(iron.BasisQuadratureSchemes.DEFAULT, gid, 3)
            if variable in [FittingVariableTypes.CAUCHY_STRESS, FittingVariableTypes.PRINCIPAL_STRESS,
                            FittingVariableTypes.AVERAGE_STRESS]:
                tensor = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS, eid, xi,
                                                          (3, 3))
            elif variable in [FittingVariableTypes.GREEN_LAGRANGE_STRAIN, FittingVariableTypes.PRINCIPAL_STRAIN,
                              FittingVariableTypes.AVERAGE_STRAIN]:
                tensor = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.GREEN_LAGRANGE_STRAIN,
                                                          eid, xi, (3, 3))
            elif variable in [FittingVariableTypes.DEFORMATION_GRADIENT, FittingVariableTypes.PRINCIPAL_DEFORMATION,
                              FittingVariableTypes.JACOBIAN]:
                tensor = equationsSet.TensorInterpolateXi(iron.EquationsSetDerivedTensorTypes.DEFORMATION_GRADIENT, eid,
                                                          xi, (3, 3))

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
                                                            iron.FieldParameterSetTypes.VALUES, gid, eid, 1, eigs[0])
                fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, gid, eid, 2, eigs[1])
                fittingField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, gid, eid, 3, eigs[2])
            elif numComponents == 1:
                eigs = np.linalg.eigvals(tensor)
                if FittingVariableTypes.JACOBIAN:
                    value = np.prod(eigs)
                else:
                    value = np.mean(eigs)
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
    fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber, region, geometricField,
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
