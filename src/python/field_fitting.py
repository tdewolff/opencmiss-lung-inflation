from opencmiss.iron import iron


def FitField(region, decomposition, fibreField, variable, tau, kappa):
    fittedFieldUserNumber = 100
    dependentFieldUserNumber = 101
    equationsSetFieldUserNumber = 102
    equationsSetUserNumber = 103
    materialsFieldUserNumber = 104
    problemUserNumber = 105

    fittedField = iron.Field()
    fittedField.CreateStart(fittedFieldUserNumber, region)
    fittedField.TypeSet(iron.FieldTypes.GENERAL)
    fittedField.MeshDecompositionSet(decomposition)
    fittedField.GeometricFieldSet(fibreField)  # geometricField? Or leave out?
    fittedField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
    fittedField.NumberOfVariablesSet(2)
    fittedField.VariableTypesSet([iron.FieldVariableTypes.U, iron.FieldVariableTypes.V])
    fittedField.VariableLabelSet(iron.FieldVariableTypes.U, "GaussStress")
    fittedField.VariableLabelSet(iron.FieldVariableTypes.V, "GaussWeight")
    fittedField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 6)
    fittedField.NumberOfComponentsSet(iron.FieldVariableTypes.V, 6)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 5, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 6, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 1, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 2, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 3, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 4, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 5, 1)
    fittedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 6, 1)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 3,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 5,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 6,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 1,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 2,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 3,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 4,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 5,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 6,
                                          iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittedField.CreateFinish()

    # Initialise Gauss point weight field to 1.0
    fittedField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                            iron.FieldParameterSetTypes.VALUES, 1, 1.0)
    fittedField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                            iron.FieldParameterSetTypes.VALUES, 2, 1.0)
    fittedField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                            iron.FieldParameterSetTypes.VALUES, 3, 1.0)
    fittedField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                            iron.FieldParameterSetTypes.VALUES, 4, 1.0)
    fittedField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                            iron.FieldParameterSetTypes.VALUES, 5, 1.0)
    fittedField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                            iron.FieldParameterSetTypes.VALUES, 6, 1.0)

    equationsSetField = iron.Field()
    equationsSet = iron.EquationsSet()
    equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                                 iron.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                                 iron.EquationsSetSubtypes.GAUSS_POINT_FITTING,
                                 iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
    equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, equationsSetSpecification,
                             equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()

    # self.equationsSet.DerivedCreateStart(fittedFieldUserNumber, fittedField)
    # self.equationsSet.DerivedVariableSet(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,
    #                                     iron.FieldVariableTypes.U)
    # self.equationsSet.DerivedCreateFinish()

    dependentField = iron.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "FittedStress")
    # Set the number of components to 2
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 6)
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 6)
    # Set the field variables to be triquadratic Lagrange
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 5, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 6, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 1, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 2, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 3, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 5, 1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 6, 1)
    equationsSet.DependentCreateFinish()

    equationsSet.IndependentCreateStart(fittedFieldUserNumber, fittedField)
    equationsSet.IndependentCreateFinish()

    materialField = iron.Field()
    equationsSet.MaterialsCreateStart(materialsFieldUserNumber, materialField)
    materialField.VariableLabelSet(iron.FieldVariableTypes.U, "SmoothingParameters")
    equationsSet.MaterialsCreateFinish()
    materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                              1, tau)
    materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                              2, kappa)

    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    equations.outputType = iron.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    problem = iron.Problem()
    problemSpecification = [iron.ProblemClasses.FITTING,
                            iron.ProblemTypes.DATA_FITTING,
                            iron.ProblemSubtypes.STATIC_FITTING]
    problem.CreateStart(problemUserNumber, problemSpecification)
    problem.CreateFinish()

    problem.ControlLoopCreateStart()
    problem.ControlLoopCreateFinish()

    solver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
    solver.outputType = iron.SolverOutputTypes.MONITOR
    problem.SolversCreateFinish()

    solverEquations = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    solver.SolverEquationsGet(solverEquations)
    solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
    solverEquations.BoundaryConditionsCreateFinish()

    # self.equationsSet.DerivedVariableCalculate(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS)
    problem.Solve()
    problem.Finalise()
