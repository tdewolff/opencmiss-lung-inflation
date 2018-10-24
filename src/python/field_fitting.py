from opencmiss.iron import iron

fittingUserNumber = 100

def FitField(name, variable, region, decomposition, geometricField, equationsSet, tau=0.01, kappa=0.0005):
    global fittingUserNumber

    fittingFieldUserNumber = fittingUserNumber + 0
    fittingEquationsSetUserNumber = fittingUserNumber + 1
    fittingEquationsSetFieldUserNumber = fittingUserNumber + 2
    fittingDependentFieldUserNumber = fittingUserNumber + 3
    fittingMaterialsFieldUserNumber = fittingUserNumber + 4
    fittingProblemUserNumber = fittingUserNumber + 5
    fittingUserNumber += 10

    fittingField = iron.Field()
    fittingField.CreateStart(fittingFieldUserNumber, region)
    fittingField.TypeSet(iron.FieldTypes.GENERAL)
    fittingField.MeshDecompositionSet(decomposition)
    fittingField.GeometricFieldSet(geometricField)
    fittingField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
    fittingField.NumberOfVariablesSet(2)
    fittingField.VariableTypesSet([iron.FieldVariableTypes.U, iron.FieldVariableTypes.V])
    fittingField.VariableLabelSet(iron.FieldVariableTypes.U, name + "_GaussStress")
    #fittingField.VariableLabelSet(iron.FieldVariableTypes.V, name + "_GaussWeight")
    fittingField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 6)
    fittingField.NumberOfComponentsSet(iron.FieldVariableTypes.V, 6)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 5, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 6, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 1, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 2, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 3, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 4, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 5, 1)
    fittingField.ComponentMeshComponentSet(iron.FieldVariableTypes.V, 6, 1)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 3,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 5,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 6,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 1,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 2,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 3,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 4,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 5,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.ComponentInterpolationSet(iron.FieldVariableTypes.V, 6,
                                           iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    fittingField.CreateFinish()

    equationsSet.DerivedCreateStart(fittingFieldUserNumber, fittingField)
    equationsSet.DerivedVariableSet(variable,  # WHAT DOES THIS DO?
                                    iron.FieldVariableTypes.U)
    equationsSet.DerivedCreateFinish()

    # Create the fitting equations_set
    fittingEquationsSetField = iron.Field()
    fittingEquationsSet = iron.EquationsSet()
    fittingEquationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                                        iron.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                                        iron.EquationsSetSubtypes.GAUSS_POINT_FITTING,
                                        iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
    fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber, region, geometricField,
                                    fittingEquationsSetSpecification,
                                    fittingEquationsSetFieldUserNumber, fittingEquationsSetField)
    fittingEquationsSet.CreateFinish()

    # Create the fitting dependent field
    fittingDependentField = iron.Field()
    fittingEquationsSet.DependentCreateStart(fittingDependentFieldUserNumber, fittingDependentField)
    fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.U, name)
    # Set the number of components to 2
    fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 6)
    fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 6)
    # Set the field variables to be triquadratic Lagrange
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 5, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 6, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 1, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 2, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 3, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 5, 1)
    fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 6, 1)
    fittingEquationsSet.DependentCreateFinish()

    # Create the fitting independent field
    fittingEquationsSet.IndependentCreateStart(fittingFieldUserNumber, fittingField)
    fittingEquationsSet.IndependentCreateFinish()

    # Initialise Gauss point weight field to 1.0
    fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                             iron.FieldParameterSetTypes.VALUES, 1, 1.0)
    fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                             iron.FieldParameterSetTypes.VALUES, 2, 1.0)
    fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                             iron.FieldParameterSetTypes.VALUES, 3, 1.0)
    fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                             iron.FieldParameterSetTypes.VALUES, 4, 1.0)
    fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                             iron.FieldParameterSetTypes.VALUES, 5, 1.0)
    fittingField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,
                                             iron.FieldParameterSetTypes.VALUES, 6, 1.0)

    # Create material field (Sobolev parameters)
    fittingMaterialField = iron.Field()
    fittingEquationsSet.MaterialsCreateStart(fittingMaterialsFieldUserNumber, fittingMaterialField)
    fittingMaterialField.VariableLabelSet(iron.FieldVariableTypes.U, "SmoothingParameters")
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

    # Calculate the stress field
    equationsSet.DerivedVariableCalculate(variable)  # WHAT DOES THIS DO?

    # Solve the fitting problem
    fittingProblem.Solve()
