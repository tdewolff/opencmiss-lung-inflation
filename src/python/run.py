from lung_inflation import *
from scipy import ndimage

# Cube example
# lungInflation = LungInflation()
# lungInflation.LoadHermiteMesh('./cube_hermite.exelem', './cube_hermite.exnode', 'Geometry')
# lungInflation.LoadDisplacements('./displacement/cube_dx.nii', './displacement/cube_dy.nii', './displacement/cube_dz.nii')
# lungInflation.Setup()
# lungInflation.Solve("./results/cube", 1)
# sys.exit()

transformationCube = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]
transformationLungMesh = [[1,0,0,0], [0,-1,0,512], [0,0,-1,0], [0,0,0,1]]  # flip y, invert z
transformationLungVessels = [[-1,0,0,512], [0,1,0,0], [0,0,-1,0], [0,0,0,1]]  # flip x, invert z

# Lung
path = '/hpc/tdew803/Lung/Data/Pig_PE_Study_HRC/'
subject = 'AP00157'
refPressure = '9.4cmH2O'
defPressure = '15cmH2O'
registration = '9_15'

print('Subject:', subject)
print('Registration:', registration)

displacementFiles = path + subject + '/reg_' + registration + '_d%s.nii'
displacement = Displacement(displacementFiles % 'x', displacementFiles % 'y', displacementFiles % 'z')

# Evaluate registration error AP00157 Left
refVesselsFilename = path + subject + '/' + refPressure + '/Vessel/AP00157_9_4cmH2O_ArteryCenterLowerLeft.exnode'
defVesselsFilename = path + subject + '/' + defPressure + '/Vessel/AP00157_15cmH2O_ArteryCenterLowerLeft_Reordered.exnode'
errors_9, deformation_9 = NodalRegistrationError(refVesselsFilename, defVesselsFilename, 'coordinates', transformationLungVessels, displacement,
                               nids=[1,2,3,4,6,7,8,11,13,14,15,21,22,25,27,28,36,65,66,78,79,80,88,89,103,104,113,114,
                                     115,116])

# Evaluate registration error AP00157 Right
refVesselsFilename = path + subject + '/' + refPressure + '/Vessel/AP00157_9_4cmH2O_ArteryCenterLowerRight_Reordered.exnode'
defVesselsFilename = path + subject + '/' + defPressure + '/Vessel/AP00157_15cmH2O_ArteryCenterLowerRight.exnode'
errors_15, deformation_15 = NodalRegistrationError(refVesselsFilename, defVesselsFilename, 'coordinates', transformationLungVessels, displacement,
                               nids=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,29,30,31,32,
                                     33,34,35,36,37,39,40,41,43,44,45,47,48,49,50,51,52,55,56,57,58,60,63,64,65,66,67,
                                     68,69,71,72,74,76,80,85,91,92,98,100,106,107,108,112,115,116,119,120])

print(np.mean(errors_9 + errors_15), np.std(errors_9 + errors_15))
print(np.mean(deformation_9 + deformation_15), np.std(deformation_9 + deformation_15))

# Plot Joint Density Histogram
# frcFilename = path + subject + '/' + refPressure + '/Raw/' + subject + '.nii'
# tlcFilename = path + subject + '/' + defPressure + '/Raw/' + subject + '.nii'
# frcMaskFilename = path + subject + '/' + refPressure + '/Lung/' + subject + '.Lung.mask.img'
# outputFilename = path + subject + '/out/' + subject + '_' + registration + '_JDH.nii'
# hist = DensityHistogram2D(frcFilename, frcMaskFilename, tlcFilename, displacement, (-1000, -300), outputFilename)
# mask = hist > 50
# hist = np.ma.masked_where(~mask, hist)
# print('Plotting...')
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.imshow(hist, cmap='jet', origin='lower')
# ax.set_xticklabels(['%d' % f for f in ax.get_xticks() - 1000])
# ax.set_yticklabels(['%d' % f for f in ax.get_yticks() - 1000])
# ax.set_xlabel('FRC density (cmH2O)')
# ax.set_ylabel('TLC density (cmH2O)')
# ax.set_title('Joint Density Histogram %s' % (subject))
# plt.show()

# Fit mechanical expansion
# exelemFilename = path + subject + '/' + refPressure + '/Lung/FEMesh/Left_RefittedNoVersion.exelem'
# exnodeFilename = path + subject + '/' + refPressure + '/Lung/FEMesh/Left_RefittedNoVersion.exnode'
# lungInflation = LungInflation()
# lungInflation.SetDisplacement(displacement)
# lungInflation.LoadHermiteMesh(exelemFilename, exnodeFilename, 'coordinates', transformationLungMesh)
# lungInflation.AddFittingField("CauchyStress", FittingVariableTypes.CAUCHY_STRESS)
# lungInflation.AddFittingField("GreenLagrangeStrain", FittingVariableTypes.GREEN_LAGRANGE_STRAIN)
# lungInflation.AddFittingField("AverageCauchyStress", FittingVariableTypes.AVERAGE_STRESS)
# lungInflation.AddFittingField("AverageGreenLagrangeStrain", FittingVariableTypes.AVERAGE_STRAIN)
# lungInflation.AddFittingField("Jacobian", FittingVariableTypes.JACOBIAN)
# lungInflation.Setup()
# lungInflation.Solve("./results/out", 1)
