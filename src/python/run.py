from lung_inflation import LungInflation

# Cube example
lungInflation = LungInflation()
lungInflation.LoadHermiteMesh('./cube_hermite.exelem', './cube_hermite.exnode', 'Geometry')
lungInflation.LoadDisplacements('./displacement/cube_dx.nii', './displacement/cube_dy.nii', './displacement/cube_dz.nii')
lungInflation.Setup()
lungInflation.Solve("./results/cube", 1)

# Lung
path = '/hpc/tdew803/Lung/Data/Pig_PE_Study_HRC/'
subject = 'AP00157'
refPressure = '9.4'
registration = '9_15'

print('Subject:', subject)
print('Registration:', registration)

exelemFilename = path + subject + '/' + refPressure + 'cmH2O/Lung/FEMesh/Left_RefittedNoVersion.exelem'
exnodeFilename = path + subject + '/' + refPressure + 'cmH2O/Lung/FEMesh/Left_RefittedNoVersion.exnode'
displacementFiles = path + subject + '/reg_' + registration + '_d%s.nii'

lungInflation = LungInflation()
lungInflation.LoadHermiteMesh(exelemFilename, exnodeFilename, 'coordinates', True)
lungInflation.LoadDisplacements(displacementFiles % 'x', displacementFiles % 'y', displacementFiles % 'z')
lungInflation.Setup()
lungInflation.Solve("./results/out", 1)
