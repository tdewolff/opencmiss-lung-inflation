import exfile
import numpy as np
from opencmiss.iron import iron
import scipy

def exfile_to_OpenCMISS(nodeFilename, elementFilename, coordinateField, region, basis, meshUserNumber,
                      interpolation='linear'):
    """Convert an exnode and exelem files to a morphic mesh.

    Only Linear lagrange elements supported.

    Keyword arguments:
    nodeFilename -- exnode filename
    elementFilename -- exelem filename
    coordinateField -- the field to read in
    dimension -- dimension of mesh to read in
    """

    # Load exfiles
    exnode = exfile.Exnode(nodeFilename)
    exelem = exfile.Exelem(elementFilename)

    totalNumberOfNodes = len(exnode.nodeids)
    totalNumberOfElements = len(exelem.elements)

    # Start the creation of a manually generated mesh in the region
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, exelem.num_dims)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(totalNumberOfElements)

    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, totalNumberOfNodes)
    nodes.UserNumbersAllSet(exnode.nodeids)
    nodes.CreateFinish()

    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    elemNums = []
    for elem in exelem.elements:
        elemNums.append(elem.number)
    elements.UserNumbersAllSet(elemNums)
    for elem_idx, elem in enumerate(exelem.elements):
        elements.NodesSet(elem_idx+1, elem.nodes)
    elements.CreateFinish()

    mesh.CreateFinish()

    # Add nodes
    if interpolation == 'linear':
        derivatives = [1]
    elif interpolation == 'hermite':
        derivatives = range(1,9)
    elif interpolation == 'cubicLagrange':
        derivatives = [1]
    coordinates = np.zeros((totalNumberOfNodes, 3, len(derivatives)))
    for node_idx, node_num in enumerate(exnode.nodeids):
        for component_idx, component in enumerate(range(1, 4)):
            for derivative_idx, derivative in enumerate(derivatives):
                component_name = ["x", "y", "z"][component - 1]
                value = exnode.node_value(coordinateField,
                                  component_name, node_num,
                                  derivative)
                coordinates[node_idx,component_idx,derivative_idx] = value
        print('Node added', node_num, coordinates[node_idx,:])

    return mesh, coordinates, exnode.nodeids

def morphic_to_OpenCMISS(morphicMesh, region, basis, meshUserNumber,
                      dimension=2, interpolation='linear',UsePressureBasis=False, pressureBasis=None):
    """Convert an exnode and exelem files to a morphic mesh.

    Only Linear lagrange elements supported.

    Keyword arguments:
    morphicMesh -- morphic mesh
    dimension -- dimension of mesh to read in
    """

    # Create mesh topology
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, 3)
    if (UsePressureBasis):
        mesh.NumberOfComponentsSet(2)
    else:
        mesh.NumberOfComponentsSet(1)

    node_list = morphicMesh.get_node_ids()[1]
    element_list = morphicMesh.get_element_ids()

    mesh.NumberOfElementsSet(len(element_list))
    nodes = iron.Nodes()
    nodes.CreateStart(region, len(node_list))
    nodes.UserNumbersAllSet((scipy.array(node_list)).astype('int32'))
    nodes.CreateFinish()

    MESH_COMPONENT1 = 1
    MESH_COMPONENT2 = 2
    elements = iron.MeshElements()
    elements.CreateStart(mesh, MESH_COMPONENT1, basis)
    elements.UserNumbersAllSet((scipy.array(element_list).astype('int32')))
    global_element_idx = 0
    for element_idx, element in enumerate(morphicMesh.elements):
        global_element_idx += 1
        elements.NodesSet(global_element_idx, scipy.array(element.node_ids, dtype='int32'))
    elements.CreateFinish()

    if (UsePressureBasis):
        pressure_elements = iron.MeshElements()
        pressure_elements.CreateStart(mesh, MESH_COMPONENT2, pressureBasis)
        pressure_elements.AllUserNumbersSet((scipy.array(element_list).astype('int32')))
        for element_idx, element in enumerate(morphicMesh.elements):
            pressure_elements.NodesSet(element.id, scipy.array(element.node_ids, dtype='int32'))
        pressure_elements.CreateFinish()

    mesh.CreateFinish()

    # Add nodes
    if interpolation == 'linear' or interpolation == 'cubicLagrange':
        derivatives = [1]
    elif interpolation == 'hermite':
        derivatives = range(1,9)
    coordinates = np.zeros((len(node_list), 3,len(derivatives)))
    for node_idx,morphic_node in enumerate(morphicMesh.nodes):
        print("node: ", morphic_node.id)
        for comp_idx in range(3):
            print("component: ", comp_idx + 1)
            for derivative_idx, derivative in enumerate(derivatives):
                coordinates[node_idx, comp_idx, derivative_idx] = morphic_node.values[comp_idx]

    return mesh, coordinates, node_list