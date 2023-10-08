# -*- coding: utf-8 -*-
"""
@author: Osikomaiya Folajimi
"""

import csv
import numpy as np
import time
from os import listdir
import os

# Definition of functions used in the code, all below

# Active power
def active_power(bus, theta, voltage, conductance, susceptance, neighbor_positions):
    p = 0
    for i in neighbor_positions[bus]:
        p = p + \
            voltage[bus] * voltage[i] * (
                conductance[bus, i] * np.cos(theta[bus] - theta[i]) +
                susceptance[bus, i] * np.sin(theta[bus] - theta[i])
            )
    return p

# Reactive power
def reactive_power(bus, theta, voltage, conductance, susceptance, neighbor_positions):
    q = 0
    for i in neighbor_positions[bus]:
        q = q + \
            voltage[bus] * voltage[i] * (
                conductance[bus, i] * np.sin(theta[bus] - theta[i]) -
                susceptance[bus, i] * np.cos(theta[bus] - theta[i])
            )
    return q

# Matrix H
def matrix_h(h, theta, voltage, conductance, susceptance, neighbor_positions, non_voltage_theta_pos):
    rows, cols = h.shape
    for i in range(rows):
        for d in range(cols):
            if i != d:
                h[i, d] = voltage[i] * voltage[d] * (
                    conductance[i, d] * np.sin(theta[i] - theta[d]) -
                    susceptance[i, d] * np.cos(theta[i] - theta[d])
                )
            else:
                h[i, d] = -1 * voltage[i] * voltage[i] * susceptance[i, d]
                for k in neighbor_positions[i]:
                    h[i, d] = h[i, d] - voltage[i] * voltage[k] * (
                        conductance[i, k] * np.sin(theta[i] - theta[k]) -
                        susceptance[i, k] * np.cos(theta[i] - theta[k])
                    )
    for i in range(rows):
        for d in range(cols):
            if i == non_voltage_theta_pos and d == non_voltage_theta_pos:
                h[i, d] = np.inf
    return h

# Matrix N
def matrix_n(n, theta, voltage, conductance, susceptance, neighbor_positions):
    rows, cols = n.shape
    for i in range(rows):
        for d in range(cols):
            if i != d:
                n[i, d] = voltage[i] * (
                    conductance[i, d] * np.cos(theta[i] - theta[d]) +
                    susceptance[i, d] * np.sin(theta[i] - theta[d])
                )
            else:
                n[i, d] = voltage[i] * conductance[i, d]
                for k in neighbor_positions[i]:
                    n[i, d] = n[i, d] + voltage[k] * (
                        conductance[i, k] * np.cos(theta[i] - theta[k]) +
                        susceptance[i, k] * np.sin(theta[i] - theta[k])
                    )
    return n

# Matrix M
def matrix_m(m, theta, voltage, conductance, susceptance, neighbor_positions):
    rows, cols = m.shape
    for i in range(rows):
        for d in range(cols):
            if i != d:
                m[i, d] = -1 * voltage[i] * voltage[d] * (
                    conductance[i, d] * np.cos(theta[i] - theta[d]) +
                    susceptance[i, d] * np.sin(theta[i] - theta[d])
                )
            else:
                m[i, d] = -1 * voltage[i] * voltage[i] * conductance[i, d]
                for k in neighbor_positions[i]:
                    m[i, d] = m[i, d] + voltage[i] * voltage[k] * (
                        conductance[i, k] * np.cos(theta[i] - theta[k]) +
                        susceptance[i, k] * np.sin(theta[i] - theta[k])
                    )
    return m

# Matrix L
def matrix_l(l, theta, voltage, conductance, susceptance, neighbor_positions, non_voltage_theta_pos, non_pq_voltage_pos):
    rows, cols = l.shape
    for i in range(rows):
        for d in range(cols):
            if i != d:
                l[i, d] = voltage[i] * (
                    conductance[i, d] * np.sin(theta[i] - theta[d]) -
                    susceptance[i, d] * np.cos(theta[i] - theta[d])
                )
            else:
                l[i, d] = -1 * voltage[i] * susceptance[i, d]
                for k in neighbor_positions[i]:
                    l[i, d] = l[i, d] + voltage[k] * (
                        conductance[i, k] * np.sin(theta[i] - theta[k]) -
                        susceptance[i, k] * np.cos(theta[i] - theta[k])
                    )
    for i in range(rows):
        for d in range(cols):
            if i == non_voltage_theta_pos and d == non_voltage_theta_pos:
                l[i, d] = np.inf
    for i in range(rows):
        for d in range(cols):
            if i in non_pq_voltage_pos and d in non_pq_voltage_pos:
                l[i, i] = np.inf
    return l


# Menu for printing
def menu():
    print('Completed.\n\nChoose one of the options below for printing the results or enter any other value to exit:\n')
    print('1 - Active losses in transmission')
    print('2 - Reactive losses in transmission')
    print('3 - Final voltages and angles')
    print('4 - Net active power injected at buses')
    print('5 - Net reactive power injected at buses')
    print('6 - Active power flow in transmission')
    print('7 - Reactive power flow in transmission')
    print('8 - N-1 Contingency Analysis')
    
    


opi = 0.0
Pmin = 0.0
Pmax = 0.0
opi_list = []

def program(perform_contingency = False, line_choice="", multiple=True):
    # List all the text files in the current directory
    txt_files = [x for x in listdir() if x.endswith('.txt')]
    choices = []

    # Generate a list of choices for the user to select from
    for i, filename in enumerate(txt_files):
        choices.append(f"{i+1} - {filename}")
    
    print("\n===============================================================")
    print("Choose a file from the list below by entering the integer value:")
    print("Note: The file should be in the same directory as the script.\n")

    # Display the available choices to the user
    for choice in choices:
        print(choice)
        
    possible_choices = [str(i) for i in range(1, len(txt_files) + 1)]
    chosen_option = input('\nYour Choice: ')

    # Keep asking for input until a valid choice is made
    while chosen_option not in possible_choices:
        print("Try again! Enter valid value!")
        chosen_option = input('\nYour Choice: ')

    error_input = input('\nEnter the value error or press Enter for the default error (10e-5): ')

    # Set the error to the default value if no input is provided
    if error_input == "":
        error_input = 0.00001
    else:
        error_input = float(error_input)

    start_time = time.time()
    
    # Open the selected text file
    with open(txt_files[int(chosen_option) - 1]) as source_txt_file:
        reader = csv.reader(source_txt_file)
        data = list(reader)
        
    # Find the positions where the data ends based on the '-999' indicator
    cut_points = []

    for i in range(len(data)):
        if data[i-1][0][:4] == '-999':
            cut_points.append(i-1)
    
    del(cut_points[-1])  # Remove the last element as it's not needed

    # Extract the system base value (sbase) from the first line of data
    sbase = float(data[0][0][31:36])
    
    # Remove the second column (containing strings) from the node data for conversion
    for i in range(2, cut_points[0]):
        data[i][0] = data[i][0][:5] + data[i][0][17:]
    
    # Create a list of nodes
    nodes = []
    num_nodes = 0
    for i in range(2, cut_points[0]):
        nodes.append(data[i])
        num_nodes += 1
    
    # Create a list of branches
    num_branches = 0
    branches = []
    for i in range(cut_points[0] + 2, cut_points[1]):
        branches.append(data[i])
        num_branches += 1
    
    # Use NumPy to create an array of node data
    nodes_data = np.arange(num_nodes * 17, dtype=float).reshape(num_nodes, 17)
    for i in range(num_nodes):
        nodes_data[i] = np.fromstring(nodes[i][0], dtype=float, sep=' ')

    """  
    1   Bus number (I) *
    2   Load flow area number (I) Don't use zero! *
    3   Loss zone number (I)
    4   Type (I) *
            0 - Unregulated (load, PQ)
            1 - Hold MVAR generation within voltage limits, (PQ)
            2 - Hold voltage within VAR limits (gen, PV)
            3 - Hold voltage and angle (swing, V-Theta) (must always have one)
    5   Final voltage, p.u. (F) *
    6   Final angle, degrees (F) *
    7   Load MW (F) *
    8   Load MVAR (F) *
    9   Generation MW (F) *
    10  Generation MVAR (F) *
    11  Base KV (F)
    12  Desired volts (pu) (F) (This is desired remote voltage if this bus is 
        controlling another bus.
    13  Maximum MVAR or voltage limit (F)
    14  Minimum MVAR or voltage limit (F)
    15  Shunt conductance G (per unit) (F) *
    16  Shunt susceptance B (per unit) (F) *
    17  Remote controlled bus number
    """
    # Using NumPy to create an array of branches
    branches_data = np.arange(num_branches * 21, dtype=float).reshape(num_branches, 21)
    for i in range(num_branches):
        branches_data[i] = np.fromstring(branches[i][0], dtype=float, sep=' ')

    # From this point, the data is ready to be processed.
    # The respective header for each column is as follows (21 columns):
    """
    1   Tap bus number (I) *
                    For transformers or phase shifters, the side of the model
                    the non-unity tap is on
    2   Z bus number (I) *
                    For transformers and phase shifters, the side of the model
                    the device impedance is on.
    3   Load flow area (I)
    4   Loss zone (I)
    5   Circuit (I) * (Use 1 for single lines)
    6   Type (I) *
                    0 - Transmission line
                    1 - Fixed tap
                    2 - Variable tap for voltage control (TCUL, LTC)
                    3 - Variable tap (turns ratio) for MVAR control
                    4 - Variable phase angle for MW control (phase shifter)
    7   Branch resistance R, per unit (F) *
    8   Branch reactance X, per unit (F) * No zero impedance lines
    9   Line charging B, per unit (F) * (total line charging, +B)
    10  Line MVA rating No 1 (I) Left justify!
    11  Line MVA rating No 2 (I) Left justify!
    12  Line MVA rating No 3 (I) Left justify!
    13  Control bus number
    14  Side (I)
                    0 - Controlled bus is one of the terminals
                    1 - Controlled bus is near the tap side
                    2 - Controlled bus is near the impedance side (Z bus)
    15  Transformer final turns ratio (F)
    16  Transformer (phase shifter) final angle (F)
    17  Minimum tap or phase shift (F)
    18  Maximum tap or phase shift (F)
    19  Step size (F)
    20  Minimum voltage, MVAR or MW limit (F)
    21  Maximum voltage, MVAR or MW limit (F)
    """

    # Create a list of neighbors and another with positions of each bus, excluding the
    # own bus. The syntax is, for example, neighbors of bus k -> neighbors[k-1]
    # Syntax for positions is similar, positions of neighbors of k (in branches_data):
    # -> neighbor_positions[k-1]

    neighbors = []            # Neighbors will not be used, only for reference
    neighbors_K = []
    neighbor_positions = []   # Where positions of neighboring buses will be allocated
    pos_neighbor_K = []
    neighbor_positions_branches = [] # List with line positions of each neighbor

    for i in range(len(nodes_data)):
        neighbors.append([])
        neighbors_K.append([])
        neighbor_positions.append([])
        pos_neighbor_K.append([])
        neighbor_positions_branches.append([])
        for j in range(len(branches_data)):
            if int(nodes_data[i][0]) == int(branches_data[j][0]):
                neighbor_positions_branches[i].append(j)
                neighbors[i].append(int(branches_data[j][1]))
                neighbors_K[i].append(int(branches_data[j][1]))
                neighbor_positions[i].append(nodes_data[:, 0].tolist().index(branches_data[j][1]))
                pos_neighbor_K[i].append(nodes_data[:, 0].tolist().index(branches_data[j][1]))
            elif int(nodes_data[i][0]) == int(branches_data[j][1]) and int(branches_data[j][0]) not in neighbors[i]:
                neighbor_positions_branches[i].append(j)
                neighbors[i].append(int(branches_data[j][0]))
                neighbors_K[i].append(int(branches_data[j][0]))
                neighbor_positions[i].append(nodes_data[:, 0].tolist().index(branches_data[j][0]))
                pos_neighbor_K[i].append(nodes_data[:, 0].tolist().index(branches_data[j][0]))
        neighbors_K[i].append(i + 1)
        pos_neighbor_K[i].append(nodes_data[:, 0].tolist().index(nodes_data[:, 0][i]))

    for i in neighbors:
        i.sort()
    for i in neighbors_K:
        i.sort()
    for i in neighbor_positions:
        i.sort()
    for i in pos_neighbor_K:
        i.sort()

    pos_nodes_in_nodes = nodes_data[:, 0]
    pos_nodes_in_nodes = pos_nodes_in_nodes.tolist()

    pos_nodes_in_branches = branches_data[:, 0::]
    pos_nodes_in_branches = pos_nodes_in_branches.tolist()

    for i in range(len(pos_nodes_in_branches)):
        for term in range(2):
            pos_nodes_in_branches[i][term] = float(pos_nodes_in_nodes.index(pos_nodes_in_branches[i][term]))

    # Initializing the Admittance Matrix. Remember that Y = G + jB, and the dimension is num_nodes x num_nodes
    ypu = np.arange(num_nodes * num_nodes).reshape(num_nodes, num_nodes).astype('complex')
    ypu.fill(0)  # Fill all elements with 0

    # print("branches_data=====================", branches_data)
    # print("nodes_data=====================", nodes_data)

    line_banches = []
    # Constructing elements outside the main diagonal
    for i in branches_data:
        ypu[int(i[0]) - 1, int(i[1]) - 1] = -1 * (1 / (i[7] * 1j + i[6]) * np.exp(-1 * (i[15] * 1j)))
        ypu[int(i[1]) - 1, int(i[0]) - 1] = -1 * (1 / (i[7] * 1j + i[6]) * np.exp(1 * (i[15] * 1j)))

        line_banches.append(f"{int(i[0])}" + " " + f"{int(i[1])}")

    # Constructing the main diagonal, always non-zero
    for i in range(num_nodes):
        ypu[i][i] = nodes_data[i, 14] + nodes_data[i, 15] * 1j
        for j in neighbor_positions_branches[i]:
            ypu[i][i] = ypu[i][i] + branches_data[j][8] * (1j/2) + (1 / (branches_data[j][7] * 1j + branches_data[j][6]))

    # global perform_contingency
    if perform_contingency == True:
        if line_choice == "":
            line_choice = input("\nWhich lines do you want to switch off (use whitespace): ")
        a, b = line_choice.split()
        
        ypu[int(a)-1][int(a)-1] += ypu[int(a)-1][int(b)-1]
        ypu[int(b)-1][int(b)-1] += ypu[int(b)-1][int(a)-1]

        ypu[int(a)-1][int(b)-1] = 0
        ypu[int(b)-1][int(a)-1] = 0

        # print("=====================")
        # print(ypu)
        # print("=====================")


    

    # print("=====================")
    # print(ypu)
    # print("=====================")

    gpu = np.real(ypu)  # Real part of Ypu
    bpu = np.imag(ypu)  # Imaginary part of Ypu

    # Calculating the number of PQ and PV equations for Jacobian matrix dimension
    npq, npv, nvtheta = (0, 0, 0)
    pos_npq = []
    pos_npv = []

    for i in range(len(nodes_data[:, 3])):
        if nodes_data[:, 3][i] == 0.0 or nodes_data[:, 3][i] == 1.0:
            npq = npq + 1
            pos_npq.append(i)
        elif nodes_data[:, 3][i] == 2.0:
            npv = npv + 1
            pos_npv.append(i)
        else:
            nvtheta = nvtheta + 1
            pos_nvtheta = i

    # Initializing submatrices of the Jacobian matrix
    h = np.arange(num_nodes * num_nodes, dtype=np.float64).reshape(num_nodes, num_nodes)
    h.fill(0)
    n = np.arange(num_nodes * num_nodes, dtype=np.float64).reshape(num_nodes, num_nodes)
    n.fill(0)
    m = np.arange(num_nodes * num_nodes, dtype=np.float64).reshape(num_nodes, num_nodes)
    m.fill(0)
    l = np.arange(num_nodes * num_nodes, dtype=np.float64).reshape(num_nodes, num_nodes)
    l.fill(0)

    jac = np.arange((2 * num_nodes) * (2 * num_nodes), dtype=np.float64).reshape(2 * num_nodes, 2 * num_nodes)
    jac.fill(0)

    # Data and Unknowns
    nodes_pos = []

    for i in range(len(nodes_data[:, 0])):
        nodes_pos.append(i)

    for i in range(len(nodes_data[:, 0])):
        nodes_pos.append(i)

    ppu = (nodes_data[:, 8] - nodes_data[:, 6]) / sbase
    qpu = (nodes_data[:, 9] - nodes_data[:, 7]) / sbase
    vpu = nodes_data[:, 4]
    theta = nodes_data[:, 5] * (np.pi/180.)  # Theta vector converted to radians

    # Arranging vectors, meaning, for PQ buses, V and theta are unknowns. For
    # PV buses, Q and theta are unknowns. For Vtheta buses, P and Q are unknowns
    # The initial guesses are based on the reference bus and are adjusted below

    ref_vpu = nodes_data[pos_nvtheta, 4]
    ref_theta = nodes_data[pos_nvtheta, 5]

    for i in pos_npq:
        vpu[i] = ref_vpu
        theta[i] = ref_theta

    for i in pos_npv:
        theta[i] = ref_theta

    # Creating vectors for different types of buses
    vector_pos_pa = list(pos_npq)
    for i in pos_npv:
        vector_pos_pa.append(i)
    vector_pos_pa.sort()

    vector_pos_pr = list(pos_npq)
    vector_pos_pr.sort()

    vector_pos_data = list(vector_pos_pa)
    for i in vector_pos_pr:
        vector_pos_data.append(i)

    # Defining the cut for theta values in the vector
    cut_theta = len(vector_pos_pa)

    # Creating NumPy arrays for positions
    array_pos_pa = np.array(vector_pos_pa).reshape(len(vector_pos_pa), 1)
    array_pos_pr = np.array(vector_pos_pr).reshape(len(vector_pos_pr), 1)
    array_pos_pa = array_pos_pa.astype(dtype=np.float64)
    array_pos_pr = array_pos_pr.astype(dtype=np.float64)

    # Creating lists for expected values
    expected_p = list(ppu)
    expected_q = list(qpu)
    expected = list(expected_p)
    for i in expected_q:
        expected.append(i)

    # Creating copies of position vectors
    deltap = np.copy(array_pos_pa)
    deltaq = np.copy(array_pos_pr)

    # Initialize deltap and deltaq with a value greater than e
    for bus in range(len(deltap)):
        deltap[bus][0] = error_input + 1
    for bus in range(len(deltaq)):
        deltaq[bus][0] = error_input + 1

    # Create delta vector
    delta = np.arange(2 * num_nodes).reshape(2 * num_nodes, 1)
    delta = delta.astype(dtype=np.float64)
    del_e = np.copy(delta)

    # Initialize delta and del_e with a value greater than e
    for i in range(len(delta)):
        delta[i][0] = error_input + 1

    # Initialize iteration count
    iteration = 0

    #-----------------------------SUBSYSTEM I-------------------------------------
    while (abs(del_e).max()) >= error_input:

        # Calculate delta values for each bus
        for i in range(2 * num_nodes):
            if i < num_nodes:
                delta[i][0] = expected[i] - active_power(nodes_pos[i], theta, vpu, gpu, bpu, pos_neighbor_K)
            else:
                delta[i][0] = expected[i] - reactive_power(nodes_pos[i], theta, vpu, gpu, bpu, pos_neighbor_K)

        # Initialize submatrices
        h.fill(0)
        n.fill(0)
        m.fill(0)
        l.fill(0)

        # Calculate submatrices
        matrix_h(h, theta, vpu, gpu, bpu, pos_neighbor_K, pos_nvtheta)
        matrix_n(n, theta, vpu, gpu, bpu, pos_neighbor_K)
        matrix_m(m, theta, vpu, gpu, bpu, pos_neighbor_K)
        matrix_l(l, theta, vpu, gpu, bpu, pos_neighbor_K, pos_nvtheta, pos_npv)

        # Include submatrix data into -J
        lin, col = h.shape
        for i in range(lin):
            for d in range(col):
                jac[i, d] = h[i, d]

        lin, col = n.shape
        for i in range(lin):
            for j in range(col):
                jac[i, d + num_nodes] = n[i, d]

        lin, col = m.shape
        for i in range(lin):
            for d in range(col):
                jac[i + num_nodes, d] = m[i, d]

        lin, col = l.shape
        for i in range(lin):
            for d in range(col):
                jac[i + num_nodes, d + num_nodes] = l[i, d]

        # Calculate the inverse of the Jacobian matrix
        jacinv = np.linalg.inv(jac)

        # Solve for solution vector
        sol = jacinv @ delta

        # Update theta and voltage magnitude values
        for i in range(num_nodes):
            theta[i] = theta[i] + sol[i][0]
            vpu[i] = vpu[i] + sol[i + num_nodes][0]

        # Update del_e vector
        for i in range(2 * num_nodes):
            if i not in vector_pos_data:
                del_e[i][0] = 0
            else:
                del_e[i][0] = delta[i][0]

        # Increment iteration count
        iteration = iteration + 1

    #----------------------------SUBSYSTEM II-------------------------------------
    # Power at the reference bus
    for i in pos_npv:
        ppu[i] = active_power(i, theta, vpu, gpu, bpu, pos_neighbor_K)

    ppu[pos_nvtheta] = active_power(pos_nvtheta, theta, vpu, gpu, bpu, pos_neighbor_K)
    qpu[pos_nvtheta] = reactive_power(pos_nvtheta, theta, vpu, gpu, bpu, pos_neighbor_K)
    #------------------------------------------------------------------------------
    
    #---------------------------SUBSYSTEM III-------------------------------------
    # Power flows
    active_power_flow = np.arange(num_nodes * num_nodes, dtype=np.float64).reshape(num_nodes, num_nodes)
    active_power_flow.fill(0)

    for i in pos_nodes_in_branches:
        conv = (1 / ((i[6] + i[7] * 1j)))

        active_power_flow[int(i[0]), int(i[1])] = np.real(conv) * (vpu[int(i[0])]**2) - \
            vpu[int(i[0])] * vpu[int(i[1])] * (np.real(conv) * np.cos(theta[int(i[0])] - \
            theta[int(i[1])]) + np.imag(conv) * np.sin(theta[int(i[0])] - \
                    theta[int(i[1])]))

        active_power_flow[int(i[1]), int(i[0])] = np.real(conv) * (vpu[int(i[1])]**2) - \
            vpu[int(i[0])] * vpu[int(i[1])] * (np.real(conv) * np.cos(theta[int(i[0])] - \
            theta[int(i[1])]) - np.imag(conv) * np.sin(theta[int(i[0])] - \
                    theta[int(i[1])]))

    # Suggestion: in the given problems, there are no phase-shifting transformers, so
    # phi_km is zero, but it can be included in the calculation

    reactive_power_flow = np.arange(num_nodes * num_nodes, dtype=np.float64).reshape(num_nodes, num_nodes)
    reactive_power_flow.fill(0)

    for i in pos_nodes_in_branches:
        conv = (1 / ((i[6] + i[7] * 1j)))

        reactive_power_flow[int(i[0]), int(i[1])] = -1 * (np.imag(conv) + (i[8] / 2)) * \
            (vpu[int(i[0])]**2) - vpu[int(i[0])] * vpu[int(i[1])] * (np.real(conv) * \
                np.sin(theta[int(i[0])] - theta[int(i[1])]) - np.imag(conv) * \
                np.cos(theta[int(i[0])] - theta[int(i[1])]))

        reactive_power_flow[int(i[1]), int(i[0])] = -1 * (np.imag(conv) + (i[8] / 2)) * \
            (vpu[int(i[1])]**2) + vpu[int(i[0])] * vpu[int(i[1])] * (np.real(conv) * \
                np.sin(theta[int(i[0])] - theta[int(i[1])]) + np.imag(conv) * \
                np.cos(theta[int(i[0])] - theta[int(i[1])]))

    print("==============================================")
    print('\nThe system converged in the ' + str(iteration) + \
        'th iteration\nThe execution time was ' + str(time.time() - start_time) + \
        ' seconds\n\nPress Enter to continue\n')
    input()
    #------------------------------------------------------------------------------
    #%%
    # Losses and system results

    active_losses = []
    reactive_losses = []
    final_voltages_and_angles = []
    active_powers = []
    reactive_powers = []
    active_flows = []
    reactive_flows = []

    for i in pos_nodes_in_branches:
        active_losses.append(['Active losses between bus ' + \
                                str(int(pos_nodes_in_nodes[int(i[0])])) + ' and ' + \
                                str(int(pos_nodes_in_nodes[int(i[1])])) + ' (or ' + \
                                str(int(pos_nodes_in_nodes[int(i[1])])) + ' and ' + \
                                str(int(pos_nodes_in_nodes[int(i[0])])) + '): ' + str \
                                (np.round(active_power_flow[int(i[0]), int(i[1])] + \
                                        active_power_flow[int(i[1]), int(i[0])], 5)) + ' pu'])

        reactive_losses.append(['Reactive losses between bus ' + \
                                str(int(pos_nodes_in_nodes[int(i[0])])) + ' and ' + \
                                str(int(pos_nodes_in_nodes[int(i[1])])) + ' (or ' + \
                                str(int(pos_nodes_in_nodes[int(i[1])])) + ' and ' + \
                                str(int(pos_nodes_in_nodes[int(i[0])])) + '): ' + str \
                                (np.round(reactive_power_flow[int(i[0]), int(i[1])] + \
                                        reactive_power_flow[int(i[1]), int(i[0])], 5)) + ' pu'])

        active_flows.append(['Active power flow from bus ' + \
                            str(int(pos_nodes_in_nodes[int(i[0])])) + ' to ' + \
                            str(int(pos_nodes_in_nodes[int(i[1])])) + ' and ' + \
                            str(int(pos_nodes_in_nodes[int(i[1])])) + ' to ' + \
                            str(int(pos_nodes_in_nodes[int(i[0])])) + \
                            ', respectively: ' + str \
                            (np.round(active_power_flow[int(i[0]), int(i[1])], 5)) + ' and ' + \
                                    str(np.round(active_power_flow[int(i[1]), int(i[0])], 5)) + \
                                    ' pu'])

        reactive_flows.append(['Reactive power flow from bus ' + \
                                str(int(pos_nodes_in_nodes[int(i[0])])) + ' to ' + \
                                str(int(pos_nodes_in_nodes[int(i[1])])) + ' and ' + \
                                str(int(pos_nodes_in_nodes[int(i[1])])) + ' to ' + \
                                str(int(pos_nodes_in_nodes[int(i[0])])) + \
                                ', respectively: ' + str \
                                (np.round(reactive_power_flow[int(i[0]), int(i[1])], 5)) + ' and ' + \
                                    str(np.round(reactive_power_flow[int(i[1]), int(i[0])], 5)) + \
                                    ' pu'])
    

    for i in range(num_nodes):
        final_voltages_and_angles.append(['The final voltage at bus ' + \
                                            str(int(pos_nodes_in_nodes[i])) + ' is ' + \
                                            str(np.round(vpu[i], 5)) + \
                                                ' pu, with an angle of ' + str(np.round(theta[i], 5)) + \
                                                ' radians (or ' + \
                                                str(np.round(np.degrees(theta[i]), 5)) + \
                                                '°) with respect to the reference'])

        active_powers.append(['Net active power injected at bus ' + \
                                str(int(pos_nodes_in_nodes[i])) + ' is ' + \
                                str(np.round(active_power(nodes_pos[i], \
                                        theta, vpu, gpu, bpu, pos_neighbor_K), 5)) + ' pu'])

        reactive_powers.append(['Net reactive power injected at bus ' + \
                                str(int(pos_nodes_in_nodes[i])) + ' is ' + \
                                str(np.round(reactive_power(nodes_pos[i], \
                                        theta, vpu, gpu, bpu, pos_neighbor_K), 5)) + ' pu'])
    
    Pi_values = []
    for i in pos_nodes_in_branches:
        Pi_values.append(np.round(active_power_flow[int(i[0]), int(i[1])], 5))

    Vi_values = []
    Vj_values = []
    impedance = []

    for i in branches_data:
        Vi_values.append(np.round(vpu[int(i[0]) - 1], 5))
        Vj_values.append(np.round(vpu[int(i[1]) - 1], 5))
        impedance.append(i[7])

    
    # Define the system parameters
    NL = len(impedance)  # Number of lines
    W = 1.0  # Weighting factor
    n = 1.0  # Penalty function

    # Calculate PIP
    pip = 0.0
    for i in range(NL):
        Pi = Pi_values[i]
        Pi_max = (Vi_values[i] * Vj_values[i]) / impedance[i]
        pip += (W / (2 * n)) * ((Pi / Pi_max) ** (2 * n))

    print(f"Active Power Performance Index (PIP): {pip:.5f}")

    

    #-----------------------------PIV---------------------------
    # Define the system parameters (you can adjust these values as needed)
    W = 1.0  # Weighting factor
    n = 1.0  # Penalty function

    Vi_values = []
    for i in range(num_nodes):
        Vi_values.append(np.round(vpu[i], 5)) 

    NG = len(Vi_values)
    piv = 0.0
    Vmax, Vmin = 1.05, 0.95
    for i in range(NG):
        # Calculate PIV for the given bus
        piv = (W / (2 * n)) * ((Vi_values[i] - Vmin) / (Vmax - Vmin))**(2 * n)
    
    print(f"Voltage Performance Index (PIP): {piv:.5f}\n")
    

    opi = pip + piv
    print(f"Overall Performance Index (OPI): {opi:.5f}\n")
    
    
    if perform_contingency == True and multiple == True:   
        global opi_list     
        opi_list.append(opi)
        global Pmin
        Pmin = min(opi_list)
        global Pmax
        Pmax = max(opi_list)
        return

    # Printing

    menu()
    option = input("\nChoice: ")

    while option in ['1', '2', '3', '4', '5', '6', '7', '8']:
        if option == '1':
            print('')
            for i in active_losses:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '2':
            print('')
            for i in reactive_losses:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '3':
            print('')
            for i in final_voltages_and_angles:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '4':
            print('')
            for i in active_powers:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '5':
            print('')
            for i in reactive_powers:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '6':
            print('')
            for i in active_flows:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '7':
            print('')
            for i in reactive_flows:
                print(i[0])
            print('\nPress Enter to continue')
            input()
            menu()
            option = input("\nChoice: ")
            print('')
        elif option == '8':            
            os.system('cls')
            choice = input("Do you want to perform contingency analysis for all the branch lines? y/n:  ")
            if choice == 'y':
                for branch in line_banches:
                    program(perform_contingency = True, line_choice=branch)
                    print(f"Transmission line {branch} completed\n")
                new_list = []
                i=0
                #  Normalization of OPI data
                for branch in line_banches:                    
                    Pn = (0.8 * (opi_list[i] - Pmin) / (Pmax - Pmin)) + 0.1                        
                    new_list.append({branch : round(Pn, 4)})

                    i += 1
                # Sort the list of dictionaries by their values in descending order
                sorted_data = sorted(new_list, key=lambda x: list(x.values())[0], reverse=True)


                print("sorted_opi_list: ", sorted_data)
                print('\nPress Enter to continue')
                input()
                menu()
                option = input("\nChoice: ")
                print('')
            else:
                program(perform_contingency = True, multiple=False)
        
    print('\nRestart? (y/n)\n')
    return

restart = ''
program()
while restart != 'y' or restart != 'n':
    restart = input()
    if restart == 'y':
        program()
    else:
        break




