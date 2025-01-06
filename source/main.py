from slirp import SLIRP
from sage.all import sin, pi
import numpy as np
import os
def main():

    #make demo for a 2 link arm
    # print("Starting SLIRP Demo for RR robot arm:")
    # rr = SLIRP(num_links=2)
    # rr.compute_lagrangian()
    # rr.compute_mass_matrix()
    # rr.compute_coriolis_terms()
    # rr.compute_gravity_vector()
    # rr.compute_eom()

    # rr.print_matrix_as_python(rr.M, "Mq")
    # input("Press Enter to continue...")
    # os.system('clear')
    # rr.print_matrix_as_python(rr.C, "Cq")
    # input("Press Enter to continue...")
    # os.system('clear')
    # rr.print_matrix_as_python(rr.G, "Gq")
    # input("Press Enter to continue...")
    # os.system('clear')
    # rr.print_matrix_as_python(rr.tau_matrix, "tau")
    # #make user press enter to contiinue with demo
    # input("Press Enter to continue...")
    # os.system('clear')

    print("Starting SLIRP Demo for RRR robot arm:")
    rrr = SLIRP(num_links=3)
    rrr.compute_lagrangian()
    rrr.compute_mass_matrix()
    rrr.compute_coriolis_terms()
    rrr.compute_gravity_vector()
    rrr.compute_eom()

    rrr.print_matrix_as_python(rrr.M, "Mq")
    print("Mq Dimensions:", rrr.M.dimensions())
    input("Press Enter to continue...")
    os.system('clear')
    rrr.print_matrix_as_python(rrr.C, "Cq")
    print("Cq Dimensions:", rrr.C.dimensions())
    input("Press Enter to continue...")
    os.system('clear')
    rrr.print_matrix_as_python(rrr.G, "Gq")
    print("Gq Dimensions:", rrr.G.dimensions())
    input("Press Enter to continue...")
    os.system('clear')
    rrr.print_matrix_as_python(rrr.tau_matrix, "tau")

    #make user press enter to contiinue with demo
    input("Press Enter to continue...")
    os.system('clear')
    # Define a trajectory for the joints
    q_traj = [sin(0.75*pi*rrr.t), 0.1*rrr.t, (pi/4)*sin(pi*rrr.t)]
    rrr.set_trajectory(q_traj)

    param_values = {
        rrr.link_lengths[0]: 1.0, rrr.link_lengths[1]: 0.8, rrr.link_lengths[2]: 0.5,
        rrr.link_masses[0]: 1.0,  rrr.link_masses[1]: 1.2,  rrr.link_masses[2]: 0.8,
        rrr.link_inertias[0]: 0.1, rrr.link_inertias[1]: 0.05, rrr.link_inertias[2]: 0.025,
        rrr.g: 9.81,
        rrr.w: 0.1
    }

    # Evaluate tau at over 10 seconds at 100 points
    trange = np.linspace(0, 10, 100)
    for i in range(trange.size):
        tau_at_t = rrr.evaluate_tau_at_time(trange[i], param_values)
        print(f"Torques at t={trange[i]:.2f}s:", tau_at_t)


    print("Done!")
    print("Starting SLIRP Demo for UR5 robot arm:")
    input("Press Enter to continue...")
    os.system('clear')
    ur5 = SLIRP(num_links=5)
    ur5.compute_lagrangian()
    ur5.compute_mass_matrix()
    ur5.compute_coriolis_terms()
    ur5.compute_gravity_vector()
    ur5.compute_eom()

    ur5.print_matrix_as_python(ur5.M, "Mq")
    print("Mq Dimensions:", ur5.M.dimensions())
    input("Press Enter to continue...")
    os.system('clear')
    ur5.print_matrix_as_python(ur5.C, "Cq")
    print("Cq Dimensions:", ur5.C.dimensions())
    input("Press Enter to continue...")
    os.system('clear')
    ur5.print_matrix_as_python(ur5.G, "Gq")
    print("Gq Dimensions:", ur5.G.dimensions())
    input("Press Enter to continue...")
    os.system('clear')
    ur5.print_matrix_as_python(ur5.tau_matrix, "tau")
    print("tau Dimensions:", ur5.tau_matrix.dimensions())
    input("Press Enter to continue...")
    os.system('clear')

    # Define a trajectory for the joints
    q_traj = [sin(0.75*pi*ur5.t), 0.1*ur5.t, (pi/4)*sin(pi*ur5.t), 0.2*ur5.t, 0.1*sin(0.5*pi*ur5.t)]
    ur5.set_trajectory(q_traj)

    param_values = {
        ur5.link_lengths[0]: 1.0, ur5.link_lengths[1]: 0.8, ur5.link_lengths[2]: 0.5, ur5.link_lengths[3]: 0.3, ur5.link_lengths[4]: 0.2,
        ur5.link_masses[0]: 1.0,  ur5.link_masses[1]: 1.2,  ur5.link_masses[2]: 0.8, ur5.link_masses[3]: 0.6, ur5.link_masses[4]: 0.4,
        ur5.link_inertias[0]: 0.1, ur5.link_inertias[1]: 0.05, ur5.link_inertias[2]: 0.02, ur5.link_inertias[3]: 0.01, ur5.link_inertias[4]: 0.005,
        ur5.g: 9.81,
        ur5.w: 0.1
    }

    # Evaluate tau at over 10 seconds at 100 points
    trange = np.linspace(0, 10, 100)
    for i in range(trange.size):
        tau_at_t = ur5.evaluate_tau_at_time(trange[i], param_values)
        print(f"Torques at t={trange[i]:.2f}s:", tau_at_t)
    print("Done!")

if __name__ == "__main__":
    main()