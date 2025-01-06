from sage.all import *
'''
Class Name: SLIRP - Serial Link Integrated Robot Physics
Inputs: num_links - Number of links in the robot arm
Description: This class automates the symbolic calculations for a generalized serial link robot with revolute joints.
The joints all rotate about the Z axis, and the end-effector is defined as the final link's end point.
The class computes the Lagrangian, Mass Matrix, Coriolis Matrix, Gravity Vector, and Equations of Motion for the robot using Sage Math.
The user can define joint trajectories q_i(t) and evaluate the torques required to move the robot along the trajectory.
'''
class SLIRP():
    def __init__(self, num_links):
        self.t = var('t')
        self.num_links = num_links
        self.g = var('g')  # Gravity
        self.w = var('W')  

        self.link_lengths = []
        self.link_masses = []
        self.link_inertias = []
        self.thetas = []    # theta_i(t)
        self.qs = []        # q_i (symbolic)
        self.qdots = []     # q_i_dot (symbolic)
        self.qddots = []    # q_i_ddot (symbolic) - new addition
        self.link_COMs = []
        self.end_effector = []

        self.theta_to_q_only = {}
        
        self.T = None #kinetic Energy 
        self.U = None #Potential Energy
        self.L = None #Langrangian
        self.tau = None #Torque
        self.tau_matrix = None #Torque Matrix
        self.J = None #Jacobian


        #list for derivatives
        self.DtDqdotT = [] #partial derivatives of dT/dqdot_i
        self.DqT = [] #partial derivatives of dT/dq_i

        self.M = None #Mass Matrix
        self.C1 = None #Coriolis Matrix
        self.C2 = None #Coriolis Matrix

        # Trajectories (to be set by user)
        self.q_trajectories = []  # joint trajectories q_i(t)
        self.ee_trajectory = []   # end-effector trajectory if needed

        self.establish_link_variables()
        self.calculate_COM_positions()
        self.define_end_effector()
        self.compute_jacobian()
    '''
    Function Name: establish_link_variables
    Inputs: None
    Outputs: None
    Description: This function initializes the symbolic variables for the link lengths, masses, inertias, joint angles, and joint velocities.
    '''
    def establish_link_variables(self):
        # 1-based indexing for convenience
        for i in range(1, self.num_links+1):
            L_i = var(f'L{i}')
            m_i = var(f'm{i}')
            Iz_i = m_i/SR(12)*(L_i**2 + self.w**2) 
            self.link_lengths.append(L_i)
            self.link_masses.append(m_i)
            self.link_inertias.append(Iz_i)

            theta_i = function(f'theta{i}')(self.t)
            self.thetas.append(theta_i)

            # Now define q_i, q_i_dot, q_i_ddot as symbolic variables
            q_i = function(f'q{i}')(self.t)
            # For partial derivatives, also define symbolic variables to substitute later:
            q_sym = var(f'q{i}')
            qdot_sym = var(f'q{i}dot')
            qddot_sym = var(f'q{i}ddot')

            # store q_syms and qdot_syms for final evaluation
            # For forming equations, use q(t), diff(q(t), t), diff(q(t), t,2)
            self.qs.append(q_sym)
            self.qdots.append(qdot_sym)
            self.qddots.append(qddot_sym)

        for i in range(self.num_links):
            self.theta_to_q_only[self.thetas[i]] = self.qs[i]
    '''
    Function Name: calculate_COM_positions
    Inputs: None
    Outputs: None
    Description: This function calculates the position of the center of mass for each link in the robot arm.
    '''
    def calculate_COM_positions(self):
        x_prev, y_prev = 0, 0
        theta_sum = 0
        for i in range(self.num_links):
            theta_sum += self.thetas[i]
            L_i = self.link_lengths[i]
            x_COM = x_prev + (L_i/2)*cos(theta_sum)
            y_COM = y_prev + (L_i/2)*sin(theta_sum)
            self.link_COMs.append([x_COM, y_COM])
            #Stre previous link to add to next links COM
            x_prev += L_i*cos(theta_sum)
            y_prev += L_i*sin(theta_sum)
    '''
    Function Name: define_end_effector
    Inputs: None
    Outputs: None
    Description: This function calculates the position of the end effector of the robot arm.
    '''
    def define_end_effector(self):
        theta_sum = 0
        x_end, y_end = 0, 0
        for i in range(self.num_links):
            L_i = self.link_lengths[i]
            theta_sum += self.thetas[i]
            x_end += L_i*cos(theta_sum)
            y_end += L_i*sin(theta_sum)
        self.end_effector = [x_end, y_end, theta_sum]
    '''
    Function Name: compute_jacobian
    Inputs: None
    Outputs: None
    Description: This function computes the Jacobian matrix for the robot arm.
    '''
    def compute_jacobian(self):
        X_q = [expr.subs(self.theta_to_q_only) for expr in self.end_effector]
        J = Matrix([[diff(expr, q) for q in self.qs] for expr in X_q])
        self.J = J.simplify()
    '''
    Function Name: compute_lagrangian
    Inputs: None
    Outputs: None
    Description: This function computes the Lagrangian for the robot arm.
    '''
    def compute_lagrangian(self):
        T = 0
        U = 0
        qdot_subs = {}
        qddot_subs = {}
        # Map diff(theta_i(t), t) to q_i_dot symbolic

        for i in range(self.num_links):
            qdot_subs[diff(self.thetas[i], self.t)] = self.qdots[i]
            qddot_subs[diff(self.thetas[i], self.t, 2)] = self.qddots[i]

        for i in range(self.num_links):
            m_i = self.link_masses[i]
            I_i = self.link_inertias[i]
            x_com_i, y_com_i = self.link_COMs[i]

            dx_com_i = diff(x_com_i, self.t)
            dy_com_i = diff(y_com_i, self.t)

            dx_ci_q = dx_com_i.subs(self.theta_to_q_only).subs(qdot_subs)
            dy_ci_q = dy_com_i.subs(self.theta_to_q_only).subs(qdot_subs)

            theta_sum_i = sum(self.thetas[:i+1])
            dtheta_sum_i = diff(theta_sum_i, self.t).subs(qdot_subs)
            # Kinetic Energy - Translational and Rotational
            T += (1/SR(2))*m_i*(dx_ci_q**2 + dy_ci_q**2) + (1/SR(2))*I_i*(dtheta_sum_i**2)

            y_com_i_q = y_com_i.subs(self.theta_to_q_only)
            U += m_i*self.g*y_com_i_q

        self.T = T
        self.U = U
        self.L = (T - U).simplify()


        self.compute_DtDqdotT()
        self.compute_DqT()
        return self.L

    '''
    Function Name: compute_DtDqdotT
    Inputs: None
    Outputs: None
    Description: This function computes the partial derivatives of dT/dqdot_i for each i.
    '''
    def compute_DtDqdotT(self):
        self.DtDqdotT = []
        for i, qdot_i in enumerate(self.qdots):
            # Compute dT/dqdot_i
            dT_dqdot_i = diff(self.T, qdot_i).simplify()
            
            # Compute the partial derivatives
            d2T_dqj_dqdot_i = [diff(dT_dqdot_i, qj).simplify() for qj in self.qs]
            d2T_dqdotj_dqdot_i = [diff(dT_dqdot_i, qdot_j).simplify() for qdot_j in self.qdots]
            
            # Store the partial derivatives
            self.DtDqdotT.append((d2T_dqj_dqdot_i, d2T_dqdotj_dqdot_i))
            
            # # Debugging: Print partial derivatives
            # print(f"DtDqdotT[{i}] - d2T/dqj_dqdot_i:")
            # for j, term in enumerate(d2T_dqj_dqdot_i):
            #     print(f"  d2T/dq{j+1}_dqdot_{i+1} = {term}")
            # print(f"DtDqdotT[{i}] - d2T/dqdotj_dqdot_i:")
            # for j, term in enumerate(d2T_dqdotj_dqdot_i):
            #     print(f"  d2T/dqdot{j+1}_dqdot_{i+1} = {term}")
    '''
    Function Name: compute_DqT
    Inputs: None
    Outputs: None
    Description: This function computes the partial derivatives of dT/dq_i for each i.
    '''
    def compute_DqT(self):
        self.DqT = []
        for i, q_i in enumerate(self.qs):
            # Compute dT/dq_i
            dT_dq_i = diff(self.T, q_i).simplify()
            
            # Compute the partial derivatives
            d2T_dqj_dqi = [diff(dT_dq_i, qj).simplify() for qj in self.qs]
            d2T_dqdotj_dqi = [diff(dT_dq_i, qdot_j).simplify() for qdot_j in self.qdots]
            
            # Store the partial derivatives
            self.DqT.append((d2T_dqj_dqi, d2T_dqdotj_dqi))
            
            # # Debugging: Print partial derivatives
            # print(f"DqT[{i}] - d2T/dqj_dqi:")
            # for j, term in enumerate(d2T_dqj_dqi):
            #     print(f"  d2T/dq{j+1}_dq{self.num_links+1} = {term}")
            # print(f"DqT[{i}] - d2T/dqdotj_dqi:")
            # for j, term in enumerate(d2T_dqdotj_dqi):
            #     print(f"  d2T/dqdot{j+1}_dq{self.num_links+1} = {term}")


    '''
    Function Name: compute_mass_matrix
    Inputs: None
    Outputs: None
    Description: This function computes the mass matrix for the robot arm.
                    This function must be called prior to calculating tau
    '''
    def compute_mass_matrix(self):
        if not hasattr(self, 'T'):
            if not hasattr(self, 'U'):
                raise ValueError("U is not defined. Make sure to compute_lagrangian() and store T and U.")
            self.T = (self.L + self.U).simplify()

        qdot = self.qdots
        n = self.num_links
        M = [[diff(diff(self.T, qdot[i]), qdot[j]).simplify() for j in range(n)] for i in range(n)]
        self.M = Matrix(SR, M)
        return self.M

    '''
    Function Name: compute_coriolis_terms
    Inputs: None
    Outputs: None
    Description: This function computes the Coriolis matrix for the robot arm.
                 C_ij = sum_k (d2T/dqj_dqdot_i * qdot_j + d2T/dqdotj_dqdot_i * qddot_j - d2T/dqdotj_dqi * qdot_j - d2T/dqj_dqi * qdot_j)
    '''
    def compute_coriolis_terms(self):
        if not self.DtDqdotT or not self.DqT:
            raise ValueError("DtDqdotT and DqT must be computed before computing Coriolis matrix.")
        
        n = self.num_links
        # Initialize matrix over the Symbolic Ring [DO NOT CHANGE, this is required for symbolic matrix operations]
        C = Matrix(SR, n, n, 0)
        
        for i in range(n):
            for j in range(n):
                # Sum over k
                C_ij = 0
                for k in range(n):
                    # From d/dt(dT/dqdot_i)
                     # NOTE: This was the only way I could get the symbolic math to work with the coriolus terms [DO NOT CHANGE WITHOUT EXTENSIVE TESTING]
                    term1 = self.DtDqdotT[i][0][j] * self.qdots[k]  # d2T/dqj_dqdot_i * qdot_j
                    term2 = self.DtDqdotT[i][1][j] * self.qddots[k]  # d2T/dqdotj_dqdot_i * qddot_j
                    
                    # From d/dt(dT/dq_i)
                    term3 = self.DqT[i][1][j] * self.qdots[k]  # d2T/dqdotj_dqi * qdot_j
                    term4 = self.DqT[i][0][j] * self.qdots[k]  # d2T/dqj_dqi * qdot_j
                    C_ij += term1 + term2 - term3 - term4
                C[i, j] = C_ij.simplify()
        
        self.C = C
        
        # Debugging: Print Coriolis matrix dimensions and contents
        
        #print("Coriolis Matrix (C):")
        #print(self.C)
        #print("Coriolis Matrix (C) dimensions:", self.C.dimensions())


    '''
    Function Name: compute_gravity_vector
    Inputs: None
    Outputs: None
    Description: This function computes the gravity vector for the robot arm.
                 G(q) = dU/dq
    '''
    def compute_gravity_vector(self):
        if not hasattr(self, 'U'):
            # U must be defined. 
            raise ValueError("U is not defined. Make sure to compute_lagrangian() and store U.")
        q = self.qs
        n = self.num_links
        G = [diff(self.U, q[i]).simplify() for i in range(n)]
        self.G = Matrix(SR, G).transpose().simplify()
        return self.G
    '''
    Function Name: compute_eom
    Inputs: None
    Outputs: None
    Description: This function computes the equations of motion for the robot arm.
                 tau = M*qddot + C*qdot + G. 
    NOTE: This function requires that compute_mass_matrix(), compute_coriolis_terms(), and compute_gravity_vector() have been called.
          AND that they are all properly implemented as symbolic matrices. Otherwise, this function will raise an error.
    '''
    def compute_eom(self):
        if not hasattr(self, 'M'):
            self.compute_mass_matrix()
        if not hasattr(self, 'C'):
            self.compute_coriolis_terms()
        if not hasattr(self, 'G'):
            self.compute_gravity_vector()

        # Convert qddots and qdots from lists to column vectors over SR
        qddot_vec = Matrix(SR, self.num_links, 1, self.qddots)
        qdot_vec = Matrix(SR, self.num_links, 1, self.qdots)

        # Debugging: Check matrix dimensions before multiplication
        # print("Mass Matrix (M) dimensions:", self.M.dimensions())
        # print("Coriolis Matrix (C) dimensions:", self.C.dimensions())
        # print("Gravity Vector (G) dimensions:", self.G.dimensions())
        # print("qddot_vec dimensions:", qddot_vec.dimensions())
        # print("qdot_vec dimensions:", qdot_vec.dimensions())

        # Now perform the matrix multiplications
        try:
            self.tau_matrix = self.M * qddot_vec + self.C * qdot_vec + self.G
        except Exception as e:
            print("Error during matrix multiplication:")
            print(f"M * qddot_vec:\n{self.M * qddot_vec}")
            print(f"C * qdot_vec:\n{self.C * qdot_vec}")
            print(f"G:\n{self.G}")
            raise e

        # Convert to list for easier handling downstream (Matrix form is stored in self.tau_matrix)
        self.tau = [self.tau_matrix[i,0].simplify() for i in range(self.num_links)]
        #print("Torque Vector (tau):")
        #for i, tau_i in enumerate(self.tau):
            #print(f"tau_{i+1} = {tau_i}")
        return self.tau

    
    '''
    Function Name: evaluate_tau
    Inputs: values_dict - Dictionary of values to substitute into the equations of motion (link lengths, masses, inertias, gravity, etc.)
    Outputs: tau_evaluated - List of evaluated torques at the given values
    Description: This function evaluates the torques at the given values.
    '''
    def evaluate_tau(self, values_dict):
        if self.tau is None:
            raise ValueError("Must compute equations of motion first.")
        tau_sub = [expr.subs(values_dict) for expr in self.tau]
        tau_evaluated = [ts.n() for ts in tau_sub]
        return tau_evaluated
    '''
    Function Name: set_trajectory
    Inputs: q_trajectories - List of symbolic expressions for the joint trajectories q_i(t)
    Outputs: None
    Description: This function sets the joint trajectories for the robot arm.
    '''
    def set_trajectory(self, q_trajectories):
        if len(q_trajectories) != self.num_links:
            raise ValueError("Number of trajectories must match num_links.")
        self.q_trajectories = q_trajectories
    '''
    Function Name: get_state_at_time
    Inputs: t_val - Time value to evaluate the joint trajectories at
    Outputs: q_values - List of joint positions at time t_val
             qdot_values - List of joint velocities at time t_val
             qddot_values - List of joint accelerations at time t_val
    Description: This function evaluates the joint positions, velocities, and accelerations at the given time t_val.
    '''
    def get_state_at_time(self, t_val):
        subs_dict = {self.t: t_val}
        q_values = [q_expr.subs(subs_dict) for q_expr in self.q_trajectories]
        qdot_values = [diff(q_expr, self.t).subs(subs_dict) for q_expr in self.q_trajectories]
        qddot_values = [diff(q_expr, self.t,2).subs(subs_dict) for q_expr in self.q_trajectories]
        return q_values, qdot_values, qddot_values
    '''
    Function Name: evaluate_tau_at_time
    Inputs: t_val - Time value to evaluate the torques at
            param_values - Dictionary of parameter values to substitute into the equations of motion
    Outputs: tau_at_t - List of evaluated torques at time t_val
    Description: This function evaluates the torques at the given time t_val and parameter values.
                 This function requires the helper function get_state_at_time() to evaluate the joint positions, velocities, and accelerations.
    '''
    def evaluate_tau_at_time(self, t_val, param_values):
        q_values, qdot_values, qddot_values = self.get_state_at_time(t_val)
        subs_dict = param_values.copy()
        for i in range(self.num_links):
            subs_dict[self.qs[i]] = q_values[i]
            subs_dict[self.qdots[i]] = qdot_values[i]
            subs_dict[self.qddots[i]] = qddot_values[i]
        return self.evaluate_tau(subs_dict)
    '''
    Function Name: set_end_effector_trajectory
    Inputs: ee_trajectory - List of symbolic expressions for the end-effector trajectory
    Outputs: None
    Description: This function sets the end-effector trajectory for the robot arm.
    '''
    def set_end_effector_trajectory(self, ee_trajectory):
        self.ee_trajectory = ee_trajectory
    '''
    Function Name: get_end_effector_state_at_time
    Inputs: t_val - Time value to evaluate the end-effector trajectory at
    Outputs: ee_values - List of end-effector positions at time t_val
             ee_dot_values - List of end-effector velocities at time t_val
             ee_ddot_values - List of end-effector accelerations at time t_val'''
    def compute_joint_trajectory_from_ee(self):
        #TODO: Implement IK for the robot arm
        pass
    '''
    Function Name: print_matrix_as_python
    Inputs: matrix_var - Matrix to print
            matrix_name - Name of the matrix
            Outputs: None
    Description: This function formats and prints a matrix in Python nested list format.
                 This function was made by Dr. Swenson and is used for debugging purposes and easy printing of matrices.'''
    # Function to format and print a matrix in Python nested list format (TAKEN FROM DR. SWENSONS RRR_PLANAR.ipynb)
    def print_matrix_as_python(self, matrix_var, matrix_name):
        matrix_python = [[str(matrix_var[i, j]) for j in range(matrix_var.ncols())] for i in range(matrix_var.nrows())]
        print(f"{matrix_name} = [")
        for row in matrix_python:
            print("    [" + ", ".join(row) + "],")
        print("]")
