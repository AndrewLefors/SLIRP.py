# AUTHORS: Andrew J. Lefors and Dr. John Swensen
# Contact: andrew.lefors@wsu.edu

# SLIRP.py (Serial Link Integrated Robot Physics) 
[Pronounced "SLUR-PEE"]

## Overview

SLIRP models a serial-link manipulator with an arbitrary number of revolute joints. Its primary capabilities include:

1. Symbolic definition of kinematic variables.
2. Forward kinematics for pose (position/orientation), velocity, and acceleration.
3. Computation of link center of mass (COM) positions.
4. End-effector position and orientation (pose).
5. Substitution between joint angles (θ) and generalized coordinates (q).
6. Computation of the Jacobian matrix.
7. Computation of kinetic and potential energy, and inclusion of non-conservative forces.
8. Incorporation of constraints using Lagrange multipliers (stretch goal).
9. Output of Equations of Motion (EOM) in generalized or workspace coordinates.

SLIRP serves as a foundation for advanced robot analysis, simulation, and control development.

---

## Key Design Goals

1. **Arbitrary Number of Links**: Dynamically handle any number of links.
2. **Symbolic Computation**: Use Sage Math for symbolic variables for analytical solutions.
3. **Forward Kinematics**: Compute end-effector pose, velocity, and acceleration from joint coordinates and derivatives.
4. **Energy and Forces**: Compute kinetic and potential energy.
5. **Constraints**: Symbolically incorporate constraints into the system’s equations.
6. **Multi-Coordinate Output**: Generate EOM in joint and workspace coordinates.
7. **User-Friendliness**: Provide clear documentation and extensibility.

---

## Class Structure and Responsibilities

### Initialization

**`__init__(self, num_links)`**

- **Inputs**: `num_links` (integer)
- **Actions**:
  1. Store `num_links`.
  2. Define symbolic time variable.
  3. Initialize link lengths, masses, angles, and generalized coordinates.
  4. Compute COMs, end-effector, substitution dictionaries, and Jacobian.

### Key Attributes

- `self.num_links`: Number of links.
- `self.time`: Symbolic time variable.
- `self.link_lengths, self.link_masses`: Lists of symbolic variables.
- `self.thetas, self.qs, self.qdots`: Joint angles and generalized coordinates.
- `self.link_COMs`: Positions of link COMs.
- `self.end_effector`: [x, y, θ] or [x, y].
- `self.J`: Jacobian matrix.

---

## Methods Overview

### Core Functions

- **`establish_link_variables()`**: Define symbolic variables for links.
- **`calculate_COM_positions()`**: Compute link COM positions.
- **`define_end_effector()`**: Compute the end-effector pose.
- **`compute_velocities()`**: Compute symbolic end-effector velocities.
- **`compute_accelerations()`**: Compute symbolic end-effector accelerations.
- **`compute_jacobian()`**: Compute the Jacobian matrix.
- **`compute_kinetic_energy()`**: Compute total kinetic energy.
- **`compute_potential_energy()`**: Compute gravitational potential energy.
- **`derive_EOM_in_generalized_coordinates()`**: Derive EOM in generalized coordinates.
- **`transform_EOM_to_workspace_coordinates()`**: Transform EOM to workspace coordinates.

### Advanced Features

- Constraints with Lagrange Multipliers (stretch goal).
- Substitution between θ and q variables for flexibility.
- Extensions for inertia, non-conservative forces, and redundant mechanisms.

---

## User Manual

### Installation

**Requirements**: Python environment with Sage Math.

### Basic Usage

1. **Initialize the Robot**:
   ```python
   from SLIRP import SLIRP
   robot = SLIRP(num_links=3)
   ```

2. **Set Link Lengths and Masses**:
   ```python
   robot.link_lengths = [1.0, 0.5, 0.75]
   robot.link_masses = [2.0, 1.5, 1.0]
   ```

3. **Compute End-Effector Pose**:
   ```python
   pose = robot.end_effector  # [x_end, y_end, theta_end]
   ```

4. **Compute Jacobian**:
   ```python
   J = robot.J
   ```

5. **Velocity/Acceleration**:
   ```python
   v_ee = robot.compute_velocities()
   a_ee = robot.compute_accelerations()
   ```

6. **Energy Calculations**:
   ```python
   T = robot.compute_kinetic_energy()
   U = robot.compute_potential_energy()
   ```

7. **Equations of Motion**:
   ```python
   EOM_q = robot.derive_EOM_in_generalized_coordinates()
   EOM_workspace = robot.transform_EOM_to_workspace_coordinates(EOM_q)
   ```

### Extending the Code

- **Add Inertial Properties**: Compute inertia matrices for kinetic energy.
- **Include Control Inputs**: Add joint torques or external forces.
- **Numerical Evaluation**: Substitute numeric values for simulation.

---

SLIRP empowers users to build, analyze, and simulate robotic systems symbolically and numerically.
