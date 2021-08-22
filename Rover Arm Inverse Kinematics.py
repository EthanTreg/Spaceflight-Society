import matplotlib.pyplot as plt
import numpy as np

def main():
    """
        Controls:
        x direction:
            w - Forwards
            s - Backwards
        y direction:
            e - Up
            q - Down
    """

    #Constants
    upper_length_m = 0.1
    lower_length_m = 0.2
    upper_angle_limit_rad = [-np.pi, 0] #Not implimented yet
    lower_angle_limit_rad = [0, np.pi / 2] #Not implimented yet
    lower_angle_rad = 0 #Initial angle - Not implimented yet
    upper_angle_rad = np.pi / 2 #Initial angle - Not implimented yet
    y_solutions = np.array((0., 0.))
    joint_position_m = [0, lower_length_m]
    head_target_position_m = np.array((upper_length_m, lower_length_m))
    arm_square_difference_m = lower_length_m ** 2 - upper_length_m ** 2

    plot_arm(joint_position_m, head_target_position_m, upper_length_m, lower_length_m)

    while True:
        direction = input()

        if direction.lower() == 'w':
            head_target_position_m[0] += 0.01

        elif direction.lower() == 's':
            head_target_position_m[0] -= 0.01

        elif direction.lower() == 'e':
            head_target_position_m[1] += 0.01

        elif direction.lower() == 'q':
            head_target_position_m[1] -= 0.01

        head_target_distance_m = np.sqrt(np.sum(head_target_position_m ** 2))

        #If the head target location is out of reach, set the coordinates to the maximum/minimum range in same direction
        if head_target_distance_m > upper_length_m + lower_length_m:
            head_target_position_m *= (upper_length_m + lower_length_m) / head_target_distance_m
            head_target_distance_m = upper_length_m + lower_length_m
        elif head_target_distance_m < lower_length_m - upper_length_m:
            head_target_position_m *= (lower_length_m - upper_length_m) / head_target_distance_m
            head_target_distance_m = lower_length_m - upper_length_m

        #Calculate y solutions
        y_weighted_mid_point_m = 0.5 * head_target_position_m[1] * (1 + arm_square_difference_m / (head_target_distance_m ** 2))
        y_mid_point_offset_m = np.sqrt(abs(
            - head_target_distance_m ** 6 + head_target_distance_m ** 4 * (head_target_position_m[1] ** 2 - 2 * arm_square_difference_m)
            + head_target_distance_m ** 2 * (2 * head_target_position_m[1] ** 2 * arm_square_difference_m - arm_square_difference_m ** 2 + 4 * head_target_position_m[0] ** 2 * lower_length_m ** 2)
            + head_target_position_m[1] ** 2 * arm_square_difference_m ** 2)) / (2 * head_target_distance_m ** 2)

        y_solutions[0] = y_weighted_mid_point_m + y_mid_point_offset_m
        y_solutions[1] = y_weighted_mid_point_m - y_mid_point_offset_m

        #Calculate correct x & y solutions
        x_solutions = np.sqrt(np.where((lower_length_m ** 2 - y_solutions ** 2 < 0) & (abs(lower_length_m ** 2 - y_solutions ** 2) < 1e-3), 0, lower_length_m ** 2 - y_solutions ** 2))
        x_solutions = np.vstack((x_solutions, -x_solutions))
        solution_arm_lengths_m = np.sqrt((head_target_position_m[1] - y_solutions) ** 2 + (head_target_position_m[0] - x_solutions) ** 2)
        solution_idx = np.argwhere(abs(solution_arm_lengths_m - upper_length_m) < 1e-4)
        solutions = np.transpose(np.vstack((x_solutions[solution_idx[:, 0], solution_idx[:, 1]], y_solutions[solution_idx[:, 1]])))

        #Find correct solution closest to current position
        solution_distance_m = np.sqrt((joint_position_m[0] - solutions[:, 0]) ** 2 + (joint_position_m[1] - solutions[:, 1]) ** 2)
        joint_position_m = solutions[np.argmin(solution_distance_m)]

        plot_arm(joint_position_m, head_target_position_m, upper_length_m, lower_length_m)

def plot_arm(joint_position_m, head_target_position_m, upper_length_m, lower_length_m):
    """
    Arm visualisation.
    ### Params:
        joint_position_m: x, y position of elbow pivot
        head_target_position_m: x, y position of claw
        upper_length_m: length of forearm
        lower_length_m: length of aftarm
    """
    plt.clf()
    ax = plt.gca()
    plt.plot((0, joint_position_m[0], head_target_position_m[0]), (0, joint_position_m[1], head_target_position_m[1]))
    plt.xlim(-upper_length_m - lower_length_m, upper_length_m + lower_length_m)
    plt.ylim(-upper_length_m - lower_length_m, upper_length_m + lower_length_m)
    ax.set_aspect('equal')
    plt.show(block=False)

main()
