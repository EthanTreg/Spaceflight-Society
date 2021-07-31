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
    lower_length_m = 0.1
    upper_angle_limit_rad = [-np.pi, 0] #Not implimented yet
    lower_angle_limit_rad = [0, np.pi / 2] #Not implimented yet
    lower_angle_rad = 0 #Initial angle - Not implimented yet
    upper_angle_rad = np.pi / 2 #Initial angle - Not implimented yet
    y_solutions = [0, 0]
    joint_position_m = [0, lower_length_m]
    head_position_m = head_target_position_m = np.array((upper_length_m, lower_length_m))
    arm_square_difference_m = lower_length_m ** 2 - upper_length_m ** 2

    #Arm visulisation
    plt.figure()
    plt.plot((0, joint_position_m[0], head_target_position_m[0]), (0, joint_position_m[1], head_target_position_m[1]))
    plt.xlim(0, upper_length_m + lower_length_m)
    plt.ylim(0, upper_length_m + lower_length_m)
    plt.show(block=False)

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

        #If the head target location is out of reach, set the coordinates to the maximum range in same direction
        if head_target_distance_m > upper_length_m + lower_length_m:
            head_target_position_m *= (upper_length_m + lower_length_m) / head_target_distance_m
            head_target_distance_m = upper_length_m + lower_length_m

        #Calculate first y solution
        y_solutions[0] = y_solutions[1] = 0.5 * head_target_position_m[1] * (1 + arm_square_difference_m / (head_target_distance_m ** 2))
        y_solutions[0] += np.sqrt(abs(
            - head_target_distance_m ** 6 + head_target_distance_m ** 4 * (head_target_position_m[1] ** 2 - 2 * arm_square_difference_m)
            + head_target_distance_m ** 2 * (2 * head_target_position_m[1] ** 2 * arm_square_difference_m - arm_square_difference_m ** 2 + 4 * head_target_position_m[0] ** 2 * lower_length_m ** 2)
            + head_target_position_m[1] ** 2 * arm_square_difference_m ** 2)) / (2 * head_target_distance_m ** 2)

        #Calculate second y solution
        y_solutions[1] -= np.sqrt(abs(
            - head_target_distance_m ** 6 + head_target_distance_m ** 4 * (head_target_position_m[1] ** 2 - 2 * arm_square_difference_m)
            + head_target_distance_m ** 2 * (2 * head_target_position_m[1] ** 2 * arm_square_difference_m - arm_square_difference_m ** 2 + 4 * head_target_position_m[0] ** 2 * lower_length_m ** 2)
            + head_target_position_m[1] ** 2 * arm_square_difference_m ** 2)) / (2 * head_target_distance_m ** 2)

        #Take the y solution colsest to the current y coordinate
        if abs(y_solutions[0] - joint_position_m[1]) < abs(y_solutions[1] - joint_position_m[1]):
            joint_position_m[1] = y_solutions[0]
            joint_position_m[0] = np.sqrt(lower_length_m ** 2 - y_solutions[0] ** 2)

        else:
            joint_position_m[1] = y_solutions[1]
            joint_position_m[0] = np.sqrt(lower_length_m ** 2 - y_solutions[1] ** 2)

        #Plots new arm
        plt.clf()
        plt.plot((0, joint_position_m[0], head_target_position_m[0]), (0, joint_position_m[1], head_target_position_m[1]))
        plt.xlim(0, upper_length_m + lower_length_m)
        plt.ylim(0, upper_length_m + lower_length_m)
        plt.show(block=False)

    return

main()