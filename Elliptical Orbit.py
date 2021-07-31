import numpy as np
from numpy.core.numeric import full
import scipy.constants as const
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib


def elliptical_orbit(x):
    """Equation for satellite's orbit.
        ### Params:
            x: x-coordinate
        ### Return:
            y: y-coordinate corresponding to x-coordinate along elliptical orbit
    """
    return semiminor_axis_m * (1 - (x + semimajor_axis_m * eccentricity) ** 2 / (semimajor_axis_m ** 2)) ** (1/2)

def quadratic_formula(a, b, c):
    """Calculates the roots of a quadratic equation ax^2+bx+c
        ### Params:
            a: Coefficient of x^2 term
            b: Coefficient of x term
            c: Constant
        ### Return:
            List of two quadratic roots
    """
    return (-b + [np.sqrt(b ** 2 - 4 * a * c), -np.sqrt(b ** 2 - 4 * a * c)]) / (2 * a)

def var_orbital_radius_m(orbital_true_anomaly_rad):
    """Calculates the orbital radius of the satellite from a given orbital true anomaly.
        ### Params:
            orbital_true_anomaly_rad: Angle of orbit relative to periareion
        ### Return:
            Radius of satellite at the given true anomaly
    """
    return semi_latus_rectum_m / (1 + eccentricity * np.cos(orbital_true_anomaly_rad))

def var_angular_frequency_rads(orbital_true_anomaly_rad):
    """Calculates the angular frequency at a given orbital position
        ### Params:
            orbital_true_anomaly_rad: True anomaly position of satellite
        ### Return:
            Angular frequency at orbital position
    """
    return radius_periareion_m * velocity_periareion_ms / (var_orbital_radius_m(orbital_true_anomaly_rad) ** 2)

def var_delta_connection_period_s(orbital_true_anomaly_rad):
    """Calculates the time required to move dθ at a given orbital position.
        ### Params:
            orbital_true_anomaly_rad: True anomaly of satellite
        ### Return:
            Time dt to move dθ at position θ
    """
    return 1 / var_angular_frequency_rads(orbital_true_anomaly_rad)

def var_slant_range_m(orbital_true_anomaly_rad, rover_true_anomaly_rad):
    """Slant range of satellite from rover at a given orbital true anomaly and rover true anomaly.
        ### Params:
            orbital_true_anomaly_rad: True anomaly of satellite
            rover_true_anomaly_rad: True anomaly of rover
        ### Return:
            Distance of satellite from rover
    """
    orbit_radius_m = var_orbital_radius_m(orbital_true_anomaly_rad)
    return np.sqrt(mars_radius_m ** 2 + orbit_radius_m ** 2 - 2 * mars_radius_m * orbit_radius_m * np.cos(orbital_true_anomaly_rad - rover_true_anomaly_rad))

def var_delta_free_space_loss_db(orbital_true_anomaly_rad, rover_true_anomaly_rad):
    """Calculates the small free space loss for a small angle dθ at a given orbital true anomaly and rover true anomaly.
        ### Params:
            orbital_true_anomaly_rad: True anomaly of satellite
            rover_true_anomaly_rad: True anomaly of rover
        ### Return:
            Small free space loss at given true anomaly of orbiter and rover
    """
    return -10 * np.log10((4 * const.pi * var_slant_range_m(orbital_true_anomaly_rad, rover_true_anomaly_rad) / wavelength_m) ** 2) / var_angular_frequency_rads(orbital_true_anomaly_rad)

def var_eccentric_anomaly_rad(true_anomaly_rad):
    """Calculates the eccentric anomaly from true anomaly.
        ### Params:
            true_anomaly_rad: True anomaly to be converted to eccentric anomaly
        ### Return:
            Eccentric anomaly equivalent of true anomaly
    """
    return 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(true_anomaly_rad / 2))

def var_orbit_time_s(eccentric_anomaly_rad, initial_time=0):
    """Calculates the difference in time to reach the given eccentric anomaly from an initial time of -periareion.
        ### Params:
            eccentric_anomaly_rad: Eccentric anomaly of satellite position
            initial_time: Negative time of periareion, default = 0
        ### Return:
            Time difference of eccentric anomaly position relative to the negative time of periareion
    """
    return (eccentric_anomaly_rad - eccentricity * np.sin(eccentric_anomaly_rad)) / mean_motion_rads - initial_time

def main():
    """Main function to calculate the free space loss at different orbiter positions.
    """

    rover_true_anomalies_rad = [0, const.pi]
    labels = ['Periareion window', 'Apoareion window']
    target_window_s = 600

    #Calculate free space loss for predefined rover true anomalies
    for i in range(len(rover_true_anomalies_rad)):
        comm_time_s, comm_free_space_loss_db = comm_time_free_space_calculation(rover_true_anomalies_rad[i], labels[i], plot=True)

        #If connection time exceeds target window time, limit the comm window to target window
        if comm_time_s > target_window_s:
            center_true_anomaly_rad = const.pi
            initial_time = (var_eccentric_anomaly_rad(center_true_anomaly_rad) - \
                eccentricity * np.sin(var_eccentric_anomaly_rad(center_true_anomaly_rad))) / \
                    mean_motion_rads + target_window_s / 2

            eccentric_anomaly_rad = optimize.fsolve(
                lambda eccentric_anomaly_rad: var_orbit_time_s(eccentric_anomaly_rad, initial_time), const.pi)[0]

            orbital_true_anomaly_rad = abs(2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(eccentric_anomaly_rad / 2)))
            orbital_true_anomaly_rad = [orbital_true_anomaly_rad, 2 * const.pi - orbital_true_anomaly_rad]

            orbital_radius_m = var_orbital_radius_m(orbital_true_anomaly_rad)
            comm_x_m = orbital_radius_m * np.cos(orbital_true_anomaly_rad)
            comm_y_m = orbital_radius_m * np.sin(orbital_true_anomaly_rad)

            comm_free_space_loss_db = integrate.quad(
                lambda orbital_true_anomaly_rad: var_delta_free_space_loss_db(orbital_true_anomaly_rad, const.pi),
                orbital_true_anomaly_rad[0], orbital_true_anomaly_rad[1])[0] / \
                    integrate.quad(var_delta_connection_period_s, orbital_true_anomaly_rad[0], orbital_true_anomaly_rad[1])[0]

            plot_orbit_cartesian(-mars_radius_m, 0, comm_x_m, comm_y_m, '600 s ' + labels[i].lower())

        print('\nRover true anomaly: {} deg\nConnection duration: {:.1f} s\nFree space loss: {:.1f} db'.format(
            rover_true_anomalies_rad[i] * 180 / const.pi, comm_time_s, comm_free_space_loss_db))

    #Calculate position in orbit with window length of the target window length
    rover_true_anomaly_rad = optimize.fsolve(
        lambda rover_true_anomaly_rad: comm_time_free_space_calculation(rover_true_anomaly_rad, '', target_window_s, full_return=False), 0.5)[0]

    comm_time_s, target_free_space_loss_db = comm_time_free_space_calculation(rover_true_anomaly_rad, '600 s window', plot=True)

    print('\nTarget Orbit:\nConnection duration: {:.1f} s\nRover angle: {:.1f} deg\nFree space loss: {:.1f} db'.format(
        comm_time_s, rover_true_anomaly_rad * 180 / const.pi, target_free_space_loss_db))

def comm_time_free_space_calculation(rover_true_anomaly_rad, label, target_time_s=0, full_return=True, plot=False):
    """Calculate connection time and free space loss at a given rover true anomaly.
        ### Params:
            rover_true_anomaly_rad: True anomaly of rover
            label: Label of comm connection being plotted
            target_time_s: target window length, default = 0
            full_return: If function should return free space loss and comm angles, default = True
            plot: If the rover comm window should be plotted, default = False
        ### Return:
            full_return = True:
                Difference in window length and target window length
                Average free space loss for window
            full_return = False:
                Difference in window length and target window length
    """
    rover_true_anomaly_rad = float(rover_true_anomaly_rad)

    comm_intersection_xy_m = comm_range_calculation(rover_true_anomaly_rad, label, plot=plot)
    comm_true_anomaly_rad = np.arctan(comm_intersection_xy_m[1, :] / comm_intersection_xy_m[0, :])
    comm_true_anomaly_rad = np.where(comm_intersection_xy_m[0] < 0, const.pi + comm_true_anomaly_rad, comm_true_anomaly_rad)
    comm_time_s = abs(integrate.quad(var_delta_connection_period_s, comm_true_anomaly_rad[0], comm_true_anomaly_rad[1])[0])
    comm_free_space_loss_db = integrate.quad(
        lambda orbital_true_anomaly_rad: var_delta_free_space_loss_db(orbital_true_anomaly_rad, rover_true_anomaly_rad),
        comm_true_anomaly_rad[0], comm_true_anomaly_rad[1])[0] / integrate.quad(var_delta_connection_period_s, comm_true_anomaly_rad[0], comm_true_anomaly_rad[1])[0]

    if full_return:
        return comm_time_s - target_time_s, comm_free_space_loss_db
    else:
        return comm_time_s - target_time_s

def comm_range_calculation(rover_true_anomaly_rad, label, plot=False):
    """Calculates the cartesian coordinates for rover location and the coordinates for satellite connection.
        ### Params:
            rover_true_anomaly_rad: True anomaly of rover
            label: Label of comm connection being plotted
            plot: If graph should be plotted, default = False
        ### Return:
            Satellite connection cartesian coordinates
    """

    rover_location_x_m = mars_radius_m * np.cos(rover_true_anomaly_rad)
    rover_location_y_m = mars_radius_m * np.sin(rover_true_anomaly_rad)

    comm_x_m, comm_y_m = comm_intersection_calculation(rover_location_x_m, rover_location_y_m, rover_true_anomaly_rad)

    if plot:
        plot_orbit_cartesian(rover_location_x_m, rover_location_y_m, comm_x_m, comm_y_m, label)

    return np.array([comm_x_m, comm_y_m])

def comm_intersection_calculation(rover_location_x_m, rover_location_y_m, rover_true_anomaly_rad):
    """Calculates the coordinates when the satellite enters and leaves the rover's field of view.
        ### Params:
            rover_location_x_m: Rover x coordinate
            rover_location_y_m: Rover y coordinate
            rover_true_anomaly_rad: True anomaly of rover
        ### Return:
            x coordinates for satellite connection
            y coordinates for satellite connection
    """

    comm_gradient = np.array((min_angle_rad + rover_true_anomaly_rad, -min_angle_rad + rover_true_anomaly_rad))

    a = 1 / semimajor_axis_m ** 2 + 1 / (semiminor_axis_m * np.tan(comm_gradient)) ** 2
    b = 2 * eccentricity / semimajor_axis_m - 2 * mars_radius_m * (np.cos(rover_true_anomaly_rad) + np.sin(rover_true_anomaly_rad) * np.tan(comm_gradient)) / (semiminor_axis_m ** 2 * np.tan(comm_gradient) ** 2)
    c = (mars_radius_m * (np.cos(rover_true_anomaly_rad) + np.sin(rover_true_anomaly_rad) * np.tan(comm_gradient))) ** 2 / (semiminor_axis_m ** 2 * np.tan(comm_gradient) ** 2) + eccentricity ** 2 - 1

    x = np.array(quadratic_formula(a, b, c)).flatten()
    y = elliptical_orbit(x)

    rover_vector_m = [rover_location_x_m, rover_location_y_m]
    rover_unit_vector = rover_vector_m / np.linalg.norm(rover_vector_m)
    rover_comm_vector_m = np.array([x - rover_location_x_m, y - rover_location_y_m])
    rover_comm_unit_vector = rover_comm_vector_m / np.linalg.norm(rover_comm_vector_m, axis=0)
    angle = np.round(np.arccos(np.dot(rover_unit_vector, rover_comm_unit_vector)) * 180 / const.pi, 0)
    angle_index = np.argwhere((angle != 80) & (angle != 100)).flatten()

    if len(angle_index) != 0:
        y[angle_index] = -y[angle_index]
        rover_comm_vector_m = np.array([x - rover_location_x_m, y - rover_location_y_m])
        rover_comm_unit_vector = rover_comm_vector_m / np.linalg.norm(rover_comm_vector_m, axis=0)
        angle = np.round(np.arccos(np.dot(rover_unit_vector, rover_comm_unit_vector)) * 180 / const.pi, 0)

    comm_index = np.argwhere(angle != 80)
    x = np.delete(x, comm_index)
    y = np.delete(y, comm_index)

    if y[0] == y[1]:
        y[1] = -y[1]

    return x, y

def plot_orbit_cartesian(rover_location_x, rover_location_y, comm_x_m, comm_y_m, label):
    """Plots a graph to show the different comm connections and satellite orbit
        ### Params:
            rover_location_x: Rover x coordinate location
            rover_location_y: Rover y coordinate location
            comm_x_m: Comm connection x coordinate
            comm_y_m: Comm connection y coordinate
            label: Label of comm connection being plotted
    """
    fig = plt.figure('Orbit Cartesian', figsize=(16, 9), constrained_layout=True)

    if fig.number != 0:
        fig.number = 0
        ax = plt.gca()
        plt.get_current_fig_manager().canvas.manager.window.wm_geometry("+%d+%d" % (0, 0))

        orbit = matplotlib.patches.Ellipse([-semimajor_axis_m * eccentricity, 0], width=(2 * semimajor_axis_m), height=(2 * semiminor_axis_m), fc='None', edgecolor='b')
        mars = plt.Circle((0, 0), mars_radius_m, color='orangered')

        ax.add_patch(orbit)
        ax.add_patch(mars)
        ax.set_aspect('equal')
        ax.set_xlabel('x distance from center of Mars / m', fontsize=12)
        ax.set_ylabel('y distance from center of Mars / m', fontsize=12)

    plt.plot([comm_x_m[0], rover_location_x, comm_x_m[1]], [comm_y_m[0], rover_location_y, comm_y_m[1]], zorder=1, label=label)
    plt.scatter(rover_location_x, rover_location_y, color='b', zorder=2)
    plt.scatter(comm_x_m, comm_y_m, color='b', zorder=2)
    plt.legend()

#Initialise Constants
period_s = 4.5 * 3600
mars_mass_kg = 6.39e23
mars_radius_m = 3389.5e3
min_angle_rad = 10 * const.pi / 180
frequency_hz = 442e6

rover_angle_rad = const.pi / 2 + min_angle_rad
radius_periareion_m = 150e3 + mars_radius_m
radius_apoareion_m = 6200e3 + mars_radius_m
wavelength_m = const.c / frequency_hz
mean_motion_rads = 2 * const.pi / period_s

semimajor_axis_m = (const.G * mars_mass_kg * period_s ** 2 / (4 * const.pi ** 2)) ** (1/3)
velocity_periareion_ms = (2 * const.G * mars_mass_kg / radius_periareion_m - const.G * mars_mass_kg / semimajor_axis_m) ** (1/2)
velocity_apoareion_ms = radius_periareion_m * velocity_periareion_ms / radius_apoareion_m
semi_latus_rectum_m = (radius_periareion_m * velocity_periareion_ms) ** 2 / (const.G * mars_mass_kg)
eccentricity = 1 - radius_periareion_m / semimajor_axis_m
semiminor_axis_m = (semimajor_axis_m * semi_latus_rectum_m) ** (1/2)

main()

plt.savefig('Rover comm elliptical orbit')
plt.show()
