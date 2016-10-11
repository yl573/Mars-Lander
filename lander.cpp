// Mars lander simulator
// Version 1.8
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, October 2014

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"

vector3d calc_drag(vector3d v = velocity);
vector3d get_wind(bool stable = true);



//function for calculating the desired velocity to enter a Hohmann transfer orbit
double desired_speed(double current_r, double desired_r)
{
	double vr = velocity * position.norm();
	double temp = (1 / desired_r - 1 / current_r) * 2 * GRAVITY * MARS_MASS + vr * vr;
	return sqrt(temp / (pow(current_r / desired_r, 2) - 1));
}

//function for calculating the desired velocity to enter a Hohmann transfer orbit
//neglects radial velocity
double ini_desired_speed(double current_r, double desired_r)
{
	double temp = (1 / desired_r - 1 / current_r) * 2 * GRAVITY * MARS_MASS;
	return sqrt(temp / (pow(current_r / desired_r, 2) - 1));
}

inline double clamp(double x)
{
	return (x > 1) ? 1 : ((x < 0) ? 0 : x);
}

inline int sign(double x)
{
	return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
	stabilized_attitude = true;
	double current_mass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY;
	double altitude = position.abs() - MARS_RADIUS;
	double vr = velocity * position.norm();
	vector3d vtheta = vector3d(velocity - position.norm()*(velocity*position.norm()));
	double hor_demand = 0;//+ve in direction of travel
	double ver_demand = 0;//+ve upwards
	double max_mass = UNLOADED_LANDER_MASS + FUEL_CAPACITY * FUEL_DENSITY;
	double delta = GRAVITY * MARS_MASS * current_mass/ (position.abs2() * MAX_THRUST);

	switch (autopilot_status) {

		case LAND:
		{
			/*
			De-orbit:
			The autopilot first starts a burn that sends the lander on an elliptic orbit with perigee at half
			the EXOSPHERE/2.This causes the orbit to decay over time, saving the fuel needed to reduce
			altitude. Eventually the lander completely enters the atmosphere and an ordinary landing starts.
			*/
			static bool danger_zone;
			if (altitude > EXOSPHERE / 2)//Orbit not clipping dense atmosphere
			{
				//calculate the speed needed for the apogee to clip the atmosphere
				if (vtheta.abs() > desired_speed(position.abs(), MARS_RADIUS + EXOSPHERE / 2)) 
				{
					hor_demand = 1.0; //reduce speed
				}
			}
			else //Orbit clipping dense atmosphere, initiate landing
			{
				if (new_autopilot)
				{
					danger_zone = false;
					new_autopilot = false;
				}
				if (safe_to_deploy_parachute() && (position.abs() - MARS_RADIUS) < EXOSPHERE / 4)
					parachute_status = DEPLOYED;
				double kh = 0.12 / ((ENGINE_DELAY + ENGINE_LAG)*0.01 + WIND_SPEED + 1);
				double v_desired;
				double kp;

				/*Deal with wind:
				Wind too powerful will require too much thrust in the horizontal 
				direction to balance and there will not be enough thrust for vertical descent.
				*/
				if (parachute_status == DEPLOYED) //Make sure the speed is actually caused by wind
					hor_demand = vtheta * closeup_coords_right; 
				
				/*
				Delay and lag:
				The altitude used by the autopilot is the altitude ENGINE_DELAY + ENGINE_LAG 
				seconds into the future. This way when the engine demand is executed, the actual 
				velocity will be close to what is used in the calculation. Also, the transition 
				between greedy and mild is smoothed out when there is a delay to prevent the feedback
				loop from becoming unstable. However, landing with high delay and lag results in a lot of 
				extra fuel consumption and the 
				*/
				//altitude when the current command takes action (after the delay and lag)
				double future_altitude = altitude + velocity * position.norm() * (ENGINE_DELAY + ENGINE_LAG);
				future_altitude = future_altitude > 0 ? future_altitude : 0;

				/*
				Minimizing Fuel:
				At high altitude a “greedy” controller is used.It calculates the velocity the lander
				needs in order for safe landing to be possible from the maximum thrust.When the altitude
				is lower, a “mild” controller similar to the one on the handout is used for the final touch down.
				The change of controllers occour when the lander is 8 seconds from hitting the ground
				*/
				if (future_altitude > velocity.abs() * (8 + (ENGINE_DELAY + ENGINE_LAG) * 10) && !danger_zone) //above danger zone	
				{	
					//greedy controller
					kp = 0.5;
					double max_deceleration = (MAX_THRUST / UNLOADED_LANDER_MASS) - (GRAVITY*MARS_MASS / (MARS_RADIUS*MARS_RADIUS));
					v_desired = -sqrt(2 * max_deceleration * future_altitude);
				}
				else
				{	
					//mild controller		
					danger_zone = true;
					kp = 0.2 / ((ENGINE_DELAY + ENGINE_LAG) + 1);
					v_desired = -(0.5 + future_altitude * kh);
				}
				double ep = v_desired - vr;
				double p_out = kp * ep;
				
				ver_demand = clamp(delta + p_out);
			}
			break;
		}

		case ORBITAL_TRANSFER:
		{
			double kp = 1;
			double vtheta_desired, vr_desired;
			static int transfer_status;
			// 0 = entering Hohmann transfer orbit
			// 1 = entering desired circular orbit
			// 2 = stabilising circular orbit
			static double ini_orb_rad, ini_orb_vel;

			//If start of orbital change, set orbit
			if (new_autopilot) 
			{
				transfer_status = 0;
				ini_orb_rad = position.abs();
				ini_orb_vel = ini_desired_speed(position.abs(), desired_altitude + MARS_RADIUS);
				new_autopilot = false;
			}
			//If orbit apagee/perigee has been reached, enter circular orbit
			if (abs(altitude - desired_altitude)/desired_altitude < 0.01 && transfer_status == 0)
			{
				transfer_status = 1;
			}
			double v_orbit = sqrt(GRAVITY * MARS_MASS / (desired_altitude + MARS_RADIUS));

			//If velocity matches that for the circular orbit, stabilize the orbit
			if ((velocity.abs() - v_orbit) / v_orbit < 0.05 && transfer_status == 1)
			{
				transfer_status = 2;
			}

			switch (transfer_status)
			{
			case 0:
				{
					//Conservation of M.O.M
					vtheta_desired = ini_orb_rad * ini_orb_vel / position.abs();
					hor_demand = -(vtheta_desired - vtheta.abs()) * kp;
					ver_demand = 0;
					break;
				}
			case 1:
				{
					//Conservation of M.O.M
					ini_orb_vel = sqrt(GRAVITY * MARS_MASS / (MARS_RADIUS + desired_altitude));
					ini_orb_rad = MARS_RADIUS + desired_altitude;
					vtheta_desired = ini_orb_rad * ini_orb_vel / position.abs();
					hor_demand = -(vtheta_desired - vtheta.abs()) * kp;
					ver_demand = 0;
					break;
				}
			case 2:
				{
					double kpver = 0.1, kphor = 0.01;
					//double vr = velocity * position / position.abs();
					double PE_required = GRAVITY * MARS_MASS * (1 / position.abs() - 1 / (MARS_RADIUS + desired_altitude));
					double v_orbit = sqrt(GRAVITY * MARS_MASS / (desired_altitude + MARS_RADIUS));

					//Target circumferential velocity approaches orbital velocity as altitude approaches target altitude
					//Also slows down if altitude is too high
					double vtheta_required = v_orbit - abs(desired_altitude - altitude)*0.05;
					double ep = abs(vtheta_required) - abs(vtheta.abs()); //Error to be corrected
					hor_demand = -kphor * ep;

					//Required radial velocity calculated by the kinetic energy needed to reach orbital altitude
					double vr_required = PE_required > 0 ? sqrt(2 * PE_required) : -sqrt(-2 * PE_required);
					double centripetal = vtheta.abs2()*current_mass / (position.abs()*MAX_THRUST);
					ver_demand = clamp((vr_required - vr) * kpver + delta - centripetal);
				}
			}
			break;
		}

		/*The autopilot controls the radial and circumferential components of thrust separately.
		At each frame, the kinetic energy required to lift the lander to orbital altitude is 
		calculated and the corresponding velocity is seen as the desired radial velocity. While the 
		circumferential velocity is set to approach the orbital velocity linearly as altitude 
		approaches the desired altitude.*/
		case INJECTION:
		{
			double kpver = 0.1, kphor = 0.01;
			double vr = velocity * position / position.abs();
			double PE_required = GRAVITY * MARS_MASS * (1 / position.abs() - 1 / (MARS_RADIUS + desired_altitude));
			hor_demand = 0;

			double v_orbit = sqrt(GRAVITY * MARS_MASS / (desired_altitude + MARS_RADIUS));
				
			double vtheta_required = v_orbit - abs(desired_altitude - altitude)*0.05;
			double ep = vtheta_required - vtheta.abs();
			hor_demand = clamp(kphor * ep);
			
			double vr_required = PE_required > 0 ? sqrt(2 * PE_required) : -sqrt(-2 * PE_required);	
			double centripetal = vtheta.abs2()*current_mass / (position.abs()*MAX_THRUST);
			ver_demand = clamp((vr_required - vr) * kpver + delta - centripetal);
			
			break;
		}
	}

	//hor and ver demand are converted into attitude and thrust
	throttle = clamp(sqrt(hor_demand * hor_demand + ver_demand * ver_demand));
	if (hor_demand == 0 || throttle < 0.0001)
		stabilized_attitude_angle = 0.0;
	else if (ver_demand == 0)
		stabilized_attitude_angle = -90 * sign(hor_demand);
	else
		stabilized_attitude_angle = atan(hor_demand / -ver_demand) * 180 / M_PI;
}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
	vector3d force, pos_temp, thrust, drag, weight, wind;
	double current_mass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY;
	static vector3d previous_position;
	if (simulation_time == 0.0)
		previous_position = position - velocity * delta_t;

	thrust = thrust_wrt_world();
	weight = -(GRAVITY * MARS_MASS * current_mass * position.norm() / position.abs2());
	drag = calc_drag();
	wind = get_wind(false);
	//cout << wind << endl;
	force = thrust + weight + drag + wind;
	pos_temp = position;
	position = 2 * position - previous_position + delta_t * delta_t * force / current_mass;
	velocity = (1 / (2 * delta_t)) * (position - previous_position);
	previous_position = pos_temp;

  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled) autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();
}

vector3d calc_drag(vector3d v)
{
	double projected_area;
	if (parachute_status != DEPLOYED)
	{
		projected_area = M_PI * LANDER_SIZE * LANDER_SIZE;
		return -0.5 * atmospheric_density(position) * DRAG_COEF_LANDER * projected_area * v.abs() * v;
	}
	else
	{
		projected_area = 20 * LANDER_SIZE * LANDER_SIZE;
		return -0.5 * atmospheric_density(position) * DRAG_COEF_CHUTE * projected_area * v.abs() * v;
	}
}

vector3d get_wind(bool stable)
{
	if (WIND_SPEED == 0)
		return vector3d(0, 0, 0);
	if (position.abs() > MARS_RADIUS + EXOSPHERE)
		return vector3d();
	vector3d wind_direction;

	wind_direction = vector3d(position ^ vector3d(1,0,0)).norm(); //wind perpendicular to position
	//cout << wind_direction*position<< endl;
	if (wind_direction.abs() == 0) //special case
		wind_direction = vector3d(0, 0, 1);

	if (!RANDOM_GUSTS)//stable wind
		return calc_drag(WIND_SPEED * wind_direction);
	else //random gusts
		return calc_drag((double)(rand()%WIND_SPEED + WIND_SPEED) * wind_direction);
}

void initialize_simulation(void)
// Lander pose initialization - selects one of 10 possible scenarios
{
	new_autopilot = true;
	autopilot_status = LAND;
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen
  // autopilot_status - the operation to be performed by the autopilot, default is LAND
  // desired_altitude - the desired altitude for orbital injection and orbital transfer 

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "aerostationary orbit";
  scenario_description[7] = "surface launch to 1000km";
  scenario_description[8] = "orbital transfer up from 1000km to 2000km";
  scenario_description[9] = "orbital transfer down from 4000km to 1000km";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(-(MARS_RADIUS + EXOSPHERE), 0.0, 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
	//aerostationary orbit
	position = vector3d(cbrt(GRAVITY*MARS_MASS*MARS_DAY*MARS_DAY/(4*M_PI*M_PI)), 0.0, 0.0);
	velocity = vector3d(0.0, 2*M_PI*cbrt(GRAVITY*MARS_MASS*MARS_DAY*MARS_DAY / (4 * M_PI*M_PI))/MARS_DAY, 0.0);
	orientation = vector3d(0.0, 90.0, 0.0);
	delta_t = 0.1;
	parachute_status = NOT_DEPLOYED;
	stabilized_attitude = false;
	autopilot_enabled = false;
    break;

  case 7:
	//surface launch to 5000km
	position = vector3d(1 + MARS_RADIUS + LANDER_SIZE / 2.0, 0.0, 0.0);
	velocity = vector3d(0.0, 0.0, 0.0);
	orientation = vector3d(0.0, 90.0, 0.0);
	delta_t = 0.1;
	parachute_status = NOT_DEPLOYED;
	stabilized_attitude = false;
	autopilot_enabled = true;
	autopilot_status = INJECTION;
	desired_altitude = 5000000;
    break;

  case 8:
	//orbital transfer from 1000km to 5000km
	position = vector3d(MARS_RADIUS + 1000000, 0.0, 0.0);
	velocity = vector3d(0.0,sqrt(GRAVITY*MARS_MASS/(MARS_RADIUS + 1000000)), 0.0);
	orientation = vector3d(0.0, 90.0, 0.0);
	delta_t = 0.1;
	parachute_status = NOT_DEPLOYED;
	stabilized_attitude = false;
	autopilot_enabled = false;
	desired_altitude = 5000000;
	autopilot_status = ORBITAL_TRANSFER;
    break;

  case 9:
	//orbital transfer from 7000km to 1000km
	position = vector3d(MARS_RADIUS + 7000000, 0.0, 0.0);
	velocity = vector3d(0.0, sqrt(GRAVITY*MARS_MASS / (MARS_RADIUS + 7000000)), 0.0);
	orientation = vector3d(0.0, 90.0, 0.0);
	delta_t = 0.1;
	parachute_status = NOT_DEPLOYED;
	stabilized_attitude = false;
	autopilot_enabled = false;
	desired_altitude = 1000000;
	autopilot_status = ORBITAL_TRANSFER;
    break;

  }
}
