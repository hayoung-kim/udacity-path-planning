#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

#include "spline.h"
#include "polysolver.h"
#include "coordinate_transforms.h"
#include "checkcollision.h"

typedef struct Vehicle {
  int id;
  double x;
  double y;
  double vx;
  double vy;
  double s;
  double d;
  double speed;
  // vector<double> future_s;
} Vehicle;

typedef struct Planner {
  double target_d;
  vector<Vehicle> obstacles;
  Vehicle target_to_follow;
  int following_target_id;
  double dist_to_target;
  MatrixXd s_trajectories;
  VectorXd s_costs;
  MatrixXd d_trajectories;
  VectorXd d_costs;
  bool obstacle_following; // if false: velocity keeping
  bool feasible_traj_exist;
  int optimal_s_id;
  int optimal_d_id;
  double minimal_cost;
  int iters;
} Planner;

int getMyLane(double d0) {
  int mylane = 1;
  if (d0 > 0 && d0 <= 4) {mylane = 0;}
  else if (d0 >4 && d0 <= 8) {mylane = 1;}
  else {mylane = 2;}
  return mylane;
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

int main() {

  uWS::Hub h;

  // debuging
  // freopen( "output.txt", "w", stdout );

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  // double max_s = 6914.15;


  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // double max_s = map_waypoints_s[map_waypoints_s.size()-1];

  VectorXd optimal_s_coeff(6);
  VectorXd optimal_d_coeff(6);
  double s_cost = 999;
  double d_cost = 999;
  int step = 0;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx\
    ,&map_waypoints_dy, &optimal_s_coeff, &optimal_d_coeff, &step, &max_s](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            double pos_x;
            double pos_y;
            double pos_yaw;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            double MPH2mps = 1.0/2.23694;
            double max_speed = 44*MPH2mps;

            int prev_path_size = previous_path_x.size();
            int max_s_waypoint = map_waypoints_s[map_waypoints_s.size()-1];

            step += 1;
            cout << "------------------------------------------------" << endl;
            cout << " [-] step: " << step << endl;
            cout << " [-] max_s: " << max_s << endl;
            // // cout << "optimal_s_coeff: \n" << optimal_s_coeff << endl;
            // cout << " [-] prev_path_size: " << prev_path_size << endl;
            //
            // cout << " [-] car_s: " << car_s << endl;
            // cout << " [-] car_d: " << car_d << endl;
            // cout << " [-] speed: " << car_speed << " (MPH)" << endl;
            cout << " [-] speed: " << car_speed*MPH2mps << " (m/s)" << endl;

            // INITIALIZE PLANNER
            vector<Planner> planners;
            for (int i=0; i<3; i++) {
              double _target_d = 2.0 + 4* i - 0.02;
              Planner planner;
              MatrixXd s_trajectories(6, 0);
              VectorXd s_costs(0);
              MatrixXd d_trajectories(6, 0);
              VectorXd d_costs(0);

              planner.s_trajectories = s_trajectories;
              planner.s_costs = s_costs;
              planner.d_trajectories = d_trajectories;
              planner.d_costs = d_costs;
              planner.target_d = _target_d;
              planner.dist_to_target = 999.9;
              planner.obstacle_following = false;
              planner.feasible_traj_exist = true;
              planner.minimal_cost = 9999999.9;
              planner.optimal_s_id = 0;
              planner.optimal_d_id = 0;
              planner.iters = -1;
              planners.push_back(planner);
            }

            // --------------------------------------------------
            // GET INFORMATION OF NEARBY VEHICLES
            // --------------------------------------------------
            vector<Vehicle> NearbyVehicles;
            for (int i=0; i<sensor_fusion.size(); i++) {
              // parsing the sensor fusion data
              // id, x, y, vx, vy, s, d
              double s_other = sensor_fusion[i][5];
              double s_dist = s_other - car_s;
              double d_other = sensor_fusion[i][6];
              // NEARBY VEHICLES
              double detect_range_front = 70.0;
              double detect_range_backward = 20.0;

              if ((s_dist < detect_range_front) && (s_dist >= - detect_range_backward) && (d_other > 0)) {
                Vehicle vehicle;
                vehicle.id = sensor_fusion[i][0];
                vehicle.x  = sensor_fusion[i][1];
                vehicle.y  = sensor_fusion[i][2];
                vehicle.vx = sensor_fusion[i][3];
                vehicle.vy = sensor_fusion[i][4];
                vehicle.s  = sensor_fusion[i][5];
                vehicle.d  = sensor_fusion[i][6];
                vehicle.speed = sqrt(vehicle.vx * vehicle.vx + vehicle.vy * vehicle.vy);

                NearbyVehicles.push_back(vehicle);
              }
            }

            int n_planning_horizon = 150;
            int n_pass = n_planning_horizon - prev_path_size;
            int start_index = n_pass - 1;
            double start_time = start_index * 0.02;

            double s0 = getPosition(optimal_s_coeff, start_time);
            if (s0 > max_s) {s0 = s0 - max_s; cout << "s0 = " << s0 << endl;}
            double s0dot = getVelocity(optimal_s_coeff, start_time);
            double s0ddot = getAcceleration(optimal_s_coeff, start_time);
            double d0 = getPosition(optimal_d_coeff, start_time);
            double d0dot = getVelocity(optimal_d_coeff, start_time);
            double d0ddot = getAcceleration(optimal_d_coeff, start_time);

            cout << " [!!] s0, d0 = " << s0 << ", " << d0 << endl;
            int mylane;
            // if (prev_path_size == 0) {mylane = getMyLane(car_d);}
            // else {mylane = getMyLane(d0);}
            mylane = getMyLane(car_d);



            // ---------------------------------------
            // EXTRACT NEAREST VEHICLE FOR EACH LANE
            // ---------------------------------------
            cout << "extract" << endl;
            for (int i=0; i<NearbyVehicles.size(); i++) {
              Vehicle _vehicle = NearbyVehicles[i];
              for (int j=0; j<planners.size(); j++) {
                if ((_vehicle.d > planners[j].target_d - 2) && (_vehicle.d <= planners[j].target_d + 2)) {

                  // frontmost obstacle
                  double from_ego_to_other = _vehicle.s - car_s;
                  if (j == mylane) {
                    if (from_ego_to_other >= 3) {
                      planners[j].obstacles.push_back(_vehicle);
                    }
                  }
                  else {planners[j].obstacles.push_back(_vehicle);}

                  if (from_ego_to_other >= -2.0){
                    if (from_ego_to_other < planners[j].dist_to_target) {
                      planners[j].dist_to_target = from_ego_to_other;
                      planners[j].target_to_follow = _vehicle;
                    }
                    if (from_ego_to_other <= 65) {
                      planners[j].obstacle_following = true;
                    }
                  }
                }
              }
            }



            int _;

            // WAYPOINTS SMOOTHING
            cout << "smoothing" << endl;
            int id_map_last = map_waypoints_x.size() - 1;
            // !- test loop
            int _close_way_point_id;
            double start_s;
            if (step < 20){
              start_s = map_waypoints_s[id_map_last-5];
              double start_x = map_waypoints_x[id_map_last-5];
              double start_y = map_waypoints_y[id_map_last-5];

              _close_way_point_id = ClosestWaypoint(start_x, start_y, map_waypoints_x, map_waypoints_y) + 1;
              if (_close_way_point_id == id_map_last) {
                for (int jj=0; jj<10; jj++) {cout << " [!!!!] CLOSE WAY POINT ID = MAP'S LAST WAYPOINT ID" << endl;}
              }
            }
            else {
              _close_way_point_id = ClosestWaypoint(car_x, car_y, map_waypoints_x, map_waypoints_y);
            }

            // -! test loop

            // int _close_way_point_id = ClosestWaypoint(car_x, car_y, map_waypoints_x, map_waypoints_y);

            int id_interp_start = _close_way_point_id - 4;
            int id_interp_end   = _close_way_point_id + 9;

            // cout << "setting a range for interpolate ... " << endl;
            if (id_interp_start < 0) {id_interp_start = 0;}
            // if (id_interp_end > id_map_last) {id_interp_end = id_map_last;}
            cout << "smoothing 2" << endl;
            vector<double> map_x_to_interp, map_y_to_interp, map_s_to_interp;
            for (int map_id=id_interp_start; map_id < id_interp_end; map_id ++) {

              if (map_id == _close_way_point_id) {
                cout << " nearest way point ! " << endl;
              }

              if (map_id == id_map_last) {
                cout << " last point ! " << endl;
              }

              if (map_id > id_map_last) {
                int _map_id = map_id - id_map_last;
                cout << map_id << endl;
                map_s_to_interp.push_back(map_waypoints_s[_map_id] + max_s);
                map_x_to_interp.push_back(map_waypoints_x[_map_id]);
                map_y_to_interp.push_back(map_waypoints_y[_map_id]);

                cout << "(s,x,y) = " << map_waypoints_s[_map_id] + max_s << \
                        ", " << map_waypoints_x[_map_id] << ", " << map_waypoints_y[_map_id] << endl;
              }
              else {
                map_s_to_interp.push_back(map_waypoints_s[map_id]);
                map_x_to_interp.push_back(map_waypoints_x[map_id]);
                map_y_to_interp.push_back(map_waypoints_y[map_id]);
                cout << "(s,x,y) = " << map_waypoints_s[map_id]<< \
                        ", " << map_waypoints_x[map_id] << ", " << map_waypoints_y[map_id] << endl;
              }


            }
            cout << "set spline" << endl;
            tk::spline x_given_s;
            tk::spline y_given_s;
            x_given_s.set_points(map_s_to_interp, map_x_to_interp);
            y_given_s.set_points(map_s_to_interp, map_y_to_interp);

            // cout << "interpolating starts...." << endl;

            vector<double> map_ss, map_xs, map_ys;
            double _s = map_s_to_interp[0];

            cout << "interp" << endl;

            while (_s < map_s_to_interp[map_s_to_interp.size()-1]) {
              double _x = x_given_s(_s);
              double _y = y_given_s(_s);
              map_ss.push_back(_s);
              map_xs.push_back(_x);
              map_ys.push_back(_y);
              _s += 0.05;
            }

            // TRAJECTORY PLANNING

            // SET INITIAL s0, d0 and their derivatives
            // prev_path_size : number of left over of previous planned trajectory
            // n_planning_horizon : n_use_previous_path + n_newly_planned_path

            // int n_pass = n_planning_horizon - prev_path_size;
            // int start_index = n_pass - 1;
            // double start_time = start_index * 0.02;
            //
            // double s0 = getPosition(optimal_s_coeff, start_time);
            // if (s0 > max_s) {s0 = s0 - max_s; cout << "s0 = " << s0 << endl;}
            // double s0dot = getVelocity(optimal_s_coeff, start_time);
            // double s0ddot = getAcceleration(optimal_s_coeff, start_time);
            // double d0 = getPosition(optimal_d_coeff, start_time);
            // double d0dot = getVelocity(optimal_d_coeff, start_time);
            // double d0ddot = getAcceleration(optimal_d_coeff, start_time);
            //
            // cout << " [!!] s0, d0 = " << s0 << ", " << d0 << endl;

            cout << "planning" << endl;

            if (prev_path_size == 0) {
              double target_s1dot = 20 * MPH2mps;
              bool in_mylane = true;
              // _ = VelocityKeepingTrajectories(car_s, car_speed, 0, target_s1dot, max_speed, \
              //                                 planners[1].s_trajectories, planners[1].s_costs);
              // for test loop
              _ = VelocityKeepingTrajectories(start_s, car_speed, 0, target_s1dot, max_speed, \
                                              planners[1].s_trajectories, planners[1].s_costs);
              _ = lateralTrajectories(car_d, 0, 0, 6.0, in_mylane, planners[1].d_trajectories, planners[1].d_costs);
              vector<int> opt_idx = optimalCombination(planners[1].s_costs, planners[1].d_costs);
              optimal_s_coeff = planners[1].s_trajectories.col(opt_idx[0]);
              optimal_d_coeff = planners[1].d_trajectories.col(opt_idx[1]);

              cout << "pushing current position ... " << endl;
              next_x_vals.push_back(car_x);
              next_y_vals.push_back(car_y);

              for (int hrz=0; hrz<n_planning_horizon; hrz++) {
                double s = getPosition(optimal_s_coeff, hrz*0.02);
                double d = getPosition(optimal_d_coeff, hrz*0.02);
                // vector<double> xy = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                vector<double> xy = getXY(s, d, map_ss, map_xs, map_ys);
                next_x_vals.push_back(xy[0]);
                next_y_vals.push_back(xy[1]);

                cout << "s=" << s << endl;
                cout << "x=" << xy[0] << endl;
              }
            }
            else {
              // GENERATE TRAJECTORY CANDIDATES AND COSTS
              double target_s1dot = (car_speed) * MPH2mps;
              if (target_s1dot > max_speed) {target_s1dot = max_speed;}

              for (int i=0; i<3; i++) {
                // LONGITUDINAL TRAJECTORY GENERATION
                if (i == mylane) {
                  bool in_mylane = true;
                  // cout << " [*] planning trajecotries to my lane ..." << endl;
                  if (planners[i].obstacle_following){
                    Vehicle target_obstacle = planners[i].target_to_follow;
                    _ = FollowingTrajectories(s0, s0dot, s0ddot, \
                          target_obstacle.s, target_obstacle.speed - 0.2, max_speed, \
                          planners[i].s_trajectories, planners[i].s_costs);

                    _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, \
                                                    target_obstacle.speed, max_speed, \
                                                    planners[i].s_trajectories, planners[i].s_costs);
                  }
                  else {
                    _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, \
                                                    target_s1dot, max_speed, \
                                                    planners[i].s_trajectories, planners[i].s_costs);
                  }

                  // stop
                  _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, \
                                                  0.0, max_speed, \
                                                  planners[i].s_trajectories, planners[i].s_costs);
                  // LATERAL TRAJECTORY GENERATION
                  _ = lateralTrajectories(d0, d0dot, d0ddot, \
                                          planners[i].target_d, in_mylane, \
                                          planners[i].d_trajectories, planners[i].d_costs);

                }
                else if ((abs(i - mylane) <= 1) && (car_speed*MPH2mps > 9.0)) {
                  // cout << " [*] planning trajectories to nearby lane ... " << endl;
                  bool in_mylane = false;
                  if (planners[i].obstacle_following){
                    Vehicle target_obstacle = planners[i].target_to_follow;
                    _ = FollowingTrajectories(s0, s0dot, s0ddot, \
                          target_obstacle.s, target_obstacle.speed - 0.2, max_speed, \
                          planners[i].s_trajectories, planners[i].s_costs);
                    _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, \
                                                    target_s1dot - 10.0, max_speed, \
                                                    planners[i].s_trajectories, planners[i].s_costs);
                  }
                  else {
                    // keeping
                    _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, \
                                                    target_s1dot, max_speed, \
                                                    planners[i].s_trajectories, planners[i].s_costs);
                  }

                  // stop
                  _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, \
                                                  0.0, max_speed, \
                                                  planners[i].s_trajectories, planners[i].s_costs);


                  // LATERAL TRAJECTORY GENERATION
                  _ = lateralTrajectories(d0, d0dot, d0ddot, \
                                          planners[i].target_d, in_mylane, \
                                          planners[i].d_trajectories, planners[i].d_costs);
                }
                else {
                  planners[i].feasible_traj_exist = false;
                  cout << " [*] infeasible because this lane is not adjacent to my lane (" << i << ")" << endl;
                }

                // OPTIMAL TRAJECTORY SELECTION
                double klon = 1.0;
                double klat = 0.5;
                int ns = planners[i].s_costs.size();
                int nd = planners[i].d_costs.size();
                int ntraj = ns * nd;

                // build cost matrix
                if (ntraj == 0) {
                  planners[i].feasible_traj_exist = false;
                  cout << " [*] infeasible due to ntraj = 0 (" << i << ")" << endl;
                }
                else {
                  MatrixXd sd_costs(ns, nd);
                  cout << " [*] build cost matrix " << i << " ..." << endl;
                  for (int ss=0; ss<ns; ss++){
                    for (int dd=0; dd<nd; dd++){
                      sd_costs(ss,dd) = klon * planners[i].s_costs[ss] \
                                      + klat * planners[i].d_costs[dd];
                    }
                  }

                  cout << " [-] sd_costs size: (" << ns << ", " << nd << ")" << endl;

                  // collision check
                  // cout << " [*] checking collision " << i+1 << " ..." << endl;
                  int max_iter = 100;
                  int iters = -1;
                  if (max_iter >= ntraj) {max_iter = ntraj;}
                  for (int k=0; k<max_iter; k++) {
                    // cout << " [*] iters " << k << endl;
                    bool crash_predicted = false;
                    int min_s_idx, min_d_idx;
                    double minCost = sd_costs.minCoeff(&min_s_idx, &min_d_idx);
                    optimal_s_coeff = planners[i].s_trajectories.col(min_s_idx);
                    optimal_d_coeff = planners[i].d_trajectories.col(min_d_idx);

                    for (int t=0; t<n_planning_horizon; t++) {
                      // my position
                      double _s = getPosition(optimal_s_coeff, t*0.02);
                      double _d = getPosition(optimal_d_coeff, t*0.02);
                      // double _sdot = getVelocity(optimal_s_coeff, t*0.02);
                      // double _ddot = getVelocity(optimal_d_coeff, t*0.02);
                      // double _heading = atan2( _ddot, _sdot );

                      // other car's position
                      for (int n=0; n<planners[i].obstacles.size(); n++) {
                        Vehicle other_car = planners[i].obstacles[n];
                        double _s_other = other_car.s + t * 0.02 * (other_car.speed - 0.2);
                        double _d_other = other_car.d;
                        int crash = checkCollision(_s, _d, 0, _s_other, _d_other, 0.0);
                        if (crash == 1) {crash_predicted = true; break;}
                      }
                      if (crash_predicted) {sd_costs(min_s_idx, min_d_idx) = 9999999.9; break;}
                    }
                    iters = k;
                    if (!crash_predicted) {
                      // cout << " [*] LANE " << i+1 << " is FREE!" << endl;
                      planners[i].optimal_s_id = min_s_idx;
                      planners[i].optimal_d_id = min_d_idx;
                      planners[i].minimal_cost = sd_costs(min_s_idx, min_d_idx);
                      break;
                    }

                  }
                  planners[i].iters = iters;
                  if (iters == max_iter-1) {
                    planners[i].feasible_traj_exist = false;
                    cout << " [*] infeasible due to iters == max_iter-1 (" << i << ")" << endl;
                  }
                }
              }

              // FIND OPTIMAL ONE
              // cout << " [*] finding optimal lane ..." << endl;
              double minimal_cost = 9999999.9;
              int opt = 1;
              for (int i=0; i<3; i++) {

                if (planners[i].feasible_traj_exist) {
                  if (planners[i].minimal_cost < minimal_cost) {
                    opt = i;
                    minimal_cost = planners[i].minimal_cost;
                  }
                }

              }

              cout << " [-] NearbyVehicles = " << NearbyVehicles.size() << endl;
              cout << " [-] mylane = " << mylane << endl;
              cout << " [-] obstacles : " << planners[0].obstacles.size() << ", " << \
                      planners[1].obstacles.size() << ", " << planners[2].obstacles.size() << endl;
              cout << " [-] feasible : " << planners[0].feasible_traj_exist << ", " << \
                      planners[1].feasible_traj_exist << ", " << planners[2].feasible_traj_exist << endl;
              cout << " [-] following : " << planners[0].obstacle_following << ", " << \
                      planners[1].obstacle_following << ", " << planners[2].obstacle_following << endl;
              cout << " [-] cost : " << planners[0].minimal_cost << ", " << \
                      planners[1].minimal_cost << ", " << planners[2].minimal_cost << endl;
              cout << " [-] iters : " << planners[0].iters << ", " << \
                      planners[1].iters << ", " << planners[2].iters << endl;
              cout << " [-] optimal lane: " << opt << endl;

              int opt_s = planners[opt].optimal_s_id;
              int opt_d = planners[opt].optimal_d_id;

              // cout << " [-] opt_s = " << opt_s << ", opt_d = " << opt_d << endl;
              optimal_s_coeff = planners[opt].s_trajectories.col(opt_s);
              optimal_d_coeff = planners[opt].d_trajectories.col(opt_d);

              // RUN!!
              for (int hrz=0; hrz<n_planning_horizon + 1; hrz++) {
                double s = getPosition(optimal_s_coeff, hrz*0.02);
                double d = getPosition(optimal_d_coeff, hrz*0.02);


                // loop solution
                // if ((s < map_s_to_interp[0]) && (map_s_to_interp[map_s_to_interp.size()-1] >= max_s)) {
                //   s = s + max_s;
                // }
                if (s > map_s_to_interp[map_s_to_interp.size()-1]) {
                  s = s - max_s;
                }
                else if ( s < map_s_to_interp[0]) {
                  s = s + max_s;
                }




                vector<double> xy = getXY(s, d, map_ss, map_xs, map_ys);
                next_x_vals.push_back(xy[0]);
                next_y_vals.push_back(xy[1]);
                if (hrz % 10 == 0){cout << " s to (x,y) = " << s << ", " <<  xy[0] << ", " << xy[1]<< endl;}
              }
            }

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;
            // msgJson["next_x"] = interpolated_x_vals;
            // msgJson["next_y"] = interpolated_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
