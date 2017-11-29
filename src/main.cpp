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

typedef struct Vehicle {
  int id;
  double x;
  double y;
  double vx;
  double vy;
  double s;
  double d;
} Vehicle;

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

  VectorXd optimal_s_coeff(6);
  VectorXd optimal_d_coeff(6);
  double s_cost = 999;
  double d_cost = 999;
  int step = 0;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx\
    ,&map_waypoints_dy, &optimal_s_coeff, &optimal_d_coeff, &step](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

            int prev_path_size = previous_path_x.size();
            step += 1;
            cout << "------------------------------------------------" << endl;
            cout << " [-] step: " << step << endl;
            // cout << "optimal_s_coeff: \n" << optimal_s_coeff << endl;
            cout << " [-] prev_path_size: " << prev_path_size << endl;

            cout << " [-] car_s: " << car_s << endl;
            cout << " [-] car_d: " << car_d << endl;
            cout << " [-] speed: " << car_speed << endl;

            // --------------------------------------------------
            // GET INFORMATION OF NEARBY VEHICLES
            // -------------------------------------------------
            vector<Vehicle> NearbyVehicles;
            for (int i=0; i<sensor_fusion.size(); i++) {
              // parsing the sensor fusion data
              // id, x, y, vx, vy, s, d
              double s_other = sensor_fusion[i][5];
              double s_dist = s_other - car_s;
              double d_other = sensor_fusion[i][6];
              // NEARBY VEHICLES
              double detect_range_front = 80.0;
              double detect_range_backward = 0.0;

              if ((s_dist < detect_range_front) && (s_dist >= - detect_range_backward) && (d_other > 0)) {
                Vehicle vehicle;
                vehicle.id = sensor_fusion[i][0];
                vehicle.x  = sensor_fusion[i][1];
                vehicle.y  = sensor_fusion[i][2];
                vehicle.vx = sensor_fusion[i][3];
                vehicle.vy = sensor_fusion[i][4];
                vehicle.s  = sensor_fusion[i][5];
                vehicle.d  = sensor_fusion[i][6];

                NearbyVehicles.push_back(vehicle);
              }
            }


            cout << " [-] # NearbyVehicles = " << NearbyVehicles.size() << endl;




            // cout << " [-] Car speed= " << car_speed << endl;
            MatrixXd s_trajectories(6, 0);
            VectorXd s_costs(0);
            MatrixXd d_trajectories(6, 0);
            VectorXd d_costs(0);

            int n_planning_horizon = 100;
            int n_use_previous_path = 0;

            int _;

            // WAYPOINTS SMOOTHING
            int _close_way_point_id = ClosestWaypoint(car_x, car_y, map_waypoints_x, map_waypoints_y);
            int id_interp_start = _close_way_point_id - 5;
            int id_interp_end   = _close_way_point_id + 5;
            int id_map_last = map_waypoints_x.size();

            cout << "setting a range for interpolate ... " << endl;
            if (id_interp_start < 0) {id_interp_start = 0;}
            if (id_interp_end > id_map_last) {id_interp_end = id_map_last;}

            vector<double> map_x_to_interp, map_y_to_interp, map_s_to_interp;
            for (int map_id=id_interp_start; map_id < id_interp_end; map_id ++) {
              map_x_to_interp.push_back(map_waypoints_x[map_id]);
              map_y_to_interp.push_back(map_waypoints_y[map_id]);
              map_s_to_interp.push_back(map_waypoints_s[map_id]);
            }

            tk::spline x_given_s;
            tk::spline y_given_s;
            x_given_s.set_points(map_s_to_interp, map_x_to_interp);
            y_given_s.set_points(map_s_to_interp, map_y_to_interp);

            // cout << "interpolating starts...." << endl;

            vector<double> map_ss, map_xs, map_ys;
            double _s = map_s_to_interp[0];
            // cout << " [-] s_start = " << _s << endl;
            // cout << " [-] s_end = " << map_s_to_interp[map_s_to_interp.size()-1] << endl;

            while (_s < map_s_to_interp[map_s_to_interp.size()-1]) {
              double _x = x_given_s(_s);
              double _y = y_given_s(_s);
              map_ss.push_back(_s);
              map_xs.push_back(_x);
              map_ys.push_back(_y);
              _s += 0.2;
            }

            // TRAJECTORY PLANNING

            if (prev_path_size == 0) {
              double target_s1dot = 20 / 2.23694;
              _ = VelocityKeepingTrajectories(car_s, car_speed, 0, target_s1dot, s_trajectories, s_costs);
              _ = lateralTrajectories(car_d, 0, 0, 6.0, d_trajectories, d_costs);
              vector<int> opt_idx = optimalCombination(s_costs, d_costs);
              optimal_s_coeff = s_trajectories.col(opt_idx[0]);
              optimal_d_coeff = d_trajectories.col(opt_idx[1]);

              for (int hrz=0; hrz<n_planning_horizon; hrz++) {
                double s = getPosition(optimal_s_coeff, hrz*0.02);
                double d = getPosition(optimal_d_coeff, hrz*0.02);
                // vector<double> xy = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                vector<double> xy = getXY(s, d, map_ss, map_xs, map_ys);
                next_x_vals.push_back(xy[0]);
                next_y_vals.push_back(xy[1]);
              }
            }
            else {
              // SEND SUBSET OF PREVIOUS PATH
              int n_subset = n_use_previous_path;
              if (n_subset >= prev_path_size) {n_subset = prev_path_size;}

              cout << "n_subset: " << n_subset << endl;

              for (int h=0; h<n_subset; h++){
                next_x_vals.push_back(previous_path_x[h]);
                next_y_vals.push_back(previous_path_y[h]);
              }

              // SET INITIAL s0, d0 and their derivatives
              // prev_path_size : number of left over of previous planned trajectory
              // n_use_previous_path : how many use previous path
              // n_planning_horizon : n_use_previous_path + n_newly_planned_path

              int n_pass = n_planning_horizon - prev_path_size;
              int start_index = n_subset + n_pass - 1;
              double start_time = start_index * 0.02;

              double s0 = getPosition(optimal_s_coeff, start_time);
              double s0dot = getVelocity(optimal_s_coeff, start_time);
              double s0ddot = getAcceleration(optimal_s_coeff, start_time);
              double d0 = getPosition(optimal_d_coeff, start_time);
              double d0dot = getVelocity(optimal_d_coeff, start_time);
              double d0ddot = getAcceleration(optimal_d_coeff, start_time);

              // CALCULATE TRAJECTORY CANDIDATES AND COSTS OF THEOREM
              double target_s1dot = (car_speed + 20) / 2.23694;
              if (target_s1dot > 45 / 2.23694) {target_s1dot = 45 / 2.23694;}

              _ = VelocityKeepingTrajectories(s0, s0dot, s0ddot, target_s1dot, s_trajectories, s_costs);
              _ = lateralTrajectories(d0, d0dot, d0ddot, 2.0, d_trajectories, d_costs);
              // _ = lateralTrajectories(d0, d0dot, d0ddot, 6.0, d_trajectories, d_costs);
              // _ = lateralTrajectories(d0, d0dot, d0ddot, 10.0, d_trajectories, d_costs);

              vector<int> opt_idx = optimalCombination(s_costs, d_costs);
              optimal_s_coeff = s_trajectories.col(opt_idx[0]);
              optimal_d_coeff = d_trajectories.col(opt_idx[1]);

              for (int hrz=0; hrz<n_planning_horizon - n_subset + 1; hrz++) {
                double s = getPosition(optimal_s_coeff, hrz*0.02);
                double d = getPosition(optimal_d_coeff, hrz*0.02);
                vector<double> xy = getXY(s, d, map_ss, map_xs, map_ys);
                next_x_vals.push_back(xy[0]);
                next_y_vals.push_back(xy[1]);
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
