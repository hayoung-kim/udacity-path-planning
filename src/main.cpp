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

            // cout << " [-] Car speed= " << car_speed << endl;
            MatrixXd s_trajectories(6, 0);
            VectorXd s_costs(0);
            MatrixXd d_trajectories(6, 0);
            VectorXd d_costs(0);

            int n_planning_horizon = 300;
            int n_use_previous_path = 30;

            int _;


            if (prev_path_size == 0) {

              double target_s1dot = 20 / 2.23694;
              _ = VelocityKeepingTrajectories(car_s, car_speed, 0, target_s1dot, s_trajectories, s_costs);
              _ = lateralTrajectories(car_d, 0, 0, 6.0, d_trajectories, d_costs);
              vector<int> opt_idx = optimalCombination(s_costs, d_costs);
              optimal_s_coeff = s_trajectories.col(opt_idx[0]);
              optimal_d_coeff = d_trajectories.col(opt_idx[1]);

              // double vel = 20 / 2.23694;
              // vector<double> ssd = getFrenet(car_x, car_y, deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);
              for (int hrz=0; hrz<n_planning_horizon; hrz++) {
                double s = getPosition(optimal_s_coeff, hrz*0.02);
                double d = getPosition(optimal_d_coeff, hrz*0.02);
                // double s = ssd[0] + vel * 0.02 * hrz;
                // double d = 6.0;
                vector<double> xy = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                next_x_vals.push_back(xy[0]);
                next_y_vals.push_back(xy[1]);
              }
            }
            else {
              for (int hrz=0; hrz<prev_path_size; hrz++) {
                next_x_vals.push_back(previous_path_x[hrz]);
                next_y_vals.push_back(previous_path_y[hrz]);
              }
            }

            //
            // if (prev_path_size == 0) {
            //   pos_x = car_x;
            //   pos_y = car_y;
            //   pos_yaw = deg2rad(car_yaw);
            //
            //   double target_s1dot = 3.0;
            //   // cout << " [*] Planning starts !!!" << endl;
            //   // cout << " [*] Generating VelocityKeepingTrajectories ..." << endl;
            //   _ = VelocityKeepingTrajectories(car_s, 0, 0, target_s1dot, s_trajectories, s_costs);
            //   // cout << " [*] Generating lateralTrajectories ..." << endl;
            //   _ = lateralTrajectories(car_d, 0, 0, car_d, d_trajectories, d_costs);
            //   vector<int> opt_idx = optimalCombination(s_costs, d_costs);
            //   // cout << " [-] index for optimal (s,d) combination: " << opt_idx[0] << ", " \
            //   //      << opt_idx[1] << endl;
            //   // cout << " [-] Trajectories list for s: \n" << s_trajectories << endl;
            //   // cout << " [-] Trajectories list for d: \n" << d_trajectories << endl;
            //   // cout << " [*] Calculating optimal coeffs for s and d ..." << endl;
            //   optimal_s_coeff = s_trajectories.col(opt_idx[0]);
            //   optimal_d_coeff = d_trajectories.col(opt_idx[1]);
            //   cout << " [-] Optimal Trajectory for s: \n" << optimal_s_coeff << endl;
            //   cout << " [-] Optimal Trajectory for d: \n" << optimal_d_coeff << endl;
            //   for (int hrz=0; hrz<n_planning_horizon; hrz++){
            //     double s = getPosition(optimal_s_coeff, hrz*0.02);
            //     double d = getPosition(optimal_d_coeff, hrz*0.02);
            //     vector<double> xy = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            //     next_x_vals.push_back(xy[0]);
            //     next_y_vals.push_back(xy[1]);
            //   }
            // }
            // else {
            //
            // }

            // /* --------------------------------------------------------- */
            // // INTERPOLATING TRAJECTORY: USING SPLINE
            // /* --------------------------------------------------------- */
            // // TRANSFORM TRAJECTORY TO VEHICLE LOCAL COORDINATES
            // for(int i=0; i<next_x_vals.size(); i++) {
            //   double x_local = next_x_vals[i] - pos_x;
            //   double y_local = next_y_vals[i] - pos_y;
            //   next_x_vals[i] = x_local * cos(0 - pos_yaw) - y_local * sin(0 - pos_yaw);
            //   next_y_vals[i] = x_local * sin(0 - pos_yaw) + y_local * cos(0 - pos_yaw);
            // }
            //
            // // GET speed
            // double dx_total = next_x_vals[next_x_vals.size()] - next_x_vals[0];
            // double dy_total = next_y_vals[next_x_vals.size()] - next_y_vals[0];
            // double avg_speed = sqrt(dx_total*dx_total + dy_total*dy_total) \
            //                     / (3.0);
            // if (avg_speed < 1) {avg_speed = 1;}
            //
            // cout << " [*] AVG SPEED CALCULATED ! : " << avg_speed << endl;

            // // fit spline
            // // for (int i=0; i<next_x_vals.size(); i++) {
            // //   cout << " [-] next_x_vals= " << next_x_vals[i] << endl;
            // // }
            //
            // tk::spline s;
            // s.set_points(next_x_vals, next_y_vals);
            // cout << " [*] SPLINE FITTED !" << endl;
            //
            // // interpolated points
            // vector<double> interpolated_x_vals;
            // vector<double> interpolated_y_vals;
            // double target_x = 30.0;
            // double target_y = s(target_x);
            // double target_dist = sqrt(target_x*target_x + target_y*target_y);
            //
            // double x_reference_prev = 0;
            //
            // cout << " [*] INTERPOLATING POINTS ..." << endl;
            //
            // for(int i=0; i<n_planning_horizon-prev_path_size; i++) {
            //   double N = target_dist/(0.02*avg_speed); // MPH
            //   double x_reference = x_reference_prev + target_x/N;
            //   double y_reference = s(x_reference);
            //   x_reference_prev = x_reference;
            //
            //   // TRANSFORM TO GLOBAL COORDINATES
            //   double x_reference_global;
            //   double y_reference_global;
            //   x_reference_global = x_reference * cos(pos_yaw) - y_reference * sin(pos_yaw);
            //   y_reference_global = x_reference * sin(pos_yaw) + y_reference * cos(pos_yaw);
            //   x_reference_global += pos_x;
            //   y_reference_global += pos_y;
            //
            //   interpolated_x_vals.push_back(x_reference_global);
            //   interpolated_y_vals.push_back(y_reference_global);
            // }
            // cout << "next_x_vals.size(): " << next_x_vals.size() << endl;
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
