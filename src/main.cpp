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
#include "spline.h"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

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

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// MY PERSONAL FUNCTIONS
VectorXd getPolynomialCoeffs(\
  double s0, double s0_dot, double s0_ddot, \
  double st, double st_dot, double st_ddot, \
  double T_terminal) {
  // return polynomial coeffs p = a0 + a1 * t + ... a5 * t^5
  VectorXd coeffs(6);

  double a0 = s0;
  double a1 = s0_dot;
  double a2 = s0_ddot / 2.0;

  double T = T_terminal;

  Matrix3d c345;
  c345 << pow(T, 3.), pow(T, 4.), pow(T, 5.),
          3*pow(T, 2.), 4*pow(T, 3.), 5*pow(T,4.),
          6*T, 12*pow(T, 2.), 20*pow(T, 3.);

  Vector3d q345;
  q345 << st - (a0 + a1 * T + a2 * pow(T, 2.)),
          st_dot - (a1 + 2 * a2 * T),
          st_ddot - 2 * a2;

  Vector3d a345;
  a345 = c345.inverse() * q345;

  coeffs << a0, a1, a2, a345(0), a345(1), a345(2);

  return coeffs;
}

int solvePolynomialsFullTerminalCond(double s0, double s0dot, double s0ddot, \
double s1, double s1dot, double s1ddot, vector<double> Tjset, vector<double> ds1set, \
MatrixXd &Trajectories, VectorXd &Costs) {
  Matrix3d M1;
  M1 << 1, 0, 0,
        0, 1, 0,
        0, 0, 2;

  Vector3d zeta0(s0, s0dot, s0ddot);
  Vector3d c012;
  c012 = M1.inverse() * zeta0;

  MatrixXd coeffs(6, ds1set.size() * Tjset.size());
  VectorXd costs(ds1set.size() * Tjset.size());

  double kj = 1.0; // weight for jerk
  double kt = 0; // weight for time
  double ks = 0; // weight for terminal state

  double acc_thres = 6.0; // acceleration threshold < 6m/s^2

  int N_ACCELERATION_CONDITION_SATISFIED = 0;

  for (int i=0; i<ds1set.size(); i++) {
    for (int j=0; j<Tjset.size(); j++) {
      // set [s1+ds1_i, s1dot, s1ddot, Tj]
      double ds1 = ds1set[i];
      double Tj = Tjset[j];
      double s1_target = s1 + ds1;

      // solve c3, c4, c5
      Matrix3d M1_T;
      Matrix3d M2_T;

      M1_T << 1, Tj, Tj*Tj,
                0, 1, 2* Tj,
                0, 0, 2;
      M2_T << Tj*Tj*Tj, Tj*Tj*Tj*Tj, Tj*Tj*Tj*Tj*Tj,
                3*Tj*Tj, 4*Tj*Tj*Tj, 5*Tj*Tj*Tj*Tj,
                6*Tj, 12*Tj*Tj, 20*Tj*Tj*Tj;
      Vector3d zeta1(s1_target, s1dot, s1ddot);
      Vector3d c345;
      c345 = M2_T.inverse() * (zeta1 - M1_T * c012);

      // calculate accleration outsized
      double acc0 = abs(6 * c345(0));
      double acc1 = abs(6 * c345(0) + 24* c345(1) * Tj + 60 * c345(2) * Tj *Tj);
      double maxp = - c345(1) / (5.0 * c345(2));
      double acc2 = abs(6 * c345(0) + 24 * c345(1) * maxp + 60 * c345(2) * maxp*maxp);
      bool ACCELERATION_OUTSIZED = (acc0 >= acc_thres) || (acc1 >= acc_thres) || (acc2 >= acc_thres);

      cout << "acceleration: " <<  max(acc0, acc1) << "\n outsized?" << ACCELERATION_OUTSIZED << endl;

      if (ACCELERATION_OUTSIZED == false){

        // calculate cost
        double cost_Jerk = 36 * c345(0) * c345(0) * Tj \
                          + 144 * c345(0) * c345(1) * Tj * Tj \
                          + (192 * c345(1) * c345(1) + 240 * c345(0) *c345(2)) * Tj * Tj * Tj \
                          + 720 * c345(1) * c345(2) * Tj*Tj*Tj*Tj \
                          + 720 * c345(2) * c345(2) * Tj*Tj*Tj*Tj*Tj;
        double cost_terminal = ds1 * ds1;
        double cost_Total = kj * cost_Jerk + kt * Tj + ks * cost_terminal;

        // cout << "c34: " << endl;
        // cout << c34 << "\n" << endl;
        //
        // cout << "cost_Jerk:  " << endl;
        // cout << cost_Jerk << "\n" << endl;
        //
        // cout << "cost: " << cost_Total <<  "\n" << endl;

        // append coeffs and cost to trajectory sets
        VectorXd c012345(6);
        c012345 << c012(0), c012(1), c012(2), c345(0), c345(1), c345(2);
        coeffs.col(N_ACCELERATION_CONDITION_SATISFIED) = c012345;
        costs(N_ACCELERATION_CONDITION_SATISFIED) = cost_Total;

        Trajectories.conservativeResize(6, Trajectories.cols()+1);
        Trajectories.col(Trajectories.cols()-1) = c012345;
        Costs.conservativeResize(Costs.size()+1);
        Costs(Costs.size()-1) = cost_Total;


        N_ACCELERATION_CONDITION_SATISFIED += 1;
      }

    }
  }

  return 0;


}

int solvePolynomialsTwoTerminalCond(double s0, double s0dot, double s0ddot, \
double s1dot, double s1ddot, vector<double> Tjset, vector<double> ds1dotset, \
MatrixXd &Trajectories, VectorXd &Costs) {
  cout << "ds1set size: " << ds1dotset.size() << endl;
  cout << "Tjset size: " << Tjset.size() << "\n" << endl;

  Matrix3d M1;
  M1 << 1, 0, 0,
        0, 1, 0,
        0, 0, 2;

  Vector3d zeta0(s0, s0dot, s0ddot);
  Vector3d c012;
  c012 = M1.inverse() * zeta0;
  Vector2d c12(c012(1), c012(2));

  MatrixXd coeffs(6, ds1dotset.size() * Tjset.size());
  VectorXd costs(ds1dotset.size() * Tjset.size());
  double kj = 1.0; // weight for jerk
  double kt = 0; // weight for time
  double ksdot = 0; // weight for terminal state

  double acc_thres = 6.0; // acceleration threshold < 6m/s^2

  int N_ACCELERATION_CONDITION_SATISFIED = 0;

  for (int i=0; i<ds1dotset.size(); i++) {
    for (int j=0; j<Tjset.size(); j++) {
      // set target state and terminal time
      double ds1dot = ds1dotset[i];
      double Tj = Tjset[j];
      double s1dot_target = s1dot + ds1dot;

      // solve c3, c4 (c5 = 0 by transversality condition)
      Vector2d zeta1(s1dot_target, s1ddot);
      Matrix2d M2;
      Matrix2d C2;
      M2 << 3 * Tj*Tj, 4 * Tj*Tj*Tj,
            6 * Tj, 12 * Tj*Tj;
      C2 << 1, 2*Tj,
            0, 2;

      Vector2d c34;
      c34 = M2.inverse() * (zeta1 - C2*c12);

      // check acceleration outsized
      double acc0 = abs(6 * c34(0));
      double acc1 = abs(6 * c34(0) + 24* c34(1) * Tj);
      bool ACCELERATION_OUTSIZED = (acc0 >= acc_thres) || (acc1 >= acc_thres);

      cout << "acceleration: " <<  max(acc0, acc1) << "\n outsized?" << ACCELERATION_OUTSIZED << endl;

      if (ACCELERATION_OUTSIZED == false){

        // calculate cost
        double cost_Jerk = 36 * c34(0) * c34(0) * Tj \
                          + 144 * c34(0) * c34(1) * Tj * Tj \
                          + 192 * c34(1) * c34(1) * Tj * Tj * Tj;
        double cost_terminal = ds1dot * ds1dot;
        double cost_Total = kj * cost_Jerk + kt * Tj + ksdot * cost_terminal;

        // cout << "c34: " << endl;
        // cout << c34 << "\n" << endl;
        //
        // cout << "cost_Jerk:  " << endl;
        // cout << cost_Jerk << "\n" << endl;
        //
        // cout << "cost: " << cost_Total <<  "\n" << endl;

        // append coeffs and cost to trajectory sets
        VectorXd c012345(6);
        c012345 << c012(0), c012(1), c012(2), c34(0), c34(1), 0;
        coeffs.col(N_ACCELERATION_CONDITION_SATISFIED) = c012345;
        costs(N_ACCELERATION_CONDITION_SATISFIED) = cost_Total;

        Trajectories.conservativeResize(6, Trajectories.cols()+1);
        Trajectories.col(Trajectories.cols()-1) = c012345;
        Costs.conservativeResize(Costs.size()+1);
        Costs(Costs.size()-1) = cost_Total;


        N_ACCELERATION_CONDITION_SATISFIED += 1;
      }//ifend
    }//for j (Tj) end
  }
  cout << "\nN_ACCELERATION_CONDITION_SATISFIED: " << N_ACCELERATION_CONDITION_SATISFIED << endl;

  cout << "*************\ncoeffs:\n" << coeffs << endl;
  cout << "*************\ncosts:\n" << costs << endl;

  // resize trajectory coeffs and costs
  // Map<VectorXd, 0> trajectory_costs(costs.data(), N_ACCELERATION_CONDITION_SATISFIED);
  // Map<MatrixXd, 0> trajectory_coeffs(coeffs.data(), coeffs.rows(), N_ACCELERATION_CONDITION_SATISFIED);
  // cout << "trajectory_costs: \n" << trajectory_costs << endl;
  // cout << "trajectory_coeffs: \n" << trajectory_coeffs << endl;



  // // append to trajectories and costs
  // MatrixXd newTrajectories(Trajectories.rows(), Trajectories.cols()+N_ACCELERATION_CONDITION_SATISFIED);
  //
  // Trajectories.conservativeResize(Trajectories.rows(), Trajectories.cols()+N_ACCELERATION_CONDITION_SATISFIED);
  // Trajectories.col()

  return 0;
}


int VelocityKeepingTrajectories() {
  vector<double> target_speed_sets = {22.5, 20, 18.5};
  cout << target_speed_sets[1] << endl;



  return 0;
}

int main() {

  // VectorXd test_ = getPolynomialCoeffs(1, 0.1, 0.7, 20, 5, 0, 5);
  // cout << test_ << endl;
  vector<double> Tjset = {1,2,3};
  vector<double> ds1set = {0, 1, 2, 3};
  MatrixXd Trajectoies(6,0);
  VectorXd Costs(0);
  // double i_ = solvePolynomialsTwoTerminalCond(0, 0, 0, 2, 0, Tjset, ds1set, Trajectoies, Costs);
  int asdf = solvePolynomialsFullTerminalCond(0, 0, 0, 3, 3, 0, Tjset, ds1set, Trajectoies, Costs);
  cout << "\n ~~~~~~~~~~~~~~ \n Trajectoies: \n" <<  Trajectoies << endl;
  cout << "\n ~~~~~~~~~~~~~~ \n Costs: \n" <<  Costs << endl;

  uWS::Hub h;

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

  double target_speed = 1.0; // m/s
  double lane = 1.0;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &target_speed, &lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

            int prev_path_size = previous_path_x.size();



            // sensor Fusion
            bool FOLLOWING = false;
            bool VELOCITY_KEEPING = true;
            double safe_distance = 50.0;

            for (int i=0; i<sensor_fusion.size(); i++) {
              float d = sensor_fusion[i][6];
              bool is_car_in_my_lane = ((d < 2+4*lane+2) && (d >2+4*lane-2));
              if (is_car_in_my_lane) {
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double speed_other = sqrt(vx*vx + vy*vy);
                double s_other = sensor_fusion[i][5];

                s_other += (double)prev_path_size*0.02*speed_other;

                bool is_car_close_front_of_me = (s_other > car_s) \
                              && (s_other - car_s < safe_distance + 15.0);

                if (is_car_close_front_of_me) {
                  FOLLOWING = true;
                  VELOCITY_KEEPING = false;
                }

              }
            }
            if (FOLLOWING) {
              target_speed -= 0.224;

            }
            else if (target_speed < 48) {
              target_speed += 0.224;
            }

            // path planning
            double pos_x;
            double pos_y;
            double pos_yaw;

            vector<double> points_for_ref_x;
            vector<double> points_for_ref_y;

            // set previous path
            for(int i=0; i<prev_path_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // get the remaining trajectory starting point
            if(prev_path_size == 0) {
              pos_x = car_x;
              pos_y = car_y;
              pos_yaw = deg2rad(car_yaw);
              points_for_ref_x.push_back(pos_x);
              points_for_ref_y.push_back(pos_y);
            }
            else {
              pos_x = previous_path_x[prev_path_size-1];
              pos_y = previous_path_y[prev_path_size-1];

              double pos_x_prev = previous_path_x[prev_path_size-2];
              double pos_y_prev = previous_path_y[prev_path_size-2];
              pos_yaw = atan2(pos_y-pos_y_prev, pos_x-pos_x_prev);

              points_for_ref_x.push_back(pos_x_prev);
              points_for_ref_x.push_back(pos_x);
              points_for_ref_y.push_back(pos_y_prev);
              points_for_ref_y.push_back(pos_y);
            }

            // get points for making reference trajectory
            vector<double> next_waypoints_1 = getXY(car_s + 40, 2+4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_waypoints_2 = getXY(car_s + 80, 2+4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_waypoints_3 = getXY(car_s + 120, 2+4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            points_for_ref_x.push_back(next_waypoints_1[0]);
            points_for_ref_x.push_back(next_waypoints_2[0]);
            points_for_ref_x.push_back(next_waypoints_3[0]);

            points_for_ref_y.push_back(next_waypoints_1[1]);
            points_for_ref_y.push_back(next_waypoints_2[1]);
            points_for_ref_y.push_back(next_waypoints_3[1]);

            // TRANSFORM TO VEHICLE LOCAL COORDINATES
            for(int i=0; i<points_for_ref_x.size(); i++) {
              double x_local = points_for_ref_x[i] - pos_x;
              double y_local = points_for_ref_y[i] - pos_y;
              points_for_ref_x[i] = x_local * cos(0 - pos_yaw) - y_local * sin(0 - pos_yaw);
              points_for_ref_y[i] = x_local * sin(0 - pos_yaw) + y_local * cos(0 - pos_yaw);
            }

            // fit spline
            tk::spline s;
            s.set_points(points_for_ref_x, points_for_ref_y);

            // get trajectory points

            double target_x = 30.0;
            double target_y  = s(target_x);
            double target_dist = sqrt(target_x*target_x + target_y*target_y);

            double x_reference_prev = 0;

            // cout << " ------------------------- " << endl;

            for(int i=0; i<50-prev_path_size; i++) {
              double N = target_dist/(0.02*target_speed/2.24); // MPH
              double x_reference = x_reference_prev + target_x/N;
              double y_reference = s(x_reference);

              x_reference_prev = x_reference;
              // cout << "(x,y):  " << x_reference << ", " <<  y_reference << endl;


              // TRANSFORM TO GLOBAL COORDINATES
              double x_reference_global;
              double y_reference_global;
              x_reference_global = x_reference * cos(pos_yaw) - y_reference * sin(pos_yaw);
              y_reference_global = x_reference * sin(pos_yaw) + y_reference * cos(pos_yaw);
              x_reference_global += pos_x;
              y_reference_global += pos_y;

              next_x_vals.push_back(x_reference_global);
              next_y_vals.push_back(y_reference_global);
            }

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

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
