#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "Eigen-3.3/Eigen/LU"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Global constants
const double MAX_S = 6945.554; // meters
const double LANE_WIDTH = 4.0; // meters

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
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
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

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
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

// lane tracking in Fernet coordinates
vector<vector<double>> laneFollowing(double car_s, vector<double> previous_path_x, vector<double> previous_path_y, 
                                      vector<double> map_waypoints_x, vector<double> map_waypoints_y, vector<double> map_waypoints_s, bool verbose)
{
  
  double target_lane = 1.0;
  double target_d = LANE_WIDTH * (0.5 + target_lane);
  double dist_inc = 0.5;

  int path_size = previous_path_x.size();

  if (verbose) {cout << "Adding " << path_size << " elements from previous path." << endl << "Previous elements: " << endl;}


  vector<double> next_x_vals;
  vector<double> next_y_vals;

  for(int i = 0; i < path_size; i++) {

      if (verbose) { cout << "#" << i << "(" << previous_path_x[i] << "," << previous_path_y[i] << ")" << endl; }
      next_x_vals.push_back(previous_path_x[i]);
      next_y_vals.push_back(previous_path_y[i]);
  }

  double prev_s;
  if (path_size == 0) {
    prev_s = car_s;
  } else {
    double y1 = previous_path_y[path_size-1];
    double y2 = previous_path_y[path_size-2];
    double x1 = previous_path_y[path_size-1];
    double x2 = previous_path_y[path_size-2];
    double prev_yaw = atan2(y1-y2, x1-x2); 
    vector<double> sd = getFrenet(previous_path_x[path_size-1], previous_path_y[path_size-1], prev_yaw, map_waypoints_x, map_waypoints_y);
    prev_s = sd[0];
  }

  // new points along the centerine
  if (verbose) { cout << "New points:" << endl; }
  for(int i = 1; i < 50 - path_size; i++)
  {
      double target_s = fmod(prev_s + dist_inc*i, MAX_S);
      vector<double> xy = getXY(target_s, target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
      if (verbose) {
        cout << "target_s = " << target_s << "m | target_d = " << target_d << "m  | ";
        cout << "xy[0] = " << xy[0] << " | xy[1] = " << xy[1] << endl;
      }
      next_x_vals.push_back(xy[0]);
      next_y_vals.push_back(xy[1]);
  }

  if (verbose) { cout << "------------" << endl << endl; }
  return {next_x_vals, next_y_vals};
} 



// this functions checks is a lane change is safe
bool change_lane_safe(double car_s, double car_speed, int targetLane, 
  vector<vector<double>> sensor_fusion, double safe_dist) {
  
  bool change_safe = true;
  for (int i=0; i<sensor_fusion.size(); i++) {
    float d = sensor_fusion[i][6];

    // check for vehicles in our lane
    bool in_target_lane = d < LANE_WIDTH*(targetLane+1) && d > LANE_WIDTH*(targetLane);
    double check_car_s = sensor_fusion[i][5];
    double rel_range = car_s-check_car_s;
    double obj_speed = sqrt(sensor_fusion[i][3]*sensor_fusion[i][3] + sensor_fusion[i][4]*sensor_fusion[i][4]);
    double obj_rel_speed = obj_speed - car_speed;
    double time_left = rel_range / obj_rel_speed;

    bool collision_time = fabs(time_left) < 5; // seconds

    if (in_target_lane && (fabs(rel_range) < safe_dist)) { // || collision_time)) {
		// only change lane if it is safe
      change_safe = false;
    } 
  }  

  return change_safe;
}


int main() {
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

  int lane = 1;

  double ref_vel = 0.0; // mph
  const double speed_limit = 49.5;

  h.onMessage([&ref_vel, &speed_limit, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

            int prev_size = previous_path_x.size();

            if (prev_size > 0) {
              car_s = end_path_s;
            }

            double ref_mph = 2.24;
            bool obstacle_detected = false;


            double speed_constraint = speed_limit; // mph
            for (int i=0; i<sensor_fusion.size(); i++) {
              float d = sensor_fusion[i][6];
              bool in_my_lane = d<(2+4*lane+2) && d > (2+4*lane-2);
              if (in_my_lane) {
                double v_x = sensor_fusion[i][3];
                double v_y = sensor_fusion[i][4];
                double check_speed = sqrt(v_x*v_x + v_y*v_y);
                double check_car_s = sensor_fusion[i][5];

                check_car_s += ((double) prev_size)*0.02*check_speed;
                bool in_front = check_car_s > car_s;
                double rel_range = check_car_s-car_s;
				// track aheads car speed if it is within 30 m distance
				// slow down if it is driving slower then us
				// if the car is driving too slow we try to pass it
                if (in_front && rel_range < 30) {
                  speed_constraint = min(speed_constraint, check_speed * ref_mph);
                  
                  if (speed_limit - speed_constraint > 3) {
                    obstacle_detected = true;
                  }
                }
              }
            }
// code to switch lanes
            if (obstacle_detected) {

              bool can_pass_left = lane > 0;
              bool can_pass_right = lane < 2;
              bool safe_left_pass;

              // check if we can safely change to the left lane 
              if (can_pass_left) {
                safe_left_pass = change_lane_safe(car_s, car_speed, lane-1, sensor_fusion, 20);
              } else {
                safe_left_pass = false;
              }

              // check if we can safely change to the right lane
              bool safe_right_pass;
              if (can_pass_right) {
                safe_right_pass = change_lane_safe(car_s, car_speed, lane+1, sensor_fusion, 20);
              } else {
                safe_right_pass = false;
              }

			// decide if it is safe to pass left, if not try right, if not maintain safe distance
              if (can_pass_left && safe_left_pass) {
                lane -= 1;
              } else if (can_pass_right && safe_right_pass) {
                lane += 1;
              } else {
                if (ref_vel >= speed_constraint) {
                  ref_vel -= 0.244;
                }

              }
            } else {
              if (ref_vel < speed_constraint) {
                cout << "  Accelerating to speed.";
                ref_vel += 0.224;
              }


              // Drive on the right lane if possible
              bool close_to_limit = (speed_limit - ref_vel) < 5; // mph
              bool safeToGetOver = change_lane_safe(car_s, car_speed, lane+1, sensor_fusion, 50);
              if (lane < 2 && close_to_limit && safeToGetOver) {
                lane += 1;
                cout << " Changing lane to " << lane <<".";
              }
              cout << endl;

            }


			// the following code is inspired by the walkthrough video
            vector<double> ptsx, ptsy;
            double ref_x = car_x;
            double ref_y = car_d;
            double ref_yaw = (car_yaw);
          	
			// for cold start
            if (prev_size < 2) {
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);

              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);
            } else {
              // for warm Start
              ref_x = previous_path_x[prev_size-1];
              ref_y = previous_path_y[prev_size-1];

              double ref_x_prev = previous_path_x[prev_size-2];
              double ref_y_prev = previous_path_y[prev_size-2];

              ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

              ptsx.push_back(ref_x_prev);
              ptsx.push_back(ref_x);

              ptsy.push_back(ref_y_prev);
              ptsy.push_back(ref_y);
            }

            // Sample distant waypoints 
            vector<double> next_wp0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          
            // Add waypoints to anchor set.
            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);


            // transform the anchor points
            for (int i=0; i<ptsx.size(); i++) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
              ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
            }

            // spline
            tk::spline s;

            s.set_points(ptsx, ptsy);

            vector<double> next_x_vals;
            vector<double> next_y_vals;


            // connect to the previous points
            for (int i=0; i<previous_path_x.size(); i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // add new points with spline and anchor points
            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt(target_x*target_x + target_y*target_y);

            double x_add_on = 0.0;
            for (int i=1; i<=50 - previous_path_x.size(); i++) {
              double N = (target_dist) / (0.02 * ref_vel /2.24);
              double x_point = x_add_on + target_x / N;
              double y_point = s(x_point);
              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              // transform coordinates
              x_point = x_ref * cos(ref_yaw) - y_ref*sin(ref_yaw);
              y_point = x_ref * sin(ref_yaw) + y_ref*cos(ref_yaw);

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
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
















































































