#ifndef _POLYSOLVER_H_
#define _POLYSOLVER_H_

//
// #include "polysolver.h"
double huber_loss(double speed, double target);

VectorXd getPolynomialCoeffs(\
  double s0, double s0_dot, double s0_ddot, \
  double st, double st_dot, double st_ddot, \
  double T_terminal);

int solvePolynomialsFullTerminalCond(double s0, double s0dot, double s0ddot, \
  double s1, double s1dot, double s1ddot, double kspeed, double max_speed, vector<double> Tjset, vector<double> ds1set, \
  MatrixXd &Trajectories, VectorXd &Costs);

int solvePolynomialsTwoTerminalCond(double s0, double s0dot, double s0ddot, \
  double s1dot, double s1ddot, double kspeed, double max_speed, vector<double> Tjset, vector<double> ds1dotset, \
  MatrixXd &Trajectories, VectorXd &Costs);

int VelocityKeepingTrajectories(double s0, double s0dot, double s0ddot, \
  double s1dot, double max_speed, MatrixXd &s_trajectories, VectorXd &s_costs);

int FollowingTrajectories(double s0, double s0dot, double s0ddot, \
    double s_lv0, double s_lv0dot, double max_speed, MatrixXd &s_trajectories, VectorXd &s_costs);

int lateralTrajectories(double d0, double d0dot, double d0ddot, \
  double d1, MatrixXd &d_trajectories, VectorXd &d_costs);

vector<int> optimalCombination(VectorXd s_costs, VectorXd d_costs);

double getPosition(VectorXd coeffs, double t);

double getVelocity(VectorXd coeffs, double t);

double getAcceleration(VectorXd coeffs, double t);


#endif
