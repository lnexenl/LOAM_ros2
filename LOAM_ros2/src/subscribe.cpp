#include <rclcpp/node.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp/utilities.hpp>
#include <Eigen/Core>

#include "PointType.h"

int main(int argc, char **argv) {
  rclcpp::init(argc, argv);
  auto node = rclcpp::Node::make_shared("subscribe");
  Eigen::Vector3f a = Eigen::Matrix<float, 3, 1>::Random().block(0, 0, 3, 1) - Eigen::Matrix<float, 3, 1>::Random().block(0, 0, 3, 1);
  return 0;
}