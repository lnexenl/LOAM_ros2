#include <Eigen/Core>
#include <rclcpp/node.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp/utilities.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include "PointType.h"
void lidar_callback(const sensor_msgs::msg::PointCloud2::SharedPtr msg) {
  std::cout << "Got point cloud with " << msg->width << " points\n";
}

void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg) { std::cout << "Got imu message\n"; }

int main(int argc, char **argv) {
  rclcpp::init(argc, argv);
  auto node = rclcpp::Node::make_shared("subscribe");
  auto lidar_sub = node->create_subscription<sensor_msgs::msg::PointCloud2>("/os_cloud_node/points", 1000, lidar_callback);
  auto imu_sub = node->create_subscription<sensor_msgs::msg::Imu>("/os_cloud_node/imu", 1000, imu_callback);
  return 0;
}