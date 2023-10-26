/*
 * OpenVINS: An Open Platform for Visual-Inertial Research
 * Copyright (C) 2018-2023 Patrick Geneva
 * Copyright (C) 2018-2023 Guoquan Huang
 * Copyright (C) 2018-2023 OpenVINS Contributors
 * Copyright (C) 2018-2019 Kevin Eckenhoff
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef OV_MSCKF_STATE_OPTIONS_H
#define OV_MSCKF_STATE_OPTIONS_H

#include <memory>

#include "types/LandmarkRepresentation.h"
#include "utils/logger.h"
#include "utils/opencv_yaml_parse.h"
#include "utils/sensor_data.h"

namespace ov_msckf {

/**
 * @brief Struct which stores all our filter options
 */
struct StateOptions {
  /// Bool to determine whether or not to do first estimate Jacobians
  bool do_fej = true;

  /// Numerical integration methods
  enum IntegrationMethod { DISCRETE, RK4, ANALYTICAL };

  /// What type of numerical integration is used during propagation
  IntegrationMethod integration_method = IntegrationMethod::RK4;

  /// Bool to determine whether or not to calibrate imu-to-lidar pose
  bool do_calib_lidar_pose = false;

  /// Bool to determine whether or not to calibrate the IMU intrinsics
  bool do_calib_imu_intrinsics = false;

  /// Bool to determine whether or not to calibrate the Gravity sensitivity
  bool do_calib_imu_g_sensitivity = false;

  /// IMU intrinsic models
  enum ImuModel { KALIBR, RPNG };

  /// What model our IMU intrinsics are
  ImuModel imu_model = ImuModel::KALIBR;

  std::shared_ptr<Logger> logger = std::make_shared<Logger>("StateOptions");

  /// Nice print function of what parameters we have loaded
  void print(const std::shared_ptr<ov_core::YamlParser> &parser = nullptr) {
    if (parser != nullptr) {
      parser->parse_config("use_fej", do_fej);

      // Integration method
      std::string integration_str = "rk4";
      parser->parse_config("integration", integration_str);
      if (integration_str == "discrete") {
        integration_method = IntegrationMethod::DISCRETE;
      } else if (integration_str == "rk4") {
        integration_method = IntegrationMethod::RK4;
      } else if (integration_str == "analytical") {
        integration_method = IntegrationMethod::ANALYTICAL;
      } else {
        logger->error("invalid imu integration model: {}", integration_str);
        logger->error("please select a valid model: discrete, rk4, analytical");
        std::exit(EXIT_FAILURE);
      }

      // Calibration booleans
      parser->parse_config("calib_lidar_extrinsics", do_calib_lidar_pose);
      parser->parse_config("calib_imu_intrinsics", do_calib_imu_intrinsics);
      parser->parse_config("calib_imu_g_sensitivity", do_calib_imu_g_sensitivity);

      // IMU model
      std::string imu_model_str = "kalibr";
      parser->parse_external("relative_config_imu", "imu0", "model", imu_model_str);
      if (imu_model_str == "kalibr" || imu_model_str == "calibrated") {
        imu_model = ImuModel::KALIBR;
      } else if (imu_model_str == "rpng") {
        imu_model = ImuModel::RPNG;
      } else {
        logger->error("invalid imu model: {}", imu_model_str);
        logger->error("please select a valid model: kalibr, rpng");
        std::exit(EXIT_FAILURE);
      }
      if (imu_model_str == "calibrated" && (do_calib_imu_intrinsics || do_calib_imu_g_sensitivity)) {
        logger->error("calibrated IMU model selected, but requested calibration!");
        logger->error("please select what model you have: kalibr, rpng");
        std::exit(EXIT_FAILURE);
      }
    }
    logger->debug("  - use_fej: {}", do_fej);
    logger->debug("  - integration: {}", integration_method);
    logger->debug("  - calib_lidar_extrinsics: {}", do_calib_lidar_pose);
    logger->debug("  - calib_imu_intrinsics: {}", do_calib_imu_intrinsics);
    logger->debug("  - calib_imu_g_sensitivity: {}", do_calib_imu_g_sensitivity);
    logger->debug("  - imu_model: {}", imu_model);
  }
};

} // namespace ov_msckf

#endif // OV_MSCKF_STATE_OPTIONS_H