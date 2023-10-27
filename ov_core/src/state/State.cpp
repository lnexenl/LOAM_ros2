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

#include "State.h"

#include "types/HMQuat.h"

using namespace ov_core;
using namespace ov_type;
using namespace ov_msckf;

State::State(StateOptions &options) {
  // Save our options
  _options = options;

  // Append the imu to the state and covariance
  int current_id = 0;
  _imu = std::make_shared<IMU>();
  _imu->set_local_id(current_id);
  _variables.push_back(_imu);
  current_id += _imu->size();

  // Append lidar extrinsic to the state and covariance
  _calib_IMUtoLid = std::make_shared<PoseHM>();
  _calib_IMUtoLid->set_local_id(current_id);
  _variables.push_back(_calib_IMUtoLid);
  current_id += _calib_IMUtoLid->size();

  // Append the imu intrinsics to the state and covariance
  // NOTE: these need to be right "next" to the IMU state in the covariance
  // NOTE: since if calibrating these will evolve / be correlated during
  // propagation
  _calib_imu_dw = std::make_shared<Vec>(6);
  _calib_imu_da = std::make_shared<Vec>(6);
  if (options.imu_model == StateOptions::ImuModel::KALIBR) {
    // lower triangular of the matrix (column-wise)
    Eigen::Matrix<double, 6, 1> _imu_default = Eigen::Matrix<double, 6, 1>::Zero();
    _imu_default << 1.0, 0.0, 0.0, 1.0, 0.0, 1.0;
    _calib_imu_dw->set_value(_imu_default);
    _calib_imu_dw->set_fej(_imu_default);
    _calib_imu_da->set_value(_imu_default);
    _calib_imu_da->set_fej(_imu_default);
  } else {
    // upper triangular of the matrix (column-wise)
    Eigen::Matrix<double, 6, 1> _imu_default = Eigen::Matrix<double, 6, 1>::Zero();
    _imu_default << 1.0, 0.0, 0.0, 1.0, 0.0, 1.0;
    _calib_imu_dw->set_value(_imu_default);
    _calib_imu_dw->set_fej(_imu_default);
    _calib_imu_da->set_value(_imu_default);
    _calib_imu_da->set_fej(_imu_default);
  }
  _calib_imu_tg = std::make_shared<Vec>(9);
  _calib_imu_GYROtoIMU = std::make_shared<HMQuat>();
  _calib_imu_ACCtoIMU = std::make_shared<HMQuat>();
  if (options.do_calib_imu_intrinsics) {
    // Gyroscope dw
    _calib_imu_dw->set_local_id(current_id);
    _variables.push_back(_calib_imu_dw);
    current_id += _calib_imu_dw->size();

    // Accelerometer da
    _calib_imu_da->set_local_id(current_id);
    _variables.push_back(_calib_imu_da);
    current_id += _calib_imu_da->size();

    // Gyroscope gravity sensitivity
    if (options.do_calib_imu_g_sensitivity) {
      _calib_imu_tg->set_local_id(current_id);
      _variables.push_back(_calib_imu_tg);
      current_id += _calib_imu_tg->size();
    }

    // If kalibr model, R_GYROtoIMU is calibrated
    // If rpng model, R_ACCtoIMU is calibrated
    if (options.imu_model == StateOptions::ImuModel::KALIBR) {
      _calib_imu_GYROtoIMU->set_local_id(current_id);
      _variables.push_back(_calib_imu_GYROtoIMU);
      current_id += _calib_imu_GYROtoIMU->size();
    } else {
      _calib_imu_ACCtoIMU->set_local_id(current_id);
      _variables.push_back(_calib_imu_ACCtoIMU);
      current_id += _calib_imu_ACCtoIMU->size();
    }
  }
}
