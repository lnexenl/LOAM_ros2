#ifndef POINTTYPE_H
#define POINTTYPE_H

#include <initializer_list>

#define PCL_NO_PRECOMPILE
#include <pcl/memory.h>
#include <pcl/pcl_macros.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

using namespace pcl;

struct EIGEN_ALIGN16 PointType {
public:
  PCL_ADD_POINT4D;
  float intensity = 0;
  float curvature = 0;
  double t = 0;
  uint8_t ring = 0;
  PCL_MAKE_ALIGNED_OPERATOR_NEW
  explicit PointType() = default;

  template <typename T> explicit PointType(T v) {
    x = static_cast<float>(v(0));
    y = static_cast<float>(v(1));
    z = static_cast<float>(v(2));
  }

  template <typename T> explicit PointType(T v0, T v1, T v2) {
    x = static_cast<float>(v0);
    y = static_cast<float>(v1);
    z = static_cast<float>(v2);
  }

  template <typename T> explicit PointType(std::initializer_list<T> list) {
    auto it = list.begin();
    x = static_cast<float>(*it++);
    y = static_cast<float>(*it++);
    z = static_cast<float>(*it);
  }
};

struct EIGEN_ALIGN16 OusterPoint {
  PCL_ADD_POINT4D;
  float intensity;
  uint32_t t;
  uint16_t reflectivity;
  uint8_t ring;
  uint16_t ambient;
  uint32_t range;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct EIGEN_ALIGN16 VelodynePoint {
  PCL_ADD_POINT4D;
  float intensity;
  float time;
  uint16_t ring;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct EIGEN_ALIGN16 PandarPoint {
  PCL_ADD_POINT4D;
  float intensity;
  double timestamp;
  uint16_t ring;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct EIGEN_ALIGN16 PointIRT {
public:
  PCL_ADD_POINT4D;
  float intensity;
  uint16_t ring;
  float time;
  PCL_MAKE_ALIGNED_OPERATOR_NEW
};

// clang-format off
///NOTE:This part must be put in global namespace
POINT_CLOUD_REGISTER_POINT_STRUCT(
  PointType,
  (float, x, x)
  (float, y, y)
  (float, z, z)
  (float, intensity, intensity)
  (float, curvature, curvature)
  (double, t, t)
  (std::uint8_t, ring, ring)
)

POINT_CLOUD_REGISTER_POINT_STRUCT(
    OusterPoint,
    (float, x, x)
    (float, y, y)
    (float, z, z)
    (float, intensity, intensity)
    // use std::uint32_t to avoid conflicting with pcl::uint32_t
    (std::uint32_t, t, t)
    (std::uint16_t, reflectivity, reflectivity)
    (std::uint8_t, ring, ring)
    (std::uint16_t, ambient, ambient)
    (std::uint32_t, range, range)
)

POINT_CLOUD_REGISTER_POINT_STRUCT(
    VelodynePoint,
    (float, x, x)
    (float, y, y)
    (float, z, z)
    (float, intensity, intensity)
    (float, time, time)
    (uint16_t, ring, ring)
)

POINT_CLOUD_REGISTER_POINT_STRUCT(
    PandarPoint,
    (float, x, x)
    (float, y, y)
    (float, z, z)
    (float, intensity, intensity)
    (double, timestamp, timestamp)
    (uint16_t, ring, ring)
)

POINT_CLOUD_REGISTER_POINT_STRUCT(
    PointIRT,
    (float, x, x)
    (float, y, y)
    (float, z, z)
    (float, intensity, intensity)
    (uint16_t, ring, ring)
    (float, time, time)
)
//clang-format on

#endif
