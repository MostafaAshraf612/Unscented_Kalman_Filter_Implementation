# 🚗 Unscented Kalman Filter (UKF)

An implementation of the **Unscented Kalman Filter** for sensor fusion, capable of accurately tracking a moving object using **Lidar** and **Radar** measurements.  

This project was developed as part of the **Udacity Self-Driving Car Engineer Nanodegree** and demonstrates robust state estimation in non-linear systems.

---

## ✨ Key Features
- **Sensor Fusion**: Supports both Lidar (px, py) and Radar (ρ, φ, ρ̇).  
- **Non-linear Handling**: Uses unscented transform with augmented sigma points.  
- **Noise Modeling**: Incorporates process and measurement noise.  
- **Angle Normalization**: Ensures numerical stability in angle computations.  
- **Modular Design**: Clean separation of prediction and update steps.  

---

## 🛠️ Build Instructions

### Prerequisites
- C++11 or higher  
- CMake ≥ 3.5  
- Eigen3 library  

### Build & Run
```bash
git clone https://github.com/MostafaAshraf612/Unscented-Kalman-Filter_Implementation.git
cd unscented-kalman-filter
mkdir build && cd build
cmake ..
make
./UnscentedKF
