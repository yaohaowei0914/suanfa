#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "third_party/eigen/Eigen/Dense"

// 样条保序回归
std::vector<double> monotone_spline_regression(const std::vector<double>& x, const std::vector<double>& y, double s = 0) {
    int n = x.size();
    Eigen::VectorXd X(n), Y(n);
    for (int i = 0; i < n; ++i) {
        X(i) = x[i];
        Y(i) = y[i];
    }

    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = std::pow(X(i) - X(j), 3);
        }
    }

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n, 4);
    for (int i = 0; i < n; ++i) {
        P(i, 0) = 1;
        P(i, 1) = X(i);
        P(i, 2) = X(i) * X(i);
        P(i, 3) = X(i) * X(i) * X(i);
    }

    Eigen::MatrixXd PT = P.transpose();
    Eigen::MatrixXd PT_P = PT * P;
    Eigen::MatrixXd PT_Y = PT * Y;

    Eigen::MatrixXd A_inv = A.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::VectorXd spline = A_inv * (Y - P * PT_P.completeOrthogonalDecomposition().pseudoInverse() * PT_Y);

    std::vector<double> y_spline(n);
    for (int i = 0; i < n; ++i) {
        y_spline[i] = spline(i);
    }

    // 保序约束
    for (int i = 1; i < n; ++i) {
        y_spline[i] = std::max(y_spline[i], y_spline[i - 1]);
    }

    return y_spline;
}

int main() {
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();

    // 读取输入文件
    std::ifstream infile("input2.txt");
    std::vector<double> x, y;
    double value_x, value_y;
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> value_x >> value_y) {
                x.push_back(value_x);
                y.push_back(value_y);
            }
        }
        infile.close();
    } else {
        std::cerr << "无法打开输入文件" << std::endl;
        return 1;
    }

    // 应用样条保序回归算法
    std::vector<double> y_spline;
    try {
        y_spline = monotone_spline_regression(x, y);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // 将结果写入输出文件
    std::ofstream outfile("output.txt");
    if (outfile.is_open()) {
        for (size_t i = 0; i < y_spline.size(); ++i) {
            outfile << x[i] << " " << y_spline[i] << "\n";
        }
        outfile.close();
    } else {
        std::cerr << "无法打开输出文件" << std::endl;
        return 1;
    }

    // 记录结束时间
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // 输出运行时间
    std::cout << "程序运行时间: " << duration.count() << " 秒" << std::endl;

    return 0;
}
