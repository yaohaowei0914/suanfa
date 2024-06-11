#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>

// 高斯核函数
double gaussian_kernel(double u) {
    return std::exp(-0.5 * u * u) / std::sqrt(2 * M_PI);
}

// 盒型核函数
double box_kernel(double u) {
    return std::abs(u) <= 1 ? 0.5 : 0.0;
}

// 三角核函数
double triangle_kernel(double u) {
    return std::abs(u) <= 1 ? 1 - std::abs(u) : 0.0;
}

// Epanechnikov核函数
double epanechnikov_kernel(double u) {
    return std::abs(u) <= 1 ? 0.75 * (1 - u * u) : 0.0;
}

// 核回归
std::vector<double> kernel_regression(const std::vector<double>& x, const std::vector<double>& y, double bandwidth, const std::vector<double>& x_new, std::function<double(double)> kernel) {
    int n = x.size();
    int m = x_new.size();
    std::vector<double> y_new(m);

    for (int i = 0; i < m; ++i) {
        double x0 = x_new[i];
        std::vector<double> weights(n);
        for (int j = 0; j < n; ++j) {
            weights[j] = kernel((x0 - x[j]) / bandwidth);
        }
        double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        double weighted_sum = std::inner_product(weights.begin(), weights.end(), y.begin(), 0.0);
        y_new[i] = weighted_sum / weight_sum;
    }

    return y_new;
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

    // 带宽
    double bandwidth = 1.0;

    // 生成新数据点
    std::vector<double> x_new(100);
    double step = (x.back() - x.front()) / (x_new.size() - 1);
    std::generate(x_new.begin(), x_new.end(), [n = 0, &step]() mutable { return n++ * step; });

    // 应用核回归算法
    std::vector<double> y_new_gaussian, y_new_box, y_new_triangle, y_new_epanechnikov;
    try {
        y_new_gaussian = kernel_regression(x, y, bandwidth, x_new, gaussian_kernel);
        y_new_box = kernel_regression(x, y, bandwidth, x_new, box_kernel);
        y_new_triangle = kernel_regression(x, y, bandwidth, x_new, triangle_kernel);
        y_new_epanechnikov = kernel_regression(x, y, bandwidth, x_new, epanechnikov_kernel);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // 将结果写入输出文件
    std::ofstream outfile("output2.txt");
    if (outfile.is_open()) {
        outfile << "Gaussian Kernel Regression:\n";
        for (size_t i = 0; i < y_new_gaussian.size(); ++i) {
            outfile << x_new[i] << " " << y_new_gaussian[i] << "\n";
        }
        outfile << "\nBox Kernel Regression:\n";
        for (size_t i = 0; i < y_new_box.size(); ++i) {
            outfile << x_new[i] << " " << y_new_box[i] << "\n";
        }
        outfile << "\nTriangle Kernel Regression:\n";
        for (size_t i = 0; i < y_new_triangle.size(); ++i) {
            outfile << x_new[i] << " " << y_new_triangle[i] << "\n";
        }
        outfile << "\nEpanechnikov Kernel Regression:\n";
        for (size_t i = 0; i < y_new_epanechnikov.size(); ++i) {
            outfile << x_new[i] << " " << y_new_epanechnikov[i] << "\n";
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
