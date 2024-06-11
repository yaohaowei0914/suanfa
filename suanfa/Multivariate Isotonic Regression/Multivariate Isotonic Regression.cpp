#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <functional>

// 一维保序回归 (使用池化相邻违背者算法)
std::vector<double> isotonic_regression_1d(const std::vector<double>& y) {
    int n = y.size();
    std::vector<double> ghat(y);
    std::vector<int> stack(n);
    int stack_size = 0;

    for (int i = 0; i < n; ++i) {
        stack[stack_size] = i;
        stack_size++;
        while (stack_size > 1 && ghat[stack[stack_size - 1]] < ghat[stack[stack_size - 2]]) {
            stack_size--;
            int k = stack[stack_size - 1];
            int j = stack[stack_size];
            double mean = std::accumulate(ghat.begin() + k, ghat.begin() + j + 1, 0.0) / (j - k + 1);
            std::fill(ghat.begin() + k, ghat.begin() + j + 1, mean);
        }
    }

    return ghat;
}

// 使用动态规划实现多元保序回归
std::vector<double> multivariate_isotonic_regression_dp(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
    int n_samples = X.size();
    int n_features = X[0].size();

    std::vector<double> f_hat(n_samples, 0.0);

    // 按维度递归求解保序回归
    std::function<std::vector<double>(int, std::vector<int>&)> dp_recursive = [&](int dim, std::vector<int>& indices) {
        if (dim == n_features) {
            std::vector<double> result(indices.size());
            for (size_t i = 0; i < indices.size(); ++i) {
                result[i] = y[indices[i]];
            }
            return result;
        }

        // 对当前维度排序
        std::sort(indices.begin(), indices.end(), [&](int a, int b) { return X[a][dim] < X[b][dim]; });

        // 递归求解下一维度的保序回归
        std::vector<double> f_next_dim = dp_recursive(dim + 1, indices);

        // 当前维度上的一维保序回归
        std::vector<double> f_current_dim = isotonic_regression_1d(f_next_dim);

        // 恢复原顺序
        std::vector<double> restored_order(indices.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            restored_order[indices[i]] = f_current_dim[i];
        }

        return restored_order;
    };

    // 计算最终的保序回归结果
    std::vector<int> indices(n_samples);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> result = dp_recursive(0, indices);

    return result;
}

// 使用投影梯度法实现多元保序回归
std::vector<double> multivariate_isotonic_regression_pgd(const std::vector<std::vector<double>>& X, const std::vector<double>& y, int max_iter = 100, double learning_rate = 0.01) {
    int n_samples = X.size();
    int n_features = X[0].size();
    std::vector<double> f_hat(n_samples, std::accumulate(y.begin(), y.end(), 0.0) / n_samples);

    for (int iteration = 0; iteration < max_iter; ++iteration) {
        // 计算梯度
        std::vector<double> gradient(n_samples);
        for (int i = 0; i < n_samples; ++i) {
            gradient[i] = -2 * (y[i] - f_hat[i]);
        }

        // 沿负梯度方向更新估计值
        for (int i = 0; i < n_samples; ++i) {
            f_hat[i] -= learning_rate * gradient[i];
        }

        // 投影操作
        for (int feature = 0; feature < n_features; ++feature) {
            std::vector<int> order(n_samples);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return X[a][feature] < X[b][feature]; });

            std::vector<double> f_hat_sorted(n_samples);
            for (int i = 0; i < n_samples; ++i) {
                f_hat_sorted[i] = f_hat[order[i]];
            }

            f_hat_sorted = isotonic_regression_1d(f_hat_sorted);

            for (int i = 0; i < n_samples; ++i) {
                f_hat[order[i]] = f_hat_sorted[i];
            }
        }
    }

    return f_hat;
}

int main() {
    // 读取输入文件
    std::ifstream infile("input1.txt");
    std::vector<std::vector<double>> X;
    std::vector<double> y;
    double value;
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::vector<double> row;
            while (iss >> value) {
                row.push_back(value);
            }
            y.push_back(row.back());
            row.pop_back();
            X.push_back(row);
        }
        infile.close();
    } else {
        std::cerr << "无法打开输入文件" << std::endl;
        return 1;
    }

    // 记录动态规划法运行时间
    auto start_dp = std::chrono::high_resolution_clock::now();
    std::vector<double> f_hat_dp = multivariate_isotonic_regression_dp(X, y);
    auto end_dp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_dp = end_dp - start_dp;
    std::cout << "动态规划法程序运行时间: " << duration_dp.count() << " 秒" << std::endl;

    // 记录投影梯度法运行时间
    auto start_pgd = std::chrono::high_resolution_clock::now();
    std::vector<double> f_hat_pgd = multivariate_isotonic_regression_pgd(X, y);
    auto end_pgd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pgd = end_pgd - start_pgd;
    std::cout << "投影梯度法程序运行时间: " << duration_pgd.count() << " 秒" << std::endl;

    // 将结果写入输出文件
    std::ofstream outfile("output1.txt");
    if (outfile.is_open()) {
        outfile << "动态规划法保序回归结果: ";
        for (double val : f_hat_dp) {
            outfile << val << " ";
        }
        outfile << std::endl;

        outfile << "投影梯度法保序回归结果: ";
        for (double val : f_hat_pgd) {
            outfile << val << " ";
        }
        outfile.close();
    } else {
        std::cerr << "无法打开输出文件" << std::endl;
        return 1;
    }

    return 0;
}
