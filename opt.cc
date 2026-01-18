// opt.cc
// Cross 场景 jerk 优化器 —— 完全修复版
#include <iostream>
#include <cmath>

bool computeOptimalJerksCross(
    double d_ego, double v_ego, double a_ego,
    double d_obs,  double v_obs,  double a_obs,
    double t_c,
    double t_hw,
    double w_agent,
    double& jerk_ego,
    double& jerk_obs
) {
    const double eps = 1e-9;
    if (t_c <= eps || w_agent < 0) return false;

    const double t = t_c;
    const double hw = t_hw;

    double t2 = t * t;
    double t3 = t2 * t;

    double K = t3 / 6.0;      // jerk → traveled 距离增益
    double Kv = 0.5 * t2;     // jerk → 速度增益

    // 常数项（不含 jerk）
    double trav_ego_const = v_ego * t + 0.5 * a_ego * t2;
    double trav_obs_const = v_obs * t + 0.5 * a_obs * t2;
    double v_obs_const = v_obs + a_obs * t;

    // ✅ 正确的安全约束（推导见注释）：
    // d_ego - trav_ego <= d_obs - trav_obs - v_obs(t) * t_hw
    // 整理为等式：A * j_ego + B * j_obs = C_rhs
    /* 

    traveled_ego(t_c) = v_ego * t_c + 0.5 * a_ego * t_c^2 + (1.0/6.0) * j_ego * t_c^3
    traveled_obs(t_c) = v_obs * t_c + 0.5 * a_obs * t_c^2 + (1.0/6.0) * j_obs * t_c^3
    v_obs(t_c) = v_obs + a_obs * t_c + 0.5 * j_obs * t_c^2

    d_ego - traveled_ego(t_c) <= d_obs - traveled_obs(t_c) - v_obs(t_c) * t_hw
    -(1.0/6.0) * j_ego * t_c^3 <= -(1.0/6.0) * j_obs * t_c^3 - 0.5 * a_obs * t_c^2 - v_obs * t_c + d_obs - d_ego + v_ego * t_c + 0.5 * a_ego * t_c^2 + (v_obs + a_obs * t_c + 0.5 * j_obs * t_c^2) * t_hw
    -(1.0/6.0) * t_c^3 * j_ego + ((1.0/6.0) * t_c^3 +0.5 * t_c^2 * t_hw) * j_obs <= - 0.5 * a_obs * t_c^2 - v_obs * t_c + d_obs - d_ego + v_ego * t_c + 0.5 * a_ego * t_c^2 + (v_obs + a_obs * t_c) * t_hw
    A = -(1.0/6.0) * t_c^3
    B = (1.0/6.0) * t_c^3 +0.5 * t_c^2 * t_hw
    C_rhs = - 0.5 * a_obs * t_c^2 - v_obs * t_c + d_obs - d_ego + v_ego * t_c + 0.5 * a_ego * t_c^2 + (v_obs + a_obs * t_c) * t_hw

    C_rhs = d_obs - d_ego - (v_obs * t_c + 0.5 * a_obs * t_c^2) + (v_ego * t_c + 0.5 * a_ego * t_c^2) + (v_obs + a_obs * t_c) * t_hw
    */

    double A = -K;                                      // j_ego 系数
    double B = K + Kv * hw;                             // j_obs 系数
    double C_rhs = (d_obs - d_ego)                      // 右侧常数项
                 - (trav_obs_const - trav_ego_const)
                 - v_obs_const * hw;                    // ⚠️ 关键：此处是 MINUS！

    if (w_agent < eps) {
        // w_agent ≈ 0: 障碍物 jerk 无成本 → 让它承担全部避让
        if (std::abs(B) < eps) return false;
        jerk_ego = 0.0;
        jerk_obs = C_rhs / B;
        return true;
    }

    // 通用情况：min j_ego^2 + w_agent * j_obs^2
    double denom = A * A + (B * B) / w_agent;
    if (std::abs(denom) < eps) return false;

    double lambda = 2.0 * C_rhs / denom;
    jerk_ego = lambda * A / 2.0;
    jerk_obs = lambda * B / (2.0 * w_agent);

    return true;
}

void testScenario(
    double d_ego, double v_ego, double a_ego,
    double d_obs,  double v_obs,  double a_obs,
    double t_c, double t_hw,
    double w_agent
) {
    double j_ego, j_obs;
    bool ok = computeOptimalJerksCross(d_ego, v_ego, a_ego,
                                       d_obs,  v_obs,  a_obs,
                                       t_c, t_hw, w_agent,
                                       j_ego, j_obs);
    std::cout << "w_agent = " << w_agent << "\n";
    if (ok) {
        std::cout << "  jerk_ego = " << j_ego << " m/s³\n";
        std::cout << "  jerk_obs  = " << j_obs << " m/s³\n";

        // 验证轨迹
        double t = t_c;
        double trav_ego = v_ego*t + 0.5*a_ego*t*t + (1.0/6.0)*j_ego*t*t*t;
        double trav_obs = v_obs*t + 0.5*a_obs*t*t + (1.0/6.0)*j_obs*t*t*t;
        double v_obs_t = v_obs + a_obs*t + 0.5*j_obs*t*t;

        double ego_remaining = d_ego - trav_ego;           // >0 表示未到
        double obs_safe_margin = d_obs - trav_obs - v_obs_t * t_hw; // >0 表示安全

        std::cout << "  Ego remaining to CP: " << ego_remaining << " m\n";
        std::cout << "  Obs safe margin:      " << obs_safe_margin << " m\n";

        // 安全条件：ego 剩余距离 ≤ obs 安全裕度
        bool satisfied = (ego_remaining <= obs_safe_margin + 1e-3);
        std::cout << "  Constraint satisfied? " << (satisfied ? "YES" : "NO") << "\n";
    } else {
        std::cout << "  Failed to solve.\n";
    }
    std::cout << "-------------------\n";
}

int main() {
    // 🚦 测试场景：真实冲突风险
    // 自车：20m 外，8 m/s → 到达时间 2.5s
    // 障碍物：18m 外，9 m/s → 到达时间 2.0s
    // 安全要求：时间差 ≥ 1.0s → 必须协调！
    double d_ego = 30.0;
    double v_ego = 8.0;
    double a_ego = 0.0;

    double d_obs = 30.0;
    double v_obs = 8.0;
    double a_obs = 0.0;

    double t_ego0 = d_ego / std::max(v_ego, 0.1);
    double t_obs0 = d_obs / std::max(v_obs, 0.1);
    double t_c = (t_ego0 + t_obs0) / 2.0; // ≈ 2.25 s

    double t_hw = 1.0;

    std::cout << "=== Cross Scenario Jerk Optimization (Fully Fixed) ===\n";
    std::cout << "Ego: d=" << d_ego << " m, v=" << v_ego << " m/s\n";
    std::cout << "Obs: d=" << d_obs << " m, v=" << v_obs << " m/s\n";
    std::cout << "t_c = " << t_c << " s, t_hw = " << t_hw << " s\n\n";

    for (double w : {0.0, 0.25, 0.5, 1.0}) {
        testScenario(d_ego, v_ego, a_ego,
                     d_obs,  v_obs,  a_obs,
                     t_c, t_hw, w);
    }

    return 0;
}