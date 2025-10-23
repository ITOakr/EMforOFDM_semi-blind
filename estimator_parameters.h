#ifndef ESTIMATOR_PARAMETERS_H
#define ESTIMATOR_PARAMETERS_H

// 推定器のパラメータを定義する構造体
struct EstimatorParameters
{
    const int Q_est = 2; // ★推定器が仮定するパス数
};

#endif // ESTIMATOR_PARAMETERS_H