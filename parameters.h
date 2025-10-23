#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <cmath>

// シミュレーションのパラメータをまとめた構造体
struct SimulationParameters
{
    const int K_ = 4;                          // サブキャリア数 K52
    const int L_ = 10;                           // 1フレームのシンボル数 L
    const int Q_ = 2;                           // 伝送路のインパルス応答のパス数 Q
    const double T_ = 3.2 * std::pow(10, -6);   // 有効シンボル長 T
    const double Tgi_ = 0.8 * std::pow(10, -6); // ガードインターバル長 Tgi
    const double Ts_ = T_ + Tgi_;               // シンボル全体の長さ
    const int NUMBER_OF_FFT = 64;               // FFTポイント数(IEEE802.11a)
    const int NUMBER_OF_PILOT = 1;              // パイロットシンボル個数
    int NUMBER_OF_BIT;                          // 変調方式のビット数 (入力で設定)
    int NUMBER_OF_SYMBOLS;                      // シンボル数 (2^NUMBER_OF_BIT)
    int seed = 100;
};

#endif /* PARAMETERS_H */