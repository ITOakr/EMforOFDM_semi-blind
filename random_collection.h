/* 
 * File:   random_collection.h
 * Author: fujii
 *
 * Created on 2021/2/12, 16:12
*/

#ifndef RANDOM_COLLECTION_H
#define RANDOM_COLLECTION_H

#include <random>       // 標準ライブラリの乱数生成用ヘッダー
#include <functional>   // std::functionとstd::bindの使用
#include <complex>      // 複素数のサポート
#include <iostream>     // デバッグや出力用

// 基本クラス: distribution
// 各分布クラスの基底クラスとして使用され、乱数生成機能を提供します。
template <typename T>
class distribution {
public:
    // ランダムな値を生成するための関数オペレーター
    T operator()() const {
        return dist_();
    }
protected:
    // 乱数生成関数オブジェクトを保持
    std::function<T()> dist_;
};

// 整数の一様分布
// 指定した範囲 [min, max] からランダムな整数を生成します。
template <typename T = int, typename E = std::mt19937>
class uniform_int_distribution : public distribution<T> {
public:
    // 分布を初期化するメソッド
    // min: 最小値, max: 最大値, seed: シード値
    void init(T min, T max, unsigned long int seed) {
        distribution<T>::dist_ = std::bind(std::uniform_int_distribution<T>(min, max), E(seed));
    }
};

// 実数の一様分布
// 指定した範囲 [min, max] からランダムな実数を生成します。
template <typename T = double, typename E = std::mt19937>
class uniform_real_distribution : public distribution<T> {
public:
    // 分布を初期化するメソッド
    // min: 最小値, max: 最大値, seed: シード値
    void init(T min, T max, unsigned long int seed) {
        distribution<T>::dist_ = std::bind(std::uniform_real_distribution<T>(min, max), E(seed));
    }
};

// 正規分布（ガウス分布）
// 指定した平均と標準偏差に基づいてランダムな実数を生成します。
template <typename T = double, typename E = std::mt19937>
class normal_distribution : public distribution<T> {
public:
    // 分布を初期化するメソッド
    // mean: 平均値, sd: 標準偏差, seed: シード値
    void init(T mean, T sd, unsigned long int seed) {
        distribution<T>::dist_ = std::bind(std::normal_distribution<T>(mean, sd), E(seed));
    }
};

// 複素数正規分布
// 各次元で独立した正規分布に基づいて複素数を生成します。
template <typename T = double, typename E = std::mt19937>
class cnormal_distribution : public distribution<T> {
public:
    // 分布を初期化するメソッド
    // mean: 平均値, sd: 標準偏差, seed: シード値
    void init(T mean, T sd, unsigned long int seed) {
        distribution<T>::dist_ = std::bind(std::normal_distribution<T>(mean, sd), E(seed));
    }

    // ランダムな複素数を生成する関数オペレーター
    std::complex<T> operator()() const {
        return std::complex<T>(distribution<T>::dist_(), distribution<T>::dist_());
    }
};

#endif /* RANDOM_COLLECTION_H */
