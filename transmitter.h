#ifndef TRANSMITTER_H
#define TRANSMITTER_H

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "parameters.h"
#include "random_collection.h"

class Transmitter
{
public:
    Transmitter(const SimulationParameters &params) : params_(params)
    {
        txData_.resize(params_.L_, params_.K_);
        X_.resize(params_.L_, params_.K_);
        symbol_.resize(params_.NUMBER_OF_SYMBOLS);
        grayNum_.resize(params_.NUMBER_OF_SYMBOLS);

        unitIntUniformRand_.init(0, params_.NUMBER_OF_SYMBOLS - 1, params_.seed);

        // グレイ符号のテーブルを作成
        setGrayNum();
        // シンボル設計とDFT行列設定
        setSymbol();
    }

    /**
     * 送信信号生成
     */
    void setX_()
    {
        for (int l = 0; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                if (l < params_.NUMBER_OF_PILOT)
                {
                    txData_(l, k) = 0;
                    X_(l, k) = params_.PILOT_SYMBOL_;
                }
                else
                {
                    txData_(l, k) = unitIntUniformRand_();
                    X_(l, k) = symbol_(txData_(l, k));
                }
            }
        }
    }

    // ゲッター群
    const Eigen::MatrixXcd &getX() const { return X_; }
    const Eigen::MatrixXi &getTxData() const { return txData_; }
    const Eigen::VectorXcd &getSymbol() const { return symbol_; }
    const std::vector<int> &getGrayNum() const { return grayNum_; }

    /**
     * [Mode 23用] 指定したサブキャリアの送信信号(X)の時間変動を出力
     */
    void exportTxSymbolTrace(int target_k, const std::string& filename)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        ofs << "l,Tx_I,Tx_Q,Tx_Abs,Tx_Phase" << std::endl;

        for (int l = 0; l < params_.L_; l++)
        {
            std::complex<double> val = X_(l, target_k);

            ofs << l << ","
                << val.real() << ","
                << val.imag() << ","
                << std::abs(val) << ","
                << std::arg(val)
                << std::endl;
        }
        ofs.close();
        std::cout << "Exported Tx Trace (k=" << target_k << ") to " << filename << std::endl;
    }

    /**
     * [Mode 24用] 指定したサブキャリアの「周波数応答の影響を受けた送信信号(HX)」の時間変動を出力
     */
    void exportFadedSymbolTrace(int target_k, const std::string& filename, const Eigen::MatrixXcd& H_current)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        ofs << "l,Faded_I,Faded_Q,Faded_Abs,Faded_Phase" << std::endl;

        for (int l = 0; l < params_.L_; l++)
        {
            std::complex<double> x_val = X_(l, target_k);
            std::complex<double> h_val = H_current(l, target_k);
            std::complex<double> val = h_val * x_val;

            ofs << l << ","
                << val.real() << ","
                << val.imag() << ","
                << std::abs(val) << ","
                << std::arg(val)
                << std::endl;
        }
        ofs.close();
        std::cout << "Exported Faded Trace (k=" << target_k << ") to " << filename << std::endl;
    }

    /**
     * [Mode 25用] 指定したサブキャリアのチャネル応答の絶対値(|H|)の時間変動を出力
     */
    void exportChannelMagnitudeTrace(int target_k, const std::string& filename, const Eigen::MatrixXcd& H_current)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        ofs << "l,H_Abs" << std::endl;

        for (int l = 0; l < params_.L_; l++)
        {
            double h_mag = std::abs(H_current(l, target_k));
            ofs << l << "," << h_mag << std::endl;
        }
        ofs.close();
        std::cout << "Exported Channel Magnitude Trace (k=" << target_k << ") to " << filename << std::endl;
    }

private:
    const SimulationParameters &params_;
    std::vector<int> grayNum_;
    Eigen::MatrixXi txData_;
    Eigen::MatrixXcd X_;
    Eigen::VectorXcd symbol_;
    uniform_int_distribution<> unitIntUniformRand_;

    /**
     * シンボル生成
     */
    void setSymbol() {
        int M = params_.NUMBER_OF_SYMBOLS;               // シンボル数 (M = 2^NUMBER_OF_BIT)
        int sqrtM = sqrt(M);                    // 実部/虚部のレベル数 (例: 16QAMならsqrtM=4)
        double P = 1.0 / (2.0 * (M - 1) / 3.0);

        // シンボル設計
        int i = 0;
        for (int v1 = 0; v1 < sqrtM; v1++) {
            for (int v2 = 0; v2 < sqrtM; v2++) {
                symbol_(i).real((2 * v1 - (sqrtM - 1)) * sqrt(P));  // 実部
                if (v1 % 2 == 0) {  // v1が偶数のときは通常の配置
                    symbol_(i).imag((2 * v2 - (sqrtM - 1)) * sqrt(P));  // 虚部
                } else {  // v1が奇数のとき、虚部の値を逆順にする
                    symbol_(i).imag(((sqrtM - 1) - 2 * v2) * sqrt(P));
                }
                i++;
            }
        }
    }

    /**
     * グレイ符号の生成
     */
    int grayCode(int num) {
        return num ^ (num >> 1);
    }

    /**
     * グレイ符号のテーブルを作成
     */
    void setGrayNum() {
        for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++) {
            grayNum_[i] = grayCode(i);
        }
    }
};

#endif // TRANSMITTER_H
