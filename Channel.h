#ifndef CHANNEL_H
#define CHANNEL_H

#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Eigen>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>
#include "parameters.h"
#include "random_collection.h"

class Channel
{
public:
    Channel(const SimulationParameters &params) : params_(params)
    {
        W_.resize(params_.K_, params_.Q_);
        setW_();

        xi_.resize(params_.Q_);
        h_.resize(params_.L_, params_.Q_);
        Cmat_.resize(params_.L_, params_.L_);
        U_.resize(params_.L_, params_.L_);
        lambda_.resize(params_.L_, params_.L_);
        H_.resize(params_.L_, params_.K_);

        // 乱数設定
        unitCNormalRand_.init(0.0, 1.0 / sqrt(2.0), params_.seed + 1);

        // 遅延プロファイル生成
        setChannelProfile();
    }

    // 周波数応答 H_ を生成する public メソッド
    void generateFrequencyResponse(double f_d)
    {
        f_d_ = f_d;
        // 伝送路のインパルス応答の生成
        seth_();
        // 伝送路の周波数応答の生成:式(19)
        for (int l = 0; l < params_.L_; l++)
        {
            H_.row(l) = (W_ * h_.row(l).transpose()).transpose();
        }
    }

    // H_ を外部から参照するための getter メソッド
    const Eigen::MatrixXcd &getH() const
    {
        return H_;
    }

private:
    const SimulationParameters &params_;
    double f_d_;

    Eigen::VectorXd xi_;   // 伝送路のインパルス応答の遅延プロファイル
    Eigen::MatrixXcd h_;   // インパルス応答
    Eigen::MatrixXd Cmat_; // 共分散行列
    Eigen::MatrixXcd U_;
    Eigen::MatrixXcd lambda_;
    Eigen::MatrixXcd H_; // 周波数応答
    Eigen::MatrixXcd W_; // DFT行列（Channelクラス内で閉じるため、Simulatorから移動）

    cnormal_distribution<> unitCNormalRand_;

    void setW_()
    {
        for (int q = 0; q < params_.Q_; ++q)
        {
            for (int k = 0; k < params_.K_ / 2; ++k)
            {
                // -26から-1番目のキャリヤ
                W_(k, q) = std::polar(1.0, -2.0 * M_PI * ((double)k - (double)params_.K_ / 2.0) * (double)q / (double)params_.NUMBER_OF_FFT);
                // 1から26番目のキャリヤ
                W_(k + params_.K_ / 2, q) = std::polar(1.0, -2.0 * M_PI * ((double)k + 1.0) * (double)q / (double)params_.NUMBER_OF_FFT);
            }
        }
    }

    void setChannelProfile()
    {
        for (int q = 0; q < params_.Q_; q++)
        {
            if (q == 0)
            {
                xi_(q) = 1.0;
            }
            else
            {
                xi_(q) = xi_(q - 1) * std::pow(10.0, -0.1);
            }
        }

        double tmp = xi_.sum();

        for (int q = 0; q < params_.Q_; q++)
        {
            xi_(q) = xi_(q) / tmp;
        }
    }

    /**
     * 共分散行列生成
     */
    void setCmat_()
    {
        for (auto l_1 = 0; l_1 < params_.L_; l_1++)
        {
            for (auto l_2 = 0; l_2 < params_.L_; l_2++)
            {
                Cmat_(l_1, l_2) = std::cyl_bessel_j(0, 2.0 * M_PI * (double)abs((l_1 - l_2) * params_.Ts_) * f_d_); // Jakesモデル
            }
        }
    }

    /**
     * インパルス応答行列生成
     */
    void seth_()
    {
        Eigen::VectorXcd h_q(params_.L_);
        Eigen::MatrixXcd A(params_.L_, params_.L_);
        Eigen::VectorXcd x(params_.L_);

        for (auto q = 0; q < params_.Q_; q++)
        {
            setCmat_();
            Cmat_ = xi_(q) * Cmat_;

            computeUnitaryAndDiagonal();
            // computeUnitaryAndDiagonal_test();

            A = U_ * lambda_.array().sqrt().matrix();

            for (auto l = 0; l < params_.L_; l++)
            {
                x(l) = unitCNormalRand_();
            }

            h_q = A * x;

            for (auto l = 0; l < params_.L_; l++)
            {
                h_(l, q) = h_q(l);
            }
        }
    }

    // 固有値分解によりユニタリ行列と対角行列を生成する関数
    void computeUnitaryAndDiagonal()
    {
        // 固有値分解を実行
        Eigen::EigenSolver<Eigen::MatrixXd> solver(Cmat_);

        // ユニタリ行列（固有ベクトルを列に持つ行列）
        U_ = solver.eigenvectors();

        // 対角行列（固有値を対角要素に持つ行列）
        lambda_ = solver.eigenvalues().asDiagonal();
    }
};
#endif