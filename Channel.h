#ifndef CHANNEL_H
#define CHANNEL_H

#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Eigen>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>
#include "parameters.h"
#include "random_collection.h"

class Channel
{
public:
    Channel(const SimulationParameters &params, const Eigen::MatrixXcd &W) : params_(params), W_(W)
    {
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
    void generateFrequencyResponse(double fd_Ts)
    {
        fd_Ts_ = fd_Ts;
        // 伝送路のインパルス応答の生成
        seth_();
        // 伝送路の周波数応答の生成:式(19)
        for (int l = 0; l < params_.L_; l++)
        {
            H_.row(l) = (W_ * h_.row(l).transpose()).transpose();
        }
        // std::cout << H_ << std::endl;
    }

    /**
     * @brief 1フレーム分のH_の平均電力を計算する
     * @return double 平均電力
     */
    double getAveragePower()
    {
        return H_.squaredNorm() / H_.size();
    }

    // H_ を外部から参照するための getter メソッド
    const Eigen::MatrixXcd &getH() const
    {
        return H_;
    }

    const Eigen::MatrixXcd &get_h() const
    {
        return h_;
    }

private:
    const SimulationParameters &params_;
    const Eigen::MatrixXcd &W_;
    double fd_Ts_;

    Eigen::VectorXd xi_;   // 伝送路のインパルス応答の遅延プロファイル
    Eigen::MatrixXcd h_;   // インパルス応答
    Eigen::MatrixXd Cmat_; // 共分散行列
    Eigen::MatrixXcd U_;
    Eigen::MatrixXcd lambda_;
    Eigen::MatrixXcd H_; // 周波数応答

    cnormal_distribution<> unitCNormalRand_;

    void setChannelProfile()
    {
        double totalPower = 0.0;
        for (int q = 0; q < params_.Q_; q++)
        {
            // if (q == 0)
            // {
            //     xi_(q) = 1.0;
            // }
            // else
            // {
            //     xi_(q) = xi_(q - 1) * std::pow(10.0, -0.1);
            // }
            xi_(q) = (double)params_.pathMask[q] * std::pow(10.0, -0.1 * q);
            totalPower += xi_(q);
        }

        if (totalPower > 0){
            for (int q = 0; q < params_.Q_; q++)
                {
                    xi_(q) = xi_(q) / totalPower;
                }
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
                Cmat_(l_1, l_2) = std::cyl_bessel_j(0, 2.0 * M_PI * std::abs(l_1 - l_2) * fd_Ts_); // Jakesモデル
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
                // std::cout << l << "=" << h_ << std::endl;
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