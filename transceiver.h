#ifndef TRANSCEIVER_H
#define TRANSCEIVER_H

#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Eigen>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>
#include "parameters.h"
#include "random_collection.h"

class Transceiver
{
public:
    Transceiver(const SimulationParameters &params) : params_(params)
    {
        // リサイズ
        W_.resize(params_.K_, params_.Q_);
        H_est_.resize(params_.K_);
        txData_.resize(params_.L_, params_.K_);
        rxData_.resize(params_.L_, params_.K_);
        Y_.resize(params_.L_, params_.K_);
        X_.resize(params_.L_, params_.K_);
        X_l.resize(params_.K_, params_.K_);
        R_.resize(params_.L_, params_.K_);
        symbol_.resize(params_.NUMBER_OF_SYMBOLS);
        xPro.resize(params_.NUMBER_OF_SYMBOLS);
        X_bar.setZero(params_.K_, params_.K_);
        R_moment.setZero(params_.K_, params_.K_);
        h_l.resize(params_.Q_);

        unitIntUniformRand_.init(0, params_.NUMBER_OF_SYMBOLS - 1, params_.seed);
        unitCNormalRand_.init(0.0, 1.0 / sqrt(2.0), params_.seed);

        // シンボル設計とDFT行列設定
        setSymbol();
        setW_();
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
                txData_(l, k) = unitIntUniformRand_();
                X_(l, k) = symbol_(txData_(l, k));
            }
        }
    }

    /**
     * 受信信号生成
     */
    void setY_(const Eigen::MatrixXcd& H, double noiseSD)
    {
        for (int l = 0; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                Y_(l, k) = H(l, k) * X_(l, k) + noiseSD * unitCNormalRand_();
            }
        }
    }

    // EMアルゴリズムによる等化と復調
    void equalizeAndDemodulate(double noiseSD)
    {
        equalizeChannelWithEM(noiseSD);
        setRxDataByML();
    }

    /**
     * ビット誤り数のカウント
     * @return 全ての誤りビット数
     */
    int getBitErrorCount()
    {
        int count = 0;
        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                count += hammingDistance(txData_(l, k), rxData_(l, k));
            }
        }
        return count;
    }

private:
    const SimulationParameters &params_;

    Eigen::MatrixXcd W_;
    Eigen::MatrixXi txData_;
    Eigen::MatrixXi rxData_;
    Eigen::VectorXcd symbol_;
    Eigen::MatrixXcd Y_;
    Eigen::MatrixXcd X_;
    Eigen::MatrixXcd X_l;
    Eigen::MatrixXcd R_;
    Eigen::VectorXcd H_est_;
    Eigen::VectorXd xPro;
    Eigen::MatrixXcd X_bar;
    Eigen::MatrixXd R_moment;
    Eigen::VectorXcd h_l;
    uniform_int_distribution<> unitIntUniformRand_;
    cnormal_distribution<> unitCNormalRand_;

    // DFT行列Wの生成:式(17)
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
        // std::cout << "W_=" << W_ << std::endl;
    }

    /**
     * シンボル生成
     */
    void setSymbol()
    {
        symbol_(0) = -1.0;
        symbol_(1) = 1.0;
    }

        // パイロットシンボルからｈの初期値を得る
    void seth_l_byPilot()
    {
        X_l = X_.row(0).asDiagonal();
        // std::cout << "Y_.row(0).transpose()" << Y_.row(0).transpose() << std::endl;
        h_l = (W_.adjoint() * X_l.adjoint() * X_l * W_).inverse() * W_.adjoint() * X_l.adjoint() * Y_.row(0).transpose();
        // std::cout << "h_l_by_p=" << h_l << std::endl;
    }

    void equalizeChannelWithEM(double noiseSD)
    {
        seth_l_byPilot();

        const int MAX_ITER = 10;
        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int iter = 0; iter < MAX_ITER; iter++)
            {
                // Eステップ
                Estep(l, noiseSD);
                // Mステップ
                Mstep(l);
            }
            // std::cout << "L=" << l << ":" << h_l << std::endl;
            for (int k = 0; k < params_.K_; k++)
            {
                R_(l, k) = Y_(l, k) / (W_.row(k) * h_l)(0);
            }
        }
        // std::cout << "R_=" << R_ << std::endl;
    }

    void Estep(int l, double noiseSD)
    {
        Eigen::VectorXcd H_current = W_ * h_l;

        for (int k = 0; k < params_.K_; k++)
        {
            double sumxP = 0.0;
            Eigen::VectorXd posterior_prob(params_.NUMBER_OF_SYMBOLS);
            for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++)
            {
                std::complex<double> s = symbol_(i);
                double norm = std::norm(Y_(l, k) - H_current(k) * s);
                double variance = noiseSD * noiseSD;
                posterior_prob(i) = std::exp(-norm / variance);
                sumxP += posterior_prob(i);
            }

            posterior_prob /= sumxP;

            std::complex<double> expected_X = 0.0;
            double expected_X_norm_sq = 0.0;

            for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++)
            {
                expected_X += posterior_prob(i) * symbol_(i);
                expected_X_norm_sq += posterior_prob(i) * std::norm(symbol_(i));
            }
            X_bar(k, k) = expected_X;
            R_moment(k, k) = expected_X_norm_sq;
        }
    }

    void Mstep(int l)
    {
        h_l = (W_.adjoint() * R_moment * W_).inverse() * W_.adjoint() * X_bar.adjoint() * Y_.row(l).transpose();
    }

    /**
     * 最尤復調
     */
    void setRxDataByML()
    {
        Eigen::VectorXd obj(params_.NUMBER_OF_SYMBOLS); // 最小化の目的関数

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++)
                {
                    // 最尤復調の周波数応答は推定値を使う？
                    obj(i) = std::norm((R_(l, k) - symbol_(i)));
                }
                Eigen::VectorXd::Index minColumn; // ノルムが最小な index（つまり受信データ）
                obj.minCoeff(&minColumn);
                rxData_(l, k) = minColumn;
            }
        }
    }

    /**
     * ハミング距離計算
     * @param 整数1，整数2
     * @return ハミング距離
     */
    int hammingDistance(int num1, int num2)
    {
        int ham = 0;
        int xorResult;
        int bitMask = 1;

        xorResult = num1 ^ num2;

        for (int i = 0; i < params_.NUMBER_OF_BIT; i++)
        {
            ham += (xorResult & bitMask) >> i;
            bitMask <<= 1;
        }

        return ham;
    }
};

#endif