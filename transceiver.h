#ifndef TRANSCEIVER_H
#define TRANSCEIVER_H

#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Eigen>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>
#include "parameters.h"
#include "random_collection.h"
#include "estimator_parameters.h"

class Transceiver
{
public:
    Transceiver(const SimulationParameters &params) : params_(params)
    {
        EstimatorParameters est_params_;
        // リサイズ
        W_est_.resize(params_.K_, est_params_.Q_est);
        H_true_.resize(params_.L_, params_.K_);
        H_est_.resize(params_.L_, params_.K_);
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
        h_l.resize(est_params_.Q_est);
        grayNum_.resize(params_.NUMBER_OF_SYMBOLS);

        unitIntUniformRand_.init(0, params_.NUMBER_OF_SYMBOLS - 1, params_.seed);
        unitCNormalRand_.init(0.0, 1.0 / sqrt(2.0), params_.seed);

        // グレイ符号のテーブルを作成
        setGrayNum();
        // シンボル設計とDFT行列設定
        setSymbol();
        setW_est_();
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
        H_true_ = H;
        for (int l = 0; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                Y_(l, k) = H(l, k) * X_(l, k) + noiseSD * unitCNormalRand_();
            }
        }
        // std::cout << "H_true_=" << H_true_ << std::endl;
    }

    // pilot信号による等化と復調
    void equalizeByPilotAndDemodulate(double noiseSD)
    {
        equalizeChannelWithPilot(noiseSD);
        setRxDataByML();
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
                count += hammingDistance(grayNum_[txData_(l, k)], grayNum_[rxData_(l, k)]);
            }
        }
        return count;
    }

    /**
     * チャネル推定のMSEを計算する
     * @return 1試行あたりの二乗誤差の合計
     */
    double getMSE()
    {
        double mse = 0.0;
        // データシンボル区間（パイロットを除く）のMSEを計算
        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            mse += (H_true_.row(l) - H_est_.row(l)).squaredNorm();
        }
        return mse;
    }

private:
    const SimulationParameters &params_;
    EstimatorParameters est_params_;
    std::vector<int> grayNum_;

    Eigen::MatrixXcd W_est_;
    Eigen::MatrixXi txData_;
    Eigen::MatrixXi rxData_;
    Eigen::VectorXcd symbol_;
    Eigen::MatrixXcd Y_;
    Eigen::MatrixXcd X_;
    Eigen::MatrixXcd X_l;
    Eigen::MatrixXcd R_;
    Eigen::MatrixXcd H_est_;
    Eigen::MatrixXcd H_true_;
    Eigen::VectorXd xPro;
    Eigen::MatrixXcd X_bar;
    Eigen::MatrixXd R_moment;
    Eigen::VectorXcd h_l;
    uniform_int_distribution<> unitIntUniformRand_;
    cnormal_distribution<> unitCNormalRand_;

    // DFT行列Wの生成:式(17)
    void setW_est_()
    {
        for (int q = 0; q < est_params_.Q_est; ++q)
        {
            for (int k = 0; k < params_.K_ / 2; ++k)
            {
                // -26から-1番目のキャリヤ
                W_est_(k, q) = std::polar(1.0, -2.0 * M_PI * ((double)k - (double)params_.K_ / 2.0) * (double)q / (double)params_.NUMBER_OF_FFT);
                // 1から26番目のキャリヤ
                W_est_(k + params_.K_ / 2, q) = std::polar(1.0, -2.0 * M_PI * ((double)k + 1.0) * (double)q / (double)params_.NUMBER_OF_FFT);
            }
        }
        // std::cout << "W_=" << W_ << std::endl;
    }

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
     * グレイ符号の生成 (追加)
     */
    int grayCode(int num) {
        return num ^ (num >> 1);
    }

    /**
     * グレイ符号のテーブルを作成 (追加)
     */
    void setGrayNum() {
        for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++) {
            grayNum_[i] = grayCode(i);
        }
    }

        // パイロットシンボルからｈの初期値を得る
    void seth_l_byPilot()
    {
        X_l = X_.row(0).asDiagonal();
        // // ▼▼▼ 正則化を追加 ▼▼▼
        // Eigen::MatrixXcd A = W_est_.adjoint() * X_l.adjoint() * X_l * W_est_;
        // double delta = 1e-10; // 正則化のための微小な値
        // Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(A.rows(), A.cols());
        // h_l = (A + delta * I).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_.row(0).transpose();
        // // ▲▲▲ 正則化を追加 ▲▲▲

        // std::cout << "Y_.row(0).transpose()" << Y_.row(0).transpose() << std::endl;
        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_.row(0).transpose();
        // std::cout << "h_l_by_p=" << h_l << std::endl;
    }

    void equalizeChannelWithPilot(double noiseSD)
    {
        seth_l_byPilot();
        H_est_.row(0) = (W_est_ * h_l).transpose();

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            H_est_.row(l) = (W_est_ * h_l).transpose();
            for (int k = 0; k < params_.K_; k++)
            {
                R_(l, k) = Y_(l, k) / (W_est_.row(k) * h_l)(0);
            }
        }
    }

    void equalizeChannelWithEM(double noiseSD)
    {
        seth_l_byPilot();

        H_est_.row(0) = (W_est_ * h_l).transpose();

        // std::cout << H_est_ << std::endl;

        const int MAX_ITER = 5;
        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int iter = 0; iter < MAX_ITER; iter++)
            {
                // Eステップ
                Estep(l, noiseSD);
                // Mステップ
                Mstep(l);
            }

            H_est_.row(l) = (W_est_ * h_l).transpose();
            // std::cout << "L=" << l << ":" << h_l << std::endl;
            for (int k = 0; k < params_.K_; k++)
            {
                R_(l, k) = Y_(l, k) / (W_est_.row(k) * h_l)(0);
            }
        }
        // std::cout << "H_est_=" << H_est_ << std::endl;
    }

    void Estep(int l, double noiseSD)
    {
        // std::cout << "h_l=" << h_l << std::endl;
        Eigen::VectorXcd H_current = W_est_ * h_l;
        double variance = noiseSD * noiseSD;
        // std::cout << "H_current=" << H_current << std::endl;
        // varianceが非常に小さい場合のアンダーフロー対策（ゼロ除算を避ける）
        if (variance < std::numeric_limits<double>::epsilon()) {
            variance = std::numeric_limits<double>::epsilon();
        }
        for (int k = 0; k < params_.K_; k++)
        {
            Eigen::VectorXd log_likelihoods(params_.NUMBER_OF_SYMBOLS);

            for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++)
            {
                std::complex<double> s = symbol_(i);
                double norm = std::norm(Y_(l, k) - H_current(k) * s);
                log_likelihoods(i) = -norm / variance;
            }
            double log_max = log_likelihoods.maxCoeff();
            double sumxP_shifted = 0.0;
            Eigen::VectorXd posterior_prob(params_.NUMBER_OF_SYMBOLS);
            for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++)
            {
                // exp(対数尤度 - 最大対数尤度) を計算
                // exp の引数は最大でも0なので、アンダーフローしにくい
                posterior_prob(i) = std::exp(log_likelihoods(i) - log_max);
                sumxP_shifted += posterior_prob(i);
            }

            if (sumxP_shifted <= 0.0) {
                 // 稀にすべてのexpの結果が0になった場合 (非常に考えにくいが念のため)
                 // 等確率を割り当てる
                 posterior_prob.fill(1.0 / static_cast<double>(params_.NUMBER_OF_SYMBOLS));
                 if (k == 0) { // デバッグ用にメッセージ表示
                      std::cerr << "Warning: sumxP_shifted is zero or negative at l=" << l << ", k=" << k << ". Assigning equal probabilities." << std::endl;
                 }
            } else {
                 // 通常の正規化
                 posterior_prob /= sumxP_shifted;
            }

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
        // // ▼▼▼ 正則化を追加 ▼▼▼
        Eigen::MatrixXcd A = W_est_.adjoint() * R_moment * W_est_;
        // double delta = 1e-10; // 正則化のための微小な値
        // Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(A.rows(), A.cols());
        // h_l = (A + delta * I).inverse() * W_est_.adjoint() * X_bar.adjoint() * Y_.row(l).transpose();
        // // ▲▲▲ 正則化を追加 ▲▲▲
        h_l = (W_est_.adjoint() * R_moment * W_est_).inverse() * W_est_.adjoint() * X_bar.adjoint() * Y_.row(l).transpose();
        // std::cout << "h_l=" << h_l << std::endl;
        // std::cout << "A=" << A << std::endl;
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