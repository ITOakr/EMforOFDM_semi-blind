#ifndef TRANSCEIVER_H
#define TRANSCEIVER_H

#include <stdexcept> // (ファイルの先頭に追加)
#include <string>    // (エラーメッセージ構築のため)
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Eigen>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>
#include "parameters.h"
#include "random_collection.h"
#include "estimator_parameters.h"

class Transceiver
{
public:
    Transceiver(const SimulationParameters &params, const Eigen::MatrixXcd &W) : params_(params), W_est_(W)
    {
        EstimatorParameters est_params_;
        // リサイズ
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
        h_l.resize(params_.Q_);
        grayNum_.resize(params_.NUMBER_OF_SYMBOLS);

        noiseVariance_ = 0.0;

        unitIntUniformRand_.init(0, params_.NUMBER_OF_SYMBOLS - 1, params_.seed);
        unitCNormalRand_.init(0.0, 1.0 / sqrt(2.0), params_.seed);

        // グレイ符号のテーブルを作成
        setGrayNum();
        // シンボル設計とDFT行列設定
        setSymbol();
        //setW_est_();
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
    void equalizeByPilotAndDemodulate()
    {
        equalizeChannelWithPilot();
        setRxDataByML();
    }

    // EMアルゴリズムによる等化と復調
    double equalizeAndDemodulate()
    {
        double avg_iter = equalizeChannelWithEM();
        setRxDataByML();
        return avg_iter;
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

    /**
     * チャネル推定のMSEを計算する
     * @return 1試行あたりの二乗誤差の合計
     */
    double getMSE_during_pilot()
    {
        double mse = 0.0;
        // データシンボル区間（パイロットを除く）のMSEを計算
        mse = (H_true_.row(0) - H_est_.row(0)).squaredNorm();
        return mse;
    }

    /**
     * 推定された雑音分散を取得する
     * @return 推定された雑音分散 noiseVariance_ の値
     */
    double getEstimatedNoiseVariance() const
    {
        return noiseVariance_;
    }

    // パイロットシンボルからhを推定し，Hを得る
    void est_H_by_initial_h(){
        set_initial_params_by_pilot();
        H_est_.row(0) = (W_est_ * h_l).transpose();
        // std::cout << "OK6" << std::endl;
    }

    // パイロットシンボルから直接Hを推定する
    void est_H_by_pilot(){
        H_est_.row(0) = (X_.row(0).asDiagonal()).inverse() * Y_.row(0).transpose();
    }


private:
    const SimulationParameters &params_;
    const Eigen::MatrixXcd &W_est_;
    EstimatorParameters est_params_;
    std::vector<int> grayNum_;

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
    double noiseVariance_;
    uniform_int_distribution<> unitIntUniformRand_;
    cnormal_distribution<> unitCNormalRand_;
    std::vector<int> activePathIndices_;

    // // DFT行列Wの生成:式(17)
    // void setW_est_()
    // {
    //     for (int q = 0; q < est_params_.Q_est; ++q)
    //     {
    //         for (int k = 0; k < params_.K_ / 2; ++k)
    //         {
    //             // -26から-1番目のキャリヤ
    //             W_est_(k, q) = std::polar(1.0, -2.0 * M_PI * ((double)k - (double)params_.K_ / 2.0) * (double)q / (double)params_.NUMBER_OF_FFT);
    //             // 1から26番目のキャリヤ
    //             W_est_(k + params_.K_ / 2, q) = std::polar(1.0, -2.0 * M_PI * ((double)k + 1.0) * (double)q / (double)params_.NUMBER_OF_FFT);
    //         }
    //     }
    //     // std::cout << "W_=" << W_ << std::endl;
    // }

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
    void set_initial_params_by_pilot()
    {
        X_l = X_.row(0).asDiagonal();
        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_.row(0).transpose();

        // 電力(norm)とインデックスのペアを作成
        std::vector<std::pair<double, int>> pathRank;
        for (int q = 0; q < params_.Q_; ++q) {
            pathRank.push_back({ std::norm(h_l(q)), q }); // 電力と元のインデックス
        }

        // 電力が大きい順にソート
        std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

        // for (const auto& p : pathRank) {
        //     // p.first が電力、p.second が元のインデックス
        //     std::cout << "Power: " << p.first << ", Index: " << p.second << std::endl;
        // }

        std::vector<double> aic_list(params_.Q_);
        std::vector<Eigen::VectorXcd> h_tilde_list(params_.Q_); // 各パスモデルにおける推定値
        std::vector<double> beta_list(params_.Q_);             // 各パスモデルにおける雑音精度

        for (int Q_tilde = 1; Q_tilde <= params_.Q_; ++Q_tilde) {
            Eigen::MatrixXcd W_tilde(params_.K_, Q_tilde);
            for (int i = 0; i < Q_tilde; ++i) {
                int original_idx = pathRank[i].second; // ソート済みの上位インデックスを取得
                W_tilde.col(i) = W_est_.col(original_idx); // 対応するDFT行列の列をコピー
            }

            // 各パスモデルにおける推定値を計算
            h_tilde_list[Q_tilde - 1] = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * W_tilde.adjoint() * X_l.adjoint() * Y_.row(0).transpose();
            // std::cout << "h_l >> " << h_tilde_list[Q_tilde - 1] << std::endl;
            double residual = (Y_.row(0).transpose() - X_l * W_tilde * h_tilde_list[Q_tilde - 1]).squaredNorm();
            beta_list[Q_tilde - 1] = (double)params_.K_ / residual;

            // AIC の計算
            double logL = params_.K_ * std::log(beta_list[Q_tilde - 1]) - params_.K_ * std::log(M_PI) - params_.K_;
            aic_list[Q_tilde - 1] = -2.0 * (logL - 2.0 * Q_tilde);
            // std::cout << "aic >> " << aic_list[Q_tilde - 1] << std::endl;
        }

        // std::cout << "OK3" << std::endl;

        // AIC が最小となるインデックスを特定
        int best_idx = std::distance(aic_list.begin(), std::min_element(aic_list.begin(), aic_list.end()));

        // std::cout << "OK4" << std::endl;
        // 確定した推定結果をクラスのメンバ変数へ格納
        int best_q_idx = std::distance(aic_list.begin(), std::max_element(aic_list.begin(), aic_list.end()));
        int best_Q_tilde = best_q_idx + 1;

        // 5. 選択されたパスの情報を整理 (0埋め処理)
        Eigen::VectorXcd h_full = Eigen::VectorXcd::Zero(params_.Q_);
        activePathIndices_.clear();
        
        for (int i = 0; i < best_Q_tilde; ++i) {
            int original_idx = pathRank[i].second;
            activePathIndices_.push_back(original_idx); // 採用したインデックスを保存
            h_full(original_idx) = h_tilde_list[best_q_idx](i); // 推定値を元の位置へ配置
        }
        this->h_l = h_full;
        this->noiseVariance_ = 1.0 / beta_list[best_idx];
        // std::cout << h_l << std::endl;
    }

    void equalizeChannelWithPilot()
    {
        set_initial_params_by_pilot();
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

    double equalizeChannelWithEM()
    {
        set_initial_params_by_pilot();

        H_est_.row(0) = (W_est_ * h_l).transpose();

        const int MAX_ITER = 100;
        const int MIN_ITER = 3;

        Eigen::VectorXi symbol_prev2(params_.K_);
        Eigen::VectorXi symbol_prev1(params_.K_);
        Eigen::VectorXi symbol_current(params_.K_);

        Eigen::VectorXd obj(params_.NUMBER_OF_SYMBOLS);

        int dataSymbolCount = params_.L_ - params_.NUMBER_OF_PILOT;
        Eigen::VectorXd iter_counts(dataSymbolCount);

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            symbol_prev2.setConstant(-1);
            symbol_prev1.setConstant(-1);

            int current_iter_count = 0;

            for (int iter = 0; iter < MAX_ITER; iter++)
            {
                current_iter_count = iter + 1;
                // Eステップ
                Estep(l);
                // Mステップ
                Mstep(l);

                for (int k = 0; k < params_.K_; k++)
                {
                    std::complex<double> H_current_k = (W_est_.row(k) *h_l)(0);
                    std::complex<double> R_current_k = Y_(l, k) / H_current_k;

                    for(int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++){
                        obj(i) = std::norm(R_current_k - symbol_(i));
                    }
                    Eigen::VectorXd::Index minColumn;
                    obj.minCoeff(&minColumn);
                    symbol_current(k) = minColumn;
                }

                if (iter >= MIN_ITER - 1)
                {
                    bool converged = (symbol_current == symbol_prev1) && (symbol_prev1 == symbol_prev2);
                    if (converged){
                        break;
                    }
                }
                symbol_prev2 = symbol_prev1;
                symbol_prev1 = symbol_current;
            }

            iter_counts(l - params_.NUMBER_OF_PILOT) = static_cast<double>(current_iter_count);

            H_est_.row(l) = (W_est_ * h_l).transpose();
            for (int k = 0; k < params_.K_; k++)
            {
                R_(l, k) = Y_(l, k) / (W_est_.row(k) * h_l)(0);
            }
        }
        if(iter_counts.size() == 0){
            return 0.0;
        }

        return iter_counts.mean();
    }

    void Estep(int l)
    {
        // std::cout << "h_l=" << h_l << std::endl;
        Eigen::VectorXcd H_current = W_est_ * h_l;
        double variance = noiseVariance_;
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
                //  posterior_prob.fill(1.0 / static_cast<double>(params_.NUMBER_OF_SYMBOLS));
                //  if (k == 0) { // デバッグ用にメッセージ表示
                //       std::cerr << "Warning: sumxP_shifted is zero or negative at l=" << l << ", k=" << k << ". Assigning equal probabilities." << std::endl;
                    // エラーメッセージを構築
                std::string error_msg = "FATAL ERROR: sumxP_shifted is zero or negative at l=" 
                          + std::to_string(l) + ", k=" + std::to_string(k);

                // どの場所でエラーが起きたかコンソールに表示
                std::cerr << error_msg << std::endl;
                
                // 例外を投げてプログラムを停止させます
                throw std::runtime_error(error_msg);   
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
        Eigen::MatrixXcd A = X_bar * W_est_;
        h_l = (W_est_.adjoint() * R_moment * W_est_).inverse() * (X_bar * W_est_).adjoint() * Y_.row(l).transpose();
        // std::cout << "h_l=" << h_l << std::endl;
        // std::cout << "A=" << A << std::endl;
        noiseVariance_ = ((Y_.row(l).transpose() - X_bar * W_est_ * h_l).squaredNorm() + ((W_est_ * h_l).adjoint() * (R_moment - X_bar.adjoint() * X_bar) * W_est_ * h_l).value().real()) / (double)params_.K_;
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