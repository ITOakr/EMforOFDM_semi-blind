#ifndef RECEIVER_H
#define RECEIVER_H

#include <stdexcept>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <Eigen/Dense>
#include "parameters.h"
#include "random_collection.h"
#include "estimator_parameters.h"

class Receiver
{
public:
    Receiver(const SimulationParameters &params, const Eigen::MatrixXcd &W) : params_(params), W_est_(W)
    {
        // リサイズ
        H_true_.resize(params_.L_, params_.K_);
        H_est_.resize(params_.L_, params_.K_);
        rxData_.resize(params_.L_, params_.K_);
        Y_.resize(params_.L_, params_.K_);
        X_l.resize(params_.K_, params_.K_);
        R_.resize(params_.L_, params_.K_);
        symbol_.resize(params_.NUMBER_OF_SYMBOLS);
        xPro.resize(params_.NUMBER_OF_SYMBOLS);
        X_bar.setZero(params_.K_, params_.K_);
        R_moment.setZero(params_.K_, params_.K_);
        h_l.resize(params_.Q_);
        grayNum_.resize(params_.NUMBER_OF_SYMBOLS);

        noiseVariance_ = 0.0;

        unitCNormalRand_.init(0.0, 1.0 / sqrt(2.0), params_.seed);

        // グレイ符号のテーブルを作成
        setGrayNum();
        // シンボル設計とDFT行列設定
        setSymbol();
    }

    /**
     * 受信信号生成
     */
    void setY_(const Eigen::MatrixXcd& H, const Eigen::MatrixXcd& X, double noiseSD)
    {
        H_true_ = H;
        for (int l = 0; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                Y_(l, k) = H(l, k) * X(l, k) + noiseSD * unitCNormalRand_();
            }
        }
    }

    // pilot信号による等化と復調
    void equalizeByPilotAndDemodulate(const Eigen::MatrixXcd& X)
    {
        equalizeChannelWithPilot(X);
        setRxDataByML();
    }

    // EMアルゴリズムによる等化と復調
    double equalizeAndDemodulate(const Eigen::MatrixXcd& X)
    {
        double avg_iter = equalizeChannelWithEM(X);
        setRxDataByML();
        return avg_iter;
    }

    /**
     * ビット誤り数のカウント
     * @return 全ての誤りビット数
     */
    int getBitErrorCount(const Eigen::MatrixXi& txData)
    {
        int count = 0;
        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                if (params_.enableDataPilots && (k == 5 || k == 19 || k == 32 || k == 46))
                {
                    continue; // Skip pilot carriers
                }
                count += hammingDistance(grayNum_[txData(l, k)], grayNum_[rxData_(l, k)]);
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
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            mse += (H_true_.row(l) - H_est_.row(l)).squaredNorm();
        }
        mse /= static_cast<double>(params_.NUMBER_OF_PILOT);
        return mse;
    }

    /**
     * フレーム先頭 (l=0) のチャネル推定MSEを計算する
     * @return l=0における二乗誤差の合計
     */
    double getMSE_at_l0()
    {
        return (H_true_.row(0) - H_est_.row(0)).squaredNorm();
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
    void est_H_by_initial_h(const Eigen::MatrixXcd& X){
        set_initial_params_by_pilot(X);
        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    // パイロットシンボルからRaghavendra GAICを用いてhを推定し，Hを得る
    void est_H_by_initial_h_RaghavendraGAIC(const Eigen::MatrixXcd& X){
        set_initial_params_by_pilot_power_sort_RaghavendraGAIC(X);
        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    void est_H_by_RaghavendraAIC(const Eigen::MatrixXcd& X){
        Eigen::RowVectorXcd H_init = set_initial_params_by_RaghavendraAIC_Update(X);
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    void est_H_by_RaghavendraAIC2(const Eigen::MatrixXcd& X){
        Eigen::RowVectorXcd H_init = set_initial_params_by_RaghavendraAIC_Update2(X);
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    // 真のモデルと既知の雑音分散を使って、パイロットからML推定を行う
    void est_H_by_known_model_and_noise(const Eigen::MatrixXcd& X, double knownNoiseVariance)
    {
        activePathIndices_.clear();
        for (int q = 0; q < params_.Q_; ++q)
        {
            if (params_.pathMask[q] == 1)
            {
                activePathIndices_.push_back(q);
            }
        }

        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();

        Eigen::MatrixXcd W_active(params_.K_, static_cast<int>(activePathIndices_.size()));
        for (int i = 0; i < static_cast<int>(activePathIndices_.size()); ++i)
        {
            W_active.col(i) = W_est_.col(activePathIndices_[i]);
        }

        Eigen::VectorXcd h_active = (W_active.adjoint() * X_l.adjoint() * X_l * W_active).inverse() * W_active.adjoint() * X_l.adjoint() * Y_avg.transpose();

        h_l = Eigen::VectorXcd::Zero(params_.Q_);
        for (int i = 0; i < static_cast<int>(activePathIndices_.size()); ++i)
        {
            h_l(activePathIndices_[i]) = h_active(i);
        }

        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }

        noiseVariance_ = knownNoiseVariance;
    }

    // 16パスあると仮定してパイロットシンボルからインパルス応答をML推定する
    void est_H_by_16paths(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();

        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_avg.transpose();

        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    // パイロットシンボルから直接Hを推定する
    void est_H_by_pilot(const Eigen::MatrixXcd& X){
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd H_init = ((X_avg.asDiagonal()).inverse() * Y_avg.transpose()).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    /**
     * 指定されたマスクを用いてインパルス応答を推定し、周波数応答を保存する
     */
    void est_H_by_given_mask(const Eigen::MatrixXcd& X, const std::vector<int>& mask)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        X_l = X_avg.asDiagonal();

        std::vector<int> active_indices;
        for (int q = 0; q < mask.size(); ++q) {
            if (mask[q] == 1) active_indices.push_back(q);
        }

        int Q_tilde = active_indices.size();
        if (Q_tilde == 0) {
            h_l = Eigen::VectorXcd::Zero(params_.Q_);
            for (int l = 0; l < params_.NUMBER_OF_PILOT; l++) {
                H_est_.row(l).setZero();
            }
            return;
        }

        Eigen::MatrixXcd W_tilde(params_.K_, Q_tilde);
        for (int i = 0; i < Q_tilde; ++i) {
            W_tilde.col(i) = W_est_.col(active_indices[i]);
        }

        // LS推定
        Eigen::VectorXcd h_active = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * W_tilde.adjoint() * X_l.adjoint() * Y_avg.transpose();

        // フルベクトルへ展開
        this->h_l = Eigen::VectorXcd::Zero(params_.Q_);
        for (int i = 0; i < Q_tilde; ++i) {
            this->h_l(active_indices[i]) = h_active(i);
        }

        // 周波数応答へ変換
        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }
    }

    /**
     * AICで選択されたパスモデル（マスク）を取得する
     */
    std::vector<int> getEstimatedPathMask()
    {
        std::vector<int> mask(params_.Q_, 0);
        for (int q = 0; q < params_.Q_; ++q)
        {
            if (std::norm(h_l(q)) > 1e-20) {
                mask[q] = 1;
            }
        }
        return mask;
    }

    /**
     * ステップ 1: AIC 8パス総当たりによる最良マスクの探索
     */
    std::vector<int> findBestMaskByExhaustiveAIC_8paths(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        X_l = X_avg.asDiagonal();

        double min_aic = 1e18; // 十分大きな値
        int best_mask_int = 0;

        for (int i = 1; i < 256; ++i) {
            std::vector<int> active_indices;
            for (int bit = 0; bit < 8; ++bit) {
                if ((i >> bit) & 1) {
                    active_indices.push_back(bit);
                }
            }
            int Q_active = static_cast<int>(active_indices.size());

            Eigen::MatrixXcd W_active(params_.K_, Q_active);
            for (int j = 0; j < Q_active; ++j) {
                W_active.col(j) = W_est_.col(active_indices[j]);
            }

            Eigen::VectorXcd h_active = (W_active.adjoint() * X_l.adjoint() * X_l * W_active).inverse() * W_active.adjoint() * X_l.adjoint() * Y_avg.transpose();
            double residual = (Y_avg.transpose() - X_l * W_active * h_active).squaredNorm();

            double beta = (double)params_.K_ / residual;
            double logL = params_.K_ * std::log(beta) - params_.K_ * std::log(M_PI) - params_.K_;
            double aic = -2.0 * (logL - 2.0 * Q_active);

            if (aic < min_aic) {
                min_aic = aic;
                best_mask_int = i;
            }
        }

        std::vector<int> result(8, 0);
        for (int bit = 0; bit < 8; ++bit) {
            if ((best_mask_int >> bit) & 1) {
                result[bit] = 1;
            }
        }
        return result;
    }

    /**
     * ステップ 1 (Raghavendra版): Raghavendra AIC 8パス総当たりによる最良マスクの探索
     */
    std::vector<int> findBestMaskByExhaustiveRaghavendraAIC_8paths(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        X_l = X_avg.asDiagonal();

        double min_gaic = 1e300;
        int best_mask_int = 0;

        for (int i = 1; i < 256; ++i) {
            std::vector<int> active_indices;
            for (int bit = 0; bit < 8; ++bit) {
                if ((i >> bit) & 1) active_indices.push_back(bit);
            }
            int Q_active = static_cast<int>(active_indices.size());

            Eigen::MatrixXcd W_active(params_.K_, Q_active);
            for (int j = 0; j < Q_active; ++j) W_active.col(j) = W_est_.col(active_indices[j]);

            Eigen::VectorXcd h_active = (W_active.adjoint() * X_l.adjoint() * X_l * W_active).inverse() * W_active.adjoint() * X_l.adjoint() * Y_avg.transpose();
            double residual = (Y_avg.transpose() - X_l * W_active * h_active).squaredNorm();
            double sigma_est_sq = residual / (double)params_.K_;

            double V_l = (double)params_.K_ * std::log(sigma_est_sq + 1e-18) / 2.0;
            double penalty = 2.0 * std::log(std::log((double)params_.K_ + 1e-12)) * ((double)Q_active + 1);
            double gaic = V_l + penalty;

            if (gaic < min_gaic) {
                min_gaic = gaic;
                best_mask_int = i;
            }
        }

        std::vector<int> result(8, 0);
        for (int bit = 0; bit < 8; ++bit) {
            if ((best_mask_int >> bit) & 1) result[bit] = 1;
        }
        return result;
    }

    /**
     * Exhaustive search (全探索) を用いて Raghavendra の GAIC を評価し、最良マスクを求めて h_l と noiseVariance_ を設定する。
     */
    void set_initial_params_by_exhaustive_RaghavendraAIC(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();

        double min_gaic = 1e300;
        int best_mask_int = 0;
        double best_sigma_est = 0.0;
        Eigen::VectorXcd best_h_active;

        const int max_mask = (1 << params_.Q_) - 1;
        for (int mask = 1; mask <= max_mask; ++mask) {
            std::vector<int> active_indices;
            for (int bit = 0; bit < params_.Q_; ++bit) {
                if ((mask >> bit) & 1) active_indices.push_back(bit);
            }

            int Q_active = static_cast<int>(active_indices.size());
            if (Q_active == 0) continue;

            Eigen::MatrixXcd W_active(params_.K_, Q_active);
            for (int j = 0; j < Q_active; ++j) W_active.col(j) = W_est_.col(active_indices[j]);

            Eigen::VectorXcd h_active = (W_active.adjoint() * X_l.adjoint() * X_l * W_active).inverse() * W_active.adjoint() * X_l.adjoint() * Y_avg.transpose();

            double residual = (Y_avg.transpose() - X_l * W_active * h_active).squaredNorm();
            double sigma_est = residual / (double)params_.K_;

            double V_l = (double)params_.K_ * std::log(sigma_est + 1e-18) / 2.0;
            double penalty = 2.0 * std::log(std::log((double)params_.K_ + 1e-12));
            double gaic = V_l + penalty * ((double)Q_active + 1);

            if (gaic < min_gaic) {
                min_gaic = gaic;
                best_mask_int = mask;
                best_sigma_est = sigma_est;
                best_h_active = h_active;
            }
        }

        this->h_l = Eigen::VectorXcd::Zero(params_.Q_);
        if (best_mask_int != 0) {
            int idx = 0;
            for (int bit = 0; bit < params_.Q_; ++bit) {
                if ((best_mask_int >> bit) & 1) {
                    this->h_l(bit) = best_h_active(idx);
                    idx++;
                }
            }
        }

        this->noiseVariance_ = best_sigma_est;
    }

    // Wrapper法によるAICモデル選択付き等化
    double equalizeWithWrapperAIC(const Eigen::MatrixXcd& X)
    {
        set_initial_params_by_pilot(X);
        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int p = 0; p < params_.NUMBER_OF_PILOT; ++p) H_est_.row(p) = H_init;

        double total_iterations_sum = 0.0;
        int dataSymbolCount = params_.L_ - params_.NUMBER_OF_PILOT;

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            activePathIndices_.clear();
            for(int q=0; q<params_.Q_; ++q) activePathIndices_.push_back(q);

            runEMLoop(l, X); // フルモデルで収束

            Eigen::VectorXcd h_full = h_l;
            double noise_full = noiseVariance_;

            std::vector<std::pair<double, int>> pathRank;
            for (int q = 0; q < params_.Q_; ++q) {
                pathRank.push_back({ std::norm(h_full(q)), q });
            }
            std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

            double min_aic = 1e18; // 十分大きな値
            Eigen::VectorXcd best_h_l = h_full;
            double best_noise = noise_full;
            int best_iter_count = 0;

            for (int q = 1; q <= params_.Q_; ++q) {
                activePathIndices_.clear();
                for(int i=0; i<q; ++i) activePathIndices_.push_back(pathRank[i].second);
                std::sort(activePathIndices_.begin(), activePathIndices_.end());

                h_l.setZero();
                for(int idx : activePathIndices_) h_l(idx) = h_full(idx);
                noiseVariance_ = noise_full;

                int iter = runEMLoop(l, X);

                double beta = 1.0 / noiseVariance_;
                Eigen::VectorXcd Y_vec = Y_.row(l).transpose();
                Eigen::VectorXcd H_est = W_est_ * h_l;
                Eigen::VectorXcd XH = X_bar * H_est;

                double term1 = Y_vec.squaredNorm();
                double term2 = (Y_vec.adjoint() * XH).value().real();
                double term3 = (XH.adjoint() * Y_vec).value().real();
                double term4 = (H_est.adjoint() * R_moment * H_est).value().real();
                double logL = params_.K_ * std::log(beta) - beta * (term1 - term2 - term3 + term4);
                double aic = -2.0 * (logL - 2.0 * q);

                if (aic < min_aic) {
                    min_aic = aic;
                    best_h_l = h_l;
                    best_noise = noiseVariance_;
                    best_iter_count = iter;
                }
            }

            h_l = best_h_l;
            noiseVariance_ = best_noise;
            total_iterations_sum += best_iter_count;

            H_est_.row(l) = (W_est_ * h_l).transpose();
            for (int k = 0; k < params_.K_; k++)
            {
                R_(l, k) = Y_(l, k) / (W_est_.row(k) * h_l)(0);
            }
        }

        return total_iterations_sum / static_cast<double>(dataSymbolCount);
    }

    // 埋め込み法（Embedded Method）による等化
    double equalizeWithEmbeddedAIC(const Eigen::MatrixXcd& X)
    {
        set_initial_params_by_pilot(X);
        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int p = 0; p < params_.NUMBER_OF_PILOT; ++p) H_est_.row(p) = H_init;

        double total_iterations_sum = 0.0;
        int dataSymbolCount = params_.L_ - params_.NUMBER_OF_PILOT;
        const int MAX_ITER = 50;
        const int MIN_ITER = 3;

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            Eigen::VectorXi symbol_prev2(params_.K_); symbol_prev2.setConstant(-1);
            Eigen::VectorXi symbol_prev1(params_.K_); symbol_prev1.setConstant(-1);
            Eigen::VectorXi symbol_current(params_.K_);
            Eigen::VectorXd obj(params_.NUMBER_OF_SYMBOLS);

            int current_iter_count = 0;

            for (int iter = 0; iter < MAX_ITER; iter++)
            {
                current_iter_count = iter + 1;

                Estep(l, X);

                Eigen::MatrixXcd G = W_est_.adjoint() * R_moment * W_est_;
                Eigen::VectorXcd b = W_est_.adjoint() * X_bar.adjoint() * Y_.row(l).transpose();
                Eigen::VectorXcd h_full = G.inverse() * b;

                std::vector<std::pair<double, int>> pathRank;
                for (int q = 0; q < params_.Q_; ++q) {
                    pathRank.push_back({ std::norm(h_full(q)), q });
                }
                std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

                double min_aic = 1e18;
                Eigen::VectorXcd best_h_l = h_full;
                double best_noise = 0.0;

                for (int q = 1; q <= params_.Q_; ++q) {
                    std::vector<int> current_indices;
                    for(int i=0; i<q; ++i) current_indices.push_back(pathRank[i].second);
                    
                    Eigen::VectorXcd h_sub = Eigen::VectorXcd::Zero(params_.Q_);
                    for(int idx : current_indices) h_sub(idx) = h_full(idx);

                    Eigen::VectorXcd Y_est = X_bar * W_est_ * h_sub;
                    double residual = (Y_.row(l).transpose() - Y_est).squaredNorm();
                    double beta_est = (double)params_.K_ / residual;
                    
                    double logL = params_.K_ * std::log(beta_est) - params_.K_ * std::log(M_PI) - params_.K_;
                    double aic = -2.0 * (logL - 2.0 * q);

                    if (aic < min_aic) {
                        min_aic = aic;
                        best_h_l = h_sub;
                        best_noise = 1.0 / beta_est;
                    }
                }

                h_l = best_h_l;
                noiseVariance_ = best_noise;

                for (int k = 0; k < params_.K_; k++)
                {
                    std::complex<double> H_current_k = (W_est_.row(k) * h_l)(0);
                    if (std::abs(H_current_k) < 1e-9) H_current_k = 1e-9;
                    
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
                    if ((symbol_current == symbol_prev1) && (symbol_prev1 == symbol_prev2)){
                        break;
                    }
                }
                symbol_prev2 = symbol_prev1;
                symbol_prev1 = symbol_current;
            }

            total_iterations_sum += current_iter_count;

            H_est_.row(l) = (W_est_ * h_l).transpose();
            for (int k = 0; k < params_.K_; k++)
            {
                std::complex<double> val = (W_est_.row(k) * h_l)(0);
                if(std::abs(val) < 1e-9) val = 1e-9;
                R_(l, k) = Y_(l, k) / val;
            }
        }

        return total_iterations_sum / static_cast<double>(dataSymbolCount);
    }

    /**
     * AICで選択されたパスモデルが真のモデルと一致しているか判定する
     */
    bool checkAICAccuracy(const std::vector<int>& trueMask)
    {
        for (int q = 0; q < params_.Q_; ++q)
        {
            bool isSelected = (std::norm(h_l(q)) > 1e-20); 
            bool isTrue = (trueMask[q] != 0);
            if (isSelected != isTrue) {
                return false;
            }
        }
        return true;
    }

    /**
     * 文献の式(61)に基づく伝送路推定評価値（SNR劣化比）を計算する
     */
    double getSNRDegradationMetric(double noiseSD)
    {
        double sum_metric = 0.0;
        int count = 0;
        double beta = 1.0 / (noiseSD * noiseSD);
        double sigma_X_sq = 1.0; 

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                double norm_H_sq = std::norm(H_true_(l, k));
                if (norm_H_sq < 1e-20) continue;

                double term1 = ((std::norm(std::conj(H_est_(l, k)) * H_true_(l, k)) / norm_H_sq )- 1) * sigma_X_sq * beta;
                double term2 = std::norm(H_est_(l, k)) / norm_H_sq;
                
                sum_metric += (term1 + term2);
                count++;
            }
        }
        
        if (count == 0) return 1.0;
        return sum_metric / count;
    }

    /**
     * パイロットAIC固定パス法
     */
    double equalizeWithPilotAICFixedPath(const Eigen::MatrixXcd& X)
    {
        set_initial_params_by_pilot(X);

        activePathIndices_.clear();
        for(int q = 0; q < params_.Q_; ++q) {
            if(std::norm(h_l(q)) > 1e-20) { 
                activePathIndices_.push_back(q);
            }
        }
        std::sort(activePathIndices_.begin(), activePathIndices_.end());

        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int p = 0; p < params_.NUMBER_OF_PILOT; ++p) H_est_.row(p) = H_init;

        double total_iterations_sum = 0.0;
        int dataSymbolCount = params_.L_ - params_.NUMBER_OF_PILOT;

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            int iter = runEMLoop(l, X);
            total_iterations_sum += iter;

            H_est_.row(l) = (W_est_ * h_l).transpose();
            for (int k = 0; k < params_.K_; k++)
            {
                std::complex<double> val = (W_est_.row(k) * h_l)(0);
                if(std::abs(val) < 1e-12) val = 1e-12;
                R_(l, k) = Y_(l, k) / val;
            }
        }

        return total_iterations_sum / static_cast<double>(dataSymbolCount);
    }

    // ゲッター群
    Eigen::VectorXcd getEstimatedPathCoefficients() const {
        return h_l;
    }

    const Eigen::MatrixXcd &getH_est() const { return H_est_; }
    const Eigen::MatrixXcd &getH_true() const { return H_true_; }
    const Eigen::MatrixXcd &getY() const { return Y_; }

    /**
     * ΔH を計算して返すヘルパー
     */
    Eigen::RowVectorXcd computeDeltaHRow(int l) const
    {
        if (l < 0 || l >= params_.L_) {
            throw std::out_of_range("l is out of range in computeDeltaHRow");
        }
        return H_est_.row(l) - H_true_.row(l);
    }

    /**
     * 添付画像の x = -H/A を使って delta 行ベクトルを作る
     */
    Eigen::RowVectorXcd computeDeltaRowFromPhotoX(int l,
                                                  double noiseSD,
                                                  double Px = 1.0,
                                                  double beta_override = -1.0) const
    {
        if (l < 0 || l >= params_.L_) {
            throw std::out_of_range("l is out of range in computeDeltaRowFromPhotoX");
        }

        double beta = (beta_override > 0.0) ? beta_override : 1.0 / (noiseSD * noiseSD);
        const double eps = 1e-12;

        Eigen::RowVectorXcd delta = Eigen::RowVectorXcd::Zero(params_.K_);
        for (int k = 0; k < params_.K_; ++k)
        {
            const std::complex<double> H = H_true_(l, k);
            const double A = std::norm(H) * Px * beta + 1.0;
            if (A < eps) {
                continue;
            }
            delta(k) = -H / A;
        }
        return delta;
    }

    // ΔH（RowVector）を与えて各サブキャリアの γ_k を返す。
    Eigen::VectorXd computeGammaFromDeltaH(const Eigen::RowVectorXcd &DeltaH_row,
                                           int l,
                                           double noiseSD,
                                           double Px = 1.0,
                                           double beta_override = -1.0,
                                           bool returnPerSubcarrier = true) const
    {
        if (DeltaH_row.size() != params_.K_) {
            throw std::invalid_argument("DeltaH_row size does not match params_.K_");
        }
        if (l < 0 || l >= params_.L_) {
            throw std::out_of_range("l is out of range");
        }

        double beta = (beta_override > 0.0) ? beta_override : 1.0 / (noiseSD * noiseSD);
        const double eps = 1e-12;
        Eigen::VectorXd gamma(params_.K_);

        for (int k = 0; k < params_.K_; ++k) {
            double H_abs_sq = std::norm(H_true_(l, k));
            if (H_abs_sq < 1e-20) {
                gamma(k) = 0.0;
                continue;
            }
            double DH_abs_sq = std::norm(DeltaH_row(k));
            double cross_abs_sq = std::norm(DeltaH_row(k) * H_true_(l, k));
            double plas_abs_sq = std::norm(H_true_(l, k) + DeltaH_row(k));

            double numerator = H_abs_sq * H_abs_sq * Px * beta;
            double denominator = cross_abs_sq * Px * beta
                                 + plas_abs_sq
                                 + eps;

            gamma(k) = numerator / denominator;
        }

        if (!returnPerSubcarrier) {
            Eigen::VectorXd out(1);
            out(0) = gamma.mean();
            return out;
        }

        Eigen::VectorXd gamma_dB = 10.0 * gamma.array().log10();
        
        return gamma_dB;
    }

    Eigen::VectorXd computeGamma_78(int l,
                                    double noiseSD,
                                    double Px = 1.0,
                                    double beta_override = -1.0,
                                    bool returnPerSubcarrier = true) const
    {
        if (l < 0 || l >= params_.L_) {
            throw std::out_of_range("l is out of range");
        }

        double beta = (beta_override > 0.0) ? beta_override : 1.0 / (noiseSD * noiseSD);
        Eigen::VectorXd gamma(params_.K_);

        for (int k = 0; k < params_.K_; ++k) {
            double H_abs_sq = std::norm(H_true_(l, k));
            if (H_abs_sq < 1e-20) {
                gamma(k) = 0.0;
                continue;
            }

            double numerator = H_abs_sq * Px * beta + 1;
            gamma(k) = numerator;
        }

        if (!returnPerSubcarrier) {
            Eigen::VectorXd out(1);
            out(0) = gamma.mean();
            return out;
        }

        Eigen::VectorXd gamma_dB = 10.0 * gamma.array().log10();
        
        return gamma_dB;
    }

    /**
     * @brief 仮定するパス数 q = 1 ... Q_ に対する通常AICおよびRaghavendra GAICを計算する
     */
    std::pair<std::vector<double>, std::vector<double>> calculateAICvsQ(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        X_l = X_avg.asDiagonal();

        Eigen::VectorXcd h_full = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_avg.transpose();

        std::vector<std::pair<double, int>> pathRank;
        for (int q = 0; q < params_.Q_; ++q) {
            pathRank.push_back({ std::norm(h_full(q)), q });
        }
        std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

        std::vector<double> aic_values(params_.Q_, 0.0);
        std::vector<double> gaic_values(params_.Q_, 0.0);

        for (int q = 1; q <= params_.Q_; ++q) {
            std::vector<int> selected_indices;
            for (int i = 0; i < q; ++i) {
                selected_indices.push_back(pathRank[i].second);
            }
            std::sort(selected_indices.begin(), selected_indices.end());

            Eigen::MatrixXcd W_tilde(params_.K_, q);
            for (int i = 0; i < q; ++i) {
                W_tilde.col(i) = W_est_.col(selected_indices[i]);
            }

            Eigen::VectorXcd h_active = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * W_tilde.adjoint() * X_l.adjoint() * Y_avg.transpose();
            
            double residual = (Y_avg.transpose() - X_l * W_tilde * h_active).squaredNorm();

            double beta = (double)params_.K_ / residual;
            double logL = params_.K_ * std::log(beta) - params_.K_ * std::log(M_PI) - params_.K_;
            aic_values[q - 1] = -2.0 * (logL - 2.0 * q);

            double sigma_est = residual / (double)params_.K_;
            double V_l = (double)params_.K_ * std::log(sigma_est + 1e-18) / 2.0;
            double penalty = 2.0 * std::log(std::log((double)params_.K_ + 1e-12)) * ((double)q + 1);
            gaic_values[q - 1] = V_l + penalty;
        }

        return {aic_values, gaic_values};
    }

private:
    const SimulationParameters &params_;
    const Eigen::MatrixXcd &W_est_;
    std::vector<int> grayNum_;
    Eigen::VectorXcd symbol_;

    Eigen::MatrixXi rxData_;
    Eigen::MatrixXcd Y_;
    Eigen::MatrixXcd H_est_;
    Eigen::MatrixXcd H_true_;
    double noiseVariance_;
    cnormal_distribution<> unitCNormalRand_;
    std::vector<int> activePathIndices_;

    Eigen::MatrixXcd X_l;
    Eigen::MatrixXcd R_;
    Eigen::VectorXd xPro;
    Eigen::MatrixXcd X_bar;
    Eigen::MatrixXd R_moment;
    Eigen::VectorXcd h_l;

    void setSymbol() {
        int M = params_.NUMBER_OF_SYMBOLS;
        int sqrtM = sqrt(M);
        double P = 1.0 / (2.0 * (M - 1) / 3.0);

        int i = 0;
        for (int v1 = 0; v1 < sqrtM; v1++) {
            for (int v2 = 0; v2 < sqrtM; v2++) {
                symbol_(i).real((2 * v1 - (sqrtM - 1)) * sqrt(P));
                if (v1 % 2 == 0) {
                    symbol_(i).imag((2 * v2 - (sqrtM - 1)) * sqrt(P));
                } else {
                    symbol_(i).imag(((sqrtM - 1) - 2 * v2) * sqrt(P));
                }
                i++;
            }
        }
    }

    int grayCode(int num) {
        return num ^ (num >> 1);
    }

    void setGrayNum() {
        for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++) {
            grayNum_[i] = grayCode(i);
        }
    }

    void set_initial_params_by_pilot(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();
        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_avg.transpose();

        std::vector<std::pair<double, int>> pathRank;
        for (int q = 0; q < params_.Q_; ++q) {
            pathRank.push_back({ std::norm(h_l(q)), q });
        }

        std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

        std::vector<double> aic_list(params_.Q_);
        std::vector<double> beta_list(params_.Q_);

        for (int Q_tilde = 1; Q_tilde <= params_.Q_; ++Q_tilde) {
            std::vector<int> selected_indices;
            for (int i = 0; i < Q_tilde; ++i) {
                selected_indices.push_back(pathRank[i].second);
            }
            std::sort(selected_indices.begin(), selected_indices.end());

            Eigen::MatrixXcd W_tilde(params_.K_, Q_tilde);
            for (int i = 0; i < Q_tilde; ++i) {
                int original_idx = selected_indices[i];
                W_tilde.col(i) = W_est_.col(original_idx);
            }

            Eigen::VectorXcd h_active = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * W_tilde.adjoint() * X_l.adjoint() * Y_avg.transpose();
            double residual = (Y_avg.transpose() - X_l * W_tilde * h_active).squaredNorm();
            beta_list[Q_tilde - 1] = (double)params_.K_ / residual;

            double logL = params_.K_ * std::log(beta_list[Q_tilde - 1]) - params_.K_ * std::log(M_PI) - params_.K_;
            aic_list[Q_tilde - 1] = -2.0 * (logL - 2.0 * Q_tilde);
        }

        int best_idx = std::distance(aic_list.begin(), std::min_element(aic_list.begin(), aic_list.end()));
        int best_Q_tilde = best_idx + 1;

        std::vector<int> final_indices;
        for (int i = 0; i < best_Q_tilde; ++i) {
            final_indices.push_back(pathRank[i].second);
        }
        std::sort(final_indices.begin(), final_indices.end());

        Eigen::MatrixXcd W_final(params_.K_, best_Q_tilde);
        for (int i = 0; i < best_Q_tilde; ++i) {
            W_final.col(i) = W_est_.col(final_indices[i]);
        }

        Eigen::VectorXcd h_final_active = (W_final.adjoint() * X_l.adjoint() * X_l * W_final).inverse() * W_final.adjoint() * X_l.adjoint() * Y_avg.transpose();

        this->h_l = Eigen::VectorXcd::Zero(params_.Q_);
        for (int i = 0; i < best_Q_tilde; ++i) {
            this->h_l(final_indices[i]) = h_final_active(i);
        }

        this->noiseVariance_ = 1.0 / beta_list[best_idx];
    }

    void set_initial_params_by_pilot_power_sort_RaghavendraGAIC(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();
        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_avg.transpose();

        std::vector<std::pair<double, int>> pathRank;
        for (int q = 0; q < params_.Q_; ++q) {
            pathRank.push_back({ std::norm(h_l(q)), q });
        }

        std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

        std::vector<double> gaic_list(params_.Q_);
        std::vector<double> residual_list(params_.Q_);

        for (int Q_tilde = 1; Q_tilde <= params_.Q_; ++Q_tilde) {
            std::vector<int> selected_indices;
            for (int i = 0; i < Q_tilde; ++i) {
                selected_indices.push_back(pathRank[i].second);
            }
            std::sort(selected_indices.begin(), selected_indices.end());

            Eigen::MatrixXcd W_tilde(params_.K_, Q_tilde);
            for (int i = 0; i < Q_tilde; ++i) {
                int original_idx = selected_indices[i];
                W_tilde.col(i) = W_est_.col(original_idx);
            }

            Eigen::VectorXcd h_active = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * W_tilde.adjoint() * X_l.adjoint() * Y_avg.transpose();
            double residual = (Y_avg.transpose() - X_l * W_tilde * h_active).squaredNorm();
            residual_list[Q_tilde - 1] = residual;

            double sigma_est = residual / (double)params_.K_;
            double V_l = (double)params_.K_ * std::log(sigma_est + 1e-18) / 2.0;
            double penalty = 2.0 * std::log(std::log((double)params_.K_ + 1e-12)) * ((double)Q_tilde + 1);
            gaic_list[Q_tilde - 1] = V_l + penalty;
        }

        int best_idx = std::distance(gaic_list.begin(), std::min_element(gaic_list.begin(), gaic_list.end()));
        int best_Q_tilde = best_idx + 1;

        std::vector<int> final_indices;
        for (int i = 0; i < best_Q_tilde; ++i) {
            final_indices.push_back(pathRank[i].second);
        }
        std::sort(final_indices.begin(), final_indices.end());

        Eigen::MatrixXcd W_final(params_.K_, best_Q_tilde);
        for (int i = 0; i < best_Q_tilde; ++i) {
            W_final.col(i) = W_est_.col(final_indices[i]);
        }

        Eigen::VectorXcd h_final_active = (W_final.adjoint() * X_l.adjoint() * X_l * W_final).inverse() * W_final.adjoint() * X_l.adjoint() * Y_avg.transpose();

        this->h_l = Eigen::VectorXcd::Zero(params_.Q_);
        for (int i = 0; i < best_Q_tilde; ++i) {
            this->h_l(final_indices[i]) = h_final_active(i);
        }

        this->noiseVariance_ = residual_list[best_idx] / (double)params_.K_;
    }

    Eigen::RowVectorXcd set_initial_params_by_RaghavendraAIC_Update(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();
        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * (W_est_.adjoint() * X_l.adjoint() * Y_avg.transpose());

        int P = params_.Q_ - 1; 
        Eigen::VectorXcd current_R = Y_avg.transpose();
        std::vector<int> selected_taps; 

        while (P >= 0) {
            std::vector<double> GAIC_list(P + 1);
            for(int i = 0; i <= P; i++){
                Eigen::VectorXcd h_temp = Eigen::VectorXcd::Zero(params_.Q_);
                for(int t = 0; t <= i; t++){
                    h_temp(t) = h_l(t);
                }
                double sigma_est = ((current_R - X_l * W_est_ * h_temp).squaredNorm()) / (double)params_.K_;
                double V_l = (params_.K_ * std::log(sigma_est + 1e-12)) / 2.0;
                GAIC_list[i] = V_l + 2.0 * std::log(std::log(params_.K_)) * (i + 1); 
            }
            int best_L = std::distance(GAIC_list.begin(), std::min_element(GAIC_list.begin(), GAIC_list.end()));
            selected_taps.push_back(best_L);

            current_R = current_R - X_l * W_est_.col(best_L) * h_l(best_L); 
            
            P = best_L - 1; 
            if (P >= 0) {
                Eigen::MatrixXcd W_rem = W_est_.leftCols(P + 1);
                h_l.head(P + 1) = (W_rem.adjoint() * X_l.adjoint() * X_l * W_rem).inverse() * (W_rem.adjoint() * X_l.adjoint() * current_R);
            }
        }

        std::sort(selected_taps.begin(), selected_taps.end());
        int num_active_taps = selected_taps.size();

        Eigen::MatrixXcd W_tilde(params_.K_, num_active_taps);
        for (int j = 0; j < num_active_taps; ++j) {
            W_tilde.col(j) = W_est_.col(selected_taps[j]);
        }

        Eigen::VectorXcd h_ils = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * (W_tilde.adjoint() * X_l.adjoint() * Y_avg.transpose());

        return (W_tilde * h_ils).transpose();
    }

    Eigen::RowVectorXcd set_initial_params_by_RaghavendraAIC_Update2(const Eigen::MatrixXcd& X)
    {
        Eigen::RowVectorXcd X_avg = X.topRows(params_.NUMBER_OF_PILOT).colwise().mean();
        Eigen::RowVectorXcd Y_avg = Y_.topRows(params_.NUMBER_OF_PILOT).colwise().mean();

        X_l = X_avg.asDiagonal();
        h_l = (W_est_.adjoint() * X_l.adjoint() * X_l * W_est_).inverse() * W_est_.adjoint() * X_l.adjoint() * Y_avg.transpose();

        int P = params_.Q_ - 1; 
        Eigen::VectorXcd current_R = Y_avg.transpose();
        std::vector<int> selected_taps; 

        while (P >= 0) {
            std::vector<double> GAIC_list(P + 1);
            for(int i = 0; i <= P; i++){
                Eigen::MatrixXcd W_L = W_est_.leftCols(i + 1);
                Eigen::VectorXcd h_L_est = (W_L.adjoint() * X_l.adjoint() * X_l * W_L).inverse() * (W_L.adjoint() * X_l.adjoint() * current_R);

                Eigen::VectorXcd h_temp = Eigen::VectorXcd::Zero(params_.Q_);
                h_temp.head(i + 1) = h_L_est;

                double sigma_est = ((current_R - X_l * W_est_ * h_temp).squaredNorm()) / (double)params_.K_;
                double V_l = (params_.K_ * std::log(sigma_est + 1e-12)) / 2.0;
                GAIC_list[i] = V_l + 2.0 * std::log(std::log(params_.K_)) * (i + 1); 
            }
            int best_L = std::distance(GAIC_list.begin(), std::min_element(GAIC_list.begin(), GAIC_list.end()));
            selected_taps.push_back(best_L);

            h_l(best_L) = 0.0;                                                                                                                                                                                      
            current_R = X_l * W_est_ * h_l;

            P = best_L - 1;
        }

        std::sort(selected_taps.begin(), selected_taps.end());
        int num_active_taps = selected_taps.size();

        Eigen::MatrixXcd W_tilde(params_.K_, num_active_taps);
        for (int j = 0; j < num_active_taps; ++j) {
            W_tilde.col(j) = W_est_.col(selected_taps[j]);
        }

        Eigen::VectorXcd h_ils = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * (W_tilde.adjoint() * X_l.adjoint() * Y_avg.transpose());

        return (W_tilde * h_ils).transpose();
    }

    void equalizeChannelWithPilot(const Eigen::MatrixXcd& X)
    {
        set_initial_params_by_pilot(X);
        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            H_est_.row(l) = H_init;
            for (int k = 0; k < params_.K_; k++)
            {
                R_(l, k) = Y_(l, k) / (W_est_.row(k) * h_l)(0);
            }
        }
    }

    double equalizeChannelWithEM(const Eigen::MatrixXcd& X)
    {
        set_initial_params_by_pilot(X);

        Eigen::RowVectorXcd H_init = (W_est_ * h_l).transpose();
        for (int l = 0; l < params_.NUMBER_OF_PILOT; l++)
        {
            H_est_.row(l) = H_init;
        }

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
                Estep(l, X);
                Mstep(l);

                for (int k = 0; k < params_.K_; k++)
                {
                    std::complex<double> H_current_k = (W_est_.row(k) * h_l)(0);
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

    int runEMLoop(int l, const Eigen::MatrixXcd& X, int max_iter = 100) {
        const int MIN_ITER = 3;
        Eigen::VectorXi symbol_prev2(params_.K_);
        Eigen::VectorXi symbol_prev1(params_.K_);
        Eigen::VectorXi symbol_current(params_.K_);
        Eigen::VectorXd obj(params_.NUMBER_OF_SYMBOLS);
        
        symbol_prev2.setConstant(-1);
        symbol_prev1.setConstant(-1);

        int current_iter_count = 0;

        for (int iter = 0; iter < max_iter; iter++)
        {
            current_iter_count = iter + 1;
            
            Estep(l, X);
            Mstep(l);

            for (int k = 0; k < params_.K_; k++)
            {
                std::complex<double> H_current_k = (W_est_.row(k) * h_l)(0);
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
        return current_iter_count;
    }

    void Estep(int l, const Eigen::MatrixXcd& X)
    {
        Eigen::VectorXcd H_current = W_est_ * h_l;
        double variance = noiseVariance_;
        if (variance < std::numeric_limits<double>::epsilon()) {
            variance = std::numeric_limits<double>::epsilon();
        }
        for (int k = 0; k < params_.K_; k++)
        {
            bool isPilot = params_.enableDataPilots && (k == 5 || k == 19 || k == 32 || k == 46);
            if (isPilot)
            {
                X_bar(k, k) = X(l, k);
                R_moment(k, k) = std::norm(X(l, k));
            }
            else
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
                    posterior_prob(i) = std::exp(log_likelihoods(i) - log_max);
                    sumxP_shifted += posterior_prob(i);
                }

                if (sumxP_shifted <= 0.0) {
                    std::string error_msg = "FATAL ERROR: sumxP_shifted is zero or negative at l=" 
                              + std::to_string(l) + ", k=" + std::to_string(k);
                    std::cerr << error_msg << std::endl;
                    throw std::runtime_error(error_msg);   
                } else {
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
    }

    void Mstep(int l)
    {
        int n_active = activePathIndices_.size();
        if (n_active == 0) return;

        Eigen::MatrixXcd W_active(params_.K_, n_active);
        for(int i = 0; i < n_active; ++i) {
            W_active.col(i) = W_est_.col(activePathIndices_[i]);
        }

        Eigen::MatrixXcd A = X_bar * W_active;
        Eigen::MatrixXcd B = W_active.adjoint() * R_moment * W_active;
        Eigen::VectorXcd h_active = B.inverse() * A.adjoint() * Y_.row(l).transpose();

        h_l.setZero();
        for(int i = 0; i < n_active; ++i) {
            h_l(activePathIndices_[i]) = h_active(i);
        }

        Eigen::VectorXcd Wh = W_active * h_active;
        double term1 = (Y_.row(l).transpose() - X_bar * Wh).squaredNorm();

        Eigen::MatrixXcd Cov_X = R_moment - X_bar.adjoint() * X_bar;
        double term2 = (Wh.adjoint() * Cov_X * Wh).value().real();

        noiseVariance_ = (term1 + term2) / (double)params_.K_;
        
        if (noiseVariance_ < 1e-10) noiseVariance_ = 1e-10;
    }

    void setRxDataByML()
    {
        Eigen::VectorXd obj(params_.NUMBER_OF_SYMBOLS);

        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            for (int k = 0; k < params_.K_; k++)
            {
                if (params_.enableDataPilots && (k == 5 || k == 19 || k == 32 || k == 46))
                {
                    continue; // Skip pilot carriers
                }
                for (int i = 0; i < params_.NUMBER_OF_SYMBOLS; i++)
                {
                    obj(i) = std::norm((R_(l, k) - symbol_(i)));
                }
                Eigen::VectorXd::Index minColumn;
                obj.minCoeff(&minColumn);
                rxData_(l, k) = minColumn;
            }
        }
    }

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

#endif // RECEIVER_H
