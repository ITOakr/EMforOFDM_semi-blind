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

    public:
    // Wrapper法によるAICモデル選択付き等化
    double equalizeWithWrapperAIC()
    {
        // 初期化: パイロットから初期推定
        set_initial_params_by_pilot();
        H_est_.row(0) = (W_est_ * h_l).transpose();

        double total_iterations_sum = 0.0;
        int dataSymbolCount = params_.L_ - params_.NUMBER_OF_PILOT;

        // データシンボルごとのループ
        for (int l = params_.NUMBER_OF_PILOT; l < params_.L_; l++)
        {
            // ---------------------------------------------------------
            // Step 1: フルモデル (全パス) でEMを収束させる
            // ---------------------------------------------------------
            activePathIndices_.clear();
            for(int q=0; q<params_.Q_; ++q) activePathIndices_.push_back(q);

            // 初期値をセット (前のシンボル or パイロットの推定値を維持してスタート)
            // ここでは h_l は前回の結果が残っているのでそのまま使う(Warm Start)
            
            runEMLoop(l); // フルモデルで収束

            // フルモデルの結果を保存
            Eigen::VectorXcd h_full = h_l;
            double noise_full = noiseVariance_;

            // ---------------------------------------------------------
            // Step 2: パスの電力ランキング作成
            // ---------------------------------------------------------
            std::vector<std::pair<double, int>> pathRank;
            for (int q = 0; q < params_.Q_; ++q) {
                pathRank.push_back({ std::norm(h_full(q)), q });
            }
            // 降順ソート
            std::sort(pathRank.begin(), pathRank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

            // ---------------------------------------------------------
            // Step 3: 候補モデルごとのWrapper評価
            // ---------------------------------------------------------
            double min_aic = 1e18; // 十分大きな値
            Eigen::VectorXcd best_h_l = h_full;
            double best_noise = noise_full;
            int best_iter_count = 0;

            // 上位 q 個のパスを使うモデルを順次評価
            for (int q = 1; q <= params_.Q_; ++q) {

                // アクティブパスの設定
                activePathIndices_.clear();
                for(int i=0; i<q; ++i) activePathIndices_.push_back(pathRank[i].second);
                std::sort(activePathIndices_.begin(), activePathIndices_.end()); // インデックス順に整理

                // パラメータのリセット (フルモデルの結果から射影してスタートするのが効率的)
                h_l.setZero();
                for(int idx : activePathIndices_) h_l(idx) = h_full(idx);
                noiseVariance_ = noise_full;

                // このモデルでEMを回す
                int iter = runEMLoop(l);

                // AICの計算
                double beta = 1.0 / noiseVariance_;
                Eigen::VectorXcd Y_vec = Y_.row(l).transpose(); // 受信信号 Y
                Eigen::VectorXcd H_est = W_est_ * h_l;          // 推定チャネル H = W * h
                Eigen::VectorXcd XH = X_bar * H_est;            // X_bar * H

                // 各項をスカラ(double)として計算
                // 第1項: ||Y||^2
                double term1 = Y_vec.squaredNorm();

                // 第2項 & 第3項: - Y^H * X_bar * H - (X_bar * H)^H * Y
                // これは「-2 * Re( Y^H * (X_bar * H) )」と同じです
                double term2 = (Y_vec.adjoint() * XH).value().real();
                double term3 = (XH.adjoint() * Y_vec).value().real();

                // 第4項: H^H * R * H
                // (H^H * R * H) はエルミート形式なので必ず実数になります
                double term4 = (H_est.adjoint() * R_moment * H_est).value().real();         // H^H * (R * H)
                // 対数尤度の簡易計算 (定数項は比較において無視可能だが、ここでは元のコードに合わせる)
                double logL = params_.K_ * std::log(beta) - beta * (term1 - term2 - term3 + term4);
                // double logL = params_.K_ * std::log(beta) - params_.K_ * std::log(M_PI) - params_.K_;
                double aic = -2.0 * (logL - 2.0 * q); // k はパラメータ数 (複素数なので自由度 2k)

                // 最良モデルの更新
                if (aic < min_aic) {
                    min_aic = aic;
                    best_h_l = h_l;
                    best_noise = noiseVariance_;
                    best_iter_count = iter; // 便宜上、最良モデルの反復数を記録
                }
            }

            // ---------------------------------------------------------
            // Step 4: 最良モデルの結果を採用
            // ---------------------------------------------------------
            h_l = best_h_l;
            noiseVariance_ = best_noise;
            total_iterations_sum += best_iter_count;

            // 結果の格納
            H_est_.row(l) = (W_est_ * h_l).transpose();
            for (int k = 0; k < params_.K_; k++)
            {
                // 等化後の信号
                R_(l, k) = Y_(l, k) / (W_est_.row(k) * h_l)(0);
            }
        }

        return total_iterations_sum / static_cast<double>(dataSymbolCount);
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
        std::vector<double> beta_list(params_.Q_);             // 各パスモデルにおける雑音精度

        // std::cout << "W_est_ >> " << W_est_.transpose() << std::endl;

        for (int Q_tilde = 1; Q_tilde <= params_.Q_; ++Q_tilde) {

            // インデックスのソート
            std::vector<int> selected_indices;
            for (int i = 0; i < Q_tilde; ++i) {
                selected_indices.push_back(pathRank[i].second);
            }
            std::sort(selected_indices.begin(), selected_indices.end());

            // W_tilde の生成
            Eigen::MatrixXcd W_tilde(params_.K_, Q_tilde);
            for (int i = 0; i < Q_tilde; ++i) {
                int original_idx = selected_indices[i]; // ソート済みの上位インデックスを取得
                W_tilde.col(i) = W_est_.col(original_idx); // 対応するDFT行列の列をコピー
            }

            // std::cout << "W_tilde >> " << W_tilde << std::endl;

            // 各パスモデルにおける推定値を計算
            Eigen::VectorXcd h_active = (W_tilde.adjoint() * X_l.adjoint() * X_l * W_tilde).inverse() * W_tilde.adjoint() * X_l.adjoint() * Y_.row(0).transpose();
            // std::cout << "h_l >> " << h_tilde_list[Q_tilde - 1] << std::endl;
            double residual = (Y_.row(0).transpose() - X_l * W_tilde * h_active).squaredNorm();
            beta_list[Q_tilde - 1] = (double)params_.K_ / residual;

            // AIC の計算
            double logL = params_.K_ * std::log(beta_list[Q_tilde - 1]) - params_.K_ * std::log(M_PI) - params_.K_;
            aic_list[Q_tilde - 1] = -2.0 * (logL - 2.0 * Q_tilde);
            // std::cout << "aic >> " << aic_list[Q_tilde - 1] << std::endl;
        }

        // std::cout << "OK3" << std::endl;

        // 最良モデルの計算フェーズ
        // AIC が最小となるインデックスを特定
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

        Eigen::VectorXcd h_final_active = (W_final.adjoint() * X_l.adjoint() * X_l * W_final).inverse() * W_final.adjoint() * X_l.adjoint() * Y_.row(0).transpose();
        // std::cout << "h_final_active >> " << h_final_active << std::endl;

        // フルサイズ配列への展開（非採用パスは0埋め）
        this->h_l = Eigen::VectorXcd::Zero(params_.Q_);
        for (int i = 0; i < best_Q_tilde; ++i) {
            this->h_l(final_indices[i]) = h_final_active(i);
        }

        // 6. パラメータ更新
        this->noiseVariance_ = 1.0 / beta_list[best_idx];
        // std::cout << "h_l >> " << h_l << std::endl;
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

    // 指定されたパス(activePathIndices_)でEMアルゴリズムを回す関数
    // 戻り値: 収束までの反復回数
    int runEMLoop(int l, int max_iter = 100) {
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
            
            // Eステップ
            Estep(l);
            
            // Mステップ (activePathIndices_ に基づいて更新)
            Mstep(l);

            // 硬判定シンボルの更新（収束判定用）
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
        // アクティブなパスの数
        int n_active = activePathIndices_.size();
        if (n_active == 0) return; // 安全策

        // 1. アクティブなパスに対応する W の部分行列を作成
        Eigen::MatrixXcd W_active(params_.K_, n_active);
        for(int i = 0; i < n_active; ++i) {
            W_active.col(i) = W_est_.col(activePathIndices_[i]);
        }

        // 2. 縮小モデルでの h (サイズ: n_active x 1) の推定
        // h_active = (W_active^H * R_moment * W_active)^-1 * (X_bar * W_active)^H * Y
        Eigen::MatrixXcd A = X_bar * W_active;
        Eigen::MatrixXcd B = W_active.adjoint() * R_moment * W_active;
        Eigen::VectorXcd h_active = B.inverse() * A.adjoint() * Y_.row(l).transpose();

        // 3. 全体の h_l (サイズ: Q x 1) にマッピング（非アクティブは0にする）
        h_l.setZero();
        for(int i = 0; i < n_active; ++i) {
            h_l(activePathIndices_[i]) = h_active(i);
        }

        // 4. 雑音分散の更新
        // J = ||Y - X_bar * W * h||^2 + h^H * W^H * (R_moment - X_bar^H * X_bar) * W * h
        Eigen::VectorXcd Wh = W_active * h_active;
        
        // 第1項: 残差ノルム
        double term1 = (Y_.row(l).transpose() - X_bar * Wh).squaredNorm();

        // 第2項: 推定誤差補正
        Eigen::MatrixXcd Cov_X = R_moment - X_bar.adjoint() * X_bar;
        double term2 = (Wh.adjoint() * Cov_X * Wh).value().real();

        noiseVariance_ = (term1 + term2) / (double)params_.K_;
        
        // 数値安定性のための下限処理
        if (noiseVariance_ < 1e-10) noiseVariance_ = 1e-10;
    }

    // void Mstep(int l)
    // {
    //     Eigen::MatrixXcd A = X_bar * W_est_;
    //     h_l = (W_est_.adjoint() * R_moment * W_est_).inverse() * (X_bar * W_est_).adjoint() * Y_.row(l).transpose();
    //     // std::cout << "h_l=" << h_l << std::endl;
    //     // std::cout << "A=" << A << std::endl;
    //     noiseVariance_ = ((Y_.row(l).transpose() - X_bar * W_est_ * h_l).squaredNorm() + ((W_est_ * h_l).adjoint() * (R_moment - X_bar.adjoint() * X_bar) * W_est_ * h_l).value().real()) / (double)params_.K_;
    // }

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