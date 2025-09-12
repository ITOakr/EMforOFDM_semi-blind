/*
 * File:   simulator.h
 * Author: Ito
 *
 * Created on 2024/12/20, 18:31
*/

#ifndef SIMULATOR_H
#define SIMULATOR_H
#define _USE_MATH_DEFINES
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Eigen>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>
#include <cmath>
#include <math.h>
#include "random_collection.h"

class Simulator {
    public:
        Simulator() {
           ////試行回数を設定
           //std::cout << "--------------------------------------------------------------------" << std::endl;
           //std::cout << "TRIAL?" << std::endl;
           //std::cout << "--------------------------------------------------------------------" << std::endl;
           //std::cin >> NUMBER_OF_TRIAL;

            // リサイズ
            W_.resize(K_ , Q_);
            xi_.resize(Q_);
            h_.resize(L_ , Q_);
            Cmat_.resize(L_ , L_);
            U_.resize(L_, L_);
            lambda_.resize(L_, L_);
            H_.resize(L_ , K_);
            H_est_.resize(K_);

            txData_.resize(L_ , K_);
            rxData_.resize(L_ , K_);
            Y_.resize(L_ , K_);
            X_.resize(L_ , K_);
            X_l.resize(K_, K_);
            R_.resize(L_ , K_);
            symbol_.resize(NUMBER_OF_SYMBOLS);
            xPro.resize(NUMBER_OF_SYMBOLS);
            X_bar.Zero(K_, K_);
            R_moment.Zero(K_, K_);
            h_l.resize(Q_);

            // DFT行列設定
            setW_();

            // シンボル設計
            setSymbol();

            // 乱数設定
            unitIntUniformRand_.init(0, NUMBER_OF_SYMBOLS - 1, seed);
            unitCNormalRand_.init(0.0, 1 / sqrt(2), seed);

            // 遅延プロファイル生成
            setChannelProfile();
        }

        virtual ~Simulator() {
        }


        /**
         * 雑音の分散設定
         * @param EbN0dB EbN0 [dB]
         */
        void setNoiseSD(double EbN0dB) {
            noiseSD_ = std::sqrt(std::pow(10.0, -0.1 * EbN0dB) / (double)NUMBER_OF_BIT);
        }

        //試行回数の設定
        // void setTrialNum(double EbN0dB) {
        //     if(EbN0dB == 0 || EbN0dB == 5 || EbN0dB == 10) {
        //         NUMBER_OF_TRIAL = 10000;
        //     }
        //     else if(EbN0dB == 15 || EbN0dB == 20) {
        //         NUMBER_OF_TRIAL = 100000;
        //     }
        //     else {
        //         NUMBER_OF_TRIAL = 1000000;
        //     }
        // }

        void setTrialNum(double EbN0dB) {
            NUMBER_OF_TRIAL = 100;
        }
        
        /**
         * ドップラー周波数設定
         * @param f_d ドップラー周波数
         */
        void setDopplerFrequence(double f_d) {
            f_d_ = f_d;
        }


        /**
         * 数値計算実験
         * @return ビット誤り率のシミュレーション値
         */
        double getBER_EM_Simulation() {
            int count = 0;
            for(int tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                setX_();
                setH_();
                setY_();

                equalizeChannelWithEM();
                setRxDataByML();

                count += getBitErrorCount();
            }
            return (double)count / ((double)NUMBER_OF_TRIAL * (double)NUMBER_OF_BIT * (double)K_ * ((double)L_ - NUMBER_OF_PILOT));
        }

        /**
         * 数値計算実験
         * @return 平均二乗誤差のシミュレーション値
         */
        double getMSESimulation() {
            double mse = 0.0;
            for(int tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                setX_();
                setH_();
                setY_();
                estimateChannelByML();
                mse += getMeanSquaredError();
            }
            return mse / NUMBER_OF_TRIAL;
        }

        //Hのノルム（1になる）
        double Hcheck() {
            double Hcheck = 0.0;
            for(auto tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                setH_();
                for(auto l = 0; l < L_; l++) {
                    for(auto k = 0; k < K_; k++) {
                        Hcheck += std::norm(H_(l, k));
                    }
                }
            }
            return Hcheck / (double)(L_ * K_) / (double)NUMBER_OF_TRIAL;
        }

        //hのノルム（1になる）
        double hcheck() {
            double hcheck = 0.0;
            for(auto tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                seth_();
                for(auto l = 0; l < L_; l++) {
                    for(auto q = 0; q < Q_; q++) {
                        hcheck += std::norm(h_(l, q));
                    }
                }
            }
            return hcheck / (double)L_ /(double)NUMBER_OF_TRIAL;
        }

        void Hset_H() {
            setX_();
            setH_();
            setY_();
            estimateChannelByML();
        }

        Eigen::MatrixXcd hh_l() {
            Eigen::MatrixXcd hh_l;
            Eigen::VectorXcd h_l;
            hh_l.resize(Q_, Q_);
            h_l.resize(Q_);
            for(auto tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                seth_();
                for(auto q = 0; q < Q_; q++) {
                    h_l(q) = h_(0, q);
                }
                hh_l += h_l * h_l.adjoint();
                std::cout << h_l << std::endl;
                std::cout << h_l.adjoint() <<std::endl;
            }
            return hh_l / (double)NUMBER_OF_TRIAL;
        }

    private:
        double noiseSD_;                            // 雑音の標準偏差
        int NUMBER_OF_TRIAL;                        // 試行回数
        const int K_ = 52;                          // サブキャリア数 K
        const int L_ = 3;                          // 1フレームのシンボル数 L
        const int Q_ = 2;                          // 伝送路のインパルス応答のパス数 Q
        const double T_ = 3.2 * std::pow(10, -6);                      // 有効シンボル長 T
        const double Tgi_ = 0.8 * std::pow(10, -6);                    // ガードインターバル長 Tgi
        const double Ts_ = T_ + Tgi_;               // シンボル全体の長さ
        const int NUMBER_OF_FFT = 64;               // FFTポイント数(IEEE802.11a)
        const int NUMBER_OF_PILOT = 2;              // パイロットシンボル個数
        const int NUMBER_OF_SYMBOLS = 2;            // 変調方式に合わせて(BPSKの例)
        const int NUMBER_OF_BIT = 1;                // 変調方式に合わせて(BPSKの例)

        double f_d_;                                // ドップラー周波数

        Eigen::MatrixXcd W_;            // DFT行列:式(17)
        Eigen::VectorXd xi_;        // 伝送路のインパルス応答の遅延プロファイル
        Eigen::MatrixXcd h_;            // インパルス応答
        Eigen::MatrixXd Cmat_;         // 共分散行列
        Eigen::MatrixXcd U_;
        Eigen::MatrixXcd lambda_;         
        Eigen::MatrixXcd H_;            // 周波数応答
        Eigen::MatrixXi txData_;        // 送信アルファベット
        Eigen::MatrixXi rxData_;        // 受信アルファベット
        Eigen::VectorXcd symbol_;       // 送信可能なシンボルベクトル
        Eigen::MatrixXcd Y_;            // 受信信号
        Eigen::MatrixXcd X_;            // 送信信号
        Eigen::MatrixXcd X_l;           // ある時刻の送信信号の対角行列
        Eigen::MatrixXcd R_;            // 等化後の受信信号
        Eigen::VectorXcd H_est_;        // 伝送路の周波数応答の推定値
        Eigen::VectorXd xPro;           // 送信信号の事後確率
        Eigen::MatrixXcd X_bar;         // Xの期待値
        Eigen::MatrixXd R_moment;       // 原点周りの2次モーメント
        Eigen::VectorXcd h_l;           // ある時刻のインパルス応答
    


        // 乱数用
        int seed = 100;
        uniform_int_distribution<> unitIntUniformRand_;     // int型一様乱数
        cnormal_distribution<> unitCNormalRand_;            // 平均0，分散1の複素正規分布（実部，虚部それぞれ平均0，分散0.5の正規分布）

        // DFT行列Wの生成:式(17)
        void setW_() {
            for (int q = 0; q < Q_; ++q) {
                for (int k = 0; k < K_ / 2; ++k) {
                    // -26から-1番目のキャリヤ
                    W_(k, q) = std::polar(1.0, -2.0 * M_PI * ((double)k - (double)K_ / 2.0) * (double)q / (double)NUMBER_OF_FFT);
                    // 1から26番目のキャリヤ
                    W_(k + K_ / 2, q) = std::polar(1.0, -2.0 * M_PI * ((double)k + 1.0) * (double)q / (double)NUMBER_OF_FFT);
                }
            }
        }

        void setChannelProfile(){
            for(int q = 0; q < Q_; q++){
                if(q == 0){
                    xi_(q) = 1.0;
                }
                else{
                    xi_(q) = xi_(q - 1) * std::pow(10.0, -0.1);
                }
            }

            double tmp = xi_.sum();
            
            for(int q = 0; q < Q_; q++){
                xi_(q) = xi_(q) / tmp;
            }
        }

        /**
         * 共分散行列生成
         */
        void setCmat_() {
            for(auto l_1 = 0; l_1 < L_; l_1++) {
                for(auto l_2 = 0; l_2 < L_; l_2++) {
                    Cmat_(l_1, l_2) = std::cyl_bessel_j(0, 2.0 * M_PI * (double)abs((l_1 - l_2) * Ts_) * f_d_); //Jakesモデル
                }
            } 
        }

        /**
         * インパルス応答行列生成
         */
        void seth_() {
            Eigen::VectorXcd h_q(L_);
            Eigen::MatrixXcd A(L_, L_);
            Eigen::VectorXcd x(L_);

            for(auto q = 0; q < Q_; q++) {
                setCmat_();
                Cmat_ = xi_(q) * Cmat_;

                computeUnitaryAndDiagonal();
                //computeUnitaryAndDiagonal_test();
                
                A = U_ * lambda_.array().sqrt().matrix();

                for(auto l = 0; l < L_; l++){
                    x(l) = unitCNormalRand_();
                }

                h_q = A * x;

                for(auto l = 0; l < L_; l++) {
                    h_(l, q) = h_q(l);
                }
            }
        }

        // 固有値分解によりユニタリ行列と対角行列を生成する関数
        void computeUnitaryAndDiagonal() {
            // 固有値分解を実行
            Eigen::EigenSolver<Eigen::MatrixXd> solver(Cmat_);

            // ユニタリ行列（固有ベクトルを列に持つ行列）
            U_ = solver.eigenvectors();

            // 対角行列（固有値を対角要素に持つ行列）
            lambda_ = solver.eigenvalues().asDiagonal();
        }

        // 固有値分解によりユニタリ行列と対角行列を生成する関数
        void computeUnitaryAndDiagonal_test() {
            // 固有値分解を実行
            Eigen::EigenSolver<Eigen::MatrixXd> solver(Cmat_);

            // ユニタリ行列（固有ベクトルを列に持つ行列）
            U_ = solver.eigenvectors();

            // 対角行列（固有値を対角要素に持つ行列）
            Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

            // NaN をゼロに置き換え
            for (int i = 0; i < eigenvalues.size(); ++i) {
                if (std::isnan(eigenvalues[i])) {
                    eigenvalues[i] = 0.0;
                }
            }

            // 対角行列を生成
            lambda_ = eigenvalues.asDiagonal();
        }


        /**
         * 送信信号生成
         */
        void setX_() {
            for(int l = 0; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    txData_(l, k) = unitIntUniformRand_();
                    X_(l, k) = symbol_(txData_(l, k));
                }
            }
        }

        /**
         * 周波数応答生成
         */
        void setH_() {
            // 伝送路のインパルス応答の生成
            seth_();
            // 伝送路の周波数応答の生成:式(19)
            for(int l = 0; l < L_; l++){
                H_.row(l) = W_ * h_.row(l);
            }
        }

        /**
         * 受信信号生成
         */
        void setY_() {
            for (int l = 0; l < L_; l++) {
                for (int k = 0; k < K_; k++) {
                    Y_(l, k) = H_(l, k) * X_(l, k) + noiseSD_ * unitCNormalRand_();
                }
            }
        }

        //パイロットシンボルからｈの初期値を得る
        void seth_l_byPilot(){
            X_l = X_.row(0).asDiagonal();
            h_l = (W_.adjoint()*X_l.adjoint()*X_l*W_).inverse()*W_.adjoint()*X_l.adjoint()*Y_.row(0).transpose();
        }

        void equalizeChannelWithEM(){
            seth_l_byPilot();

            const int MAX_ITER = 10;
            for(int l = NUMBER_OF_PILOT; l < L_; l++) {
                for(int iter = 0; iter < MAX_ITER; iter++) {
                    //Eステップ
                    Estep(l);
                    //Mステップ
                    Mstep(l);
                }
                for(int k = 0; k < K_; k++) {
                    R_(l, k) = Y_(l, k) / (W_.row(k) * h_l);
                }
            }
        }

        void Estep(int l) {
            Eigen::VectorXd H_current = W_ * h_l;

            for(int k = 0; k < K_; k++) {
                double sumxP = 0.0;
                Eigen::VectorXd posterior_prob(NUMBER_OF_SYMBOLS);
                for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                    std::complex<double> s = symbol_(i);
                    double norm = std::norm(Y_(l, k) - H_current(k) * s);
                    double variance = noiseSD_ * noiseSD_;
                    posterior_prob(i) = std::exp(-norm / variance);
                    sumxP += posterior_prob(i);
                }

                posterior_prob /= sumxP;

                std::complex<double> expected_X = 0.0;
                double expected_X_norm_sq = 0.0;

                for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                    expected_X += posterior_prob(i) * symbol_(i);
                    expected_X_norm_sq += posterior_prob(i) * std::norm(symbol_(i));
                }
                X_bar(k, k) = expected_X;
                R_moment(k, k) = expected_X_norm_sq;
            }
        }

        void Mstep(int l) {
            h_l = (W_.adjoint()*R_moment*W_).inverse()*W_.adjoint()*X_bar.adjoint()*Y_.row(l).transpose();
        }

        /**
         * シンボル生成
         */
        void setSymbol() {
            symbol_(0) = -1.0;
            symbol_(1) = 1.0;
        }
        
        /**
         * 最尤復調
         */
        void setRxDataByML() {
            Eigen::VectorXd obj(NUMBER_OF_SYMBOLS);     // 最小化の目的関数

            for(int l = NUMBER_OF_PILOT; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                        // 最尤復調の周波数応答は推定値を使う？
                        obj(i) = std::norm((R_(l, k) - symbol_(i)));
                    }
                    Eigen::VectorXd::Index minColumn;       // ノルムが最小な index（つまり受信データ）
                    obj.minCoeff(&minColumn);
                    rxData_(l, k) = minColumn;
                }
            }
        }

        /**
         * ビット誤り数のカウント
         * @return 全ての誤りビット数
         */
        int getBitErrorCount() {
            int count = 0;
            for(int l = NUMBER_OF_PILOT; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    count += hammingDistance(txData_(l, k), rxData_(l, k));
                }
            }
            return count;
        }

        /**
         * ハミング距離計算
         * @param 整数1，整数2
         * @return ハミング距離
         */
        int hammingDistance(int num1, int num2) {
            int ham = 0;
            int xorResult;
            int bitMask = 1;

            xorResult = num1 ^ num2;

            for(int i = 0; i < NUMBER_OF_BIT; i++) {
                ham += (xorResult & bitMask) >> i;
                bitMask <<= 1;
            }
            
            return ham;
        }
};

#endif /* SIMULATOR_H */