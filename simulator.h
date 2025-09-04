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
            W_.resize(Q_ , K_);
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
        double getBERSimulation() {
            int count = 0;
            for(int tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                setX_();
                setH_();
                setY_();
                estimateChannelByML();
                equalizeByEstimatedChannel();
                setRxDataByML();
                //equalizeByEstimatedChannel_test();
                //setRxDataByML_test();
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

    //private:
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
                    W_(q, k) = std::polar(1.0, -2.0 * M_PI * ((double)k - (double)K_ / 2.0) * (double)q / (double)NUMBER_OF_FFT);
                    // 1から26番目のキャリヤ
                    W_(q, k + K_ / 2) = std::polar(1.0, -2.0 * M_PI * ((double)k + 1.0) * (double)q / (double)NUMBER_OF_FFT);
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
            H_ = h_ * W_;
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

        /**
         * 平均二乗誤差（L2ノルム）(RMSE)
         * @return 誤差のL2ノルム
         */
        double getMeanSquaredError() {
            double sum = 0;
            for(int k = 0; k < K_; k++) {
                for(int l = 0; l < L_; l++) {
                    sum += std::norm(H_(l, k) - H_est_(k));
                }
            }
            return sqrt(sum / (double)L_ / (double)K_);
        }

        /**
         * MLによる伝送路推定
         * @return 誤差のL2ノルム
         */
        void estimateChannelByML() {
            Eigen::VectorXcd sum(K_);
            sum.resize(K_);
            sum.setZero();
            for(int k = 0; k < K_; k++) {
                for(int i = 0; i < NUMBER_OF_PILOT; i++) {
                    // パイロットシンボル区間での推定
                    sum(k) += Y_(i, k) / X_(i, k);
                }
            }
            H_est_ = sum / (double)NUMBER_OF_PILOT;
        }

        /**
         * 等化
         */
        void equalizeByEstimatedChannel() {
            for(int l = NUMBER_OF_PILOT; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    R_(l, k) = Y_(l, k) / H_est_(k);
                }
            }
        }

         void equalizeByEstimatedChannel_test() {
            for(int l = 0; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    R_(l, k) = Y_(l, k) / H_(l,k);
                }
            }
        }

        /**
         * シンボル生成
         */
        void setSymbol() {
            symbol_(0) = -1.0;
            symbol_(1) = 1.0;
        }
        
        //Eステップ
        // void Estep(int l) {
        //     for(int k = 0; k < K_; k++) {
        //         double sumxP = 0.0;

        //         for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
        //             auto s = symbol_(i);
        //             double expArg;
        //             expArg = -std::norm(Y_.col(l) - W_ * h_l * s);
        //             xPro[i] = std::exp(expArg);
        //             sumxP += xPro[i];
        //         }

        //         for(int idx = 0; idx < NUMBER_OF_SYMBOLS; idx++) {
        //             xPro[idx] = xPro[idx]/sumxP;
        //         }

        //         std::complex<double> xExp = 0.0;
        //         double rVar = 0.0;

        //         for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
        //             xExp += xPro[i] * symbol_(i);
        //             rVar += xPro[i] * std::norm(symbol_(i));
        //         }
        //         X_bar(k, k) = xExp;
        //         R_moment(k, k) = rVar; 
        //     }
        // }

        void Mstep(int l){

        }
        
        //パイロットシンボルからｈの初期値を得る
        void setH_byPilot(){
            h_l = (W_.adjoint()*X_.col(0).adjoint()*X_.col(0)*W_).inverse()*W_.adjoint()*X_.col(0).adjoint()*Y_.col(0);
        }

        // EMアルゴリズムでHを作る
        // void setH_byEM(){
        //     for(auto l = 1; l < L_; l++){
        //         for( int iter = 0; iter < 10; iter++){
        //             Estep(l);
        //             Mstep(l);
        //         }
        //     }
        // }

        /**
         * 最尤復調
         */
        void setRxDataByML() {
            Eigen::VectorXd obj(NUMBER_OF_SYMBOLS);     // 最小化の目的関数

            for(int l = NUMBER_OF_PILOT; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                        // 最尤復調の周波数応答は推定値を使う？
                        obj(i) = std::norm(H_est_(k)) * std::norm((R_(l, k) - symbol_(i)));
                    }
                    Eigen::VectorXd::Index minColumn;       // ノルムが最小な index（つまり受信データ）
                    obj.minCoeff(&minColumn);
                    rxData_(l, k) = minColumn;
                }
            }
        }

        void setRxDataByML_test() {
            Eigen::VectorXd obj(NUMBER_OF_SYMBOLS);     // 最小化の目的関数

            for(int l = 0; l < L_; l++) {
                for(int k = 0; k < K_; k++) {
                    for(int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                        // 最尤復調の周波数応答は推定値を使う？
                        obj(i) = std::norm(H_(l, k)) * std::norm((R_(l, k) - symbol_(i)));
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