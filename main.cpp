/*
 * File:   power_ber.cpp
 * Author: Ito
 *
 * Created on 2024/12/20, 19:36
 */
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <chrono>
#include "simulator.h"
#include "parameters.h"

// パラメータ
static const int EbN0dBmin = 0;
static const int EbN0dBmax = 30;
static const int EbN0dBstp = 5;

static const int dopplerMin = 0;
static const double dopplerMax = 0.02;
static const double dopplerStep = 0.0002;

// ファイル
std::string fileName;       // ファイル名
std::ofstream ofs;        // 出力ファイル  
const std::string outputDir = "C:/Users/Akira Ito/code2025/results/Sim_20260204~/";

// BER
double ber;

// MSE
double mse;

double avgPower;

// ドップラー周波数
double dopplerFrequency;

// 試行回数
int numberOfTrials;

std::string getCurrentTimeString() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now;
    
    // Windows(MSVC)の場合は localtime_s、Linux(GCC)等の場合は localtime_r や localtime を使用
#if defined(_MSC_VER)
    localtime_s(&tm_now, &time_t_now);
#else
    tm_now = *std::localtime(&time_t_now);
#endif

    std::ostringstream oss;
    // フォーマット: YYYYMMDD_HHMMSS (例: 20241220_193600)
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

std::string getModulationSchemeName(int numberOfBits) {
    switch (numberOfBits) {
        case 2: return "QPSK";
        case 4: return "16QAM";
        case 6: return "64QAM";
        case 8: return "256QAM";
        default: return "Unknown";
    }
}

int main()
{
	SimulationParameters params;

    std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << "Number of Bit? (QPSK:2, 16QAM:4, 64QAM:6, 256QAM:8)" << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cin >> params.NUMBER_OF_BIT;
	params.NUMBER_OF_SYMBOLS = std::pow(2, params.NUMBER_OF_BIT);

	Simulator sim(params);

	int mode_select = 0;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << "Select Simulation Mode" << std::endl;
	std::cout << "1: Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "2: Doppler sweep (fixed Eb/N0)" << std::endl;
    std::cout << "3: MSE Simulation" << std::endl;
    std::cout << "4: Average Power Simulation" << std::endl;
	std::cout << "5: only pilot Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "6: MSE Simulation pilot" << std::endl;
	std::cout << "7: only pilot Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "8: MSE Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "9: MSE Doppler sweep only pilot (fixed Eb/N0)" << std::endl;
	std::cout << "10: Channel Magnitude Response (|H(k, l)|) CSV Output" << std::endl;
	std::cout << "11: Noise Variance MSE Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "12: H MSE by initial h_est Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "13: H MSE by initial H_est Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "14: MSE vs Doppler sweep (parallel) (fixed Eb/N0)" << std::endl;
	std::cout << "15: AIC Model Selection Accuracy vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "16: AIC Path Selection Accuracy vs Doppler (fixed Eb/N0)" << std::endl;
    std::cout << "17: AIC F-Measure vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "18: Embedded AIC Method MSE vs Eb/N0" << std::endl;
	std::cout << "19: Wrapper AIC Method SNR Degradation Ratio vs Doppler" << std::endl;
	std::cout << "20: Pilot Only SNR Degradation Ratio vs Doppler" << std::endl;
	std::cout << "21: Pilot AIC Fixed Path MSE vs Doppler" << std::endl;
	std::cout << "23: Export Tx Waveform (k=10) vs Time" << std::endl;
    std::cout << "24: Export Faded Tx Waveform (k=10) vs Time" << std::endl;
	std::cout << "25: Export Channel Magnitude (|H|) (k=10) vs Time" << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cin >> mode_select;

    std::cout << "Enter number of trials:" ;
	std::cin >> numberOfTrials;
	sim.setTrialNum(numberOfTrials);

    std::string modulationName = getModulationSchemeName(params.NUMBER_OF_BIT);
	std::string timeStr = getCurrentTimeString();

	if (mode_select == 1)
	{
		// --- モード1: Eb/N0スイープ ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_BER_vs_EbN0_2path.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			ber = sim.getBER_EM_Simulation();
			double avg_iter = sim.getAverageIterations();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", BER = " << ber << std::endl;
			std::cout << "Avg_Iteration = " << avg_iter << std::endl;
			ofs << EbN0dB << "," << ber << std::endl;
		}
	}
	else if (mode_select == 2)
	{
		// --- モード2: ドップラー周波数スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_BER_vs_Doppler_2path.csv";
		ofs.open(fileName);

		sim.setNoiseSD(fixedEbN0dB);

		for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
			sim.setDopplerFrequency(dopplerFrequency);

			ber = sim.getBER_EM_Simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "f_dT_s = " << dopplerFrequency << ", BER = " << ber << std::endl;
			ofs << dopplerFrequency << "," << ber << std::endl;
		}
	}
    else if (mode_select == 3)
	{
		// --- モード3: Eb/N0スイープ ---
		double dopplerFrequency;
        double mse;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "MSE_2path.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			mse = sim.getMSE_simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}
	}
    else if (mode_select == 4)
    {
        // --- モード4: 平均電力計算 ---
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:" ;
        std::cin >> dopplerFrequency;

        sim.setDopplerFrequency(dopplerFrequency);
        
        avgPower = sim.getAveragePower_simulation();
        
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Average power of the true channel response H over " << numberOfTrials << " trials." << std::endl;
        std::cout << "f_d*T_s = " << dopplerFrequency << ", Average Power = " << avgPower << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;

        // このモードではファイル出力はせず、コンソール表示のみとします。
        return 0; // 結果を表示して終了
    }
	else if (mode_select == 5)
	{
		// --- モード1: Eb/N0スイープ ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_BER_vs_EbN0_2path_pilot.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			ber = sim.getBER_Simulation_only_pilot();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", BER = " << ber << std::endl;
			ofs << EbN0dB << "," << ber << std::endl;
		}
	}
	else if (mode_select == 6)
	{
		// --- モード3: Eb/N0スイープ ---
		double dopplerFrequency;
        double mse;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "MSE_pilot.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			mse = sim.getMSE_simulation_only_pilot();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}
	}
	else if (mode_select == 7) 
	{
		// --- モード7: パイロットのみのドップラー周波数スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		// ファイル名の設定
		fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_BER_vs_Doppler_2path_only_pilot.csv";
		ofs.open(fileName);

		sim.setNoiseSD(fixedEbN0dB);

		for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
			sim.setDopplerFrequency(dopplerFrequency);

			ber = ber = sim.getBER_Simulation_only_pilot();
			
			std::cout << "-----------" << std::endl;
			std::cout << "f_dT_s = " << dopplerFrequency << ", BER = " << ber << std::endl;
			ofs << dopplerFrequency << "," << ber << std::endl;
		}
	}
	else if (mode_select == 8)
	{
		// --- モード8: MSE vs ドップラー周波数スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_MSE_vs_Doppler_2path.csv";
		ofs.open(fileName);

		sim.setNoiseSD(fixedEbN0dB);
		double mse;

		for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
			sim.setDopplerFrequency(dopplerFrequency);

			mse = sim.getMSE_simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
			ofs << dopplerFrequency << "," << mse << std::endl;

			// --- ここで次回の値を決定 ---
            if (dopplerFrequency == 0.0) {
                // 現在が0なら、次は最小ステップ幅(0.0002)にする
                dopplerFrequency = dopplerStep;
            } else {
                // それ以外なら2倍にする (0.0002 -> 0.0004 -> 0.0008...)
                dopplerFrequency *= 2.0;
            }
		}
	}
	else if (mode_select == 9)
	{
		// --- モード9: MSE vs ドップラー周波数スイープ (パイロットのみ) ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_MSE_vs_Doppler_2path_only_pilot.csv";
		ofs.open(fileName);

		sim.setNoiseSD(fixedEbN0dB);
		double mse;

		for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
			sim.setDopplerFrequency(dopplerFrequency);

			mse = sim.getMSE_simulation_only_pilot();
			
			std::cout << "-----------" << std::endl;
			std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
			ofs << dopplerFrequency << "," << mse << std::endl;

            // ★追加: ドップラー周波数の更新ロジック (0 -> step -> 2倍...)
            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
		}
	}
	else if (mode_select == 10) // ★ モード10の追加
	{
		// --- モード10: 周波数応答の大きさ (|H(k, l)|) のCSV出力 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s for single trial:" ;
		std::cin >> dopplerFrequency;

		// ファイル名の設定: 固定の f_d*T_s に基づくファイル名
		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Channel_Magnitude_Response.csv";
		ofs.open(fileName);

		sim.setDopplerFrequency(dopplerFrequency); // 周波数応答生成のために設定
		
		// 新しいシミュレーション関数を呼び出す
		sim.saveChannelMagnitudeResponseToCSV(ofs, dopplerFrequency);
		
		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Channel magnitude response saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}
	else if (mode_select == 11)
	{
		// --- モード11: 雑音分散 MSE vs ドップラー周波数スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_NoiseVar_MSE_vs_Doppler_2path.csv";
		ofs.open(fileName);

		// シミュレータにノイズSDを設定（データ生成用）
		sim.setNoiseSD(fixedEbN0dB);
		double noise_var_mse;
        double dopplerFrequency;

		for (dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
			sim.setDopplerFrequency(dopplerFrequency);
            
            // 新しいシミュレーション関数を呼び出す
			noise_var_mse = sim.getNoiseVarianceMSE_simulation(fixedEbN0dB);
            
            double fdTs = dopplerFrequency; // main.cppの他のモードに合わせて、dopplerFrequencyがf_d*T_sの値を表すと仮定
            
			std::cout << "-----------" << std::endl;
			std::cout << "f_dT_s = " << fdTs << ", Noise Var MSE = " << noise_var_mse << std::endl;
			ofs << fdTs << "," << noise_var_mse << std::endl;
		}
	}
	else if (mode_select == 12)
	{
		// --- モード12: Eb/N0スイープ,h_est_MSE 式67の確認---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_h_est_MODE12.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			mse = sim.get_h_MSE_Simulation_during_pilot();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}
	}
	else if (mode_select == 13)
	{
		// --- モード13: Eb/N0スイープ,H_est_MSE 式69の確認---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_H_est_MODE13.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			mse = sim.get_H_est_MSE_Simulation_during_pilot();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}
	}
	else if (mode_select == 14)
	{
		// --- モード8: MSE vs ドップラー周波数スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_MSE_vs_Doppler_2path_parallel.csv";
		ofs.open(fileName);

		sim.setNoiseSD(fixedEbN0dB);
		double mse;

		// ★合計時間の計測開始
        auto start_total = std::chrono::high_resolution_clock::now();

		for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
			sim.setDopplerFrequency(dopplerFrequency);

			mse = sim.getMSE_simulation_parallel();

			std::cout << "-----------" << std::endl;
			std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
			ofs << dopplerFrequency << "," << mse << std::endl;

			// --- ここで次回の値を決定 ---
            if (dopplerFrequency == 0.0) {
                // 現在が0なら、次は最小ステップ幅(0.0002)にする
                dopplerFrequency = dopplerStep;
            } else {
                // それ以外なら2倍にする (0.0002 -> 0.0004 -> 0.0008...)
                dopplerFrequency *= 2.0;
            }
		}
		// ★合計時間の計測終了と表示
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        
        std::cout << "========================================" << std::endl;
        std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
        std::cout << "========================================" << std::endl;
	}
	else if (mode_select == 15)
    {
        // --- モード15: AICモデル選択正答率 vs ドップラー周波数 ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_AIC_Model_Accuracy_vs_Doppler.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);
        double accuracy;

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            accuracy = sim.getAICAccuracy_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", Accuracy = " << accuracy << std::endl;
            ofs << dopplerFrequency << "," << accuracy << std::endl;

             // ドップラー周波数の更新 (0 -> step -> 2*step ...)
            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
        }
    }
	else if (mode_select == 16)
    {
        // --- モード16: 正答率 (Accuracy) のみ出力 ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        // ファイル名に Accuracy を明記
        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_AIC_Path_Accuracy.csv";
        ofs.open(fileName);
        
        ofs << "f_dT_s,Accuracy" << std::endl;

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            // 計算実行 (両方計算されるが、Accuracyだけ使う)
            std::pair<double, double> result = sim.getAIC_Metrics_pilot();
            double accuracy = result.second;
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", Accuracy = " << accuracy << std::endl;
            
            ofs << dopplerFrequency << "," << accuracy << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
    }
    else if (mode_select == 17)
    {
        // --- モード17: F値 (F-Measure) のみ出力 ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        // ファイル名に F_Measure を明記
        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_AIC_F_Measure.csv";
        ofs.open(fileName);
        
        ofs << "f_dT_s,F_Measure" << std::endl;

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            // 計算実行 (両方計算されるが、F-Measureだけ使う)
            std::pair<double, double> result = sim.getAIC_Metrics_pilot();
            double f_measure = result.first;
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", F-Measure = " << f_measure << std::endl;
            
            ofs << dopplerFrequency << "," << f_measure << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
    }
	else if (mode_select == 18)
    {
        // --- モード18: 埋め込み法 MSE vs Doppler (Eb/N0固定) ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]: ";
        std::cin >> fixedEbN0dB;

        // ファイル名を変更
        fileName = outputDir + timeStr + "_" + modulationName + "_EmbeddedAIC_MSE_vs_Doppler_EbN0_" + std::to_string(fixedEbN0dB) + ".csv";
        ofs.open(fileName);

        sim.setNoiseSD((double)fixedEbN0dB);

		// ★合計時間の計測開始
        auto start_total = std::chrono::high_resolution_clock::now();

        // ドップラー周波数をスイープ
        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            
            // 埋め込み法シミュレーション (MSE)
            double mse = sim.getMSE_EmbeddedAIC_Simulation();

            std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
            ofs << dopplerFrequency << "," << mse << std::endl;

            // ループ更新
            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
        }
		// ★合計時間の計測終了と表示
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        
        std::cout << "========================================" << std::endl;
        std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
        std::cout << "========================================" << std::endl;
    }
	else if (mode_select == 19)
    {
        // --- モード19: Wrapper法 SNR劣化比 (Degradation Ratio) vs Doppler ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]: ";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_WrapperAIC_SNRDegradation.csv";
        ofs.open(fileName);
        ofs << "f_dT_s,DegradationRatio" << std::endl;

        sim.setNoiseSD((double)fixedEbN0dB);

        // 合計時間の計測開始
        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            
            std::cout << "Target: f_dT_s = " << dopplerFrequency << " processing..." << std::endl;
            
            // Wrapper法でのSNR劣化比シミュレーション
            double degradation = sim.getSNRDegradation_WrapperAIC_Simulation();

            std::cout << " Result: Ratio = " << degradation << " (1.0 is Ideal)" << std::endl;
            ofs << dopplerFrequency << "," << degradation << std::endl;

            // ドップラー周波数の更新
            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }

        // 合計時間の計測終了と表示
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        
        std::cout << "========================================" << std::endl;
        std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
        std::cout << "========================================" << std::endl;
    }
	else if (mode_select == 20)
    {
        // --- モード20: Pilot Only SNR劣化比 (Degradation Ratio) vs Doppler ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]: ";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_PilotOnly_SNRDegradation.csv";
        ofs.open(fileName);
        ofs << "f_dT_s,DegradationRatio" << std::endl;

        sim.setNoiseSD((double)fixedEbN0dB);

        // 合計時間の計測開始
        auto start_total = std::chrono::high_resolution_clock::now();

        // ドップラー周波数をスイープ (0 -> step -> *2 -> *2 ...)
        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            
            std::cout << "Target: f_dT_s = " << dopplerFrequency << " processing..." << std::endl;
            
            // パイロットのみでのSNR劣化比シミュレーション
            double degradation = sim.getSNRDegradation_PilotOnly_Simulation();

            std::cout << " Result: Ratio = " << degradation << " (1.0 is Ideal)" << std::endl;
            ofs << dopplerFrequency << "," << degradation << std::endl;

            // ドップラー周波数の更新
            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
        }

        // 合計時間の計測終了と表示
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        
        std::cout << "========================================" << std::endl;
        std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
        std::cout << "========================================" << std::endl;
    }
	else if (mode_select == 21)
    {
        // --- モード21: パイロットAIC固定パス法 MSE ---
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]: ";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_PilotAICFixed_MSE.csv";
        ofs.open(fileName);

        sim.setNoiseSD((double)fixedEbN0dB);
        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            std::cout << "Target: f_dT_s = " << dopplerFrequency << " processing..." << std::endl;

            double mse = sim.getMSE_PilotAICFixedPath_Simulation();

            std::cout << " Result: MSE = " << mse << std::endl;
            ofs << dopplerFrequency << "," << mse << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        std::cout << "Total Time: " << elapsed_total.count() << "s" << std::endl;
    }
	else if (mode_select == 23)
    {
        // --- モード23: 送信信号 X の出力 ---
        // ※送信信号そのものはドップラー周波数の影響を受けませんが、
        //   シミュレーションの都合上、パラメータ設定などはそのまま通します。
        
        fileName = outputDir + timeStr + "_" + modulationName + "_TxWaveform_k10.csv";
        
        // ターゲットサブキャリア
        int target_k = 10;
        
        std::cout << "Exporting Tx Signal for k=" << target_k << " ..." << std::endl;
        sim.runExportTxWaveform(target_k, fileName);
    }
    else if (mode_select == 24)
    {
        // --- モード24: フェージングを受けた信号 HX の出力 ---
        
        // 動き（ドップラー）を設定させたい場合は入力を受け付ける
        double inputDoppler;
        std::cout << "Enter Normalized Doppler Frequency (f_d T_s): ";
        std::cin >> inputDoppler;
        
        sim.setDopplerFrequency(inputDoppler);

        fileName = outputDir + timeStr + "_" + modulationName + "_FadedWaveform_k10.csv";
        
        // ターゲットサブキャリア
        int target_k = 10;

        std::cout << "Exporting Faded Signal for k=" << target_k << " ..." << std::endl;
        sim.runExportFadedWaveform(target_k, fileName);
    }
	else if (mode_select == 25)
    {
        // --- モード25: チャネル応答 |H| の出力 ---
        double inputDoppler;
        std::cout << "Enter Normalized Doppler Frequency (f_d T_s): ";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int target_k;
        std::cout << "Enter target subcarrier index (k): ";
        std::cin >> target_k; // ターミナルからkを入力

        fileName = outputDir + timeStr + "_" + modulationName + "_ChannelMagnitude_k" + std::to_string(target_k) + ".csv";
        
        std::cout << "Exporting Channel Magnitude for k=" << target_k << " ..." << std::endl;
        sim.runExportChannelMagnitude(target_k, fileName);
    }
	else
	{
		std::cout << "Invalid mode selected." << std::endl;
		return 1; // エラー終了
	}

	ofs.close();
	
	return 0;
}