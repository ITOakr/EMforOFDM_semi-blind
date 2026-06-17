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
	std::cout << "1: BER vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "2: BER vs Doppler sweep (fixed Eb/N0)" << std::endl;
    std::cout << "3: MSE vs Eb/N0 sweep (fixed Doppler)" << std::endl;
    std::cout << "4: Average power of the true channel response H" << std::endl;
	std::cout << "5: BER vs Eb/N0 sweep for pilot-only mode (fixed Doppler)" << std::endl;
	std::cout << "6: MSE vs Eb/N0 sweep for pilot-only mode (fixed Doppler)" << std::endl;
	std::cout << "7: BER vs Doppler sweep for pilot-only mode (fixed Eb/N0)" << std::endl;
	std::cout << "8: MSE vs Doppler sweep (fixed Eb/N0, 0 -> step -> double)" << std::endl;
	std::cout << "9: MSE vs Doppler sweep for pilot-only mode (fixed Eb/N0, 0 -> step -> double)" << std::endl;
	std::cout << "10: CSV output of channel magnitude response |H(k, l)|" << std::endl;
	std::cout << "11: MSE of estimated noise variance vs Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "12: MSE of h by initial h estimation vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "13: MSE of H by pilot estimation vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "14: Parallel MSE vs Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "15: AIC model selection accuracy vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "16: AIC path selection accuracy vs Doppler (fixed Eb/N0)" << std::endl;
    std::cout << "17: AIC path selection F-measure vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "18: Embedded AIC method MSE vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "19: Wrapper AIC method SNR degradation ratio vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "20: Pilot-only SNR degradation ratio vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "21: Pilot AIC fixed-path MSE vs Doppler (fixed Eb/N0)" << std::endl;
	std::cout << "22: Unused" << std::endl;
	std::cout << "23: Export transmit waveform X over time (k=10)" << std::endl;
    std::cout << "24: Export faded waveform HX over time (k=10)" << std::endl;
	std::cout << "25: Export channel magnitude |H| over time (k=10)" << std::endl;
	std::cout << "26: Export frequency response |H(k, l)| along k for fixed l" << std::endl;
	std::cout << "27: Export average impulse response along q for fixed l" << std::endl;
	std::cout << "28: Export estimated impulse response (l=0, Q=16) to CSV" << std::endl;
	std::cout << "29: Sweep frame length L and evaluate MSE (fixed Eb/N0 and Doppler)" << std::endl;
	std::cout << "30: Sweep CRLB-based MSE vs Eb/N0 (fixed Doppler)" << std::endl;
	std::cout << "31: Simulate impulse-response h MSE vs Eb/N0" << std::endl;
	std::cout << "32: Known-model, known-noise ML validation (H MSE vs Eb/N0)" << std::endl;
	std::cout << "33: Known-model, known-noise ML validation (H MSE at l=0 vs Eb/N0)" << std::endl;
	std::cout << "34: Instantaneous SNR (gamma) mean vs Eb/N0 at (l,k)=(0,0)" << std::endl;
	std::cout << "35: Instantaneous SNR (gamma) vs Eb/N0 at (l,k)=(0,0)" << std::endl;
	std::cout << "36: Instantaneous SNR (gamma) vs Doppler at (l,k)=(0,0)" << std::endl;
	std::cout << "37: MSE of h 16 paths by initial h estimation vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "38: MSE of Pilot Raghavendra AIC h vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "39: MSE of Pilot Raghavendra AIC H vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "40: MSE of Random Path Model AIC vs Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "41: AIC Exhaustive Search (8 paths) Accuracy vs Eb/N0 sweep" << std::endl;
	std::cout << "42: MSE of Random Path Model with Exhaustive Raghavendra AIC vs Eb/N0" << std::endl;
	std::cout << "43: AIC Exhaustive Search (8 paths) MSE vs Eb/N0 sweep (Fixed Mask)" << std::endl;
	std::cout << "44: Raghavendra AIC Exhaustive Search (8 paths) MSE vs Eb/N0 sweep (Fixed Mask)" << std::endl;
	std::cout << "45: MSE of Raghavendra AIC update vs Eb/N0 (Fixed Mask)" << std::endl;
	std::cout << "46: MSE of Random Path Model with Known Mask vs Eb/N0" << std::endl;
	std::cout << "47: AIC vs Raghavendra GAIC comparison (Single Trial)" << std::endl;
	std::cout << "48: AIC vs Raghavendra GAIC comparison (Average over Trials)" << std::endl;
	std::cout << "49: MSE of h by initial h estimation with Raghavendra GAIC vs Eb/N0 sweep (fixed Doppler)" << std::endl;
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
		std::cout << "Enter fixed Eb/N0 [dB]:";
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
		std::cout << "Enter fixed Eb/N0 [dB]:";
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
		std::cout << "Enter fixed Eb/N0 [dB]:";
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
		std::cout << "Enter fixed Eb/N0 [dB]:";
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
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
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
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int target_k;
		std::cout << "Enter target subcarrier index (k): ";
        std::cin >> target_k; // ターミナルからkを入力

        fileName = outputDir + timeStr + "_" + modulationName + "_ChannelMagnitude_k" + std::to_string(target_k) + ".csv";
        
		std::cout << "Exporting Channel Magnitude for k=" << target_k << " ..." << std::endl;
        sim.runExportChannelMagnitude(target_k, fileName);
    }
	else if (mode_select == 26)
	{
		// --- モード26: 横軸kの周波数応答 ---
		double inputDoppler;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> inputDoppler;
		sim.setDopplerFrequency(inputDoppler);

		int target_l;
		std::cout << "Enter target symbol index (l): ";
		std::cin >> target_l;

		fileName = outputDir + timeStr + "_" + modulationName + "_FreqResp_k_at_l" + std::to_string(target_l) + ".csv";
		
		std::cout << "Exporting Frequency Response vs k..." << std::endl;
		sim.saveFrequencyResponseByK(target_l, fileName);
		std::cout << "Saved to: " << fileName << std::endl;
	}
	else if (mode_select == 27)
	{
		// --- モード27: 平均インパルス応答 vs q ---
		double inputDoppler;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> inputDoppler;
		sim.setDopplerFrequency(inputDoppler);

		// 試行回数は main で最初に入力した numberOfTrials が使われます
		int target_l = 0; // 特定のシンボル時刻（通常は0でOK）

		fileName = outputDir + timeStr + "_" + modulationName + "_AverageImpulseResp_q.csv";
		
		sim.saveAverageImpulseResponseByQ(target_l, fileName);
		
		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Average impulse response saved to: " << fileName << std::endl;
		std::cout << "Check 'Power_dB' column to see the -1dB/path decay." << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}
	else if (mode_select == 28)
	{
		// --- モード28: 推定インパルス応答の出力 ---
		double inputDoppler;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> inputDoppler;
		sim.setDopplerFrequency(inputDoppler);

		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:";
		std::cin >> fixedEbN0dB;
		sim.setNoiseSD(fixedEbN0dB);

		fileName = outputDir + timeStr + "_" + modulationName + "_EbN0_" + std::to_string(fixedEbN0dB) + "_EstImpulseResponse_l0.csv";
		ofs.open(fileName);

		std::cout << "Exporting Estimated Impulse Response (l=0, Q=16)..." << std::endl;
		sim.saveEstimatedImpulseResponseToCSV(ofs, inputDoppler);

		std::cout << "Saved to: " << fileName << std::endl;
	}
	else if (mode_select == 29)
	{
		// --- モード29: フレーム長 L スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:";
		std::cin >> fixedEbN0dB;

		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		int L_start, L_end, L_step;
		std::cout << "Enter L start (min 2): "; std::cin >> L_start;
		std::cout << "Enter L end: "; std::cin >> L_end;
		std::cout << "Enter L step: "; std::cin >> L_step;

		fileName = outputDir + timeStr + "_" + modulationName + "_EbN0_" + std::to_string(fixedEbN0dB) + "_fdTs_" + std::to_string(dopplerFrequency) + "_MSE_vs_L.csv";
		ofs.open(fileName);
		ofs << "L,MSE" << std::endl;

		std::cout << "Starting L sweep simulation..." << std::endl;
		auto start_total = std::chrono::high_resolution_clock::now();

		for (int L = L_start; L <= L_end; L += L_step) {
			params.L_ = L;
			
			// Lが変更されるたびに、行列のサイズをリサイズするためSimulatorを再生成します
			Simulator sim_L(params);
			sim_L.setTrialNum(numberOfTrials);
			sim_L.setDopplerFrequency(dopplerFrequency);
			sim_L.setNoiseSD(fixedEbN0dB);

			// モード14（並列化版MSE）のロジックを使用してMSEを計算
			double mse_val = sim_L.getMSE_simulation_parallel();

			std::cout << "L = " << L << ", MSE = " << mse_val << std::endl;
			ofs << L << "," << mse_val << std::endl;
		}

		auto end_total = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_total = end_total - start_total;
		std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
		std::cout << "Results saved to: " << fileName << std::endl;
	}
	else if (mode_select == 30)
	{
		// --- モード30: CRLB 理論線出力 ---
		fileName = outputDir + timeStr + "_" + modulationName + "_CRLB_vs_EbN0.csv";
		ofs.open(fileName);
        
		ofs << "EbN0dB,CRLB_MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			// シミュレータにノイズSDを設定 (ここで noiseSD_ が計算・保持される)
			sim.setNoiseSD(EbN0dB);
			
			// 保持されている noiseSD_ を使って理論値を計算
			double crlb_mse = sim.getTheoreticalCRLB_H_MSE_FinalForm();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", CRLB MSE = " << crlb_mse << std::endl;
			ofs << EbN0dB << "," << crlb_mse << std::endl;
		}
	}
	else if (mode_select == 31)
	{
		// --- モード31: インパルス応答 (h) のMSEシミュレーション ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		// ファイル名の設定
		fileName = outputDir + timeStr + "_" + modulationName + "_fdTs_" + std::to_string(dopplerFrequency) + "_ImpulseResponse_MSE_vs_EbN0.csv";
		ofs.open(fileName);
        
		// CSVのヘッダーを出力
		ofs << "EbN0dB,Impulse_MSE" << std::endl;

		// Eb/N0をスイープさせてシミュレーションを実行
		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			mse = sim.getImpulseResponseMSE_simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", Impulse MSE = " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}
	}
	else if (mode_select == 32)
	{
		// --- モード32: 真のパスモデルと既知雑音分散での ML 推定の検証 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_KnownModelKnownNoise_HMSE_vs_EbN0.csv";
		ofs.open(fileName);
		ofs << "EbN0dB,H_MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double H_mse = sim.get_H_MSE_with_known_model_and_noise_during_pilot();

			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
			ofs << EbN0dB << "," << H_mse << std::endl;
		}
	}
	else if (mode_select == 33)
	{
		// --- モード33: 真のパスモデルと既知雑音分散での ML 推定の検証 先頭シンボル分 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_KnownModelKnownNoise_HMSE_vs_EbN0.csv";
		ofs.open(fileName);
		ofs << "EbN0dB,H_MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double H_mse = sim.get_H_MSE_with_known_model_and_noise_at_l0();

			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
			ofs << EbN0dB << "," << H_mse << std::endl;
		}
	}
	else if (mode_select == 34)
	{
		// --- モード34: 瞬時信号対雑音電力比 γ(ΔH) の平均を Eb/N0 スイープで出力 ---
		double dopplerFrequency;
		double gamma_mean;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Gamma_vs_EbN0_MODE34.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			gamma_mean = sim.getInstantaneousSNR_Mode34_simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", Mean Gamma(l=0,k=0) = " << gamma_mean << std::endl;
			ofs << EbN0dB << "," << gamma_mean << std::endl;
		}
	}
	else if (mode_select == 35)
	{
		// --- モード35: 瞬時信号対雑音電力比 γ(ΔH) の平均を Eb/N0 スイープで出力 ---
		double dopplerFrequency;
		double gamma_mean;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Gamma_vs_EbN0_MODE35.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			gamma_mean = sim.getInstantaneousSNR_Mode35_simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", Mean Gamma(l=0,k=0) = " << gamma_mean << std::endl;
			ofs << EbN0dB << "," << gamma_mean << std::endl;
		}
	}
	else if (mode_select == 36)
	{
		// --- モード36: 瞬時信号対雑音電力比 γ(ΔH) の平均を Eb/N0 スイープで出力 ---
		double dopplerFrequency;
		double gamma_mean;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Gamma_vs_EbN0_MODE36.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			gamma_mean = sim.getInstantaneousSNR_Mode36_simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", Mean Gamma(l=0,k=0) = " << gamma_mean << std::endl;
			ofs << EbN0dB << "," << gamma_mean << std::endl;
		}
	}
	else if (mode_select == 37)
	{
		// --- モード37: 16パスあると仮定してインパルス応答を推定 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_16Paths_HMSE_vs_EbN0.csv";
		ofs.open(fileName);
		ofs << "EbN0dB,H_MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double H_mse = sim.get_H_MSE_with_16paths_during_pilot();

			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
			ofs << EbN0dB << "," << H_mse << std::endl;
		}
	}
	else if (mode_select == 38)
	{
		// --- モード38: Raghavendraの提案するAICによるモデル選択 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RaghavendraAIC_MSE_vs_EbN0.csv";
		ofs.open(fileName);
		ofs << "EbN0dB,H_MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double H_mse = sim.get_H_MSE_with_RaghavendraAIC_during_pilot();

			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
			ofs << EbN0dB << "," << H_mse << std::endl;
		}
	}
	else if (mode_select == 39)
	{
		// --- モード39: 16パス推定時の平均インパルス応答電力 (l=0) ---
		double inputDoppler;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> inputDoppler;
		sim.setDopplerFrequency(inputDoppler);

		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:";
		std::cin >> fixedEbN0dB;
		sim.setNoiseSD(fixedEbN0dB);

		fileName = outputDir + timeStr + "_" + modulationName + "_EbN0_" + std::to_string(fixedEbN0dB) + "_16Paths_AverageImpulseResp_q.csv";
		
		std::cout << "Exporting Average Estimated Impulse Response Power (16 paths, l=0)..." << std::endl;
		sim.saveAverageEstimatedImpulseResponseByQ_16paths(fileName);
		
		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 40)
	{
		// --- モード40: ランダムパスモデルによる平均MSE (Mode 12ベース) ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE40.csv";
		ofs.open(fileName);

		std::cout << "Starting Random Path MSE Simulation (Mode 40)..." << std::endl;
		std::cout << "Eb/N0 [dB], MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double mse = sim.getMSE_RandomPath_Mode12_Simulation();

			std::cout << EbN0dB << ", " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 41)
	{
		// --- モード41: AIC 8パス総当たり正答率の検証 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_ExhaustiveAIC8_Accuracy_EbN0.csv";
		ofs.open(fileName);

		std::cout << "Starting Exhaustive AIC (8 paths) Accuracy Simulation (Mode 41)..." << std::endl;
		std::cout << "Eb/N0 [dB], Accuracy" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double accuracy = sim.getExhaustiveAICAccuracy_8paths_Simulation();

			std::cout << EbN0dB << ", " << accuracy << std::endl;
			ofs << EbN0dB << "," << accuracy << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 42)
	{
		// --- モード42: ランダムパスモデルによる平均MSE (Raghavendra AIC 全探索) ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE42_Raghavendra.csv";
		ofs.open(fileName);

		std::cout << "Starting Random Path MSE Simulation with Raghavendra AIC (Mode 42)..." << std::endl;
		std::cout << "Eb/N0 [dB], MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double mse = sim.getMSE_RandomPath_RaghavendraAIC_Simulation();

			std::cout << EbN0dB << ", " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 43)
	{
		// --- モード43: AIC 8パス総当たりによる固定マスクMSEの検証 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_ExhaustiveAIC8_FixedMask_MSE_EbN0.csv";
		ofs.open(fileName);

		std::cout << "Starting Exhaustive AIC (8 paths) Fixed Mask MSE Simulation (Mode 43)..." << std::endl;
		std::cout << "Eb/N0 [dB], MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double mse = sim.getMSE_ExhaustiveAIC_8paths_fixedMask_Simulation();

			std::cout << EbN0dB << ", " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 44)
	{
		// --- モード44: Raghavendra AIC 8パス総当たりによる固定マスクMSEの検証 ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_ExhaustiveRaghavendraAIC8_FixedMask_MSE_EbN0.csv";
		ofs.open(fileName);

		std::cout << "Starting Exhaustive Raghavendra AIC (8 paths) Fixed Mask MSE Simulation (Mode 44)..." << std::endl;
		std::cout << "Eb/N0 [dB], MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double mse = sim.getMSE_ExhaustiveRaghavendraAIC_8paths_fixedMask_Simulation();

			std::cout << EbN0dB << ", " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 45)
	{
		// --- モード45: ランダムパスモデルによる平均MSE (Raghavendra AIC 改訂版) ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE45_Raghavendra.csv";
		ofs.open(fileName);

		std::cout << "Starting Random Path MSE Simulation with Raghavendra AIC (Mode 45)..." << std::endl;
		std::cout << "Eb/N0 [dB], MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double mse = sim.getMSE_RandomPath_RaghavendraAIC_Simulation2();

			std::cout << EbN0dB << ", " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 46)
	{
		// --- モード46: ランダムパスモデルによる平均MSE (真のパスマスク既知) ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE46_KnownMask.csv";
		ofs.open(fileName);

		std::cout << "Starting Random Path MSE Simulation with Known Mask (Mode 46)..." << std::endl;
		std::cout << "Eb/N0 [dB], MSE" << std::endl;

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);

			double mse = sim.getMSE_RandomPath_KnownMask_Simulation();

			std::cout << EbN0dB << ", " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 47)
	{
		// --- モード47: AIC vs Raghavendra GAIC (単一試行) ---
		double dopplerFrequency, fixedEbN0dB;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;
		std::cout << "Enter fixed Eb/N0 [dB] for this trial:";
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_EbN0=" + std::to_string(fixedEbN0dB) + "_AIC_vs_GAIC_SingleTrial_MODE47.csv";
		ofs.open(fileName);

		std::cout << "Running AIC vs Raghavendra GAIC Single Trial Simulation (Mode 47)..." << std::endl;
		std::cout << "Q, AIC, Raghavendra_GAIC" << std::endl;
		ofs << "Q,AIC,Raghavendra_GAIC" << std::endl;

		auto [aic, gaic] = sim.getAICvsQ_SingleTrial_Simulation(dopplerFrequency, fixedEbN0dB);

		for (size_t q = 0; q < aic.size(); ++q) {
			std::cout << (q + 1) << ", " << aic[q] << ", " << gaic[q] << std::endl;
			ofs << (q + 1) << "," << aic[q] << "," << gaic[q] << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 48)
	{
		// --- モード48: AIC vs Raghavendra GAIC (複数回平均) ---
		double dopplerFrequency, fixedEbN0dB;
		std::cout << "Enter normalized Doppler frequency (f_d T_s):";
		std::cin >> dopplerFrequency;
		std::cout << "Enter fixed Eb/N0 [dB] for average simulation:";
		std::cin >> fixedEbN0dB;

		fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_EbN0=" + std::to_string(fixedEbN0dB) + "_AIC_vs_GAIC_Average_MODE48.csv";
		ofs.open(fileName);

		std::cout << "Running AIC vs Raghavendra GAIC Average Simulation (" << numberOfTrials << " trials) (Mode 48)..." << std::endl;
		std::cout << "Q, Average_AIC, Average_Raghavendra_GAIC" << std::endl;
		ofs << "Q,Average_AIC,Average_Raghavendra_GAIC" << std::endl;

		auto [avg_aic, avg_gaic] = sim.getAICvsQ_Average_Simulation(dopplerFrequency, fixedEbN0dB);

		for (size_t q = 0; q < avg_aic.size(); ++q) {
			std::cout << (q + 1) << ", " << avg_aic[q] << ", " << avg_gaic[q] << std::endl;
			ofs << (q + 1) << "," << avg_aic[q] << "," << avg_gaic[q] << std::endl;
		}

		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
	}

	else if (mode_select == 49)
	{
		// --- モード49: Eb/N0スイープ, Raghavendra GAICを用いたh_est_MSEの確認 ---
		// 電力ソート法
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_h_est_MODE49_RaghavendraGAIC.csv";
		ofs.open(fileName);

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			mse = sim.get_h_MSE_Simulation_during_pilot_RaghavendraGAIC();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
			ofs << EbN0dB << "," << mse << std::endl;
		}
	}

	else
	{
		std::cout << "Invalid mode selected." << std::endl;
		return 1; // エラー終了
	}

	ofs.close();
	
	return 0;
}