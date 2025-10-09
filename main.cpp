/*
 * File:   power_ber.cpp
 * Author: Ito
 *
 * Created on 2024/12/20, 19:36
 */
#include <iostream>
#include <fstream>
#include <string>
#include "simulator.h"
#include "parameters.h"

// パラメータ
static const int EbN0dBmin = 0;
static const int EbN0dBmax = 30;
static const int EbN0dBstp = 5;

static const int dopplerMin = 0;
static const int dopplerMax = 25000;
static const int dopplerStep = 5000;

// ファイル
std::string fileName;       // ファイル名
std::ofstream ofs;        // 出力ファイル  

// BER
double ber;

// ドップラー周波数
double dopplerFrequency;

// 試行回数
int numberOfTrials;

int main()
{
	SimulationParameters params;
	Simulator sim(params);

	int mode_select = 0;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << "Select Simulation Mode" << std::endl;
	std::cout << "1: Eb/N0 sweep (fixed Doppler)" << std::endl;
	std::cout << "2: Doppler sweep (fixed Eb/N0)" << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cin >> mode_select;

    std::cout << "Enter number of trials:" ;
	std::cin >> numberOfTrials;
	sim.setTrialNum(numberOfTrials);

	if (mode_select == 1)
	{
		// --- モード1: Eb/N0スイープ ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" << std::endl;
		std::cin >> dopplerFrequency;

		fileName = "f_d*T_s =" + std::to_string(dopplerFrequency) + "_BER_vs_EbN0.csv";
		ofs.open(fileName);
		ofs << "EbN0[dB],BER" << std::endl; // CSVヘッダを書き込み

		for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
			sim.setDopplerFrequency(dopplerFrequency);
			sim.setNoiseSD(EbN0dB);
			
			ber = sim.getBER_EM_Simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "EbN0dB = " << EbN0dB << ", BER = " << ber << std::endl;
			ofs << EbN0dB << "," << ber << std::endl;
		}
	}
	else if (mode_select == 2)
	{
		// --- モード2: ドップラー周波数スイープ ---
		int fixedEbN0dB;
		std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
		std::cin >> fixedEbN0dB;

		fileName = "EbN0_" + std::to_string(fixedEbN0dB) + "_BER_vs_Doppler.csv";
		ofs.open(fileName);
		ofs << "Doppler_Frequency[Hz],BER" << std::endl; // CSVヘッダを書き込み

		sim.setNoiseSD(fixedEbN0dB);

		for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
			sim.setDopplerFrequency(dopplerFrequency);

			ber = sim.getBER_EM_Simulation();
			
			std::cout << "-----------" << std::endl;
			std::cout << "f_d = " << dopplerFrequency << ", BER = " << ber << std::endl;
			ofs << dopplerFrequency << "," << ber << std::endl;
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