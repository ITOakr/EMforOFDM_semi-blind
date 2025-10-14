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
static const double dopplerMax = 0.1;
static const double dopplerStep = 0.01;

// ファイル
std::string fileName;       // ファイル名
std::ofstream ofs;        // 出力ファイル  

// BER
double ber;

// ドップラー周波数
double dopplerFrequency;

// 試行回数
int numberOfTrials;

std::string getModulationSchemeName(int numberOfBits) {
    switch (numberOfBits) {
        case 1: return "BPSK";
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
	std::cout << "Number of Bit? (BPSK:1, QPSK:2, 16QAM:4, 64QAM:6, 256QAM:8)" << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cin >> params.NUMBER_OF_BIT;
	params.NUMBER_OF_SYMBOLS = std::pow(2, params.NUMBER_OF_BIT);

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

    std::string modulationName = getModulationSchemeName(params.NUMBER_OF_BIT);

	if (mode_select == 1)
	{
		// --- モード1: Eb/N0スイープ ---
		double dopplerFrequency;
		std::cout << "Enter normalized Doppler f_d*T_s:" ;
		std::cin >> dopplerFrequency;

		fileName = modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_BER_vs_EbN0.csv";
		ofs.open(fileName);

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

		fileName = modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_BER_vs_Doppler.csv";
		ofs.open(fileName);

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