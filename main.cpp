/*
 * File:   power_ber.cpp
 * Author: Ito
 *
 * Created on 2024/12/20, 19:36
 */
#include <iostream>
#include <fstream>
#include "simulator.h"
#include "parameters.h"

// パラメータ
static const int EbN0dBmin = 0;
static const int EbN0dBmax = 30;
static const int EbN0dBstp = 5;

// ファイル
std::string fileName;       // ファイル名
std::ofstream ofs;        // 出力ファイル  

// BER
double ber;

// ドップラー周波数
double dopplerFrequence;

int main()
{
	fileName = "BER.csv";
	ofs.open(fileName);

    SimulationParameters params;
    Simulator sim(params);

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "f_d?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> dopplerFrequence;
    sim.setDopplerFrequence(dopplerFrequence);

    for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
        //試行回数設定
        sim.setTrialNum(EbN0dB);
        
        // SN設定
        sim.setNoiseSD(EbN0dB);
        
        // シミュレーション
        ber = sim.getBER_EM_Simulation();
        // 標準出力
        std::cout << "-----------" << std::endl;
        std::cout << EbN0dB << "," << ber << std::endl;
        // CSV出力
		ofs << EbN0dB << "," << ber << std::endl;
    }
    ofs.close();
    
    return 0;
}