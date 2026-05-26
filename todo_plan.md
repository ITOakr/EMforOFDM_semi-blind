# 実装計画詳細: ランダムパスモデルによるMSEシミュレーション (Mode 40)

## 概要
各試行ごとにチャネルのパス有無をランダム（50%）に決定し、それらあらゆるモデルを平均したMSEを算出するMode 40を実装します。既存のクラス設計（カプセル化）を尊重しつつ、数学的正確さ（SNRの一定性）を担保する実装を行います。

## ステップ 1: [x] `Channel.h` への生成メソッドの追加
`Channel` クラスに以下の公開メソッドを追加しました。

*   **メソッド名**: `void generateRandomPathFrequencyResponse(double fd_Ts, uniform_int_distribution<>& dist)`
*   **実装内容**:
    ```cpp
    void generateRandomPathFrequencyResponse(double fd_Ts, uniform_int_distribution<>& dist) {
        double totalPower = 0.0;
        // 1. パスの有無をランダムに決定し、未正規化の電力を計算
        do {
            totalPower = 0.0;
            for (int q = 0; q < params_.Q_; q++) {
                if (dist() == 1) { // 50%の確率
                    xi_(q) = std::pow(10.0, -0.1 * q);
                    totalPower += xi_(q);
                } else {
                    xi_(q) = 0.0;
                }
            }
        } while (totalPower == 0); // 万が一全パス0になった場合は再生成

        // 2. 合計電力を1.0に正規化
        for (int q = 0; q < params_.Q_; q++) {
            xi_(q) /= totalPower;
        }

        // 3. 既存の生成ロジックを実行
        fd_Ts_ = fd_Ts;
        seth_(); // 内部で xi_ を使用してインパルス応答 h_ を生成
        for (int l = 0; l < params_.L_; l++) {
            H_.row(l) = (W_ * h_.row(l).transpose()).transpose(); // 周波数応答 H_ を生成
        }
    }
    ```

## ステップ 2: [x] `simulator.h` へのシミュレーションメソッドの実装
`Simulator` クラスに、試行を制御するメソッドを追加しました。

*   **メソッド名**: `double getMSE_RandomPath_Mode12_Simulation()`
*   **実装内容**:
    ```cpp
    double getMSE_RandomPath_Mode12_Simulation() {
        double totalSquaredError = 0.0;
        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed); // 0 or 1 を生成

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
            transceiver_.setX_();
            // ステップ1で作ったランダム生成を呼び出す
            channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
            transceiver_.setY_(channel_.getH(), noiseSD_);

            // 推定（Mode 12ベース: AICによる初期推定）
            transceiver_.est_H_by_initial_h();

            // MSEの累積（パイロット区間のMSEを使用）
            totalSquaredError += transceiver_.getMSE_during_pilot();
        }
        // 平均化
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }
    ```

## ステップ 3: [x] `main.cpp` への組み込み
CLIの条件分岐を追加し、結果をCSV保存するようにしました。

*   **追加箇所**: `main.cpp` 内の `else if (mode_select == 40)`
*   **処理の流れ**:
    1.  ユーザーから `f_d*T_s` を入力。
    2.  CSVファイル名を決定 (`..._RandomPath_MSE_MODE40.csv`)。
    3.  `EbN0dB` スイープループ:
        *   `sim.setDopplerFrequency(dopplerFrequency);`
        *   `sim.setNoiseSD(EbN0dB);`
        *   `mse = sim.getMSE_RandomPath_Mode12_Simulation();`
        *   `ofs << EbN0dB << "," << mse << std::endl;`

## ステップ 4: [x] 最終確認
*   **SNRの安定性**: パス数が変わっても `totalPower = 1.0` に保たれているため、正しい Eb/N0 条件下での比較になっていることを確認済み。
*   **副作用の排除**: `generateRandomPathFrequencyResponse` は `pathMask` メンバ変数を書き換えない（`xi_` を直接操作する）ため、既存モードの固定パス構成を壊さないことを確認済み。
*   **統計的妥当性**: `NUMBER_OF_TRIAL` を十分に大きく取ることで、様々なパスモデルの期待値としてのMSEが得られる設計。
