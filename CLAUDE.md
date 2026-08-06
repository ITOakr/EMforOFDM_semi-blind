# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

（ユーザーは日本語のみを使用するため、このファイルおよびClaude Codeからの応答は日本語で記述しています。）

## これは何か

OFDMのチャネル推定を研究するためのC++/Eigen製シミュレータです。レイリーフェージングチャネル（Jakesモデル）上で、チャネル推定・等化アルゴリズム（LS推定、AIC/GAICによるモデル次数選択、EMベースの反復トラッキング、判定帰還トラッキングなど）を数値的に評価し、BER/MSEの結果をCSVに出力してオフライン分析（表計算ソフトやプロット用ツール等）に使います。自動テストは存在せず、正しさはシミュレーション結果（BER/MSE曲線）を理論値と比較することで判断します。推定アルゴリズムの数式を変更した場合は、「テストが通る」ことではなく、既知のモード（例: AICはMode 12、Raghavendra GAICはMode 38）と比較して数値的な妥当性を確認してください。

## ビルド

Makefile/CMakeは無く、各エントリポイントを g++ + Eigen で直接ビルドします（`.vscode/tasks.json` と同じ設定）。

```powershell
g++ -fopenmp -g -I "C:\Users\Akira Ito\eigen-3.4.0" main_em.cpp -o main_em.exe
g++ -fopenmp -g -I "C:\Users\Akira Ito\eigen-3.4.0" main_pilot.cpp -o main_pilot.exe
g++ -fopenmp -g -I "C:\Users\Akira Ito\eigen-3.4.0" main_misc.cpp -o main_misc.exe
```

`-fopenmp` が必要な理由は、receiver/simulator の一部の処理でモンテカルロ試行の並列化にOpenMPを使っているためです。Eigenのインクルードパスは環境依存なので、`eigen-3.4.0` を展開した実際の場所に合わせて `-I` を変更してください。

ビルド済みの `.exe` はソースと一緒にgit管理されています（リビルド後に `git status` で `.exe` が変更扱いになるのはこのため）。mainをリビルドすると対応する `.exe` が差分として現れるのは想定内の挙動であり、問題ではありません。

## 実行方法

各実行ファイルは対話式CLIです。標準入力で変調多値数（bit数）、数値で指定するシミュレーション「モード」（後述）、モードごとのパラメータ（Eb/N0、正規化ドップラー `f_d*T_s`、試行回数など）を順に尋ねます。結果はEb/N0やドップラーなどでスイープしながらCSVに書き出されます。

**出力先ディレクトリは各 `main_*.cpp` の先頭でハードコードされています**（`outputDir`、現状は `C:/Users/Akira Ito/code2025/results/...` 配下のパス）。別の環境で実行する前にこの定数を更新してください（他の `main_*.cpp` と揃えたい場合は同様に修正）。更新しないとファイルオープンに失敗して実行できません。

## アーキテクチャ

### 3つの並列エントリポイントと共通のシミュレーションコア

単一の `main.cpp` は存在しません。もともと1つだったmainを分割した経緯があります（ブランチ `refactor-code-split` / コミット "Update main.cpp" 参照）。現在は3つのCLIフロントエンドがあり、それぞれ異なる `Simulator*` サブクラスに対する薄いモード分岐（`switch`/`if` の連鎖）になっています。

- `main_em.cpp` → `SimulatorEM`（`simulator_em.h`）— EMベースの反復チャネルトラッキング、判定帰還（Decision-Feedback）トラッキング、フルパス/埋め込みAICによるBER・MSEモード群。
- `main_pilot.cpp` → `SimulatorPilot`（`simulator_pilot.h`）— パイロット/LSベースの推定、AICおよびRaghavendra-GAICによるパス次数選択の精度・MSEモード群、全探索AICによる検証モード群。
- `main_misc.cpp` → `SimulatorMisc`（`simulator_misc.h`）— 診断・エクスポート系のモード（波形/チャネルのトレースをCSV出力、CRLB、瞬時SNRなど）。BER/MSEスイープ系列には属さない。

3つの `Simulator*` クラスはいずれも `SimulatorBase`（`simulator_base.h`）を継承しており、共通の構成要素と、全モードで使われるEb/N0→ノイズ標準偏差変換・ドップラー・試行回数のセッターを持ちます。

```
SimulatorBase
 ├─ SimulationParameters params_      (parameters.h)
 ├─ Channel   channel_                (Channel.h)
 ├─ Transmitter transmitter_          (transmitter.h)
 └─ Receiver  receiver_               (receiver.h)
```

典型的な1試行（`Simulator*::get*_Simulation()` 系メソッド内で `NUMBER_OF_TRIAL` 回繰り返される）の流れ:

1. `transmitter_.setX_()` — ランダムなQAM/QPSKデータを生成し、パイロットシンボルを埋め込む（先頭の時間方向プリアンブル行、および `enableDataPilots` 有効時はデータ部の固定パイロットサブキャリア k=5,19,32,46）。
2. `channel_.generateFrequencyResponse(fd_Ts_)` — `params_.pathMask`（どのタップが有効か）に基づき、時間相関行列の固有値分解を用いて相関のある（Jakesモデルの）フェージングチャネルを1回分生成する。
3. `receiver_.setY_(...)` — チャネルとAWGNを適用し、受信信号を得る。
4. `receiver_.est_H_by_*()` / `equalizeWith*()` / `track_H_by_*()` のいずれか — チャネル推定・等化（EM/判定帰還モードの場合はシンボル単位のトラッキング）。
5. 試行ごとにBER（`getBitErrorCount`）またはMSE（`getMSE*`）を積算し、シミュレータのメソッドはEb/N0またはドップラー1点分の平均値を返す。

`main_*.cpp` 側ではこれをEb/N0・ドップラー・フレーム長などでスイープし、`EbN0dB,value`（またはそれに準じた形式）の行をCSVへ流し込みます。

### `receiver.h` がアルゴリズムの中核

群を抜いて最大のファイル（約1700行）です。実際の推定・等化理論はここに実装されています: フルQ次元のLS推定、標準AICと「Raghavendra GAIC」によるモデル次数/パス選択（wrapper、embedded、power-sort、簡易版、全探索/全パターン探索など複数のバリエーション）、EMのE-step/M-stepによる反復精緻化（`Estep`/`Mstep`/`track_H_by_EM`）、判定帰還によるチャネルトラッキング、理論的な参照値として使う既知モデル・既知ノイズのML推定。新しい推定モードを追加する場合、ほぼ必ずこのファイルを拡張することになります。一方 `Simulator*` クラス側は、試行ループ・`receiver_` の呼び出し・指標の平均化といったオーケストレーションが中心で、推定ロジック自体はほとんど持ちません。

### パラメータ関連

- `parameters.h` — `SimulationParameters`: OFDM/システムの定数（サブキャリア数K_、フレーム内シンボル数L_、最大タップ数Q_=16、FFTサイズ、パイロット配置）に加え、可変な `pathMask`（Q_本のタップのうちどれが実際に有効か＝推定器が復元しようとしている「真値」）とRNGの `seed` を保持。`Channel` と `Receiver` で共有されるDFT行列生成の静的メソッド `generateW()` もここにある。
- `estimator_parameters.h` — `EstimatorParameters::Q_est`: 推定器側が仮定するタップ数。チャネルの真の `Q_` とはあえて分離されており、モデル次数がミスマッチした状況をシミュレートできるようになっている。
- `random_collection.h` — `<random>` を薄くラップしたクラス群（`uniform_int_distribution`、`cnormal_distribution` など）。明示的なシード指定に対応しており、モンテカルロ試行の再現性を担保するために全体で使われている。

### モード番号について

各 `main_*.cpp` は `mode_select` を標準入力から読む前に、番号付きメニュー（`std::cout << "N: 説明"`）を表示します。モード番号は機能が追加された順に歴史的に割り振られており、グルーピングされていたり再利用されたりはしません。また3つのファイル間で番号が重複しないよう緩く配慮されてはいますが、グローバルに一元管理されているわけではありません。新しいモードを追加する際は、対象の `main_*.cpp` 内で `std::cout << "<N>:` を軽くgrepして既存メニューを確認し、既存の番号を振り直すのではなく、そのファイル内で未使用の次の番号を使ってください。

### 計画ドキュメント

`round_robin_plan.md` と `todo_plan.md` は、特定の過去/進行中のモードについての日本語の実装計画（チェックリスト付き）です。あるモードが「なぜ」そのような構造になっているかを知る設計上の参考資料として有用ですが、あくまである時点での計画であり、常に最新化されているドキュメントではありません。参照する際は実際のコードと突き合わせて確認してください。

## 開発上の慣習

- コードコメントおよびコミットメッセージは日本語で書かれています。既存ファイルを編集する際やコミットする際はこのスタイルに合わせてください。
- 新しいシミュレーション機能は通常、既存モードの挙動を変更するのではなく、新しい番号のモード＋新しい `get*_Simulation()`/`est_H_by_*()` メソッドとして追加します。既存モードは過去の結果を再現するための参照点として扱われています。
