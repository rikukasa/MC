import numpy as np
import matplotlib.pyplot as plt

# --- データファイルの読み込み ---
# C++が出力したモンテカルロ計算結果を読み込む
mc_x, mc_flux, mc_stdev = np.loadtxt('mc_results.txt', unpack=True)

# C++が出力した理論解を読み込む
an_x, an_y = np.loadtxt('analytical_results.txt', unpack=True)


# --- グラフの描画 ---
print("Plotting results using matplotlib...")

# グラフの準備
plt.figure(figsize=(10, 7))

# モンテカルロ計算結果をエラーバー付きでプロット
plt.errorbar(mc_x, mc_flux, yerr=mc_stdev,
             fmt="o", markersize=2, capsize=3,
             label="Monte Carlo Simulation")

# 理論解を線でプロット
plt.plot(an_x, an_y, linewidth=2, color='orangered',
         label="Diffusion Approximation (Analytical)")

# --- グラフの装飾 ---
# y軸を対数スケールに設定
plt.yscale("log")

# 軸の範囲を設定
plt.xlim(0.0, 10.0)
plt.ylim(0.0001, 100.0)

# ラベル、タイトル、凡例、グリッドを追加
plt.xlabel("Radius", fontsize=12)
plt.ylabel("Flux", fontsize=12)
plt.title("Monte Carlo Particle Flux Simulation", fontsize=14)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# グラフを表示
plt.show()

print("Visualization finished.")
