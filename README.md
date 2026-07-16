# heidong2 — 黑洞X射线双星光谱态演化 × MAD × 史瓦西半径

用公开数据库（MAXI/GSC, Swift/BAT, RXTE/ASM, VizieR, Bahramian & Rushton Lr–Lx 数据库）系统检验黑洞X射线双星态转变、磁阻吸积盘(MAD)指标与史瓦西半径 R_s=2GM/c² 的相关性。

## 结构
```
docs/        文献调研(literature_review.md)、方案(proposal.md)、分析日志、结果总结(results_summary.md)
todos.md     全程任务清单（随进程持续更新）
code/        分析流水线（见下）
data/raw     监测光变原始文件（MAXI/BAT ASCII、ASM FITS、VizieR 射电表）
data/processed  清洗光变 lc_*.csv、爆发目录
data/tables  黑洞样本参数表(含出处)、否决窗口表
logs/        运行日志 + 下载 checkpoint + 数据清单
output/figures  诊断与发表级图
output/results  transitions.csv, tau_table.csv, phi_jet_*.csv, hysteresis.csv, stats_summary.txt
```

## 复现步骤
```bash
pip install -r requirements.txt
cd code
python3 download_xray.py      # 下载监测数据（checkpoint 续跑）
python3 lightcurve.py         # 清洗+爆发识别
python3 transitions.py        # 态转变与转变光度
python3 ccf.py                # ICCF 单元测试
python3 radio_lags.py         # GX 339-4 radio–X-ray 延迟
python3 phi_jet.py            # MAD 磁通饱和度代理
python3 stats.py              # H1/H2/H3 统计检验与图
```

## 主要结果（详见 docs/results_summary.md）
- decay 态转变 Eddington 比中位 1.5%，随黑洞质量的斜率 −1.23(+0.67/−0.70)（边缘负相关）。
- GX 339-4 多爆发 radio–X-ray 延迟：decay 段稳定 ~+8 d，与 MAXI J1820+070 发表值一致 → MAD 形成时标不随 R_s 标度。
- 硬态峰值 φ_jet ≤ φ_MAD≈50（H1743-322 触及饱和）；φ_max 与回滞幅度显著负相关（ρ=−0.89, p=0.019）。
