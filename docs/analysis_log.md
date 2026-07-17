# 分析日志（关键决策与迭代记录）

## 数据获取
- MAXI star_data、Swift/BAT transient monitor 直接 ASCII 下载；4 处站点标识错误修正（J1348-632、IGRJ17464-3213、4U1543-47、SWIFTJ1727.8-1613）。
- RXTE/ASM：MIT 镜像不可达 → 采用 HEASARC ASMProducts definitive_1dwell FITS。**踩坑**：dwell 文件 TIME 单位是天（TIMEREF=49353），且 SUM 与 A/B/C 分带分别在 lightcurves/ 与 colors/ 两套文件——初版误按秒换算导致 1996-2005 爆发全部丢失，经与已知爆发史核对发现并修复。
- 下载器：并发=4，JSON checkpoint 断点续跑，日志实时进度（logs/download_xray.log）。

## 爆发识别迭代（3 轮）
1. 阈值 0.02 Crab：噪声段过多（J1820 出现"发现前爆发"）。
2. +峰值>0.1 Crab、点数≥10：清除噪声但丢失 V404 2015（短）与 J1118（暗、硬态-only）。
3. per-source 覆盖（J1118 峰值阈 0.025；V404 允许 BAT 触发、min_days=8）→ 51 段全部与文献爆发史一致（目检 lc_*.png 确认）。

## 态转变提取
- 硬度阈值不手选：用全样本 log-HR 分布双峰谷值（MAXI 0.40 / ASM 0.95），文献常用值（MAXI HR~0.4-0.5）相符。
- 转变要求跨阈持久性(≥67%)，光变 bootstrap(300 次) 给时间/流量误差。
- 交叉验证：GX 339-4 decay 转变 Eddington 比 1.3–4%（文献 ~2% ✓）；J1820 decay MJD 58388（文献 58380-58393 ✓）。
- 否决表（veto_windows.csv）：J1727 发现前 MAXI 段（场源污染）、4U 1543-47 2011 假段、J1820 发现前段、GRO J1655 1996 rise（ASM 早期采样差）、H1743 2003 rise（流量高估）、J1118（永不进软态，非典型转变）、Cyg X-1/GRS 1915（持续源，另行处理）。

## CCF（H2）
- ICCF 双向插值 + 0.8 峰值质心 + FR/RSS 400 次；单元测试：注入 τ=7.0 d 恢复 6.7±0.8 d ✓。
- GX 339-4 2010-11 decay：+8.0 d，与 arXiv:2605.19473 独立发表值一致 ✓（方法验证）。
- 55200_decay 的 MC 有 -20 d 混叠次模态（采样窗有限），已在误差(+2.4/−26)中如实保留。

## φ_jet（H3）
- 喷流功率标定固有 ~1 dex 不确定 → 三标定 (0.4/2/8×10³⁶ erg/s @L_R=10³⁰) 传播为误差棒；结论基于相对排序。
- Eddington 窗 [1e-4, 0.05] 限定硬态爆发段（排除宁静态不同物理与转变污染）。
- 目检 phi_vs_ledd.png：无异常轨迹；H1743 φ 最高（其独特的 outlier track 与近 MAD 状态自洽）。

## 统计（H6 阶段）
- 贝叶斯误差回归（emcee, 32×4000, 弃烧一半），留一法稳健性；J1118 剔除前后 decay 斜率从 −1.0±1.4 变为 −1.23±0.7（该源本就该剔除：从不进入软态）。
- 全部统计输出：output/results/stats_summary.txt。

## 事件级重构（v2，应外部评审意见收窄主张）
- 新增 data/tables/spin_table.csv：CF 与反射双列自旋（逐源文献出处；GRO J1655 已知两法张力如实保留）。
- code/events.py：构建 events.csv（50 事件×15 源）与 source_table.csv（r_ISCO(BPT72)、Ω_H、η_BZ、射电响度 ξ）。ξ 用 Lr–Lx 数据库对全局硬态轨迹 (斜率 0.6) 的逐源中位残差 + bootstrap 误差（单点源误差下限 0.15 dex）；避免"不同年代/频段射电峰值混拼"。
- code/hier_model.py：两层方差分解（解析边缘化随机截距）+ WAIC 模型比较。
- 关键数字：回滞 f_source=0.04（事件层主导）；ξ 各自旋参数化 WAIC 等价；GX 339-4 固定自旋事件散射 0.7 dex。
- 逐事件射电峰值仅 GX 339-4 由机读数据获得（5 事件）；其他源留空并注明扩展路径（ThunderKAT/AMI 释放），不做记忆性文献数值填充。
