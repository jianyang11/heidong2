# 文献调研：黑洞X射线双星光谱态演化 × 磁阻吸积盘(MAD) × 史瓦西半径

更新：2026-07-16（本文件随研究进程持续更新）

## 1. 研究背景

### 1.1 黑洞X射线双星(BHXRB)的光谱态演化
- BHXRB爆发遵循经典的硬度-强度图(HID)"q"形回滞轨迹：宁静态 → 硬态(LHS)上升 → 硬中间态(HIMS) → 软中间态(SIMS) → 软态(HSS) → 沿低光度回到硬态 → 宁静态 (Homan & Belloni 2005; Remillard & McClintock 2006; Belloni 2010)。
- 关键观测事实：**软→硬回转光度**（decay 端）集中在 L_trans ≈ 0.02 L_Edd 附近（Maccarone 2003; Vahdat Motlagh et al. 2019; Tetarenko et al. 2016 WATCHDOG），而硬→软转变光度（rise 端）弥散大（回滞现象）。这一"普适" 2% L_Edd 是否真正与黑洞质量（即史瓦西半径 R_s=2GM/c²）无关，尚缺乏对大样本的系统定量检验——**这是本项目的切入点之一**。

### 1.2 磁阻吸积盘(MAD)理论与数值模拟
- MAD 概念：吸积把大尺度极向磁场拖入并在视界附近堆积，当磁压足以平衡吸积 ram 压时，吸积被磁场"阻滞"，无量纲磁通量饱和于 φ_BH ≈ 50 (高斯单位, Tchekhovskoy, Narayan & McKinney 2011; Narayan et al. 2003)。MAD 状态下 Blandford-Znajek 喷流效率可超过 100%：P_jet ≈ (κ/4π) φ² (a/2)² Ṁc²。
- **Dhang, Bai & White (2023, arXiv:2309.15926)**：首个辐射双温 GRMHD 模拟覆盖 Ṁ~1e-10→1e-2 Ṁ_Edd 的爆发上升过程，证明磁通量在态转换中与 Ṁ 同等重要；磁通饱和程度决定 rise 端行为——为"回滞由磁通量演化驱动"提供理论支撑。
- **(arXiv:2309.16092)**：长时标 GRMHD 薄盘+注入模拟再现软→软中间→硬中间→硬→宁静的完整衰减序列。
- **Begelman, Scepi & Dexter (2022, MNRAS 511, 2040)**：提出 BHXRB 整个硬态本质上就是 MAD 的观点；MAD 磁通演化自然解释回滞。
- **You et al. (2023, Science 381, 961; arXiv:2309.00200)**：对 MAXI J1820+070 的 2018 年爆发的多波段分析发现射电滞后X射线 ~8 天、光学滞后 ~17 天，解释为硬态中膨胀冕放大磁场、在射电峰值时刻形成 MAD——**首个 BHXRB 中 MAD 形成的观测证据**，也是本项目"用公开监测数据求 radio/X-ray 延迟作为 MAD 形成指标"方法论的直接依据。
- **arXiv:2605.19473 (2026, GX 339-4)**：用 ICCF 测得 GX 339-4 2010-11 爆发 rise 硬态射电超前 Compton 光度 ~3 天、decay 硬态滞后 ~8 天，用内区磁场输运解释——证明该延迟可系统测量并具有磁场输运含义。**本项目将把该单源测量推广为多源样本，并检验其与 R_s 的标度关系。**

### 1.3 喷流—吸积耦合与射电/X射线相关性
- 硬态普适相关 L_R ∝ L_X^0.6（Gallo, Fender & Pooley 2003; Corbel et al. 2013），存在"outliers/steep track" (L_R ∝ L_X^~1.4)。标准盘-喷流耦合与黑洞基本面平面（Merloni et al. 2003; Falcke et al. 2004: L_R ∝ L_X^0.6 M^0.8）表明质量（R_s）应进入喷流标度。
- 公开汇编数据库：**Bahramian & Rushton, XRB Lr-Lx database (Zenodo/GitHub bersavosh/XRB-LrLx_pub)** —— 本项目已获取。
- 喷流功率-磁通诊断：φ_jet = sqrt(P_jet/(κ (a/2)² Ṁ c²))·(4π)^{1/2}，可由射电光度→喷流功率标定（如 Willott et al. 1999 校准的 P_jet–L_R 关系）与 Eddington 标度化 X 射线光度→Ṁ 估计，进而判断源在硬态峰值处是否接近 MAD 饱和 φ_max≈50。

### 1.4 黑洞质量（史瓦西半径）测量
- 动力学质量：BlackCAT (Corral-Santana et al. 2016，持续更新) 与 WATCHDOG (Tetarenko et al. 2016) 目录；近年 Gaia 距离修正（如 MAXI J1820+070: M=8.5 M_⊙, D=3.0 kpc, Torres et al. 2020）。
- R_s = 2GM/c² = 2.95 km × (M/M_⊙)。恒星级黑洞 R_s 跨度约 12–60 km（4–21 M_⊙），配合已发表的自旋测量可进一步区分 R_s 与 ISCO 半径效应。

## 2. 研究空白（创新点定位）
1. **尚无人对多源样本系统测量"射电-X射线时间延迟"（MAD 形成时标指标）并检验其随黑洞质量/史瓦西半径、Eddington 比的标度**。理论预期：MAD 形成时标 ∝ 磁通输运时标，与内区尺度 (∝R_s) 及 Ṁ 有关。
2. **软→硬转变光度 2% L_Edd 的"普适性"未与 MAD 磁通指标联合检验**：若回滞由磁通演化驱动 (Begelman+22; Dhang+23)，则 decay 转变光度应与喷流磁通代理量相关，残差可能随 R_s 系统变化。
3. **将 φ_jet（MAD 饱和度代理）作为逐源、逐爆发的可测参数并放入 (L_trans/L_Edd, φ_jet, R_s) 三维参数空间做统计推断**，是全新的观测检验框架。

## 3. 数据可得性核验（已完成）
| 数据 | 来源 | 状态 |
|---|---|---|
| MAXI/GSC 1-day 光变(2-20 keV 分band) | maxi.riken.jp star_data | ✅ 已验证可下载(每源~几百KB) |
| Swift/BAT transient monitor (15-50 keV) | swift.gsfc.nasa.gov/results/transients | 待验证 |
| RXTE/ASM (1.5-12 keV, 1996-2011) | ADS/HEASARC xte mirror | 待验证 |
| 射电-X射线数据库 | bersavosh/XRB-LrLx_pub | ✅ 已克隆(BH表~千条) |
| 黑洞质量/距离/自旋 | BlackCAT + WATCHDOG + 文献汇编 | ✅ BlackCAT 可访问 |

## 4. 关键参考文献
- Tchekhovskoy, Narayan & McKinney 2011, MNRAS 418, L79 (MAD φ_max≈50, BZ 喷流)
- Narayan, Igumenshchev & Abramowicz 2003, PASJ 55, L69 (MAD 概念)
- You et al. 2023, Science 381, 961 (MAXI J1820+070 MAD 形成观测证据)
- Begelman, Scepi & Dexter 2022, MNRAS 511, 2040 (硬态=MAD)
- Dhang, Bai & White 2023 (arXiv:2309.15926, 辐射2T GRMHD 爆发模拟)
- arXiv:2309.16092 (GRMHD 衰减序列模拟)
- arXiv:2605.19473 (GX 339-4 radio–X-ray 时间延迟)
- Corbel et al. 2013, MNRAS 428, 2500 (GX 339-4 Lr–Lx)
- Gallo et al. 2003; 2014 (Lr–Lx 普适相关)
- Merloni et al. 2003; Falcke et al. 2004 (基本面平面, 质量标度)
- Maccarone 2003, A&A 409, 697 (2% L_Edd 转变光度)
- Vahdat Motlagh et al. 2019, MNRAS 485, 2744 (RXTE 时代转变光度统计)
- Tetarenko et al. 2016, ApJS 222, 15 (WATCHDOG)
- Corral-Santana et al. 2016, A&A 587, A61 (BlackCAT)
- Homan & Belloni 2005; Belloni 2010 (态分类)
- Russell et al. 2019 (MAXI J1820 jet ejection); Bright et al. 2020 (J1820 射电监测)
