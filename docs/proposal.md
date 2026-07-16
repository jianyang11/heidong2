# 研究方案：黑洞X射线双星光谱态演化、MAD 与史瓦西半径的相关性——基于公开监测数据的系统检验

版本 v1.0（2026-07-16）。本方案为活文档，随实验进程迭代（变更记录见文末）。

## 1. 科学目标

**核心问题**：黑洞X射线双星(BHXRB)的光谱态演化（尤其是态转变光度与回滞）是否受磁阻吸积盘(MAD)物理调控？若是，其可观测指标是否随黑洞史瓦西半径 R_s = 2GM/c²（即黑洞质量）呈系统标度？

**三个可检验假设**：
- **H1（转变光度普适性）**：软→硬(decay)转变 Eddington 比 L_trans/L_Edd 与 R_s 无关（普适 ~2%）；硬→软(rise)转变光度弥散大且与磁通历史相关。对大样本做定量检验并给出对 R_s 依赖的置信区间。
- **H2（MAD 形成时标）**：以 radio–X-ray 互相关时间延迟 τ 为 MAD 形成/磁通输运时标指标（You+23; arXiv:2605.19473），预期 τ 在多源样本中可测，且检验 τ 与 M(R_s)、Eddington 比峰值、轨道周期（外盘尺度）的标度关系。**多源系统测量为本项目首创**。
- **H3（MAD 饱和度—态演化耦合）**：由硬态射电光度估计喷流功率 P_jet，结合 Ṁ ≈ L_X/(η c²) 估计无量纲磁通 φ_jet ∝ sqrt(P_jet/Ṁc²)/(a/2)，检验 (i) 各源硬态峰值 φ_jet 是否逼近 MAD 饱和值 φ_max≈50；(ii) φ_jet 与 decay 转变光度残差、回滞幅度（HID 面积/rise-decay 光度比）的相关性；(iii) 上述关系对 R_s 的偏相关。

**创新性**（目标：可发表成果，拟投 MNRAS/ApJ）：
1. 首次多源、逐爆发地系统测量 radio–X-ray 延迟并构建 τ–M–Ṁ 标度关系；
2. 首次把 φ_jet 作为逐爆发 MAD 饱和度代理放入 (L_trans/L_Edd, φ_jet, R_s) 联合统计框架；
3. 用贝叶斯层级回归同时处理测量误差、距离系统误差与上限值，给出对"硬态=MAD"（Begelman+22）与 GRMHD 预言（Dhang+23）的定量观测检验。

## 2. 样本与数据（全部公开数据库）

### 2.1 样本
以 WATCHDOG/BlackCAT 中具动力学质量测量的 BHXRB 为核心样本（金样本，~15-20 源），外加质量函数约束源（银样本）。必需参数：M、D、i、P_orb、自旋（如有）——全部取自发表文献并记录出处（`data/tables/bh_sample.csv`）。

候选金样本：GX 339-4, MAXI J1820+070, XTE J1550-564, GRO J1655-40, 4U 1543-47, XTE J1859+226, GRS 1915+105, H1743-322(银), V404 Cyg, MAXI J1348-630, MAXI J1535-571(银), XTE J1118+480, GRS 1124-684, A0620-00, Cyg X-1(持续源,对照), Swift J1727.8-1613 等。

### 2.2 数据源（体量估计 <5 GB，远低于 50 GB 限额）
| 数据 | 用途 | 获取方式 |
|---|---|---|
| MAXI/GSC 1-day/orbit 光变 (2-4/4-10/10-20 keV) | 2009 年后爆发的 HID、态分类、L_X | maxi.riken.jp `star_data` ASCII |
| Swift/BAT transient monitor (15-50 keV) | 硬X监测、硬度比、Compton 光度代理 | swift.gsfc.nasa.gov transients ASCII |
| RXTE/ASM (1.5-12 keV, A/B/C band) | 1996-2011 爆发 | HEASARC/MIT ASCII |
| bersavosh/XRB-LrLx_pub | 硬态 (L_R, L_X) 逐测光点 | GitHub（已获取） |
| 已发表射电监测光变（ATCA/VLA/AMI/MeerKAT: Corbel+13 GX339-4, Bright+20 J1820, Russell+19/20, Carotenuto+21 J1348 等） | τ 测量、P_jet | 论文表格/CDS VizieR 机读表 |
| BlackCAT/WATCHDOG/文献 | M, D, i, P_orb, spin | 网页表+文献 |

数据筛选纪律：只下载样本源的监测光变 ASCII 与所需 VizieR 表；不下载原始事件级数据（不需要谱拟合层级即可完成 H1-H3；若后续确需谱指数，改用 MAXI on-demand 每源少量谱或 Swift/XRT 产品生成器在线产品，避免本地堆积）。

## 3. 方法

### 3.1 光变与 HID 流水线
1. 下载、清洗（剔除负流量/大误差点，S/N>3），统一 MJD 时间轴；交叉标定 MAXI(2-20)/ASM(1.5-12)/BAT(15-50) 用 Crab 单位归一。
2. 爆发识别：基线+阈值法（如 >5σ 高于宁静基线持续 >10 天），人工目检确认（保存每源光变图）。
3. 硬度比 HR = F(4-10)/F(2-4)（MAXI）或 BAT/MAXI 比值；构建 HID；态分类采用文献标定阈值（HR 与光度联合；对照各源已发表态历史核验）。
4. 转变光度：识别 rise 硬→软与 decay 软→硬转变时刻（HR 快速演化的中点，误差由 bootstrap 光变重采样给出），流量→L_X 用 Crab 谱假设+band 修正，除以 L_Edd=1.26e38 (M/M_⊙) erg/s。
5. 回滞量化：HID 回路面积、L_rise/L_decay 比。

### 3.2 radio–X-ray 延迟测量（H2）
- 对每个有射电监测覆盖的爆发，用 ICCF（插值互相关，Gaskell & Peterson 1987; 实现参考 PyCCF 算法）+ DCF 交叉验证，测 τ（射电 vs 软X/硬X 各一组），误差用 flux randomization + random subset selection (FR/RSS, Peterson+98) MC。
- rise/decay 分段测量（arXiv:2605.19473 表明两段符号不同）。
- 质量控制：要求重叠基线 >3× 预期延迟、采样间隔 < τ/3；对不满足者只给上/下限。

### 3.3 MAD 饱和度代理 φ_jet（H3）
- P_jet 由硬态核喷流射电光度标定：P_jet = 4.79e35 (L_R/1e30 erg/s)^{12/17} W 之类的标定关系（采用 Fender+01/Heinz & Grimm 2005 校准，敏感性分析用多组标定）。
- Ṁ = L_X_bol/(η c²)，硬态 η 取 0.1×(L/L_Edd 依赖的 RIAF 修正)；进行系统误差传播。
- φ_jet = sqrt(4π P_jet/(κ ω_H² Ṁ c²))，κ≈0.05，ω_H 由发表自旋（无自旋测量则边缘化 a∈[0,1]）。
- 输出逐爆发 φ_jet(t) 轨迹与硬态峰值 φ_max,obs，与 φ_MAD≈50 对比。

### 3.4 统计推断
- 相关性：Spearman/Kendall + 偏相关（控制距离、Eddington 比）；上限值用生存分析（Kaplan-Meier, cenken）。
- 回归：贝叶斯层级线性回归（emcee 实现，含 x/y 测量误差与内禀散射，等价 linmix），模型比较用 WAIC/贝叶斯因子。
- 关键输出：∂log(L_trans/L_Edd)/∂log R_s 的后验；τ–M 标度指数后验；φ_jet 与回滞量的相关显著性。
- 稳健性：留一法（逐源剔除）、距离系统误差 MC、不同 P_jet 标定的敏感性矩阵。

### 3.5 与理论对照
把结果与 GRMHD 预言（Dhang+23 磁通-态转换; TNM11 φ_max）及 You+23 单源结论对比，讨论 MAD 形成时标的物理标度（粘滞/磁通输运时标 ∝ R_s/c 的倍数 vs 外盘时标 ∝ P_orb）。

## 4. 工程规范
- 目录：`code/`（模块化 python + 脚本）、`data/raw|processed|tables/`、`logs/`（每脚本一日志，实时 flush 进度）、`output/figures|results/`、`docs/`。
- 断点续跑：所有批量下载/拟合脚本以"每源每爆发"为原子单位落盘 checkpoint（JSON/parquet），启动时跳过已完成项。
- 并行：下载与逐源分析用进程池（受 I/O 与站点礼貌限制，下载并发≤4，计算并发=CPU 核数）。
- 无超时：长任务后台运行 + 日志轮询监控；对 >10 min 无日志输出的进程排查并改进。
- 质检闭环：每个阶段产出诊断图 → read 工具目检 → 记录问题 → 迭代；与已发表结果交叉验证（如 GX 339-4 各爆发态转变日期对照 Corbel+13、J1820 τ≈8d 对照 You+23）。
- 可视化清单（发表级，matplotlib，统一样式）：全样本光变+态色带图、HID 图集、转变光度 vs M 图、ICCF/DCF 图集、τ–M/τ–L 图、φ_jet 轨迹图、Lr–Lx 平面（按 φ 着色）、贝叶斯后验 corner 图、敏感性矩阵热图。

## 5. 里程碑
M1 文献+方案+样本表（今日）→ M2 监测数据管线跑通（2-3 源试点）→ M3 全样本 HID/转变光度 → M4 τ 测量 → M5 φ_jet 与统计推断 → M6 图表定稿+论文级结果文档 → M7 上传 GitHub。

## 6. 风险与失败预案（不降级核心目标）
- 射电监测覆盖不足致 τ 样本小 → 补充 OVRO/AMI/MeerKAT ThunderKAT 公开数据、扩展至银样本、用发表光变数字化（WebPlotDigitizer 级精度声明）；仍不足则以逐爆发上限做生存分析，保持 H2 的统计检验不降级。
- MAXI/ASM 站点不可达 → 换 HEASARC 镜像 / astroquery.heasarc。
- 态分类模糊源 → 用 BAT/MAXI 双波段硬度 + 已发表态历史仲裁；必要时引入 Swift/XRT 在线产品谱指数。
- φ_jet 标定系统误差大 → 报告多标定敏感性并以相对量（源间比较）为主要结论载体。

## 变更记录
- v1.0 初版。
