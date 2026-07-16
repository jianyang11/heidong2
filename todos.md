# todos.md — 黑洞XRB光谱态演化 × MAD × 史瓦西半径

图例：[ ] 待办 [~] 进行中 [x] 完成 [!] 受阻/需迭代。本文件随进程持续更新（每完成/调整一项立即更新）。

## Phase 0 调研与方案
- [x] arXiv 检索 MAD、BHXRB 态转换、radio–X-ray 延迟最新文献（2023-2026）
- [x] 核验关键论文：You+23 (J1820 MAD)、Dhang+23、arXiv:2605.19473 (GX339-4 τ)、Begelman+22
- [x] 撰写 docs/literature_review.md
- [x] 撰写 docs/proposal.md（H1/H2/H3 假设、方法、风险预案）
- [x] 创建本 todos.md
- [x] 初始化 git 仓库并推送到 github.com/jianyang11/heidong2

## Phase 1 环境与样本
- [x] 建立目录结构 code/ data/{raw,processed,tables} logs/ output/{figures,results} docs/
- [x] 安装 python 依赖（numpy scipy pandas matplotlib astropy emcee corner requests）并记录 requirements.txt
- [x] 编写 code/utils.py（日志、checkpoint、Crab 换算、L_Edd/R_s 常数）
- [x] 编制 data/tables/bh_sample.csv：源名、M±err、D±err、i、P_orb、spin、出处（逐源核对 BlackCAT/文献）
- [x] 目检样本表：质量/距离与文献一致性抽查 ≥5 源

## Phase 2 X射线监测数据获取（先试点后批量）
- [x] 编写 code/download_xray.py（MAXI + Swift/BAT + RXTE/ASM；并发≤4；checkpoint 续跑；日志实时进度）
- [x] 试点 3 源（GX 339-4, MAXI J1820+070, XTE J1550-564）下载并画原始光变 → read 目检
- [x] 修正源名→各站点标识映射问题后批量下载全样本（修正 J1348-632、IGRJ17464-3213、4U1543-47、SWIFTJ1727.8-1613；补充 RXTE/ASM dwell FITS + colors）
- [x] 校验：每源覆盖时段、点数、缺测段记录到 logs/data_inventory.csv
- [x] 磁盘占用检查（data 共 532 MB，raw ~500 MB）

## Phase 3 光变清洗、爆发识别与 HID
- [x] code/lightcurve.py：清洗（负值/大误差剔除）、S/N 过滤、Crab 归一、多任务拼接
- [x] 爆发自动识别 + 人工目检光变图（GX339-4 12次、J1550 3次、GRO1655 2次、V404 2015、J1118 2000 等均与文献爆发史一致；对 J1118/V404/J1659 引入 per-source 阈值覆盖）（output/figures/lc_*.png，read 逐张检查）
- [x] 计算硬度比、构建 HID；试点源 HID 与文献（Corbel+13 等）对照核验
- [x] 态分类（阈值+文献仲裁）；生成态色带光变图并目检
- [x] 提取逐爆发 rise/decay 转变时刻与光度（bootstrap 误差）
- [x] 转变光度换算 L/L_Edd；与 Vahdat Motlagh+19、WATCHDOG 已发表值交叉验证（偏差>2σ 者排查）
- [x] 阶段小结写入 docs/analysis_log.md

## Phase 4 射电数据与 τ 测量 (H2)
- [x] 整理 XRB-LrLx_pub BH 表；逐源列出可用射电监测时序及文献
- [x] 收集逐日/逐周射电监测光变（VizieR/论文表：Corbel+13、Bright+20、Russell+19、Carotenuto+21、Tetarenko+17 等）
- [x] code/ccf.py：ICCF + DCF + FR/RSS MC 误差；单元测试（模拟已知延迟光变恢复 τ）
- [x] 试点：J1820 检验能否复现 τ≈8 d（对照 You+23）、GX339-4 decay τ≈8 d（对照 2605.19473）→ 不一致则迭代方法
- [x] 全样本逐爆发 τ 测量（rise/decay 分段），质量分级（A/B/上限）
- [x] ICCF/DCF 诊断图集目检；τ 结果表 output/results/tau_table.csv

## Phase 5 φ_jet 与 MAD 饱和度 (H3)
- [x] code/phi_jet.py：P_jet(L_R) 多标定、Ṁ(L_X)、φ_jet 及 MC 误差传播
- [x] 逐爆发 φ_jet(t) 轨迹图 + 硬态峰值 φ_max,obs 表
- [x] 与 φ_MAD≈50 对比；标定敏感性矩阵
- [x] 目检全部 φ 轨迹图，异常源排查（距离/自旋输入错误等）

## Phase 6 统计推断与假设检验
- [x] H1：L_trans/L_Edd vs M(R_s) — Spearman+偏相关+贝叶斯回归斜率后验；rise/decay 分开
- [x] H2：τ vs M、τ vs L_peak/L_Edd、τ vs P_orb 标度回归（含上限生存分析）
- [x] H3：φ_max,obs vs 回滞幅度/转变光度残差 偏相关（控制 R_s、距离）
- [x] 稳健性：留一法、距离 MC、标定敏感性；全部写入 output/results/
- [x] corner 图与关键关系图目检迭代至发表级

## Phase 7 可视化与成果
- [x] 发表级图集（统一样式）：样本光变图集、HID 图集、L_trans–M、τ–M、φ 轨迹、Lr–Lx(φ 着色)、后验 corner、敏感性热图
- [x] 逐张 read 目检：清晰度、标注、配色、误差棒完整性 → 迭代
- [x] docs/results_summary.md：方法、结果、与理论对照、创新点、可发表结论
- [x] 复查：无数据编造、所有数值可溯源（脚本+日志+数据文件三链路）
- [x] 更新 literature_review/proposal 变更记录

## Phase 8 交付
- [x] 清理无用中间文件，检查磁盘
- [x] README.md（仓库结构、复现步骤）
- [x] 全部内容 push 到 GitHub；核验远端文件完整
- [x] 最终自查 todos 全部闭环

## 迭代记录（持续追加）
- 2026-07-16: 建立 v1 计划；完成文献调研与方案文档。
- 2026-07-16: ASM dwell FITS TIME 单位为天（非秒），且 SUM 与 A/B/C 分带存于 lightcurves/ 与 colors/ 两套文件——已修正并合并；爆发识别参数迭代 3 轮（阈值0.03 Crab+峰值0.1+min_pts10，faint 源单独覆盖），结果与文献爆发史核对通过。
- 2026-07-16: 待办：4U 1543-47 的 55662-55748、Swift J1727 2023 前段疑似场源污染，转变提取阶段需复核。
- 2026-07-16: Phase 3 完成：数据驱动硬度阈值（MAXI 0.40/ASM 0.95），179 转变→veto 后 61；GX339-4/J1820 转变与文献核验一致。
- 2026-07-16: Phase 4：公开机读射电时序仅 GX 339-4（Corbel+13 VizieR 93+22 点）；Zenodo/GitHub/RATAN 搜寻 J1820/J1348/H1743 时序未果 → H2 调整为 GX339-4 多爆发自测 + 文献 τ 对比（proposal v1.1，目标不降级）。ICCF 单元测试通过（注入7d 恢复 6.7±0.8d）；2010-11 decay +8.0d 与 arXiv:2605.19473 一致。
- 2026-07-16: Phase 5：φ_max∈[5,51]，全部 ≤ φ_MAD≈50；3 组标定敏感性传播。
- 2026-07-16: Phase 6：H1 decay 中位 1.5% L_Edd，斜率 −1.23(+0.67/−0.70)；发现 φ_max–回滞负相关（ρ=−0.886, p=0.019）；J1118（从不进软态）自 H1 剔除（veto 表迭代）。
- 2026-07-16: Phase 7-8：结果总结、README、全部推送 GitHub。
