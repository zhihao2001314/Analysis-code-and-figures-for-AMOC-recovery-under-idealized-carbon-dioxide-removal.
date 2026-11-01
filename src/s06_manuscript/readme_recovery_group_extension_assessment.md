# CAS与UK扩展恢复分组的依据

当前为分组评估，未替换正式membership或重算现有复合图。计算入口 `n08_recovery_group_extension_assessment.py`；结果位于 `output/s06_manuscript/n08_recovery_group_extension_assessment/`。

主判据延续原40°N、500m及以下AMOC：未平滑强迫年值，末期271–280年均值/初始1–10年均值，R≥0.8基本恢复，否则Non-recovery。CAS101.86%、UK64.04%，分别进入恢复与未恢复候选扩展组。26.5°N为89.84%及75.03%，分类一致。末期改20或30年、固定初始十年参考，两者分类也保持。没有使用旧piControl异常列，没有新增漂移订正。

扩展候选恢复组为GFDL、Nor、CAS；未恢复组为CESM2、CNRM、UK。这解决恢复类别问题，不等于证明六模式初始区域模拟同样合理。用初始March N²与现代IAP2005–2024年March平均，在所有模式一致的15–995m、每10m共99个深度上计算等深度RMSD：LI中CAS0.2671、UK0.4455；NS中CAS0.2528、UK0.4462，单位10⁻⁵s⁻²。CAS剖面偏差较小，但现代观测与理想初态时期不同，这只是描述性比较，不能作为新建且已验证的综合技巧阈值。

扩展组宜直接命名Substantial recovery / Non-recovery，不将六者统一称为Better-simulated。若论文需要保留原始四模式区域质量筛选，可将3+3作为扩展敏感性分析；该选择由论文最终研究设计决定。下游组均值和显著性需按新固定成员重新计算，不能沿用2+2数值。

维度：AMOC CSV为Model×Simulation year，取明确命名的Forced AMOC strength列；N² CSV为model×depth_m。没有读取原始NetCDF。源文件SHA256及范围、公式与局限保存metadata.json。GFDL起始源6年约定不在本次被重新解释或更改；只复用现行AMOC定义。
