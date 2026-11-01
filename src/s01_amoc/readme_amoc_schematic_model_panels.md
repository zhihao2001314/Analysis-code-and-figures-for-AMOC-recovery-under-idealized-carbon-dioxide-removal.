# n57 / n58：三组示意图与六个实际模式路径

用户要求三种示意图竖排，每组右侧对应两个模式，共3×3九联图，所有横轴使用CO₂。本图为通用科研工作稿，未指定期刊版心或投稿阶段，不宣称期刊合规。

## 布局和数据口径

第一行：基本恢复组示意、GFDL-ESM4、NorESM2-LM。第二行：未充分恢复组示意、CESM2、CNRM-ESM2-1。第三行：按用户要求标为Larger-bias models · Substantial recovery · n = 4，右侧仍为CanESM5、UKESM1-0-LL。第三组原成员仍为ACCESS、CanESM5、MIROC、UK四个，n=4不是只剩两个；本图不套用TOA完整三成员子集。

横轴全部为物理CO₂浓度(ppm)，不是时间或对数浓度，所有面板范围250–1200ppm。纵轴为每个模式相对自身初始1–10年原始AMOC均值的强度，0–1.2共同尺度。保留n22已有11年平滑，不新增平滑、不读原始NetCDF、不做漂移订正。蓝实线为增碳，红虚线为降碳，箭头沿实际时间；140年转折点两段共用。原始曲折、短窗口边缘和实际初值均保留，不将实际模式曲线强行改成左侧示意形状。

左列复用n55现有601节点、稀疏锚点PCHIP示意，每段归一化浓度转换为真实ppm，纵向数值完全不改。空心圆是概念初始1；红方块是末期271–280年相对初始1–10年的成员恢复率均值，放在概念C0上；第三行细灰棒/短横是全部四成员恢复率范围/个体，不是CI。填色仅表示两分支分离，不是置信带。右侧模式图不画概念端点符号，避免将平滑瞬时端点与十年恢复率混称。

图内不加总标题或Data-informed schematic长标签；左列简写Schematic用于区分右侧数据曲线。Data-informed指依据实际组平均路径抽锚点概括，不代表完整逐年轨迹，更不是平衡分支。不得从本图单独证明永久不可逆。

## 第三组选择依据

预先按用户“最符合第三种路径”解释为完整增碳和降碳形状接近现有第三组示意。将该示意插值到每个候选模式各年的实际CO₂上，分段算均方误差；两个阶段等权，再开平方得到相似度RMSE。各阶段内每年等权；140年共用。分数只用于示例选择，不是模式技能评分或独立验证，不能据此推广整个第三组。

| 模式 | 增碳RMSE | 降碳RMSE | 等阶段RMSE | 选入 |
|---|---:|---:|---:|---|
| CanESM5 |0.023883|0.097915|0.071266|是|
| UKESM1-0-LL |0.071090|0.119079|0.098065|是|
| MIROC-ES2L |0.087122|0.174622|0.137991|否|
| ACCESS-ESM1-5 |0.124795|0.313777|0.238778|否|

CanESM5末期恢复率约1.0203，UK约0.6404；第三组四成员平均约0.8056，但个体不同。Substantial recovery是本次用户指定的组显示名称，不改各成员R值、原R≥0.8判据或暗示所有成员末期均超过0.8，也不意味着单调持续恢复。全部候选分数保存，未删去未入选模式数据。

## 文件和维度

计算src/s01_amoc/n57_amoc_schematic_model_panels.py；保存output/s01_amoc/n57_amoc_schematic_model_panels/的panels.csv、coordinates.csv、arrows.csv、third_group_similarity.csv和metadata.json。coordinates为panel×phase×node；示意各601节点，模式增碳140节点/降碳141节点。selection只读既有n22/n55 output，源哈希写入metadata。

绘图plot/s01_amoc/n58_amoc_schematic_model_panels.py只读上述output；正式PNG/PDF/SVG与清单位于figure/s01_amoc/n58_amoc_schematic_model_panels/。9.2×8.4英寸、600dpi PNG5520×5040，白背景，PDF嵌入字体、SVG保留文字。独立QC逐点核对6模式全序列、3示意转换、18箭头、选择排名和源哈希，以及实际Line2D的18组坐标哈希。

scientific-visualization技能用于区分示意与数据、保持共同尺度及真实曲折、记录示例筛选和导出审阅；软件引用沿用本阶段原说明，不新增期刊合规声明。

2026-10-26第三行改名：仅显示名称由Mixed recovery改为Substantial recovery，保留Larger-bias models、n=4、原示意和所选模式。旧输出、图件及4份代码/说明归档至archive/202610262000_amoc_third_group_title，源目标数量和SHA256核验。QC额外验证坐标、箭头、筛选表文件字节哈希与18条实际绘线哈希均不变。未修改全局分组注册表和其他图件。

独立QC通过：18组实际绘制坐标、18箭头、六模式全序列和三示意均与源表一致，四候选排名复核通过。最终PNG与独立PDF渲染均已目视，无标题/轴/图例裁切或缺字。原三联图及旧QC/临时诊断已核验归档至F:/sea/Antic sea ice constraint/archive/20261026141408_amoc_nine_panel_layout/，n55科学输入保持原位且哈希未变。
