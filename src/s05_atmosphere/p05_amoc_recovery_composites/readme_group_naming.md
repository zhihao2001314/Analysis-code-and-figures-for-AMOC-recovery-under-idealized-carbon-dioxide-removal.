# n49 统一分组显示名称

2026-10-22。按用户要求沿用两种AMOC滞回路径中的Substantial recovery和Incomplete recovery。当前2/2/4图使用Better-simulated / Incomplete recovery、Better-simulated / Substantial recovery及Larger-bias models。成员、R门槛、变量、压力层、时间与季节定义、模式/观测口径均保持原样。

正式命名表、首次源代码改名计划、18目录归档清单、科学文件哈希、导出日志和本轮QC位于`output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/`。说明及旧筛选图的适用边界见`docs/group_naming_convention.md`。第三组不赋予统一AMOC恢复状态，Better-simulated也不表示全球模式能力排名或已完成严格MLD观测验证。

这是文字显示与导出修订，没有读取原始NetCDF/MAT，没有新增维序推断或科学计算。223份既有科学输出逐字节保留；CSV/NetCDF机器键和历史panel_title字段未改，以免破坏下游读取。当前主要绘图脚本从group_labels.json读显示名称，并覆盖旧显示字段。原图文件名与路径不改；其中n066原目录为空，本轮从保存的output重出配套图。

计算侧脚本`n49_group_naming.py`只生成命名表、来源清单和补丁，不重算科学场，也不直接修改源代码；补丁通过apply_patch应用。`plot/.../n49_reexport_group_labels.py`只编排已有绘图脚本，不运行任何科学计算脚本。主图n059/n061/n48另外执行，已退役的n40只更新源代码以避免将来恢复时退回旧名称，不重出其活动图。概念滞回n23名称原本一致，不改。

本轮独立QC脚本为`check_group_naming.py`：核对223份科学输入SHA256、固定成员、59张图的三种格式及旧文件名、所有SVG文本中无旧恢复标签、PNG不透明及尺寸不变、PDF尺寸/字体、主N²/MLD和九面板13个绘图区与旧图像素一致、无自编函数/类，并输出5页目视检查缩略图。不以此声称新增显著性或因果验证。旧科学QC报告及其旧代码哈希保持历史原样。

重出图件或QC前必须按archive/README.md归档对应正式文件，不能直接覆盖。首次归档、标签换行修订和最终n30换行/QC归档的精确路径见图外命名规则。所有已更新的绘图脚本均可用`D:/Download/Anaconda/envs/viper/python.exe -X utf8`运行。
