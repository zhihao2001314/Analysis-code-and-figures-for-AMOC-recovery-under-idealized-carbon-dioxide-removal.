> 当前40°N试排已由 n73/n74 的滞回形态分组图替代；见 `readme_hysteresis_morphology_groups.md`。下文为旧版生成口径，旧图已归档。

# 四行十四面板的40°N对比版

当前布局已按作者要求改为首行a占第1–2列，b占第3–4列，两者绘图区宽均为画布0.428，高0.230；与26.5°N版同步。下方12面板像素和所有科学路径及箭头保持不变。上一版及布局元数据18文件归档archive/202610292039_equal_two_column_first_row/并核对，当前导出QA与目视检查通过。

最新同步修订：k改为与c同型的平滑轮廓，端点及开口不变，与26.5°N版相同。修订前两版已归档至archive/202610292033_smooth_open_k_both_latitudes/并核对文件数和哈希。此次只改变k内部路径及箭头，其他面板像素、坐标、实际模式曲线完全不变；具体检查见output/s01_amoc/n70_four_row_pathway_40n/export_and_manual_review.json。

按作者“再画一版40°N”要求，只将n68/n69试排图的a换成40°N九模式AMOC异常，其余b–n不变。计算n70_four_row_pathway_40n.py，绘图plot/s01_amoc/n71_four_row_pathway_40n.py。图件位于figure/s01_amoc/n71_four_row_pathway_40n/amoc_pathways_four_rows_40n_preview.{png,pdf,svg}，绘图参数区明确列出完整保存路径并逐格式打印。

a读取已有l_40原始Forced AMOC strength年度列，每模式减自身未平滑初始1–10年均值，在完整已存年度序列上计算11年居中平均后显示1–280年。九模式等权均值与总体SD重算；与已有40°N初始参照图数据逐点一致。没有piControl漂移订正或原始NetCDF重读。40°N指各模式最接近目标的原生纬度，映射保存CSV。

其余7个绘图输入文件逐字节复制n68。26条路径、26个箭头、14个主轴边界与范围、两个图例边界均与n69相同；下三行及b概念区PNG逐像素完全相同。下排本来就是40°N，因此本版所有实际AMOC面板使用同一目标纬度，恢复组仍为3+3，初始偏差组为剩余三模式。

13.2×12英寸，PNG5940×5400、450dpi；PDF单页、SVG文字保留。已查看最终PNG和独立PDF栅格，无裁切；源SHA256及无自编函数结构检查通过。新目录保存，不覆盖26.5°N对比版或原Fig1；两个纬度版本按作者要求并存供比较。CSV维度为Model×Simulation year转simulation_year280×model9，路径长表维度继承n68并按具名字段读取。
