# AMOC 滞回形态分组试排

本次针对作者提供的 40°N 十四面板图调整，首行 a/b 保持原数据和各占两列的布局。

## 当前路径

- 计算：`src/s01_amoc/n73_hysteresis_morphology_groups.py`
- 输入整理及诊断：`output/s01_amoc/n73_hysteresis_morphology_groups/`
- 绘图：`plot/s01_amoc/n74_hysteresis_morphology_groups.py`
- 图件：`figure/s01_amoc/n74_hysteresis_morphology_groups/amoc_hysteresis_morphology_40n.png`，以及同名 PDF/SVG。

绘图脚本参数区明确列出以上绝对输出路径，仅读取 output，不读取原始模式数据。

## 面板顺序与试排名称

| 行 | 示意图 | 模式及面板 | 标签 |
|---|---|---|---|
| 1 | b 为面积概念图 | a 为九模式40°N异常 | 首行不变 |
| 2 | c | d ACCESS-ESM1-5；e CAS-ESM2-0 | Closed, small area |
| 3 | f | g GFDL-ESM4；h NorESM2-LM；i CanESM5 | Near-closed, larger area |
| 4 | j | k CESM2；l CNRM-ESM2-1；m UKESM1-0-LL；n MIROC-ES2L | Open, larger area |

这对应旧 f/m 交换、旧 n 进入四模式组、旧第三组移到第二行。三个示意图保持相同宽度，实际模式图按每行2/3/4个数量分配剩余宽度。

## 科学口径

九条实际模式曲线和时间箭头逐点复用 n70，未改变平滑、原初十年归一化、数据取值、经纬度或漂移处理。第一行 PNG 上部1700行像素与旧版完全相同。

分组为当前图的描述性试排，没有建立面积阈值，也未改变后续恢复分组、统计或其他图件。GFDL、CanESM 端点也接近初始值，因此没有采用会把整行标成严格开放的名称。Near-closed 表示近似返回，不表示精确闭合。两组 larger area 的面积范围有重叠，其区别同时涉及末期返回程度。

ACCESS 和 CAS 的面积不是零。将每个模式的绝对分支差在共同 CO₂ 浓度范围内积分，再除以该范围，可得约0.045和0.065；其余模式约0.155–0.356。积分使用现有11年平滑路径和原有CO₂节点，补入所有交叉点后逐段线性积分，记录在 `model_hysteresis_metrics.csv`。该量的面积单位为归一化AMOC×ppm，除以CO₂范围后为无量纲，不等同于时间滞后。

c 为作者要求的近零面积概念图：两条近直线、端点闭合、内部仅有0.016 sin(πx)的微小分离；不是实际 ACCESS/CAS 的拟合或组平均。f/j 沿用旧概念轮廓，只更新概念端点的参考值。

## 维度与验证

- 所有源数据为已有 CSV，不新增 NetCDF/MAT 读取。
- overview：`simulation_year`280行×9模式，另有多模式均值及总体标准差上下界。
- 路径：`panel/phase/node` 长表；模式每分支140或141点，概念图每分支601点。
- 所有9模式×2分支与原图数据哈希一致，18条时间箭头一致，a/b数据与布局不变。
- PNG 为5940×5400、450 dpi；PDF/SVG保留矢量路径。PNG及PDF栅格已目视检查。
- 当前为作者审阅试排，不代表已完成期刊最终制作尺寸审查。

旧40°N与26.5°N试排图在新版验证后归档。作者随后要求第一组直接标记 Closed, small area；两版均已同步此标签，未强制修改任何实际模式的端点。

26.5°N版本已纠正为整图同纬度：`src/s01_amoc/n75_hysteresis_morphology_26p5.py` → `output/s01_amoc/n75_hysteresis_morphology_26p5/` → `plot/s01_amoc/n76_hysteresis_morphology_26p5.py` → `figure/s01_amoc/n76_hysteresis_morphology_26p5/amoc_hysteresis_morphology_26p5.png`，另有同名PDF/SVG。图a和九模式实际路径现在均来自26.5°N；此前只更换图a、下方仍为40°N的混合版本已归档，不再作为当前结果。新26.5°N版详见 `readme_hysteresis_morphology_26p5.md`。
