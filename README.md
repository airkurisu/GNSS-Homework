# WHUT 卫星导航程序设计课设

C++14 实现的双频伪距单点定位（SPP）与双频多普勒测速。使用 RINEX 广播星历和观测文件，支持 GPS、北斗、Galileo、GLONASS，分别估计各系统接收机钟差与钟漂。

**定位和测速始终使用广播星历。** 精密轨道、钟差与 SINEX 仅可通过独立的离线工具分析误差，不参与 SPP 方程。

详细问题清单、修复依据、验证结果及尚未解决的限制见 [审查报告](docs/REVIEW.md)；精密产品离线参考分析见 [误差分析](docs/REFERENCE_ANALYSIS.md)。

## 构建与测试

需要 CMake 3.20+ 和支持 C++14 的编译器，无第三方 C++ 库依赖。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

Windows 使用 MinGW 时，可在配置命令中加 `-G Ninja`，并确保 MinGW 与 Ninja 在 PATH 中。Visual Studio 多配置生成器的程序位于 `build/Release/`；Ninja/Make 通常位于 `build/`。下列示例按单配置生成器书写，Windows 对应 `.exe` 文件。

## 运行

从仓库根目录运行，默认处理自带文件的全天数据：

```sh
./build/mygnss10 --mode multi --output build/multi.csv
./build/mygnss10 --mode gps --start 2 --end 6 --output build/gps.csv
./build/mygnss10 --nav data/brdc1590.24p --obs data/jfng1590.24o --systems GC --output build/gps-bds.csv
./build/mygnss10 --cold-start --start 2 --end 2.1 --output build/cold.csv
./build/mygnss10 --help
```

- `--mode gps|multi`：仅 GPS，或 GPS/BDS/Galileo/GLONASS。
- `--systems GCER`：直接指定系统字符；与 `--mode` 同时传入时，以最后一个选项为准。
- `--start` / `--end`：GPST 当天十进制小时，包含端点；默认 0–24，不支持跨午夜的倒序时段。
- `--nav` / `--obs`：输入文件路径，默认相对当前工作目录。
- `--cold-start`：每历元从零坐标开始，仍保留文件头坐标用于输出误差对照。
- `--output`：CSV 文件；省略则输出到标准输出，进度信息输出到标准错误。

退出码：`0` 至少一个历元定位成功；`1` 参数、读取或输出错误；`2` 有输入历元但全部定位失败。部分失败须查看 CSV 的 `status` 和终端汇总。

## CSV 与算法约定

| 字段 | 含义 |
| --- | --- |
| `week,tow` | GPS 周、GPST 周内秒 |
| `status` | `ok` 或具体失败原因；测速失败时有效位置仍输出 |
| `ns,nv` | 最终用于定位、测速的卫星数 |
| `pdop` | 无权几何位置精度衰减因子，包含实际参与系统的钟差未知数 |
| `x,y,z` | ECEF 位置，米 |
| `de,dn,du` | 相对观测文件头概略坐标的 ENU 差，米；不是经验证的真值误差 |
| `ve,vn,vu` | ENU 速度，米/秒 |
| `rejected` | 该历元被粗差检查剔除的伪距观测数 |

缺失或无效数值为 `nan`。系统钟差和钟漂还可通过 `SPPResult::clockBias/clockDrift` 获取，单位分别为秒和秒/秒。

使用无电离层组合 `P_IF = (f1² P1 - f2² P2)/(f1²-f2²)`，多普勒先转换为 `-D*c/f` 再组合。支持的频率组合为 GPS L1/L2、BDS B1I/B3I 或 B1I/B2I、Galileo E1/E5a 或 E1/E5b、GLONASS G1/G2。卫星钟差包含相对论改正；北斗组合单独处理群延迟。对流层采用标准气象条件下的 Hopfield + NMF，10° 高度角截止、`sin²(el)` 加权及北斗 GEO 降权。

位置用迭代加权最小二乘，速度方程用线性加权最小二乘，均由 Householder QR 求解。伪距粗差检查使用 IF 噪声放大系数和残差投影方差，至少需要两个冗余观测，以保守的 6σ 阈值逐颗剔除，最多五颗；这是基本粗差筛查，不是完整 RAIM。

## 结果统计与可选离线分析

Python 工具仅使用标准库：

```sh
python tools/summarize_results.py build/multi.csv
python tools/summarize_results.py build/multi.csv --start 2 --end 6
```

若本地另有 SP3/CLK/SINEX，可以按 [误差分析文档](docs/REFERENCE_ANALYSIS.md) 运行独立比对工具。这些参考文件不属于项目依赖，也不需要下载才能构建、测试或运行 SPP。

## 目录与边界

- `src/rinex.cpp`：RINEX 3.00–3.04 子集读取与信号选择。
- `src/gtime.cpp`、`src/orbit.cpp`：时间系统、广播轨道、速度和钟差。
- `src/spp.cpp`：两种模式共用的解算核心、坐标转换与结果输出。
- `src/tropo.cpp`、`src/matrix.cpp`：对流层与线性代数。
- `tests/`：算法、解析、异常输入、合成数据和实测回归。
- `tools/`：统计和离线分析，不向定位核心提供精密产品。

该程序仍是课程设计：未实现完整码偏差、天线相位中心、测站位移、精密坐标框架转换及完整性监测；不支持 RINEX 2/4、压缩输入、中途头更新、非单位观测缩放和非零历元接收机钟差偏移。具体边界见审查报告。
