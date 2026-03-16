# 文件功能说明（sl_interference_simulate）

> 约定（重要）：
> - “复采样点”= 一个复数采样点（I/Q 一对），在二进制里占 `4` 字节（`int16 I` + `int16 Q`，little-endian，I/Q 交织）。
> - 绝大多数后续脚本按 **0-based 复采样点索引** 计数（第 0 个点就是文件数据区的第一个复采样点）。
> - 当前工作流默认处理 **纯数据文件（无文件头）**；只有 `read_iq_file.m` 用于解析带 100B 头的原始文件。

---

## 数据文件（.iq）

- `20250912222305_part1.iq`
  - 原始采集文件（包含 **100 字节文件头** + IQ 数据）。
  - 仅当你需要研究文件头/确认格式时使用；后续处理建议先裁切成无头纯数据。

- `20250912222305_part1_cut1.iq`
  - 裁切输出文件（通常约定为 **纯数据无头**）。
  - 具体由裁切脚本/参数决定。

- `20250912222305_part1_cut2.iq`
  - 裁切输出文件（通常约定为 **纯数据无头**）。
  - 当前多数学分析脚本默认以它作为输入。

-`20250912222305_part1_cut2.iq`
303104-348159
---


## 读取/格式识别

- `read_iq_file.m`
  - 作用：读取并打印 **100 字节文件头**（HEX+ASCII），并做“数据类型探测”（例如更像 `int16` 还是 `float32`）。
  - 典型用途：首次拿到原始文件时确认格式、查看头字段的可能含义。
  - 备注：这个脚本会从文件偏移 `100` 处开始 probe 数据（即假定头长 100）。

- `iq_read_int16_le.m`
  - 作用：从 **纯数据文件** 按 `int16 little-endian`、I/Q 交织格式读取，返回复数向量 `x`。
  - 接口：`[x, meta] = iq_read_int16_le(filename, startSample, numSamples, headerBytes)`
    - `startSample`：0-based 复采样点起点
    - `numSamples`：读取复采样点数
    - `headerBytes`：默认 0；只有在你明确要读取带头文件时才用（当前工作流一般不用）
  - 依赖：无

---

## 裁切（输出纯数据）

- `cut_iq_by_samples.m`
  - 作用：按复采样点范围裁切 IQ 文件，输出到新文件。
  - 用法：编辑脚本顶部参数：
    - `startSample`（0-based，起点）
    - `endSample`（0-based，包含端点）
  - 输出：写出 `outFile`，按脚本注释约定为 **纯数据无头**。
  - 重要备注：此脚本内部使用 `fseek(fin, startSample*4, 'bof')`，因此输入文件应当是“数据区从文件起始开始”的纯数据文件；如果你拿它直接裁 `20250912222305_part1.iq`（带 100B 头），需要你自己先处理头偏移，否则会错位裁切。

---

## 时域/频域可视化

- `plot_time_domain.m`
  - 作用：读取一段 IQ 数据并绘制 **时域幅度 |IQ|**（横轴为样点索引）。
  - 特点：启用 `zoom/pan/datacursormode`，可放大与点选查看坐标。
  - 新增：可选的“时域突起(burst)区间”自动识别（短时能量+阈值），并在图上用背景色标注区间、在命令窗打印起止样点。
  - 输入：默认按纯数据读取（通过 `iq_read_int16_le`）。
  === Time-domain Burst Detection ===
range=[0..19999999], winLen=4096, hop=1024, bins=19528
thr=102329, bursts=13
burst #1: [303104..348159] (len=45056)
burst #2: [1940480..1985535] (len=45056)
burst #3: [3579904..3624959] (len=45056)
burst #4: [5218304..5262335] (len=44032)
burst #5: [6856704..6900735] (len=44032)
burst #6: [8495104..8539135] (len=44032)
burst #7: [10133504..10177535] (len=44032)
burst #8: [11770880..11816959] (len=46080)
burst #9: [13409280..13455359] (len=46080)
burst #10: [15048704..15092735] (len=44032)
burst #11: [16687104..16731135] (len=44032)
burst #12: [18325504..18369535] (len=44032)
burst #13: [19962880..19999743] (len=36864)


burst #1: [303104..348159] (len=45056)
burst #2: [1940480..1985535] (len=45056)
burst #3: [3579904..3624959] (len=45056)
burst #4: [5218304..5262335] (len=44032)
burst #5: [6856704..6900735] (len=44032)
burst #6: [8495104..8539135] (len=44032)
burst #7: [10133504..10177535] (len=44032)
burst #8: [11770880..11816959] (len=46080)
burst #9: [13409280..13455359] (len=46080)
burst #10: [15048704..15092735] (len=44032)
burst #11: [16687104..16731135] (len=44032)
burst #12: [18325504..18369535] (len=44032)
burst #13: [19962880..20007935] (len=45056)
burst #14: [21602304..21647359] (len=45056)
burst #15: [23239680..23284735] (len=45056)
burst #16: [24879104..24923135] (len=44032)
burst #17: [26516480..26561535] (len=45056)
burst #18: [28155904..28199935] (len=44032)
burst #19: [29794304..29838335] (len=44032)
burst #20: [31432704..31475711] (len=43008)
burst #21: [33071104..33115135] (len=44032)
burst #22: [34709504..34753535] (len=44032)
burst #23: [36347904..36391935] (len=44032)
burst #24: [37985280..38031359] (len=46080)
burst #25: [39623680..39668735] (len=45056)
burst #26: [41263104..41307135] (len=44032)
burst #27: [42900480..42945535] (len=45056)
burst #28: [44538880..44583935] (len=45056)
burst #29: [46178304..46222335] (len=44032)
burst #30: [47815680..47860735] (len=45056)
burst #31: [49455104..49499135] (len=44032)
burst #32: [51092480..51137535] (len=45056)
burst #33: [52730880..52775935] (len=45056)
burst #34: [54370304..54413311] (len=43008)
burst #35: [56008704..56052735] (len=44032)
burst #36: [57647104..57691135] (len=44032)
burst #37: [59284480..59329535] (len=45056)
burst #38: [60922880..60967935] (len=45056)
burst #39: [62561280..62606335] (len=45056)
burst #40: [64200704..64244735] (len=44032)
burst #41: [65838080..65883135] (len=45056)
burst #42: [67476480..67521535] (len=45056)
burst #43: [69115904..69158911] (len=43008)
burst #44: [70754304..70797311] (len=43008)
burst #45: [72392704..72436735] (len=44032)
burst #46: [74031104..74074111] (len=43008)
burst #47: [75669504..75712511] (len=43008)
burst #48: [77306880..77351935] (len=45056)
burst #49: [78945280..78990335] (len=45056)
burst #50: [80583680..80627711] (len=44032)
burst #51: [82223104..82267135] (len=44032)
burst #52: [83860480..83905535] (len=45056)
burst #53: [85499904..85543935] (len=44032)
burst #54: [87137280..87181311] (len=44032)
burst #55: [88775680..88819711] (len=44032)
burst #56: [90414080..90459135] (len=45056)
burst #57: [92052480..92097535] (len=45056)
burst #58: [93690880..93734911] (len=44032)
burst #59: [95329280..95374335] (len=45056)
burst #60: [96967680..97011711] (len=44032)
burst #61: [98606080..98650111] (len=44032)
burst #62: [100244480..100288511] (len=44032)
burst #63: [101882880..101927935] (len=45056)
burst #64: [103521280..103566335] (len=45056)
burst #65: [105159680..105204735] (len=45056)
burst #66: [106798080..106843135] (len=45056)
burst #67: [108435456..108481535] (len=46080)
burst #68: [110074880..110118911] (len=44032)
burst #69: [111713280..111757311] (len=44032)
burst #70: [113351680..113396735] (len=45056)
burst #71: [114989056..115034111] (len=45056)
burst #72: [116628480..116672511] (len=44032)
burst #73: [118265856..118310911] (len=45056)
burst #74: [119904256..119949311] (len=45056)
burst #75: [121543680..121587711] (len=44032)
burst #76: [123182080..123227135] (len=45056)
burst #77: [124819456..124864511] (len=45056)
burst #78: [126458880..126503935] (len=45056)
burst #79: [128097280..128142335] (len=45056)
burst #80: [129735680..129780735] (len=45056)
burst #81: [131374080..131418111] (len=44032)
burst #82: [133012480..133056511] (len=44032)
burst #83: [134650880..134693887] (len=43008)
burst #84: [136289280..136333311] (len=44032)
burst #85: [137927680..137971711] (len=44032)
burst #86: [139566080..139610111] (len=44032)
burst #87: [141204480..141248511] (len=44032)
burst #88: [142842880..142886911] (len=44032)
burst #89: [144481280..144526335] (len=45056)
burst #90: [146119680..146163711] (len=44032)
burst #91: [147758080..147802111] (len=44032)
burst #92: [149396480..149440511] (len=44032)
burst #93: [151033856..151078911] (len=45056)
burst #94: [152673280..152717311] (len=44032)
burst #95: [154311680..154355711] (len=44032)
burst #96: [155950080..155994111] (len=44032)
burst #97: [157588480..157632511] (len=44032)
burst #98: [159226880..159269887] (len=43008)
burst #99: [160864256..160909311] (len=45056)
burst #100: [162502656..162547711] (len=45056)
burst #101: [164142080..164186111] (len=44032)
burst #102: [165780480..165824511] (len=44032)
burst #103: [167417856..167462911] (len=45056)
burst #104: [169056256..169101311] (len=45056)
burst #105: [170695680..170739711] (len=44032)
burst #106: [172334080..172378111] (len=44032)
burst #107: [173971456..174016511] (len=45056)
burst #108: [175610880..175653887] (len=43008)
burst #109: [176157696..176200703] (len=43008)
burst #110: [177249280..177293311] (len=44032)
burst #111: [178887680..178930687] (len=43008)
burst #112: [180526080..180570111] (len=44032)
burst #113: [182164480..182207487] (len=43008)
burst #114: [183802880..183846911] (len=44032)
burst #115: [185440256..185485311] (len=45056)
burst #116: [187079680..187122687] (len=43008)
burst #117: [188717056..188762111] (len=45056)
burst #118: [190356480..190400511] (len=44032)
burst #119: [191993856..192038911] (len=45056)
burst #120: [193633280..193677311] (len=44032)
burst #121: [195270656..195314687] (len=44032)
burst #122: [196910080..196954111] (len=44032)
burst #123: [198548480..198592511] (len=44032)
burst #124: [200186880..200229887] (len=43008)
burst #125: [201825280..201868287] (len=43008)
burst #126: [203463680..203506687] (len=43008)
burst #127: [205102080..205145087] (len=43008)
burst #128: [206740480..206784511] (len=44032)
burst #129: [208377856..208422911] (len=45056)
burst #130: [210016256..210061311] (len=45056)
burst #131: [211654656..211699711] (len=45056)
burst #132: [213293056..213337087] (len=44032)
burst #133: [214384640..214436863] (len=52224)
burst #134: [214931456..214976511] (len=45056)
burst #135: [216570880..216613887] (len=43008)
burst #136: [217116672..217166847] (len=50176)
burst #137: [217662464..217706495] (len=44032)
burst #138: [218209280..218252287] (len=43008)
burst #139: [219847680..219890687] (len=43008)
burst #140: [221485056..221529087] (len=44032)
burst #141: [223123456..223168511] (len=45056)
burst #142: [224762880..224805887] (len=43008)
burst #143: [226400256..226444287] (len=44032)
burst #144: [228039680..228082687] (len=43008)
burst #145: [229131264..229211135] (len=79872)
burst #146: [229677056..229721087] (len=44032)
burst #147: [231316480..231359487] (len=43008)
burst #148: [231862272..231942143] (len=79872)
burst #149: [232407040..232480767] (len=73728)
burst #150: [232953856..232997887] (len=44032)
burst #151: [234592256..234637311] (len=45056)
burst #152: [236230656..236275711] (len=45056)
burst #153: [237868032..237913087] (len=45056)
burst #154: [239506432..239552511] (len=46080)
burst #155: [241145856..241189887] (len=44032)
burst #156: [242784256..242828287] (len=44032)
burst #157: [244422656..244466687] (len=44032)
burst #158: [246061056..246105087] (len=44032)
burst #159: [247700480..247743487] (len=43008)
burst #160: [249336832..249381887] (len=45056)
burst #161: [250976256..251020287] (len=44032)
burst #162: [252614656..252658687] (len=44032)
burst #163: [254253056..254297087] (len=44032)
burst #164: [255891456..255935487] (len=44032)
burst #165: [257529856..257573887] (len=44032)
burst #166: [259168256..259212287] (len=44032)
burst #167: [259714048..259766271] (len=52224)
burst #168: [260259840..260304895] (len=45056)
burst #169: [260806656..260850687] (len=44032)
burst #170: [262445056..262489087] (len=44032)
burst #171: [264083456..264127487] (len=44032)
burst #172: [265721856..265764863] (len=43008)
burst #173: [267360256..267404287] (len=44032)

- `autocorr_test1.m`
  - 作用：对 `test1.iq`（或你指定的任意纯数据 IQ 文件）做**自相关**并绘图，重点观察在 `lag=N`、`lag=Lsym=N+Ng` 处是否存在相关峰。
  - 实现：FFT 方式计算非负 lag 的线性自相关，输出归一化相关 `rho(lag)=r(lag)/r(0)` 的幅度(dB)（只看峰）。
  - 特点：不依赖工具箱；图上用竖线标注 `N`、`Lsym`。

- `autocorr_test2_slide128.m`
  - 作用：对 `test2.iq` 做**128 点滑动窗口自相关**（一步=1点），得到随时间变化的局部自相关。
  - 实现：对每个窗口起点 `s` 计算 `rho_s[k]=r_s[k]/r_s[0]`（`k=0..127`），并用热力图显示 `|rho|`(dB) 随 `lag` 和窗口起点变化。
  - 特点：纯自相关；不做峰检测/OFDM 标注/相位分析。

- `autocorr_test2_slide128_curve.m`
  - 作用：对 `test2.iq` 做**128 点滑动窗口自相关**（一步=1点），但以“谱线”方式显示：横轴 `lag`，纵轴为自相关值（默认 `|rho|`）。
  - 实现：对每个窗口起点 `s` 计算 `rho_s[k]=r_s[k]/r_s[0]`（`k=0..127`），逐窗向前滑动并刷新同一条曲线。
  - 特点：纯自相关；不做峰检测/OFDM 标注/相位分析；可选 `plotDb=true` 显示 dB。

- `plot_psd.m`
  - 作用：对一段 IQ 数据做 Welch 方法功率谱密度（PSD）估计并绘图。
  - 特点：Welch PSD 为脚本内自实现（Hann 窗 + 分段平均 + `fftshift`），尽量避免工具箱依赖。
  - 可选：显示绝对频率 `Fc+f` 或基带频率 `f`。

- `plot_psd_blocks.m`
  - 作用：在 `[startSample, endSample]` 范围内分块读取，每块计算一张 PSD，并按网格分页显示。
  - 适用：观察频谱随时间/样点变化。
  - 备注：脚本注释写“每 50000 点”，实际参数由 `blockLen` 控制；以脚本参数为准。

- `plot_full_signal_energy.m`
  - 作用：对整段或大范围 IQ 数据进行能量可视化，快速定位高能量 burst 区间。
  - 适用：做粗粒度“先找哪里有信号、再做精细解调”的预筛选。

- `plot_subcarrier_energy.m`
  - 作用：将指定区间信号经 DDC、CFO、Farrow 重采样后，按 1024 点 FFT 统计**平均子载波功率**，用于肉眼确认有效子载波边界。
  - 当前参数风格：与 SSS 解调脚本统一为 3 个手动参数：
    - `read_start_sample`：开始读取点
    - `read_length`：读取长度
    - `sss_decode_start_idx`：重采样后分析起点
  - 图形：支持主图+局部放大（用于观察“断崖边缘”）；可用于反推固定索引区间。

---

## SSS 解调与调试

- `sss_demodulation.m`
  - 作用：**单点解调**脚本。对指定偏移 `target_offset` 完成 SSS 频域相位补偿、QPSK 硬判决、Hex 输出。
  - 输入参数（手动填写）：
    - `read_start_sample`
    - `read_length`
    - `sss_decode_start_idx`
  - 关键流程：
    - 409.6MHz 下先做 `-63MHz` DDC
    - CFO 补偿
    - LPF + Farrow 重采样到 60MHz
    - 1024 点 FFT
    - 有效子载波检测（动态阈值 + `dc_guard`）
    - 分段 `unwrap` + 分段 `polyfit` 做 4 次方域斜率估计并补偿
  - 当前输出：
    - 星座图
    - 4 组旋转歧义候选序列（0°/90°/180°/270°）的完整 Hex（命令窗打印）
    - 支持“全部 1024 子载波”或“仅有效子载波”两种解调模式。

- `sss_demodulation_sweep.m`
  - 作用：**偏移遍历实验**脚本。对 `offset_range` 中每个偏移重复解调并输出摘要，便于挑选最佳对齐点。
  - 与单点脚本分离原因：避免“最终解调逻辑”和“扫描试验逻辑”互相污染。
  - 关键特性：
    - 使用与 `sss_demodulation.m` 同一套前处理与相位补偿策略
    - 每个 offset 输出：有效载波数/解调片段 Hex
    - 支持暂停逐帧查看（`pause_each=true`）
    - 支持全子载波解调模式开关。

- `estimate_sss_cp_60mhz.m`
  - 作用：围绕 60MHz 域的 SSS/CP 位置估计与诊断，辅助给出解调窗口的候选起点。

- `check_sss_read_location.m`
  - 作用：用于核验读取起点与窗口位置是否落在预期 burst/符号附近，减少“读偏”导致的误判。

- `SSS_Demodulation_Process.md`
  - 作用：SSS 解调流程说明文档，记录从粗同步到最终比特解调的步骤与关键参数。

---

## OFDM 符号定位（基于 CP 相关）

- `locate_ofdm_symbols_cp.m`
  - 作用：在给定 `N`（FFT 长度）与 `Ng`（CP 长度）下，通过 CP 相关度量定位 OFDM 符号边界。
  - 核心度量：
    - $P(d)=\sum_{k=0}^{Ng-1} x[d+k]\,\overline{x[d+k+N]}$
    - $R(d)=\sum_{k=0}^{Ng-1} |x[d+k+N]|^2$
    - $M(d)=\frac{|P(d)|^2}{R(d)^2+\epsilon}$
  - 输出变量：
    - `symCpStart0`：每个符号 CP 起点（0-based 复采样点索引）
    - `symDataStart0`：每个符号数据起点（0-based）= `symCpStart0 + Ng`
  - 鲁棒性增强：当前版本会尝试多个 top 峰作为锚点，按“整列符号边界的一致性（symRel 上 M 的均值）”选择最优锚点，减少真实信号里被杂峰带跑偏的概率。
  - 诊断输出：打印 `M(d)` 的统计量/分位数、`symRel` 处的 `M` 统计、全局 top 峰位置与峰间隔，用于判断“度量是否尖锐/是否呈现 Lsym 周期”。

- `estimate_ofdm_cp_locations.m`
  - 作用：将 OFDM CP 定位逻辑封装成**可复用函数**（供主入口或其他脚本调用）。
  - 输入：`inFile, N, Ng, startSample0, endSample0, Fs, opts`
  - 输出：结构体 `out`，包含：
    - `out.symCpStart0 / out.symDataStart0`
    - `out.diag`：包含 `M(d)`、锚点选择信息、top 峰、统计量等
  - 备注：`locate_ofdm_symbols_cp.m` 现在是一个薄封装脚本，内部调用本函数并负责绘图。

---

## 估计功能主入口

- `estimate_main.m`
  - 作用：将“估计/检测类功能”统一到一个主入口，通过 `action` 参数选择要运行的功能。
  - 当前支持：
    - `action='ofdm_cp'`：OFDM CP 定位（调用 `estimate_ofdm_cp_locations.m`）
    - `action='guard_intervals'`：帧间保护间隔检测（调用 `find_frame_guard_intervals.m`）
  - 明确不包含：验证类功能（例如 `segment_head_tail_similarity.m` / `run_segment_similarity.m`）。

- `run_estimate_main.m`
  - 作用：**点击运行版**参数脚本。在脚本顶部填写 `action` 与参数，直接运行即可调用 `estimate_main`。
  - 支持：
    - `action='ofdm_cp'`：运行 OFDM CP 定位，并可选绘制 `M(d)` 曲线与符号起点标记
    - `action='guard_intervals'`：运行帧间保护间隔检测，并可选绘图

---

## “切片自证”（CP 首尾相似度验证）

- `segment_head_tail_similarity.m`
  - 作用：量化两段序列的“首尾相似程度”，常用于 OFDM CP 自证：
    - head = `x[d0 : d0+L-1]`
    - tail = `x[d0+offset : d0+offset+L-1]`（OFDM 常用 `offset=N`，`L=Ng`）
  - 输出 `result`：包含 `rho`（归一化复相关系数）、`|rho|`、相位、幅度 RMS 比、NMSE/NMSE(dB) 等，并可选画图对比。
  - 备注：为兼容无工具箱环境，幅度 RMS 使用 `sqrt(mean(|x|^2))` 实现。

- `run_segment_similarity.m`
  - 作用：参数入口脚本，只负责：选择输入文件、选择 `d0/N/Ng`、设置扫描/绘图选项，然后调用 `segment_head_tail_similarity`。
  - 功能：
    - `d0` 附近扫描（`scanRadius/scanStep`），自动找 `|rho|` 最大（或 NMSE(dB) 最小）的点
    - 输出 center/best/top candidates，判断是否“只差少量采样点”
    - 可选参数扫描：`offsetCandidates` / `LCandidates`，快速排查 `N/Ng` 是否配错
  - 备注：脚本会 `clearvars -except symCpStart0 symDataStart0`，以便你先运行 `locate_ofdm_symbols_cp.m` 后复用其输出。

---

## 帧间保护间隔（Guard Interval）检测

- `find_frame_guard_intervals.m`
  - 作用：在真实信号中自动寻找“帧（burst）—保护间隔（低能量gap）—帧（burst）”结构，输出每帧区间与帧间保护间隔区间。
  - 方法：短时能量包络（滑动窗口）→ 平滑 → 自适应阈值（默认 median+K*MAD）→ 连通段提取。
  - 输入：IQ 纯数据文件路径 + 搜索范围 `[startSample0, endSample0]`（均为 0-based 复采样点索引）。
  - 输出：
    - `framesStart0/framesEnd0`：每个帧（burst）区间
    - `guardsStart0/guardsEnd0`：相邻两帧之间的保护间隔区间
    - `energyBins0/energy/threshold`：能量bin索引、能量序列、阈值（便于诊断）
  - 备注：边界定位受 `winLen` 影响（约 `winLen` 量级的时间模糊），需要更精确边界时可对结果附近二次精细化。

---

## 其他

- `.gitignore`
  - Git 忽略规则文件（不影响 MATLAB 功能）。
-
- `.git/`
  - Git 仓库元数据目录。

---

## 全量文件速查（按当前目录）

> 说明：下面按“单文件 + 批量模式”覆盖当前目录全部文件，描述保持简短。

### 1) 系统/仓库/文档文件

- `.DS_Store`：macOS 目录缓存文件（可忽略）。
- `.gitignore`：Git 忽略规则。
- `FILES.md`：本文件，项目文件功能说明。
- `SSS_Demodulation_Process.md`：SSS 解调流程文档。
- `burst_info.txt`：burst 区间或统计结果文本。
- `上行信号分析进展与疑问2.10.docx`：中文分析报告文档。
- `上行信号分析进展与疑问2.10(1).docx`：中文分析报告文档（版本副本）。
- `对上行信号的分析(1)(2).docx`：中文分析报告文档（版本副本）。

### 2) 原始/裁切 IQ 数据文件

- `20250912222305_part1.iq`：原始采集 IQ（通常含头）。
- `20250912222305_part1_cut1.iq`：裁切 IQ 数据片段。
- `20250912222305_part1_cut2.iq`：裁切 IQ 数据片段。
- `20250912222305_part1_cut3.iq`：裁切 IQ 数据片段。

### 3) 批量测试 IQ 数据文件

- `sigtest1.iq` ~ `sigtest173.iq`：批量测试 IQ 数据集（编号样本）。
- `sigtest1_test.iq`、`sigtest2_test.iq`、`sigtest3_test.iq`：特定测试版样本。
- `test1.iq`、`test2.iq`、`test3.iq`、`test4.iq`：独立测试样本文件。

### 4) 读取/裁切/基础处理脚本

- `read_iq_file.m`：带头文件读取与头部探测。
- `iq_read_int16_le.m`：纯数据 IQ 读取函数（int16 little-endian）。
- `cut_iq_by_samples.m`：按复采样点区间裁切 IQ 文件。
- `batch_cut_bursts.m`：按 burst 批量切片导出。
- `read_iq_file.m`：文件头和数据类型检查。

### 5) 时域、频域、星座图可视化

- `plot_time_domain.m`：时域波形与可选 burst 标注。
- `plot_psd.m`：单段 Welch PSD。
- `plot_psd_blocks.m`：分块 PSD（随时间观察）。
- `plot_full_signal_energy.m`：全段能量可视化。
- `plot_subcarrier_energy.m`：1024 子载波平均功率/边界观察。
- `constellation.m`：基础星座图绘制。
- `constellation_interp.m`：插值后星座分析。
- `plot_constellation_fixed_D.m`：固定参数 D 的星座显示。
- `plot_constellation_range.m`：参数范围扫描的星座显示。
- `check_ofdm_constellation.m`：OFDM 星座检查。
- `check_ofdm_psd.m`：OFDM 频谱检查。
- `check_modulation_type.m`：调制方式粗判别。

### 6) 自相关/互相关/结构搜索

- `autocorr_test1.m`：自相关基础验证。
- `autocorr_test2_slide128.m`：滑窗自相关热图。
- `autocorr_test2_slide128_curve.m`：滑窗自相关曲线动画。
- `autocorr_test2_blind_search.m`：盲搜参数下的自相关扫描。
- `autocorr_test2_visual_scan.m`：可视化批量扫描窗口长度/参数。
- `autocorr_find_structure_auto.m`：自动化结构搜索（自相关方向）。
- `cross_corr_sig1_sig2.m`：两段信号互相关定位。
- `corr_sigtest1_sigtest2.m`：sigtest1/sigtest2 相关对比。
- `cross_corr_interp.m`：插值互相关。
- `find_best_match_across_files.m`：跨文件最佳匹配搜索。
- `analyze_data_structure.m`：数据结构分析总脚本。
- `highlight_structure_sigtest1.m`：sigtest1 结构高亮查看。

### 7) OFDM 同步/估计/验证链路

- `locate_ofdm_symbols_cp.m`：基于 CP 的 OFDM 符号起点定位。
- `estimate_ofdm_cp_locations.m`：CP 定位函数封装。
- `estimate_main.m`：估计主入口（action 分发）。
- `run_estimate_main.m`：估计入口的运行参数脚本。
- `find_frame_guard_intervals.m`：帧间保护间隔检测。
- `segment_head_tail_similarity.m`：CP 首尾相似度验证。
- `run_segment_similarity.m`：相似度验证入口脚本。
- `estimate_sss_cp.m`：SSS/CP 位置估计脚本（早期/通用版）。
- `estimate_sss_cp_60mhz.m`：60MHz 域 SSS/CP 估计脚本。
- `check_sss_read_location.m`：读取窗口和 SSS 位置核验。

### 8) SSS/PSS 解调与同步脚本

- `sss_demodulation.m`：SSS 单点解调主脚本（含相位补偿、Hex 输出）。
- `sss_demodulation_sweep.m`：SSS 偏移遍历解调脚本。
- `gardner_farrow_timing_recovery.m`：Gardner + Farrow 定时恢复链路。
- `gardner_sps_demo.m`：Gardner 定时恢复演示。
- `precise_sps_resample_demo.m`：精细采样率/插值重采样演示。
- `pss.m`：PSS 相关功能脚本/函数。
- `pss_test.m`：PSS 测试主脚本。
- `pss_sdpsk_fullchain_demo.m`：PSS/SDPSK 全链路演示。
- `mseq_matched_filter_sync.m`：m 序列匹配滤波同步。
- `check_mseq.m`：m 序列检查脚本。

### 9) 其他算法实验脚本

- `can_val_align.m`：候选值/序列对齐实验脚本。
- `find_best_match_across_files.m`：跨文件匹配评估。
- `optimize_D_kmeans.m`：基于 KMeans 的参数 D 优化尝试。
- `test_farrow_length.m`：Farrow 长度敏感性测试。
- `test_farrow_length2.m`：Farrow 长度敏感性测试（第二版）。
- `test_mseq_qpsk_correlation.m`：m 序列与 QPSK 相关测试。

### 10) MATLAB 自动备份文件

- `autocorr_test2_slide128.asv`：`autocorr_test2_slide128.m` 自动保存副本。
- `pss_test.asv`：`pss_test.m` 自动保存副本。
- `sss_demodulation.asv`：`sss_demodulation.m` 自动保存副本。

