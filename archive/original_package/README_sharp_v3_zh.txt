多项解析多项式零点：最新 TeX、PDF 与附件材料包
交付日期：2026-10-08

一、最新稿
polyanalytic_sharp_v3.tex：包含截至第二轮九小时研究结束时已写入的结果。
polyanalytic_sharp_v3.pdf：由本包同名 TeX 重新编译；198 个物理页面，含 1 页封面，
正文与目录等编号页共 197 页。本次交付修正了封面信息溢到下一页的排版，
封面后全部源文、197 页的文字与公式内容均与研究结束版本一致。
这份研究记录仍有未决端点，不能称为全部结果已达到全局 sharp。

最新主要进展（精确范围和证明以正文为准）：
- 对每个偶数 n>=10，M_(2,n)=9n-10，达到全局上界。
- E_(2,7) 属于 {45,47}；已证明存在恰有 45 个非退化零点的实整数系数构造，
  并证明该临界种子邻域内最大值为 45；尚未决定全局 47 端点。
- 偶数 n>=8 时，E_(2,n) 属于 {9n-18,9n-16}；奇数 n>=9 时，
  E_(2,n) 属于 {9n-20,9n-18,9n-16}。
- 88<=E_(3,10)<=100；至少 84 个零点的未被支配次数组合为 (1,18),(2,12),(3,10)。
- M_(2,7)、一般 n=3m-1 缺口和其他尚未闭合端点仍按正文列为未决。

二、包内结构
supplement/：39 个已公开的 verify_*.py 脚本，39 份 verification.json 报告，
             7 份实际可用的精确 certificate.json 数据。脚本和数据保持同目录，
             以支持 with_name 路径及脚本之间的局部导入。
audits/：V41 局部模型与限制映射的有限精确检查脚本及报告、V42 实种子检查报告、
        V43 研究结束时的历史总审查，以及本次交付的 delivery_verification.json。
original_inputs/：用户原先提供的四份原始材料，原样保留，仅作为历史输入；
                 根目录的 v3 稿才是当前最新稿。
MISSING_HISTORICAL_ATTACHMENTS.txt：原稿曾引用但未随原始材料提供的旧附件清单。
requirements.txt：附件脚本使用的第三方 Python 依赖及本次环境版本。
SHA256.json、check_package_integrity.py：交付文件的校验清单和校验工具。

三、复现
建议 Python 3.12；本次环境为 Python 3.12.14。
解压后先在包根目录执行：
  python check_package_integrity.py
  python -m pip install -r requirements.txt
  python supplement/verify_pure_2_7_45.py
  python supplement/verify_pure_2_7_41_fixed.py
  python supplement/verify_pure_3_10_88.py
其他脚本通常可同样运行；需要自定义证书或报告路径的脚本有 --help 参数。
不要使用 python -O 或 -OO；核验依赖精确比较和断言。
多数脚本会更新同目录报告，运行后的时间字段可能变化，因此请先核对交付哈希。
历史脚本耗时不尽相同；本次交付重新运行并通过了上述三个脚本，
其余报告为研究期间保留的结果，没有声称本次重新执行全部历史检查。

编译 TeX 需要完整 TeX Live（本次使用 2024）或相应的 XeLaTeX、ctex 等宏包。
在包根目录运行以下命令，首次构建重复到交叉引用稳定（通常三次）：
  xelatex -interaction=nonstopmode -halt-on-error polyanalytic_sharp_v3.tex
源稿为独立单文件，不需要外部 .bib、图片或其他项目文件。
本次最终两遍 XeLaTeX 编译通过，无错误、警告、未解析引用或排版溢出。
已检查全稿页面概览及封面、新证明和结论的清晰页面。
内置编辑器的编译服务仍有 Windows 沙箱初始化错误；这不影响本包已成功导出的 PDF。

四、核验范围
verify_pure_2_7_45.py 验证 45 零点解析构造的精确代数恒等式、插值条件及局部展开。
它不是某个固定 epsilon 的根方框证书；全局计数依赖正文中完整的小参数证明。
verify_pure_2_7_41_fixed.py 验证另一个固定参数实例：41 个实切片零点与 8 个非物理
复交点的精确收缩证书。41 与 45 的脚本、报告和适用范围必须区分。
verify_pure_3_10_88.py 是 88 零点构造的精确代数及局部模型检查，
同样不宣称固定最终小参数的全局根方框证书。
V41/V42 审查中的有限次数检查不替代正文的任意次数证明。
旧报告中的 source_sha256 / verifier_sha256 与本次源稿的 SHA256 不是同一概念；
研究末尾审查记录的是封面排版修正之前的源稿，当前完整文件以 SHA256.json 为准。
历史材料缺口已单独列出；不以缺失旧附件冒充本次已复现的证据。
