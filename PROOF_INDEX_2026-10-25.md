# 孙智伟 Conjecture VI1–VI3：证明稿入口

更新：2026-10-26。三式的完整证明稿已获 DeepSeek 独立 AI 重建复核，
并已逐项复算其意见；修订后的三条证明链均成立，尚未经专业数学审稿认可。
本文不构成首次证明或奖金资格的认定。

最新逐条处理见 [复核回应](./gpt_review_response_2026-10-26.md)。
P1 的两种归纳差值均正确；其余意见已用于澄清记号、几何主控范围和全阶证明依据。
报告另有 Y 二阶项、解基首项等转述错误，已独立定位，不沿用错误表述。

本地审稿 ZIP：PROOF_REVIEW_BUNDLE_2026-10-25.zip。
隔离复算步骤及环境见 [PROOF_REPRODUCTION_2026-10-25.md](./PROOF_REPRODUCTION_2026-10-25.md)。
ZIP 只含清单中的证明材料与隔离核验日志，不含整个研究工作区。

| 命题 | 当前正文 | 全阶独立核验 |
|---|---|---|
| VI1（包括 B1）、VI2 | [完整证明](./VI1_VI2_PROOF_2026-10-25.md) | verify_vi12_elliptic_reduction_210.py、verify_vi12_hypergeometric_211.py |
| VI3 | [完整证明](./VI3_PROOF_2026-10-25.md) | verify_vi3_proof_220.py |

两份正文都给出原题、参数缩放、实际周期、固定参数权导数、正常收敛与端点、
外部古典定理依赖和最终常数。主证明文档 §205–214 保留发现过程。
不必先读此前两百节失败探索才能审查新证明。

## 最短复现

~~~powershell
$py = 'python'
& $py verify_vi12_elliptic_reduction_210.py
& $py verify_vi12_hypergeometric_211.py
& $py verify_vi3_proof_220.py
~~~

这些核验器读取已明确保存的证书，在有理函数域重建并检查恒等式。
它们均已实际运行、退出码 0。日志在 proof_logs_2026-10-25/。
不需要重新运行较慢的构造器 217、219。
核验器检查代数内容；正常收敛、实际周期和主控换极限必须同时审查正文。

VI3 正文 §8 另有独立复核补充及显式二阶余项界；
audit_vi3_uniform_endpoint_remainder_221.py 核对其代数和实际 Puiseux 展开，
将原先 o(t) 的匹配余项加强为 O(t^(3/2))。

VI3 正文 §9 给出不经过 R3 的直接 K、E 消去；
audit_state_and_vi3_cm_direct_222.py 同时独立推导通用状态矩阵，
并只用七次奇异模数、alpha(7) 和 Legendre 关系核对最终值。

[VI3 字面证书附录](./VI3_EXPLICIT_CERTIFICATES_2026-10-25.md)
完整展开三阶方程和固定参数导数的两组行向量。
verify_vi3_literal_certificates_223.py 不读 JSON、级数缓存或最终求值数据，
可单独复算附录的全阶恒等式；解析边界仍须结合正文审查。

## 古典定理输入

1. 前两式使用 Chan–Liaw 正式出版稿列出的两条 Ramanujan 级数，参数 2/27、4/125：
   [Pacific J. Math. 192 (2000)，219 页 (1.1)、(1.2)](https://msp.org/pjm/2000/192-2/pjm-v192-n2-p02-p.pdf)。
   1998 会议稿 (1.1) 的印刷式有误，不能代替正式版。
2. 第三式使用七次奇异模数与 alpha(7) 的精确值：
   [Borwein–Borwein (1987)，57 页 N=7](https://www.cecm.sfu.ca/~pborwein/PAPERS/P35.pdf)。
   由此得到的经典辅助式 R3 已在主证明 §8 写出代数证明。
3. Clausen、Gauss 求和、Gamma 反射和 Legendre 关系等标准特殊函数恒等式。

## 优先审稿点

VI3 的难点是原周期三维组合和非零导数边界项，而不是最后的根式整理。
独立审稿提示词见 [PROOF_REVIEW_PROMPT_2026-10-25.md](./PROOF_REVIEW_PROMPT_2026-10-25.md)。
已完成独立 AI 复核及作者方逐项回应；专业数学审稿、首次性和奖金认定仍待进行。
新增 audit_deepseek_final_review_225.py 保存本轮对复核意见的精确检查；
它不取代上述四个全阶核验器及两份正文中的分析论证。
