# Sun's type VI π-series: VI1–VI3 proof drafts

孙智伟 VI 型三条 π 级数的计算机辅助证明稿、全阶证书与复算程序。

This repository contains computer-assisted proof drafts for Zhi-Wei Sun's
three type VI π-series. It provides all-order telescoping certificates and
analytic arguments for convergence, period selection, and boundary constants.
The drafts have undergone AI-assisted cross-checking; professional mathematical
review, priority, and prize eligibility have not been established.

## 题目来源

原始出处是孙智伟的报告讲义 [Mysterious π-Series (II)](http://maths.nju.edu.cn/~zwsun/Pi-series(II).pdf#page=40)，
封面日期为 **2023 年 1 月 12 日**，所研究的三式位于 **第 40/68 页**，
页面标题为 **My conjectural series of type VI**。
该版本在三式之后提出针对“三式首次严格证明”的 300 美元悬赏。

本仓库仅处理该页顶部三条无穷级数等式，不包含下方的模 p² 同余猜想。
逐题对应、原始链接和版本核对说明见 [PROBLEM_SOURCE.md](./PROBLEM_SOURCE.md)。
另提供一份含出处的 [知乎回答草稿](./ZHIHU_ANSWER.md)。

## 阅读顺序

**优先阅读：[三式完整统一 PDF（15 页，无需运行 Python）](./paper/sun_vi_standalone_proof_2026-10-26.pdf)。**
该稿包含全部证书、有限代数核查规则、解析证明及最后求值；
[LaTeX 源文件和章节说明](./paper/README.md)一并提供。
下面的分稿与核验代码保留作补充材料。

1. [VI1、VI2 证明正文](./VI1_VI2_PROOF_2026-10-25.md)：同一参数族到两条经典 Ramanujan 级数。
2. [VI3 证明正文](./VI3_PROOF_2026-10-25.md)：新参数族、实际周期的三项匹配与非零边界常数。
3. [VI3 全部字面证书](./VI3_EXPLICIT_CERTIFICATES_2026-10-25.md)：不依赖 JSON 即可检查。
4. [逐条复核回应](./gpt_review_response_2026-10-26.md)：裁定 DeepSeek 的 P1–P4，并纠正复核报告中的转述错误。
5. [证明索引](./PROOF_INDEX_2026-10-25.md)、[复算说明](./PROOF_REPRODUCTION_2026-10-25.md)。

## 三条命题

定义 $T_n(b,c)=[x^n](x^2+bx+c)^n$。证明稿对应以下三式：

$$
\sum_{n=0}^{\infty}(66n+17)\frac{T_n(10,121)^3}{55296^n}
=\frac{540\sqrt2}{11\pi},
$$

$$
\sum_{n=0}^{\infty}(126n+31)\frac{T_n(22,441)^3}{(-512000)^n}
=\frac{880\sqrt5}{21\pi},
$$

$$
\sum_{n=0}^{\infty}(3990n+1147)\frac{T_n(62,9025)^3}{(-23887872)^n}
=\frac{432(195\sqrt{14}+94\sqrt2)}{95\pi}.
$$

题目出处及历史悬赏说明见上文“题目来源”；本仓库不据此主张奖金资格。

## 最短复算

建议 Python 3.11；本次实际复算使用 Python 3.11.15、SymPy 1.14.0。
在仓库根目录运行：

```sh
python -m pip install -r requirements-proof.txt
python -I verify_vi12_elliptic_reduction_210.py
python -I verify_vi12_hypergeometric_211.py
python -I verify_vi3_proof_220.py
python -I verify_vi3_literal_certificates_223.py
```

前三项核查已明确列出的有理函数证书、拉回和常数；第四项从字面公式重建，不读取 JSON。
所检验的证书中 n 是符号，区别于只检验前若干项的级数猜测。
正常收敛、实际周期及极限交换的理由写在正文，不能被脚本退出码或数值吻合替代。

可选的补充审计：

```sh
python -I audit_vi3_uniform_endpoint_remainder_221.py
python -I audit_state_and_vi3_cm_direct_222.py
python -I audit_deepseek_final_review_225.py
```

`package_proof_review_224.py` 从证明清单建立独立目录，用 `-I` 重跑四个核验器，
检查内部哈希后生成 [最小审稿 ZIP](./PROOF_REVIEW_BUNDLE_2026-10-25.zip)。
环境与实际退出码见 [验证记录](./VALIDATION.md) 和 `proof_bundle_validation_224.json`。

另附 [DeepSeek 独立 AI 复核材料](./review/README.md)。
其中包含已指出的问题，阅读时须同时看逐条回应；AI 复核不是期刊审稿。

## 证明路线与外部输入

VI1、VI2：由三项式母函数建立四维状态与全阶望远镜证书，识别一个 ${}_3F_2$ 拉回；
独立处理固定参数权导数后，分别化到 $z=2/27$、$4/125$ 的经典求值。

VI3：另取一参数族，以全阶证书得到三阶方程；用可求和主控匹配 $1,\sqrt t,t$
三项确定实际椭圆周期组合，再独立确定非零齐次常数 $3/(4\pi)$。
在 CM 点先求导、再使用七次奇异模数与 $\alpha(7)$ 数据，精确还原右端。

主要外部输入：

- Chan–Liaw, *Cubic modular equations and new Ramanujan-type series for 1/π*,
  Pacific J. Math. 192 (2000), p.219, (1.1)–(1.2)：[正式论文](https://msp.org/pjm/2000/192-2/pjm-v192-n2-p02-p.pdf)。
- Borwein–Borwein, *Explicit Ramanujan-type approximations to pi of high order*,
  Proc. Indian Acad. Sci. (Math. Sci.) 97 (1987), p.57, N=7：
  [论文](https://www.cecm.sfu.ca/~pborwein/PAPERS/P35.pdf)。
- Clausen、Gauss 求和、Gamma 反射和 Legendre 关系，具体用途见正文。

第三方论文仅提供链接，不在此仓库重新分发。

## 材料来源与状态

研究过程使用了 ChatGPT/Codex 和 DeepSeek，涉及推导、证书构造、程序核验和交叉复核；
不应将 AI 的参与写成专业数学家的独立审稿。详细说明见 [公开范围与来源](./RELEASE_SCOPE.md)。

公开副本对 Markdown 公式分隔符和本机解释器示例做了格式转换；数学证书保持原样。
逐文件对应关系见 [SOURCE_PROVENANCE.json](./SOURCE_PROVENANCE.json)。
本仓库不包含长期探索草稿、下载论文、个人通信或整个原工作区。
