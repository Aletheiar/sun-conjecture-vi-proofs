# 三式证明材料的隔离复算

实际核验环境：Python 3.11.15、SymPy 1.14.0。
精确核验器只使用 Python 标准库和 SymPy；数值探索和旧研究脚本的依赖不属于必需项。

解压审稿包后，先读 PROOF_INDEX_2026-10-25.md 和两份正文。
若需另建环境，可在自选的独立 Python 环境安装 requirements-proof.txt，
然后从解压目录运行：

~~~powershell
python -I verify_vi12_elliptic_reduction_210.py
python -I verify_vi12_hypergeometric_211.py
python -I verify_vi3_proof_220.py
python -I verify_vi3_literal_certificates_223.py
~~~

-I 使用 Python 的隔离模式，避免当前目录或用户级 Python 路径提供额外模块。
证书 JSON 是明确的数学输入，已列入 proof_manifest_2026-10-25.json，
与“依赖原工作区的未交付缓存”不同。
第四条核验器不读取 JSON，可直接抄算附录的字面恒等式。

前两式正文还列出 208、209 两个可选构造器，已一并收录，
但快速核验不需要重复运行构造过程。217、219 的原构造器同理。

package_proof_review_224.py 可从清单新建独立临时目录、复制明确列出的文件，
运行上述四条核验，检查核验后文件哈希仍一致，再生成本地审稿 ZIP。
它不安装软件，不联网，不发送任何材料，不改写外部审稿目录。
隔离复算的日志与环境元数据放在 ZIP 内的 ISOLATED_VALIDATION/。

隔离复算只证明交付材料能独立运行及所列代数断言成立。
解析收敛、实际周期分支、端点主控和外部古典定理仍须审查正文；
它不代表数学界认可、首次证明或奖金资格。

## 2026-10-26 独立 AI 复核后的修订

DeepSeek 的四项意见及其报告中另需纠正的表述，已逐条复算，
处理结果见 gpt_review_response_2026-10-26.md。
数学证书不变；两份正文澄清归纳差值、椭圆积分导数记号、端点主控和全阶依据。
另可运行本次精确审计：

~~~powershell
python -I audit_deepseek_final_review_225.py
~~~

225 的有限 Taylor 诊断只定位复核报告中的手算错误；
原三式的全阶依据仍是上述四个核验器、证书和正文分析论证。
旧稿及原 ZIP 已单独保存在工作区 proof_revisions_2026-10-26/pre_deepseek_response/，
不包含在当前最小证明 ZIP 中。
完整 DeepSeek 目录也不随最小包分发；本包只含其意见的自有回应与证明材料。
