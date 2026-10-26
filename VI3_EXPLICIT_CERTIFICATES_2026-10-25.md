# VI3：可直接抄算的全阶求和证书

本文是 [VI3 证明稿](./VI3_PROOF_2026-10-25.md) §2 的附录。
它把此前保存在 JSON 中的两组证书完全展开，读者无需研究日志或原构造器，
即可用有理函数运算独立核查。解析收敛、端点及实际周期仍见正文。

沿用正文 (3) 的状态矩阵 $M_n,N_n$，令
$$
D=t^2-t-1,\quad b=2(t^2+t-1),\quad
d=16t(t^2-1),\quad w=\frac1{64(t-1)(t+1)^2}.
$$
母函数微分中的两个系数简化为
$$
\sigma=\frac{2t-1}{D},\qquad \tau=-\frac{4(t^2+1)}{D}.
$$
令行向量 $e_0=(1,0,0,0)$，并递归定义
$e_{j+1}=e_j'+e_jM_n$。于是 $D_t^j(U_n^3w^n)=e_jX_n$。

## A. 三阶周期方程的四分量证书

对正文 (5) 明确给出的 $a_0,a_1,a_2$，取
$$
\begin{aligned}
R_{n,0}&=
\frac{n^3(6t^5+3t^4+10t^3-22t^2+8t+3)}
 {t^2(t-1)^2(t+1)D^3},\\
R_{n,1}&=\frac{6n^2}{t^2(t-1)(t+1)D^3}
\bigl(12nt^6-4nt^5-36nt^4-4nt^3-10nt^2-20nt-2n\\
&\hspace{35mm}{}-4t^5-15t^4-7t^3+8t^2+t+1\bigr),\\
R_{n,2}&=\frac{24n^2}{t(t-1)D^3}
\bigl(2nt^6-20nt^5+10nt^4-4nt^3+36nt^2-4nt-12n\\
&\hspace{35mm}{}+t^6-t^5+8t^4+7t^3-15t^2+4t\bigr),\\
R_{n,3}&=-\frac{64n^3(3t^5-8t^4-22t^3-10t^2+3t-6)}{D^3}.
\end{aligned}
$$
直接通分约分即可核查四个有理函数恒等式
$$
\boxed{\ R_{n+1}N_n-R_n=e_3+a_2e_2+a_1e_1+a_0e_0\ },\qquad R_0=0.
\tag{A}
$$
乘以 $X_n$ 后求和，就得到正文的三阶方程。无穷远边界由正文正常收敛证明，
不能由仅有 $R_0=0$ 省略。

## B. 固定参数权导数的二分量证书

沿用正文 (7) 的 $\mathcal A,\mathcal B,\mathcal R$，写
$$
\rho=\frac{\mathcal R'}{\mathcal R}
=\frac{2t}{t^2+1}+\frac1{t+1}+\frac1{2(t-1)}
-\frac1{2t}-\frac{2(2t-1)}D .
$$
定义
$$
Y_n=ne_0-\mathcal A e_0-\mathcal B e_1,\qquad
S_n=\left(
0,\;
-\frac{3n^2(t-1)(t+1)(t^2+1)(2t^2+7t+1)}{tD^3},\;
-\frac{12n^2(t-1)(t+1)^2(t^2+1)(t^2-7t+2)}{D^3},\;
0\right).
$$
另一组完整的有理函数恒等式为
$$
\boxed{\ S_{n+1}N_n-S_n=Y_n'+Y_nM_n-\rho Y_n\ },\qquad S_0=0.
\tag{B}
$$
故 $Y=H-\mathcal A G-\mathcal B G'$ 满足 $Y'=\rho Y$，即 $Y=C\mathcal R$。
正文从原级数端点严格确定 C=3/(4π)，而不是把齐次常数设成零。

## C. 无 JSON 依赖的复算

verify_vi3_literal_certificates_223.py 对本页字面公式重建 M、N，
通分核对 (A)、(B) 的全部分量及两个 n=0 边界。
它不读取原构造输出 JSON、有限系数缓存或最终求值数据。

~~~powershell
$py = 'python'
& $py verify_vi3_literal_certificates_223.py
~~~

与主证明一样，该证书附录不代替外部数学审稿，且不承担首证或奖金资格认定。
