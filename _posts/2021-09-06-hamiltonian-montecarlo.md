---
layout: post
title:  "ハミルトニアンモンテカルロの一般的な導出"
date:   2021-09-06 00:10:00 +0900
categories: jekyll update
---

ハミルトニアンモンテカルロ（HMC）は、既知関数$S(x)$によって$p(x) \propto e^{-S(x)}$のように書ける確率分布から変数をサンプリングする、マルコフ連鎖モンテカルロ法（MCMC）の一種です。

HMCを特殊例に含むMetropolis-Hastings法（MH法）は、現在の値$x$をもとに候補$x'$を何らかの方法で作り、確率$\alpha(x \to x')$で更新する（メトロポリステスト）というルーチンを繰り返します。
目的分布$p(x)$のサポートをまんべんなく訪問できるように候補を作り、適切な受理確率を設計することで、$\lbrace x_1, x_2, ... \rbrace$の分布が$p(x)$に収束することが保証されるというのが基本的な考え方です。

サンプル分布が早く目的分布に収束するためには連続したサンプル間の相関が小さいことが重要ですが、候補$x'$を適当に$x$から大きく離れた値にするとしばしば$\alpha(x \to x')$が非常に小さい値となり更新頻度が著しく下がってしまうという問題があります。このため、実用的なMCMC手法では、$x$から離れつつ受理確率が高くなるような候補$x'$を提案できるような工夫が行われています。

本記事では、代表的なMCMC手法の一つであるHMCを、その名前の由来となっているHamiltonian dynamicsを用いるものに限らず、なるべく一般化した形で導出することを試みます。

筆者は初学者のため誤りがあるかもしれません。メールや[twitter](https://twitter.com/quolc)でご指摘頂けると幸いです。
また、数学的な議論はかなり大雑把になっていると思います。

### 1. MCMCの条件

*本節の内容は[花田・松浦『ゼロからできるMCMC マルコフ連鎖モンテカルロ法の実践的入門』](https://www.amazon.co.jp/dp/4065201748)に基づきます。*

MCMCでは「サポートをまんべんなく訪問できるように候補を作り」「適切な受理確率を設計する」ことで目的分布$p(x)$が近似できると先述しました。これは具体的には次の4条件として書き下すことができます。

- マルコフ性

候補$x'$が、過去のサンプルの履歴${x_1, ..., x_{n-1}}$に陽に依存することなく、現在の値$x_n$だけから決まるということです。これは大抵アルゴリズムの設計によって明らかに満たされます。

- 既約性

$p(x)>0, p(x')>0$ である任意のペア$\lbrace x, x'\rbrace$について、有限回数の更新で互いに移り合うことが可能であるということです。これにより、どのような初期値$x_0$からサンプリングを始めても、分布に含まれるあらゆる$x$がサンプルできることが保証されます。
$p(x)$や候補提案の詳細に依存するので本記事では詳しく扱いません。

- 非周期性

$x \to x' \to ... \to x$という周期的な経路について、その長さ$n_s$の最大公約数が1であるという条件です。例えば、必ず偶数回のステップで回帰するようなケースでは満たされません。メトロポリステストを伴うMH法では極値以外では棄却により$n_s=1$が得られるので、ほぼ自動的に成立すると思われます。

- （詳細）釣り合い条件

定常分布$p(x)$に関して、遷移確率$T(x \to x')$が任意の$x, x'$について次を満たすのが詳細釣り合い条件です。

$$
\begin{equation}
p(x)T(x \to x') = p(x')T(x' \to x)
\end{equation}
$$

MH法における遷移確率$T$は、現在の値が$x$のときに$x'$を候補に選ぶ確率と受理確率を掛け合わせたものです。$x'=x$のときは$x$に留まる確率を意味します。

「詳細」の付かない釣り合い条件は、各$x$で正味の流出と流入がつりあっていることです。

$$
\begin{equation}
p(x) \int dx' T(x \to x') = \int dx' p(x')T(x' \to x)
\end{equation}
$$

ここで$\int dx' T(x \to x') = 1$なので次のように書き直すと定常分布の保存と理解できます。

$$
\begin{equation}
p(x) = \int dx' p(x')T(x' \to x) 
\end{equation}
$$

MCMCが成立するためには釣り合い条件だけを満たせば十分ですが、多くのアルゴリズムは詳細釣り合いを満たすように設計されています。
最近は詳細釣り合いを満たさないMCMCについても研究が盛んに行われているそうです。（[これ](https://arxiv.org/abs/1007.2262)とか[これ](https://arxiv.org/abs/1409.5191)とか）

これより、以下のような流れで説明をしていきます。

**2. self-inverse transformの詳細釣り合い**で、ある広いクラスに属する変換 $x' = \Phi(x)$ が詳細釣り合い条件を満たすことを示します。この変換はdeterministicなものですが、**3. 補助変数$r$の導入と周辺分布$p(x)$の詳細釣り合い**で、ノイズ $r$ を導入することで $x$ の確率的な遷移についても詳細釣り合いを示せることを見て、一般的なHMCを定義します。**4-5**で教科書に載っている普通のハミルトニアンモンテカルロが**3**で定義したHMCの一例になっていることを見ます。最後に**6. Non-Volume Preserving HMC**で、「普通でない」HMCの事例を少し見てみることにします。

### 2. self-inverse transformの詳細釣り合い

全域で微分可能な全単射 $\Phi: X \to X$ （ここで $X$ は $\mathbb{R}^d$ など $x$ の住んでいる空間）で、$\Phi = \Phi^{-1}$ つまり $\Phi(\Phi(x)) = x$ を満たすものを考えます。

確率分布 $p(x) \propto \exp(-S(x))$ に従う変数 $x$ に対して、受理確率 $\alpha(x \to \Phi(x))$ をうまく定義すると、$x$ から $\Phi(x)$ への更新（またはそのまま留まる）ステップの前後で確率分布が一定に保たれることを見ます。

更新後の変数を $y$ と置くと、この分布 $p_y(y)$ は次のようになります。

$$
\begin{equation}
p_y(y) = p(x)\left|\frac{\partial x}{\partial y}\right|\alpha(x \to y) + p(y)(1 - \alpha(y\to x))
\end{equation}
$$

ただし $x := \Phi^{-1}(y) = \Phi(y)$ と置きました。一項目が $x$ から $y$ へと遷移してきた寄与で、二項目がもともと $y$ にいたのが留まった分の寄与です。一項目には確率密度変数の変数変換なのでヤコビアンが現れています。

$$
\begin{equation}
\left|\frac{\partial x}{\partial y}\right|
= \left|\frac{\partial \Phi^{-1}(y)}{\partial y}\right|
= \left|\frac{\partial x}{\partial \Phi(x)}\right|
\end{equation}
$$

ここで、受理確率を次のように置いてみます。

$$
\begin{equation}
\alpha(x \to y) = \min\left\lbrace 1, \frac{p(y)}{p(x)}\left|\frac{\partial y}{\partial x}\right| \right\rbrace
\end{equation}
$$

$\min$ の条件を場合分けして $p_y(y)$ を計算します。

- $\frac{p(y)}{p(x)}\left\lvert\frac{\partial y}{\partial x}\right\rvert<1$ のとき

$$
\begin{equation}
p_y(y) = p(x)\left|\frac{\partial x}{\partial y}\right| \frac{p(y)}{p(x)}\left|\frac{\partial y}{\partial x}\right| + p(y)(1 - 1) = p(y)
\end{equation}
$$

- $\frac{p(y)}{p(x)}\left\lvert\frac{\partial y}{\partial x}\right\rvert>1$ のとき

$$
\begin{equation}
p_y(y) = p(x)\left|\frac{\partial x}{\partial y}\right| + p(y)\left(1 - \frac{p(x)}{p(y)}\left|\frac{\partial x}{\partial y}\right|\right) = p(y)
\end{equation}
$$

従って、ステップ前後の分布$p$, $p_y$が一致しており、釣り合い条件が満たされていることが分かります。

詳細釣り合い条件も同様に示せます。任意の $x, x'$ について

$$
\begin{eqnarray}
p(x)T(x \to x') &=& p(x)\delta(\Phi(x) - x')\alpha(x \to x') \\
&=& p(x)\left|\frac{\partial x}{\partial \Phi(x)}\right|\delta(x - \Phi(x'))\alpha(x \to x') \\
p(x')T(x' \to x) &=& p(x')\delta(\Phi(x') - x)\alpha(x' \to x)
\end{eqnarray}
$$

であり（ここではヤコビアンはデルタ関数の公式から出て来ました）、$x' \neq \Phi(x)$のときはデルタ関数がゼロになるので適宜 $x'$ を $y := \Phi(x)$ で置き換えてしまえば

$$
\begin{eqnarray}
p(x)T(x \to x') &=& p(x)\left| \frac{\partial x}{\partial y} \right|\delta(x - \Phi(x'))\alpha(x \to y) \\
p(x')T(x' \to x) &=& p(y)\delta(\Phi(x') - x)\alpha(y \to x)
\end{eqnarray}
$$

となります。
あとは釣り合い条件とまったく同様に場合分けをすれば

- $\frac{p(y)}{p(x)}\left\lvert\frac{\partial y}{\partial x}\right\rvert<1$ のとき

$$
\begin{eqnarray}
p(x)T(x \to x') &=& p(y)\delta(x - \Phi(x')) \\
p(x')T(x' \to x) &=& p(y)\delta(\Phi(x') - x)
\end{eqnarray}
$$

- $\frac{p(y)}{p(x)}\left\lvert\frac{\partial y}{\partial x}\right\rvert>1$ のとき

$$
\begin{eqnarray}
p(x)T(x \to x') &=& p(x)\left| \frac{\partial x}{\partial y} \right|\delta(x - \Phi(x')) \\
p(x')T(x' \to x) &=& p(x)\left| \frac{\partial x}{\partial y} \right|\delta(\Phi(x') - x)
\end{eqnarray}
$$

といずれのケースも $p(x)T(x \to x') = p(x')T(x' \to x)$ が成り立っていることが分かりました。

ところで議論の先取りになりますが、self-inverse transformは随分奇妙で人工的な変換に見えますが、ハミルトン方程式に従う運動の時間発展のような時間反転対称性を持つ変換は自然にこのような $\Phi$ と見なすことができ（詳細は後述）、これが一般化したHMCがMCMCとして最低限機能する原理になっています。
エネルギー保存などハミルトン方程式の持つ他の性質は実用上は重要ですが、MCMCの条件を満たす上では関係ありません。（ただし繰り返すように、 $\Phi(\Phi(x))$ は $x$ に戻り、エネルギーを含め全ての値が元通りになる必要があります）

逆に、時間反転対称性を持たない運動、例えば摩擦がある場合は時間発展するたびに単調にエネルギーを失っていって元の状態 $x$ に戻ることはないので、こうした $\Phi$ と見なすことはできません。変換を繰り返すと一つの静止状態に落ち着いてしまうので、サンプリングには使えないということは直観的にも納得できるかと思います。

### 3. 補助変数$r$の導入と周辺分布$p(x)$の詳細釣り合い

$\Phi$ と $\alpha$ で定まる更新ルールが $p(x)$ に関する詳細釣り合い条件を満たすことは分かりましたが、これだけでMCMCの条件を満たすわけではありません。なぜなら任意の $x$ は $\Phi(x)$ との間を往復するだけで、既約性を満たさないからです。$x$ がどこにでも行けるようにするために、ノイズの役割を果たす補助変数を導入します。

サンプリングしたい変数 $x$ とは独立の確率分布既知の変数 $r \sim p_r(r)$ を導入し、同時分布 $p(x,r) = p(x) p_r(r)$を考えます。
ここでエネルギー関数 $H(x, r)$ を次のように定義すると、

$$
\begin{eqnarray}
H(x, r) &=& -\log p(x) - \log p_r(r) \\
&=& S(x) + T(r)
\end{eqnarray}
$$

同時分布は $p(x, r) \propto \exp(-H(x, r))$ のように書くことができます。

さて、2での議論を $x$ から $\lbrace x, r \rbrace$ に置き換えれば、self-inverseな変換 $\lbrace x', r' \rbrace = \Phi(x, r)$ に対して同時分布 $p(x, r)$ について詳細釣り合いを満たすような受理確率は次のように与えられます。

$$
\begin{equation}
\alpha(\lbrace x, r \rbrace \to \lbrace x', r' \rbrace)
= \min\left\lbrace 1, \frac{p(x', r')}{p(x, r)}\left|\frac{\partial (x', r')}{\partial (x, r)}\right| \right\rbrace
\end{equation}
$$

エネルギー関数を使って書けば次のようになります。

$$
\begin{equation}
\alpha(\lbrace x, r \rbrace \to \lbrace x', r' \rbrace)
= \min\left\lbrace 1, \exp(-H(x',r') + H(x,r))\left|\frac{\partial (x', r')}{\partial (x, r)}\right| \right\rbrace
\end{equation}
$$

ここで興味のある変数 $x$ のみに注目して、「ランダムに $r\sim p_r(r)$ をサンプルして $\Phi$ と $\alpha$ で更新ステップを実行してから、 $x'$ だけを取り出す」というプロセスを考えると、 $x$ の周辺分布に詳細釣り合い条件が成立することが示せます。

2での詳細釣り合いの証明と同じ方針で、 $x, x'$ 間の流れを考えると

$$
\begin{eqnarray}
p(x) T(x \to x') 
&=& p(x) \int dr dr' p_r(r) \delta(\Phi(x,r) - \{x', r'\})\alpha(\{x, r\} \to \{x', r'\}) \\
&=& p(x) \int dr dr' p_r(r) \left|\frac{\partial (x, r)}{\partial \Phi(x, r)}\right| \delta(\{x, r\} - \Phi(\{x', r'\}))\alpha(\{x, r\} \to \{x', r'\}) \\
p(x') T(x' \to x) 
&=&  p(x') \int dr' dr p_r(r') \delta(\Phi(x', r') - \{x, r\}) \alpha(\{x', r'\} \to \{x, r\})
\end{eqnarray}
$$

となり、次に受理確率の $\min$ について条件分けすると

- $\frac{p(x', r')}{p(x, r)}\left\lvert\frac{\partial (x', r')}{\partial (x, r)}\right\rvert < 1$ のとき

$$
\begin{eqnarray}
p(x) T(x \to x')
&=& p(x) \int dr dr' p_r(r) \left|\frac{\partial (x, r)}{\partial \Phi(x, r)}\right| \delta(\{x, r\} - \Phi(\{x', r'\}))\frac{p(x', r')}{p(x, r)}\left\lvert\frac{\partial (x', r')}{\partial (x, r)}\right\rvert \\
&=& p(x') \int dr dr' p_r(r')\delta(\{x, r\} - \Phi(\{x', r'\}))\\
p(x') T(x' \to x) 
&=&  p(x') \int dr' dr p_r(r') \delta(\Phi(x', r') - \{x, r\})
\end{eqnarray}
$$

- $\frac{p(x', r')}{p(x, r)}\left\lvert\frac{\partial (x', r')}{\partial (x, r)}\right\rvert > 1$ のとき

$$
\begin{eqnarray}
p(x) T(x \to x')
&=& p(x) \int dr dr' p_r(r) \left|\frac{\partial (x, r)}{\partial \Phi(x, r)}\right| \delta(\{x, r\} - \Phi(\{x', r'\})) \\
p(x') T(x' \to x) 
&=&  p(x') \int dr' dr p_r(r') \delta(\Phi(x', r') - \{x, r\}) \frac{p(x, r)}{p(x', r')}\left\lvert\frac{\partial (x, r)}{\partial (x', r')}\right\rvert\\
&=& p(x) \int dr dr' p_r(r) \left|\frac{\partial (x, r)}{\partial \Phi(x, r)}\right| \delta(\Phi(x', r') - \{x, r\})
\end{eqnarray}
$$

で、いずれの場合も $p(x)T(x \to x') = p(x')T(x' \to x)$ が成立していることが分かりました。

これで、$p_r, \Phi, \alpha$ によって決まる $x$ の更新が、目的分布 $p(x)$ に関する詳細釣り合いを満たしていることが分かりました。**2**での詳細釣り合いとは異なり、 $x$ は確率的に色々な $x'$ に飛ぶことができるのが重要です。これで、一般化された形でのHMCのアルゴリズムを導出することができました。

#### アルゴリズム: 一般化されたHMC
> $x$ に適当な初期値を設定し、以下を繰り返す
> 1. $r \sim p_r(r)$ をサンプルする
> 2. self-inverse transform $\Phi$ によって、 $\lbrace x', r' \rbrace = \Phi(x, r)$ を計算する
> 3. $s \sim \mathcal{U}[0, 1]$ をサンプルする
> 4. $s < \alpha(\lbrace x, r \rbrace \to \lbrace x', r' \rbrace)$ ならば、 $x$ を $x'$ に更新する。そうでなければ棄却する
> 5. 現在の値 $x$ をアンサンブルに加える

ただし、このアルゴリズムがMCMCとして成立するためには既約性が必要で、これは $p_r$ と $\Phi$ の設計に依存します。例えば $\Phi(x,r) = \lbrace \Phi_x(x), \Phi_r(r)\rbrace$ のように $x$ と $r$ をバラバラに変換させてしまうと、 $r$ の導入は $x$ の更新に影響を与えませんので **2** と同じようにdeterministicな変換しか得られません。また、 $r \sim \delta(r)$ のように $r$ にランダムネスがない場合も当然 $\Phi$ によらず $x$ の更新はdeterministicになってしまいます。

以下、 $x, r$ が相互作用するような変換の具体的な実現方法を見ていきます。

### 4. 厳密なHamiltonian dynamics

**3**でエネルギー関数 $H(x, r)$ を導入しましたが、実はこの $H$ は $x, r$ を正準変数とするハミルトニアンだったのだと見なしてしまい、ハミルトン方程式に従う時間発展で $\Phi$ を作ろう、というのが教科書的なHMCです。

$x, r$ を正準変数とするハミルトン方程式は以下で与えられます。

$$
\begin{equation}
\frac{dx}{dt} = \frac{\partial H}{\partial r}\ \ \ \ 
\frac{dr}{dt} = -\frac{\partial H}{\partial x}
\end{equation}
$$

ハミルトニアンに含まれるkinetic term $T(r)$ の決め方には任意性がありますが、補助変数の確率分布としてガウス分布 $p_r(r) = \mathcal{N}(0, I)$ を採用すると、非常に見慣れた運動エネルギーが得られます。

$$
\begin{equation}
T(r) = -\log p_r(r) = \frac{1}{2}r^2 + const.
\end{equation}
$$

ここでは質量は $m=1$ として無視しましたが、rの成分ごとに異なる質量を割り当てても構いません。これはガウス分布の分散を成分ごとに変更することに対応します。

$$
\begin{equation}
T(r) = -\log p_r(r) = \sum_i \frac{r_i^2}{2m_i} + const.
\end{equation}
$$

重たい成分 $i$ では運動エネルギーの分散は $m_i$ に比例して大きくなりますが、速度 $dx_i/dt = r_i/m_i$ は概ね $1/\sqrt{m_i}$ に比例して小さくなります。

ハミルトン方程式の時間反転対称性のお陰でself-inverse transform $\Phi_H$ はシンプルに適当な時間幅 $\tau$ の時間発展で得られます。ただし、 $\Phi = \Phi^{-1}$ とするために、最後に時間を反転させる（つまり運動量ベクトルを反転させる）必要があります。
$\tau$ 秒進んだところで進行方向を逆向きにして、再び $\tau$ 秒進めてまた向きを戻すと、最初の位置と速度に戻るという具合です。

$$
\begin{equation}
\Phi_H(x, r) = \lbrace x(\tau), -r(\tau) \rbrace
\end{equation}
$$

筆者はよく混乱するので補足すると、摩擦があるなど時間反転対称性のない系では、 $\tau$ 秒後に運動量を反転させて再び $\tau$ 秒進めても最初の状態には戻ってきません。

$\Phi$ としては色々考えられる中でHamilton dynamicsを用いる利点は次の２つです。

- エネルギー保存

ハミルトニアンが時刻に陽に依存しないとき、エネルギー $H(x(t),r(t))$ は保存量です。
また、運動エネルギー項が $r \leftrightarrow -r$ について不変であれば、変換 $\Phi_H$ の最後の運動量ベクトルの反転でもエネルギーは変化しません。

$$
\begin{equation}
H(x', r') := H(x(\tau), -r(\tau)) = H(x(\tau), r(\tau)) = H(x, r)
\end{equation}
$$

これにより受理確率の式に現れる $-H(x', r') + H(x, r)$ が常にゼロになります。

- 位相体積保存（リウヴィルの定理）

受理確率の式に現れるヤコビアン行列式 $\left\lvert\frac{\partial \Phi(x,r)}{\partial (x,r)}\right\rvert$ が1になり、計算が簡単になります。

また、この２つを合わせると、受理確率が常に1になることが分かります。

$$
\begin{equation}
\alpha(\lbrace x, r \rbrace \to \lbrace x', r' \rbrace)
= \min\left\lbrace 1, \exp(-H(x',r') + H(x,r))\left|\frac{\partial (x', r')}{\partial (x, r)}\right| \right\rbrace = 1
\end{equation}
$$

発展時間 $\tau$ を大きくすれば $x$ は遠くまで移動することができますが、このときも必ず更新できるというのはMCMCの実用上極めて有用な性質と考えられます。

### 5. 近似的なHamiltonian dynamics

Hamiltonian dynamicsに従うHMCの有用性は分かりましたが、実際には単振動（これは目的分布 $p(x)$ もガウシアンということで、MCMCを使う意味がありません）のような限られた可積分系を除いて時間発展を厳密に計算することはできません。

そこで、時間発展を近似計算することで理想に近いパフォーマンスを得ることを期待します。
この近似計算がMCMCを成立させるために必要な条件は、**3**で強調したように時間反転対称性を持つということです。エネルギー保存や体積保存は望ましい性質ですが、満たす必要はありません。

こうした条件を満たす近似アルゴリズムとしてはleapfrog法がよく知られており、一般的に使われています。leapfrog法の詳細は世の中に解説がたくさんあるので本記事では説明を省略します。

leapfrog法を含むシンプレクティック法では厳密なエネルギー保存は満たされませんが、長時間の積分を行ってもエネルギーの誤差が蓄積せずに一定範囲に収まる性質を持っており、受理確率が低下しづらい恩恵があります。

また、leapfrog法は位相体積を保存する性質も持っています。これによりヤコビアンを考慮する必要がありません。このため、多くのHMCの解説では受理確率は単に次のように書かれています。

$$
\begin{equation}
\alpha(\lbrace x, r \rbrace \to \lbrace x', r' \rbrace)
= \min\left\lbrace 1, \exp(-H(x',r') + H(x,r)) \right\rbrace
\end{equation}
$$

筆者が本記事を執筆したきっかけは、位相体積を保存しない場合にはこの式をどのように補正すればよいのか気になったことでした。

### 6. Non-Volume Preserving HMC

ここまでなるべく一般化された形でHMCを導出し、その自然で効率的な実現例として近似的なHamiltonian dynamicsを用いた教科書的なHMCがあることを見てきました。
また、HMCの成立には時間反転対称性だけが必須で、他の性質はオプショナルであることも見ました。

|   | 時間反転対称性 | エネルギー保存 | 位相体積保存 |
| - | - | - | - |
| 厳密なHamiltonian dynamics  | Yes | Yes | Yes |
| Leapfrog HMC | Yes | Approximate | Yes |
| Non-Volume Preserving HMC | Yes | Approximate | No |

それでは、位相体積を保存しないような非ハミルトン的なHMCというのは考えられるのかというのは自然な疑問だと思います。
そこで簡単に調べてみたところ、こうしたNon-Volume Preserving HMCは最近になって研究が進められている分野であることを知りました。

例えば2021年のAISTATSで[Non-Volume Preserving Hamiltonian Monte Carlo and No-U-TurnSamplers (Afshar et al.)](http://proceedings.mlr.press/v130/mohasel-afshar21a.html)という論文が発表されています。
詳しく読めていませんが、目的分布 $p(x)$ がpiecewise continuousなとき、非連続面で位相体積の変化を許容する代わりにハミルトニアンの誤差を抑えるtransitionを導入することで、受理確率を向上させるテクニックが提案されているようです。
上で書いた一般の $\Phi$ に対する詳細釣り合いの導出はこの論文で与えられている証明を参考にしました。（この論文の著者らはこれまで証明がされていなかった"*a rigorous proof for its correct
convergence has been missing*"と主張していますが、下の論文などを見るに単に常識扱いされていたようにも見えます）

また、2018年にはICLRで[Generalizing Hamiltonian Monte Carlo with Neural Networks](https://openreview.net/forum?id=B1n8LexRZ)という論文が発表されています。ここではより効率的なサンプリングを実現するための $\Phi$ として、ニューラルネットワークの出力に基づいて運動量や座標を変化させるようなleapfrogを拡張した変換を導入しています。この変換では位相体積が保存しないため、本記事で見たのと同様のヤコビアン要素が導入されていて具体的な計算方法が示されています。
なお、この論文は[RealNVP (Dinh et al., ICLR 2017)](https://arxiv.org/abs/1605.08803)に着想を得たものということです。

### 参考文献

本記事の作成には以下の資料を参考にさせていただきました。

- [花田・松浦『ゼロからできるMCMC マルコフ連鎖モンテカルロ法の実践的入門』](https://www.amazon.co.jp/dp/4065201748)

  このテキストでMCMCの勉強をしていて、理解を深めるために今回の記事を書くことにしました。

- Gregory Gundersen, [Hamiltonian Monte Carlo](https://gregorygundersen.com/blog/2020/07/05/hmc/)

  この著者のブログは幅広いトピックを簡潔に解説していて大変参考になります。

- Duane et al. (1987) [Hybrid Monte Carlo](https://www.sciencedirect.com/science/article/abs/pii/037026938791197X)

  HMCを提案した論文です。周辺分布 $p(x)$ の詳細釣り合いの証明はこちらを参考にしました。

- Neal et al. (2012) [MCMC using Hamiltonian dynamics](https://arxiv.org/abs/1206.1901)

  包括的なレビューです。ほぼ読めていません。
